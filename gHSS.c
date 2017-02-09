/*************************************************************************

 gHSS - (incremental) greedy hypervolume subset selection in 2D and 3D

 ---------------------------------------------------------------------

                        Copyright (c) 2015-2017
                Andreia P. Guerreiro <apg@dei.uc.pt>
             

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License,
 version 3, as published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA

 ----------------------------------------------------------------------

 Relevant literature:

 [1] A. P. Guerreiro, C. M. Fonseca, and L. Paquete, “Greedy hypervolume subset selection in low dimensions,” Evolutionary Computation, vol. 24, pp. 521-544, Fall 2016.
 [2] A. P. Guerreiro, C. M. Fonseca, and L. Paquete, “Greedy hypervolume subset selection in the three-objective case,” in Proceedings of the 2015 on Genetic and Evolutionary Computation Conference, GECCO '15, (Madrid, Spain), pp. 671-678, ACM, 2015.

*************************************************************************/

#include "gHSS.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>


#if __GNUC__ >= 3
# define __ghss_unused    __attribute__ ((unused))
#else
# define __ghss_unused    /* no 'unused' attribute available */
#endif


FILE * infoFile;


/* ------------------------------------ Data structure ------------------------------------------*/

typedef struct dlnode {
  double x[3];          // Point
  int in;               //True or False - indicates whether the points has been selected (True) or if still left out (False)
  int updated;          //if in == False, then 'updated' indicates whether the contribution of this points was already updated
 
 //current next (for a list of 'in' points that is modified along the execution)
  struct dlnode * cnext[2]; 
  
  //current next and prev (list with some of the 'out' points)
  struct dlnode * cnextout[2]; 
  struct dlnode * cprevout[2];
  
  //global circular doubly linked list - keeps the points sorted according to coordinates 1 to 3 (all points, ie, both 'in' points and 'out' points)
  struct dlnode * next[3];
  struct dlnode *prev[3]; 
  
  //aditional info for computing contributions
  double area;          // area of 2D projections at z = lastSlicez
  double contrib;       // contribution
  double oldcontrib;    //temporary (save last contribution)
  double lastSlicez;    // up to which value of z the contribution is computed
  struct dlnode * replaced; //the point this one replaced (prev in the paper)
  int dom; //is this a dominated point?
  int id;
  
} dlnode_t;




/* ------------------------------------ Print functions ------------------------------------------*/

#if VARIANT < 2
static void printPoint(double * x, int d){
 
    int i;
    for(i = 0; i < d; i++){
        if(x[i] == -DBL_MAX){
           fprintf(infoFile, "-inf\t");
        }else{
            fprintf(infoFile, "%.2f\t", x[i]);
        }
    }
    
    fprintf(infoFile, "\n");
}




static void printPointList(dlnode_t * list, int d, int i){

    dlnode_t * p = list->next[i];
    do{
        printPoint(p->x, d);
        
        p = p->next[i];
        
    }while(p != list);
    
    
}


static void printList(dlnode_t * list, int d, int ix){

    dlnode_t * p = list->next[ix];
    do{
        
        printPoint(p, d);
        
        p = p->next[ix];
        
    }while(p != list);
    
    
}

//ix = 0 ou 1
static void printListIn(dlnode_t * p, int d, int ix){

    dlnode_t * q = p->cnext[1-ix];
    dlnode_t * stop = p->cnext[0];
    do{
        
        fprintf(infoFile, "(%d|%d|%d): ", q->id, q->in, q->updated);
        printPoint(q, d);
        
        q = q->cnext[ix];
        
    }while(q != stop);
        fprintf(infoFile, "(%d|%d|%d): ", q->id, q->in, q->updated);
    printPoint(q);
    
}



//ix = 0 or 1
static void printListOut(dlnode_t * p, int d, int ix){

    dlnode_t * q = p->cnextout[1-ix];
    dlnode_t * stop = p->cnextout[ix];
    do{
        
        fprintf(infoFile, "(%d|%d|%d): ", q->id, q->in, q->updated);
        printPoint(q);
        
        q = q->cnextout[ix];
        
    }while(q != stop);
    fprintf(infoFile, "(%d|%d|%d): ", q->id, q->in, q->updated);
    printPoint(q);
    
}


static void printContrib(dlnode_t * p){
    fprintf(infoFile, "%-16.15g\n", p->contrib);
}


static void printAllContrib(dlnode_t * list, int ix){
 
    dlnode_t * p = list->next[ix]->next[ix];
    
    while(p != list->prev[ix]){
        
        printContrib(p);
        p = p->next[ix];
        
    }
    
}
#endif



/* -------------------------------------- Setup Data ---------------------------------------------*/


static void copyPoint(double * source, double * dest, int d){
    int i;
    for(i = 0; i < d; i++)
        dest[i] = source[i];
}



static int compare_node(const void *p1, const void* p2)
{
    const double x1 = *((*(const dlnode_t **)p1)->x);
    const double x2 = *((*(const dlnode_t **)p2)->x);

    return (x1 < x2) ? -1 : (x1 > x2) ? 1 : 0;
}

static int compare_node2d(const void *p1, const void* p2)
{
    const double x1 = *((*(const dlnode_t **)p1)->x+1);
    const double x2 = *((*(const dlnode_t **)p2)->x+1);

    return (x1 < x2) ? -1 : (x1 > x2) ? 1 : 0;
}


static int compare_node3d(const void *p1, const void* p2)
{
    const double x1 = *((*(const dlnode_t **)p1)->x+2);
    const double x2 = *((*(const dlnode_t **)p2)->x+2);

    return (x1 < x2) ? -1 : (x1 > x2) ? 1 : 0;
}





/*
 * Setup circular double-linked list in each dimension (with two sentinels).
 * Initialize data.
 */

static dlnode_t *
setup_cdllist(double *data, int d, int n)
{
    dlnode_t *head;
    dlnode_t **scratch;
    int i, j;

    head  = malloc ((n+2) * sizeof(dlnode_t));
    head[0].id = -1;
    head[0].in = 1;
    head[0].area = 0; head[0].contrib = 0; head[0].oldcontrib = 0; head[0].lastSlicez = 0;
    head[0].updated = 0; head[0].dom = 0;
    
    for(i = 0; i < d; i++){
        head[0].x[i] = -1;
        head[n+1].x[i] = -1;
    }
    head[n+1].id = -2;
    head[n+1].in = 1;
    head[n+1].area = 0; head[n+1].contrib = 0; head[n+1].oldcontrib = 0; head[n+1].lastSlicez = 0;
    head[n+1].updated = 0; head[n+1].dom = 0;

    for (i = 1; i <= n; i++) {
//         head[i].x = head[i-1].x + d ;/* this will be fixed a few lines below... */
        copyPoint(&(data[(i-1)*d]), head[i].x, d);
        head[i].area = 0;
        head[i].contrib = 0;
        head[i].oldcontrib = 0;
        head[i].lastSlicez = 0;
        head[i].id = i-1;
        head[i].in = 0;
        head[i].updated = 0;
        head[i].dom = 0;
    }

    scratch = malloc(n * sizeof(dlnode_t*));
    for (i = 0; i < n; i++)
        scratch[i] = head + i + 1;

    for (j = d-1; j >= 0; j--) {
        if(j == 2) qsort(scratch, n, sizeof(dlnode_t*), compare_node3d);
        else if(j == 1) qsort(scratch, n, sizeof(dlnode_t*), compare_node2d);
        else qsort(scratch, n, sizeof(dlnode_t*), compare_node);
        head->next[j] = scratch[0];
        scratch[0]->prev[j] = head;
        for (i = 1; i < n; i++) {
            scratch[i-1]->next[j] = scratch[i];
            scratch[i]->prev[j] = scratch[i-1];
        }
//         scratch[n-1]->next[j] = head;
//         head->prev[j] = scratch[n-1];
        scratch[n-1]->next[j] = head+n+1;
        (head+n+1)->prev[j] = scratch[n-1];
        (head+n+1)->next[j] = head;
        head->prev[j] = (head+n+1);
    }

    free(scratch);

    return head;
}




/* -------------------------------------- Misc ----------------------------------------------*/



static inline double max(double a, double b){
    return (a > b) ? a : b;
}



static inline double min(double a, double b){
    return (a < b) ? a : b;
}



/* -------------------------------------- Algorithms ----------------------------------------------*/



static void updateVolume(dlnode_t * p, double z){
    
    p->contrib += p->area * (z - p->lastSlicez);
    p->lastSlicez = z;
    
}





/*
 * Find maximum contributor and, since all points are visited, also set the
 * 'updated' flag of all points to false as they may have to be updated
 */
static dlnode_t * maximumOutContributor(dlnode_t * list){

    double c = -DBL_MAX;
    dlnode_t * maxp = list;
    dlnode_t * p = list->next[0];
    dlnode_t * stop = list->prev[0];
    while(p != stop){
        p->updated = 0;
        if(!p->in && ((p->contrib > c) || (p->contrib == c && p->id < maxp->id))){
            c = p->contrib;
            maxp = p;
        }
        p = p->next[0];
    }
    return maxp;
}



/*
 * p - the newly selected point
 * xi, yi, zi - indexes of the coordinates
 * ref - reference point
 * 
 * 
 * (Step 1 of Algorithms 3 and 4 in the paper)
 * Set up the list that will be used to maintain the points already selected that are not dominated in
 * the (xi,yi)-projection at a given value of coordinate zi (here is assumed to be p->x[zi], i.e.,
 * this function sets up the list of delimiters of p at z = p->x[zi]).
 * 
 * Note: The first and last elements of such list will be stored in p->cnext[0] (the rightmost
 * delimiter and below) and in p->cnext[1] (delimiter above and to the left)
 * 
 */
static void createFloor(dlnode_t * list, dlnode_t * p, int xi, int yi, int zi, const double * ref){
    
    dlnode_t * q = list->prev[yi];
    
    //set up sentinels
    list->x[xi] = ref[xi];
    list->x[yi] = -DBL_MAX;
    list->x[zi] = -DBL_MAX;
    
    q->x[xi] = -DBL_MAX;
    q->x[yi] = ref[yi];
    q->x[zi] = -DBL_MAX;
    
    
    dlnode_t * xrightbelow = list;
    q = list->next[yi];
    
    
    //find the closest point to p according to the x-coordinate that has lower or equal yi- and zi- coordinates (xrightbelow)
    while(q->x[yi] <= p->x[yi]){
        if(q->in && q->x[zi] <= p->x[zi] && q->x[xi] <= xrightbelow->x[xi] && q->x[xi] > p->x[xi])
            xrightbelow = q;
        
        q = q->next[yi];
    }
    
    //the rightmost delimiter of p area to the right
    p->cnext[0] = xrightbelow;
    
    dlnode_t * last = xrightbelow;
    
    q = p->next[yi];
    
    //set up the list (using cnext)
    while(!q->in || q->x[xi] > p->x[xi] || q->x[zi] > p->x[zi]){
        
        if(q->in && q->x[zi] <= p->x[zi] && q->x[xi] < last->x[xi] && q->x[xi] > p->x[xi]){
            
            if(q->x[yi] == last->x[yi]){
                last = last->cnext[0];
            }
            q->cnext[0] = last;
            last->cnext[1] = q;
            last = q;
            
        }
        q = q->next[yi];
    }
    
    //the delimiter of p area above and to the left
    q->cnext[0] = last;
    last->cnext[1] = q;
    p->cnext[1] = q;
    
    
}




/* Compute the area exclusive dominated by p 
 * (the area is divided in horizontal bars and their areas are summed up)
 */
static double computeArea(dlnode_t * p, int xi, int yi){

    dlnode_t * q = p->cnext[0];
    dlnode_t * qnext = q->cnext[1];
    
    double area = (q->x[xi] - p->x[xi]) * (qnext->x[yi] - p->x[yi]);
    
    q = qnext;
    while(q != p->cnext[1]){
        qnext = q->cnext[1];
        area += (q->x[xi] - p->x[xi]) * (qnext->x[yi] - q->x[yi]);
        q = qnext;
    }
    
    return area;
    
}



/*
 * p - is the newly selected point
 * xi, yi, zi - indexes of the coordinates
 * 
 * (steps 2 and 3 of Algorithm 4)
 * cnextout and cprevout will be used to maintain the list of 'out' points that are being updated.
 * The base area for the points in that list is computed.
 * 
 * Note: The first and last elements of the list of 'out' will be stored in p->cnextout[0] and in p->cnextout[1]
 * 
 */
static void createAndInitializeBases(dlnode_t * list, dlnode_t * p, int xi, int yi, int zi){
    
    dlnode_t * q;
    dlnode_t * stop = p->cnext[1];
    double parea = p->area;
    
    
    //q is set to the first point dominated by p in the list of all points
    // (care must be taken with points with equal yi-coordinate to that of p)
    if(p->prev[yi]->x[yi] == p->x[yi]){
        q = p;
        while(q->x[yi] == p->x[yi]) q = q->prev[yi];    //deal with points with equal yi-coordinate to p but that are before p in the list
        q = q->next[yi];
    }else{
        q = p->next[yi];
    }
    
    dlnode_t * in = p->cnext[0];                        //'in' keeps track of the last 'in' point visited
    double area = 0;
    dlnode_t * out = list;                              // list is used as sentinel

    
    p->cnextout[0] = list; 
    
    
    //setup the list of 'out' points that have to be updated
    //and do the first part of the computation of their base area
    while(q != stop){
        if(q != p){                                     // if p->prev[yi]->x[yi] == p->x[yi], then p will be visited in this while loop and has to be skipped

            if(q->dom){
                q->updated = 1;
            
            }else if(q->in == 0){                       // q is out

                if(q->updated == 0){                    // q has to be updated
                    if(p->x[xi] <= q->x[xi] && p->x[yi] <= q->x[yi] && p->x[zi] <= q->x[zi]){
                        // q is dominated by p then, its contribution is reduced to 0

                        q->area = 0;
                        q->contrib = 0;
                        q->oldcontrib = 0;
                        q->updated = 1;
                        q->dom = 1;
                            
                    }else if(p->x[xi] <= q->x[xi] && p->x[yi] <= q->x[yi] 
                        && (in->x[xi] > q->x[xi] || in->x[yi] > q->x[yi])
                        && (in->cnext[1]->x[xi] > q->x[xi] || in->cnext[1]->x[yi] > q->x[yi])){//check if the contribution of q will be reduced because of p
                        
                        q->oldcontrib = q->contrib;
                        q->contrib = 0;
                        q->area = parea - (in->x[xi] - q->x[xi])*(q->x[yi] - p->x[yi]) - area;
                        q->lastSlicez = p->x[zi];
                        out->cnextout[1] = q;
                        q->cprevout[1] = out;
                        out = q;
                    }
                }
                
            }else{ // q is in

                if(q == in->cnext[1]){                  //if q is a delimiter of p
                    area += (in->x[xi] - q->x[xi]) * (q->x[yi] - p->x[yi]);
                    in = q;
                }
            }
        }
        q = q->next[yi];
    }
    
    q = list->prev[yi];
    out->cnextout[1] = q;
    q->cprevout[1] = out;
    p->cnextout[1] = q;
    out = q;

    
    
    area = 0;
    in = p->cnext[1];
    stop = p->cnext[0];
    
    if(p->prev[xi]->x[xi] == p->x[xi]){
        q = p;
        while(q->x[xi] == p->x[xi]) q = q->prev[xi];
        q = q->next[xi];
    }else{
        q = p->next[xi];
    }
    
    
    //do the second part of the computation of the base area of the 'out' points to be updated
    while(q != stop){
        if(q != p){

            if(q->dom){

                q->updated = 1;
            
            }else if(q->in == 0){           //q is out
                if(q->updated == 0){        //q has to be updated (its contribution has to be reduced)
                    if(p->x[xi] <= q->x[xi] && p->x[yi] <= q->x[yi] 
                        && (in->x[xi] > q->x[xi] || in->x[yi] > q->x[yi]) 
                        && (in->cnext[0]->x[xi] > q->x[xi] || in->cnext[0]->x[yi] > q->x[yi])){
                    
                        q->area -= (q->x[xi] - p->x[xi])*(in->x[yi] - p->x[yi]) + area;
                        out->cnextout[0] = q;
                        q->cprevout[0] = out;
                        out = q;
                    }
                }
            }else{                          // q is in

                if(q == in->cnext[0]){
                    area += (q->x[xi] - p->x[xi]) * (in->x[yi] - q->x[yi]);
                    in = q;
                }
            }
        }
        q = q->next[xi];
    }
    
    q = list;
    out->cnextout[0] = q;
    q->cprevout[0] = out;
    
    
}

/*
 * p - the newly selected point
 * cutter - a delimiter of p. The area dominated by p and the 'cutter' point is removed from p.
 * xi, yi, zi - indexes of the coordinates
 * xic - indicate the coordinate used for visiting points. Points will be visited in ascending
 *       order of coordinate xi if xic == 0 and of coordinate yi if xic == 1.
 * 
 * The area dominated by p is updated and so is the list of points that delimit the area of p at z = cutter->x[zi].
 * Moreover, the volume and areas of some of the 'out' points below p in zi are updated.  
 */
static double cutOffPartial(dlnode_t * p, dlnode_t * cutter, int xi, int yi, int zi, int xic){
    
    int yic = 1 - xic;
    
    dlnode_t * in = p->cnext[yic];
    dlnode_t * out = p->cnextout[yic];
    dlnode_t * stop;
    
    double area = 0;
    
    while(cutter->x[yi] <= in->x[yi]){
        in = in->cnext[xic];
    }
    dlnode_t * upperLeft = in;
    
    while(out->x[xi] < in->x[xi]){
        out = out->cnextout[xic];
    }
    out = out->cprevout[xic];
 
    stop = p->cnext[yic];

    
    while(in != stop){
        
        if(in->cnext[yic]->x[xi] > out->x[xi] || out->x[xi] < p->x[xi]){

            in = in->cnext[yic];
            area += (in->cnext[xic]->x[xi] - max(in->x[xi], p->x[xi])) * (in->x[yi] - cutter->x[yi]);

        }else{ 
            
            updateVolume(out, cutter->x[zi]);
            
            if(out->x[yi] >= cutter->x[yi]){        // 'out' has no more contribution above z = in->x[zi]
                out->area = 0;
                out->contrib = out->oldcontrib - out->contrib;
   
                //remove points completely updated
                out->updated = 1;
                out->cprevout[xic]->cnextout[xic] = out->cnextout[xic];
                out->cnextout[xic]->cprevout[xic] = out->cprevout[xic];
                out->cprevout[yic]->cnextout[yic] = out->cnextout[yic];
                out->cnextout[yic]->cprevout[yic] = out->cprevout[yic];
                
            }else{

                out->area -= area + (in->x[xi] - out->x[xi]) * (in->cnext[yic]->x[yi] - cutter->x[yi]);
            }
            out = out->cprevout[xic];
        }
     
    }
    
    //insert point 'cutter' as the head of the list (the in points dominated by 'cutter' in the (xi, yi)-plane are implicitly removed)
    p->cnext[yic] = cutter;
    cutter->cnext[xic] = upperLeft;
    upperLeft->cnext[yic] = cutter;
    
    p->area -= area;
    return area;
    
    
}



/*
 * p - the new point that will be added to the set of selected points
 * zi - indicates which is the third coordinate (the remaining ones are deduced in this function)
 * 
 * Note: This function corresponds to Algorithms 3 and 4 of the paper which are done together in this
 * function so as to avoid sweeping all points a second time according to the z coordinate and to
 * avoid repeating some computations.
 * 
 */
static void updateOut(dlnode_t * list, dlnode_t * p, int zi, const double * ref){

    int d = 3;
    int xi = (zi + 1) % d;  //first coordinate
    int yi = 3 - (zi + xi); //second coordinate
    dlnode_t * q = list;
    
    createFloor(list, p, xi, yi, zi, ref);
    p->area = computeArea(p, xi, yi);
    createAndInitializeBases(list, p, xi, yi, zi);
    
    dlnode_t * stop = list->prev[zi];
    stop->x[zi] = ref[zi];
    dlnode_t * domr = list;
    
    q = p->next[zi];

    while(q != stop){

        if(q->in){
            // update the area of p, update domr volume and area (Alg. 3) and do Algorithm 4 (lines 6 to 22) 

            if(q->x[xi] <= p->x[xi] && q->x[yi] <= p->x[yi]){ //q is the last delimiter of p (p has no contribution above q->x[zi])
                break;
                
            }else if(q->x[xi] <= p->x[xi] && q->x[yi] > p->x[yi] && q->x[yi] < p->cnext[1]->x[yi]){ //q is to the left of p
                updateVolume(domr, q->x[zi]);                                                       // (Alg. 3, line 16) 
                domr->area -= cutOffPartial(p, q, xi, yi, zi, 0);                             // (Alg. 3, line 10 and 17) (Alg. 4, line 12 - 14)
                
            }else if(q->x[xi] > p->x[xi] && q->x[yi] <= p->x[yi] && q->x[xi] < p->cnext[0]->x[xi]){ //q is below p
                updateVolume(domr, q->x[zi]);                                                       // (Alg. 3, line 16)
                domr->area -= cutOffPartial(p, q, yi, xi, zi, 1);                             // (Alg. 3, line 14 and 17) (Alg. 4, line 20 - 22)
            
            }
            
            
        }else{                                                                       //q is an 'out' point
            if(q->dom == 0 && q->updated == 0 && q->x[xi] <= p->x[xi] && q->x[yi] <= p->x[yi]){     //q* < p* (Alg. 3, lines 19 - 24)
                
                q->oldcontrib = q->contrib;
                q->contrib = 0;
                q->replaced = domr;
                updateVolume(domr, q->x[zi]);
                q->area = p->area;
                q->lastSlicez = q->x[zi];
                domr = q;
                
            }
            
        }
        q = q->next[zi];
        
    }
    
    //(Alg. 3, lines 25 - 30)
    updateVolume(domr, q->x[zi]);
    double vol = 0;
    while(domr != list){
        vol += domr->contrib;
        domr->contrib = domr->oldcontrib - vol;
        domr->updated = 1;
        domr = domr->replaced;
    }
    
    dlnode_t * q2;
    q2 = p->cnextout[0]->cnextout[1];
    
    while(q2 != p->cnextout[1]){
            updateVolume(q2, q->x[zi]);
            q2->contrib = q2->oldcontrib - q2->contrib;
            q2->updated = 1;
        q2 = q2->cnextout[1];
    }
    
    
    
}


static void gHSS3D(dlnode_t * list, const int k, int * selected, const double * ref){
    
    int i;
    dlnode_t * maxp = NULL;
    dlnode_t * p = list->next[0];
    dlnode_t *stop = list->prev[0];
    while(p != stop){
        if(p->dom)
            p->contrib = 0; //if p does not strongly dominate the reference point
        else
            p->contrib = (ref[0] - p->x[0]) * (ref[1] - p->x[1]) * (ref[2] - p->x[2]);
        p = p->next[0];
    }
    
    
    for(i = 0; i < k-1; i++){
        
        maxp = maximumOutContributor(list);
        if(maxp->dom == 0){
            //update contribution of the points not yet selected (out points)
            updateOut(list, maxp, 2, ref); // order (x,y,z)
            updateOut(list, maxp, 1, ref); // order (z,x,y)
            updateOut(list, maxp, 0, ref); // order (y,z,x)
        }
        
        selected[i] = maxp->id;
        maxp->in = 1;   // point 'maxp' is now part of the set of selected points
        
        
    }
    
    maxp = maximumOutContributor(list);
    selected[i] = maxp->id;
    maxp->in = 1;
    
}




static void gHSS2D(dlnode_t * list, const int k, int * selected, const double * ref){
    
    int i;
    dlnode_t * maxp = NULL;
    dlnode_t * p = list->next[0];
    dlnode_t * q = list;
    dlnode_t *stop = list->prev[0];
    dlnode_t * rightin, * upin;
    
    // set sentinels
    list->x[0] = -DBL_MAX;
    list->x[1] = ref[1];
    
    stop->x[0] = ref[0];
    stop->x[1] = -DBL_MAX;
    
    // setup list with cnext, excluding dominated points
    while(p != stop){
        // q is dominated
        if(p->x[0] == q->x[0] && q->x[1] >= p->x[1]){
            q->dom = 1;
            q->contrib = 0;
            q->cnext[1]->cnext[0] = p;
            p->cnext[1] = q->cnext[1];
            p->contrib = (ref[0] - p->x[0]) * (ref[1] - p->x[1]);
            q = p;
        //p is dominated
        }else if(p->x[1] >= q->x[1] || p->dom){
            p->dom = 1;
            p->contrib = 0;
        }else{
            p->contrib = (ref[0] - p->x[0]) * (ref[1] - p->x[1]);
            q->cnext[0] = p;
            p->cnext[1] = q;
            q = p;
        }
        p = p->next[0];
    }
    q->cnext[0] = stop;
    stop->cnext[1] = q;
    
    // greedy subset selection in 2D
    for(i = 0; i < k-1; i++){
        maxp = maximumOutContributor(list); //find the point that contributes the most to the already selected points
        if(maxp->dom == 0){
            upin = maxp->cnext[1];

            while(!upin->in) upin = upin->cnext[1];
                
            rightin = maxp->cnext[0];
            while(!rightin->in) rightin = rightin->cnext[0];
            
            p = maxp->cnext[0];
            while(p != rightin){
                p->contrib -= (rightin->x[0]-p->x[0]) * (upin->x[1]-maxp->x[1]);
                p = p->cnext[0];
            }
            
            p = maxp->cnext[1];
            while(p != upin){
                p->contrib -= (rightin->x[0]-maxp->x[0]) * (upin->x[1]-p->x[1]);
                p = p->cnext[1];
            }
            
        }else{
            maxp->contrib = 0;
        }
        
        selected[i] = maxp->id;
        maxp->in = 1;           // point 'maxp' is in now part of the set of selected points
    }
    
    // no need to update the data structure after selecting the k-th point
    maxp = maximumOutContributor(list);
    selected[i] = maxp->id;
    maxp->in = 1;
    
}


/*
 * mark and initialize the points that do not strongly dominate the reference point and return
 * how many of such points exist  
 */
int markInvalidPoints(dlnode_t * list, int d, const double * ref){
    
    int di;
    dlnode_t * p;
    dlnode_t * stop;
    
    int nmarked = 0;
    
    for(di = 0; di < d; di++){
        
        p = list->prev[di]->prev[di];
        stop = list;
        while(p != stop && p->x[di] >= ref[di]){
            if(p->dom == 0){
                p->dom = 1;
                p->contrib = 0;
                p->oldcontrib = 0;
                p->area = 0;
                nmarked++;
            }
            p = p->prev[di];
        }
    }
    
    return nmarked;
    
}


/* Input:
 * data - array containing all 3D points
 * n - number of points
 * k - subset size (select the k most promising points, one at a time)
 * ref - reference point
 * 
 * Output:
 * the total volume of the subset selected is returned
 * 'contribs' - the contribution of the selected points at the time their
 *              were selected (ex.: contribs[i] holds the contribution of
 *              the i-th selected point)
 * 'selected' - the index of the selected points regarding their order
 *              in 'data' (ex.: selected[i] holds the index of the i-th
 *              selected point. selected[i] holds a value in the range [0,...,n-1]) 
 */
double greedyhss(double *data, int d, int n, const int k, const double *ref, double * contribs, int * selected)
{
    
    
    infoFile = stdout;
    double totalhv = 0;

    dlnode_t *list;

    list = setup_cdllist(data, d, n);
    
    int nmarked = markInvalidPoints(list, d, ref);
    if(nmarked == n){
        int i;
        for(i = 0; i < n; i++){
            selected[i] = i;
            contribs[i] = 0;
        }
        free(list);
        return 0;
    }
    
    
    if (d == 2){
        gHSS2D(list, k, selected, ref);
    }else if(d == 3){
        gHSS3D(list, k, selected, ref);
    }else{
        free(list);
        return -1;
    }
    int * sel2idx = (int *) malloc(n * sizeof(int));
    dlnode_t * p = list->next[0];
    dlnode_t * stop = list->prev[0];
    int i = 0;
    
    for(i = 0; i < n; i++){
        sel2idx[i] = n;
    }
    for(i = 0; i < k; i++){
        sel2idx[selected[i]] = i;
    }
    
    
    while(p != stop){
        if(sel2idx[p->id] < k){
            contribs[sel2idx[p->id]] = p->contrib;
            totalhv += p->contrib;
        }
        p = p->next[0];
        i++;
    }
    

    free(sel2idx);
    free(list);
    
    return totalhv;
}




