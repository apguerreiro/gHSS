/*************************************************************************

 ghss.h

 ---------------------------------------------------------------------

                       Copyright (c) 2016, 2017
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



*************************************************************************/
#ifndef GHSS_H_
#define GHSS_H_

#ifdef __cplusplus
extern "C" {
#endif

double greedyhss(double *data, int d, int n, const int k, const double *ref, double * volumes, int * selected);

#ifdef __cplusplus
}
#endif

#endif
