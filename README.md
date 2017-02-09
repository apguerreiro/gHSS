gHSS
=====

This software implements algorithms for the incremental greedy approximation to the Hypervolume Subset Selection Problem (HSSP) in two and three dimensions (gHSS2D and gHSS3D [1], respectively), which provide an approximation guarantee of 1-1/e. Minimization is assumed.


**Note**: Although only *nondominated* points that strongly dominate the reference point contribute to the value of the hypervolume indicator, the code is prepared to deal with all other points, including *repeated* points. Warnings will be raised if any of the points do not strongly dominate the reference point. Moreover, different points may have (some) equal coordinates.


**Note 2**: In two dimensions, the exact solution can also be computed efficiently, see [here](https://eden.dei.uc.pt/~paquete/HSSP/) and [here](http://hpi.de/friedrich/docs/code/ssp.zip).

License
--------


Except where indicated otherwise in individual source files, this software is Copyright © 2016, 2017 Andreia P. Guerreiro.

This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Appropriate reference to this software should be made when describing research in which it played a substantive role, so that it may be replicated and verified by others. The algorithms which this software implements are described in detail in [1]. 





Building
--------


In GNU/Linux, the program can be compiled from source by invoking:

    make

It is recommended that you compile it specifically for your architecture. Depending on the compiler and version of the compiler you use, there are different ways to achieve this. For recent GCC versions, make will pick a suitable -march argument based on the processor of the build machine. This can be overridden by passing a `MARCH=` argument to make. Similarly if you use the Intel C compiler, it will pick a sensible default architecture (`-xHOST`) for you. If you want to override this, pass `XARCH=` to `make`. So, to build for an Intel *Core2* with the GCC compiler, you would use:

    make MARCH=core2

For the Intel C compiler, you should use:

    make XARCH=SSSE3

Generally, `make` will try to pick good flags for you, but you can override them by passing an `OPT_CFLAGS` argument to `make` if you wish. To build an unoptimized version of `hv4d` you could run:

    make OPT_CFLAGS="-O0 -g"

Finally, if you do not want to see the command line of each compiler invocation, pass `S=1` to make.



General Usage
-------------


**SYNOPSIS** 

    gHSS [OPTIONS] [FILE...]
    
        
**DESCRIPTION**

Compute the incremental greedy approximation to the hypervolume subset selection problem for the data set(s) in FILE(s).

With no FILE, read from the standard input.

**COMMAND LINE OPTIONS**

	 -h, --help          print this summary and exit.                          
	     --version       print version number and exit.                        
	 -v, --verbose       print some information (time, coordinate-wise maximum and minimum, etc)                                     
	 -q, --quiet         print only the results (as opposed to --verbose). 
	 -u, --union         treat all input sets within a FILE as a single set.   
	 -r, --reference=POINT use POINT as the reference point. POINT must be within quotes, e.g.,
		                 "10 10 10". If no reference point is given, it is taken as the
		                 coordinate-wise maximum of all input points.                                     
	 -s, --suffix=STRING Create an output file for each input file by appending this suffix.
		                 This is ignored when reading from stdin. If missing, output is sent
		                 to stdout.             
	 -k, --subsetsize=k  select k points (a value between 1 and n, where n is the size of the
		                 input data set. The default is n/2)   
	 -f, --format=(0|..|4) output format
		                 (0: print indices followed by the hypervolume indicator of the selected subset (default))        
		                 (1: print indices of the selected points)             
		                 (2: print the hypervolume indicator of the selected subset)    
		                 (3: print indices and the corresponding contributions to the previous subset)
		                 (4: print indices and the corresponding accumulated hypervolume)           
		                        
                               

       

Usage
-----

**Run**

The program reads sets of points from the file(s) specified in the command line:

    ./gHSS data

or standard input:

    cat data | ./gHSS

In input files, each point is given in a separate line, and point coordinates in each line are separated by whitespace. An empty line, or a line beginning with a  hash sign (#), denotes a separate set.


Sets in an input file may be treated as a single set by using option `-u`: 

    ./gHSS data -u


The reference point can be set by giving option `-r`.

    ./gHSS -r "10 10 10 10" data

 If no reference point is given, the default is the coordinate-wise maximum of all input points in all files.

For the other options available, check the output of `./gHSS --help`.
 
By default, the default subset size *k* is set to *n/2* where *n* is the data set size. The subset size can be explicitly specified using option `-k`:
 
    ./gHSS -r "10 10 10 10" data -k 10

For the other options available, check the output of `./gHSS --help`.



Examples
-------


**Input File(s)**

Empty lines and lines starting with a hash sign (#) at the beginning and/or at the end of the file are ignored.

Example of valid content of input files:

    1   1   4   4
    4   4   1   1
    
    2   2   3   3
    3   3   2   2

Another example:

    #
    6 4 9 6
    7 3 7 1
    8 2 3 2
    9 1 2 4
    #
    1 9 5 4
    2 8 3 1
    #
    3 7 8 1
    4 6 3 9
    5 5 5 8
    #
    
       
**Compilation**

Example of basic compilation: 

    make march=corei7

    
**Execution**

The code can be tested with any of the data sets provided in folder `examples`. File `test.inp` in that folder contains a 3-dimensional example with the following set of *10* points:

    0.16 0.86 0.47 
    0.66 0.37 0.29 
    0.79 0.79 0.04 
    0.28 0.99 0.29 
    0.51 0.37 0.38
    0.92 0.62 0.07 
    0.16 0.53 0.70 
    0.01 0.98 0.94 
    0.67 0.17 0.54 
    0.79 0.72 0.05 



After compilation, to compute the greedy solution for a subset of size 5 given the reference point `(1, 1, 1)`, run:

    ./gHSS examples/test.inp -r "1 1 1" -k 5

which produces the following output:

    4
    6
    8
    9
    1
    0.304494   


In the *default* output format (`-f 0`), the first *k* lines contain the indices (between *0* and *n-1*) of the selected points in the order in which they were selected. In the above example, the first point selected was the (4+1)-th point in the file, i.e., `(0.51, 0.37, 0.38)`, the second point selected was the (6+1)-th point in the file, i.e., `(0.16, 0.53, 0.70)`, and so forth. The last line contains the hypervolume indicator of the greedy solution, i.e., of the subset containing the *k* selected points. If option `-f 1` is given, the outptut contains only the indices (i.e, the first *k* lines in option `-f 0`). If option `-f 2` is given instead, only the hypervolume indicator of the selected subset is printed.

If option `-f 3`is given:


    ./gHSS examples/test.inp -r "1 1 1" -k 5 -f 3
    
the corresponding output is:


    4       0.191394       
    6       0.04935        
    8       0.03036        
    9       0.019404       
    1       0.013986



The output contains *k* lines with two columns. The first column is the same as for option `-f 1`. The second column shows the contribution of the corresponding point to the subset containing the points previously selected. In the example above, the first point selected was point `(0.51, 0.37, 0.38)`, and its contribution (to the empty set) was `0.191394`. The last point selected (i.e., point `(0.66, 0.37, 0.29)`) contributed `0.013986` to the hypervolume indicator of the previous subset of *k-1* points. Note that, since each point selected contributes the most to the previous subset, the second column is always in nonincreasing order.

Finally, if option `-f 4` is given, i.e.,

    ./gHSS examples/test.inp -r "1 1 1" -k 5 -f 4

The corresponding output is:

    4       0.191394       
    6       0.240744       
    8       0.271104       
    9       0.290508       
    1       0.304494   

In this case, the second column shows the accumulated hypervolume. Consequently, the second column of the *i*-th line shows the hypervolume indicator of the subset of the first *i* points selected.



References
----------


[1] A. P. Guerreiro, C. M. Fonseca, and L. Paquete, “Greedy hypervolume subset selection in low dimensions,” Evolutionary Computation, vol. 24, pp. 521-544, Fall 2016. [ [DOI](http://dx.doi.org/10.1162/EVCO_a_00188) | [pdf](https://eden.dei.uc.pt/~cmfonsec/GreedyHSS-ECJ2016-authorVersion.pdf) ] 

[2] A. P. Guerreiro, C. M. Fonseca, and L. Paquete, “Greedy hypervolume subset selection in the three-objective case,” in Proceedings of the 2015 on Genetic and Evolutionary Computation Conference, GECCO '15, (Madrid, Spain), pp. 671-678, ACM, 2015. EMO Track Best Paper Award. [ [DOI](http://dx.doi.org/10.1145/2739480.2754812) ]



 



