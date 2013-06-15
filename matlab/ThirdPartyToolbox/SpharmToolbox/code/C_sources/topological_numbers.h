/* 
** Spherical Harmonic Modeling and Analysis Toolkit (SPHARM-MAT) is a 3D 
** shape modeling and analysis toolkit. 
** It is a software package developed at Shenlab in Center for Neuroimaging, 
** Indiana University (email: SpharmMat@gmail.com)
** It is available to the scientific community as copyright freeware 
** under the terms of the GNU General Public Licence.
** 
** Copyright 2009, 2010, ShenLab, Center for Neuroimaging, Indiana University
** 
** This file is part of SPHARM-MAT.
** 
** SPHARM-MAT is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** SPHARM-MAT is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with SPHARM-MAT. If not, see <http://www.gnu.org/licenses/>.
*/

/* This code is based on Florent Ségonne's original implementation available 
at http://people.csail.mit.edu/fsegonne/research/Topology/topology_02.html
*/

#ifndef TOPO_INCLUDED
#define TOPO_INCLUDED
/* This file provides some code to compute the topological numbers
   introduced by Giles Bertrand and Gregoire Malandin.
   
   This code is the result of my understanding of their paper, and is
   not optimized. In fact, it is possible to generate look-up tables
   to greatly speed-up computations...

   Given a binary 3*3*3 neighborhood ( the structure
   TOPOLOGICAL_NEIGHBORHOOD or NBH ) and a pair of consistent
   connectivities ( 1=(6+,18), 2=(18,6+), 3=(6,26), 4=(26,6) ), 
   the function checkTn computes the topological number associated with
   the structure NBH and a specific connectivity. The fucntion
   checkSimple checks if the point is simple or not (see definitions below).

   The topological neighborhood NBH has to be initialized as a binary
   object 0-1, where 1 represents the Foreground Object and 0 the
   Background Object. The first connectivity in the pair of digital
   connectivities ( 1=(6+,18), 2=(18,6+), 3=(6,26), 4=(26,6) ) refers to
   the Foreground Object, the second one to the Background Object
*/

/*
//////////////////////////////////////////////////////////////////////
//     COMPATIBLE CONNECTIVITIES / TOPOLOGICAL NUMBERS 
//
//     TOPOLOGICAL CONVENTION
//     0:            No topological constraint
//     1:            (6+,18)
//     2:            (18,6+)
//     3:            (6,26)
//     4:            (26,6)   
//     default:      (6+,18)
////////////////////////////////////////////////////////////////////
*/


typedef unsigned char TOPOLOGICAL_NEIGHBORHOOD[3][3][3];

#define NBH TOPOLOGICAL_NEIGHBORHOOD 
#define TMP    100    
#define  MAXIMUM_NUMBER_OF_COMPONENTS 10 /* maximum number of components in a NBH */
#define MAX_COMP MAXIMUM_NUMBER_OF_COMPONENTS 

/* Function Declaration */

/*The connectivity number is equal to the maximum allowed distance
  (||.||1 norm) used to generate the topological neighborhood*/
int connectivityNumber(int connectivity);

/* Digital topology requires a pair of compatible connectivities */
int associatedConnectivity(int connectivity);

NBH* reverseNBH(NBH* nbh_src,NBH *nbh_dst);

NBH *N_6_1(NBH* nbh_src,NBH* nbh_dst);

NBH* N_6_2(NBH* nbh_src,NBH* nbh_dst);

NBH* N_6_3(NBH* nbh_src,NBH* nbh_dst);

NBH *N_18_1(NBH* nbh_src,NBH* nbh_dst);

NBH* N_18_2(NBH* nbh_src,NBH* nbh_dst);

NBH *N_26_1(NBH* nbh_src,NBH* nbh_dst);

NBH* Nnk(NBH* nbh_src,NBH *nbh_dst,int connectivity);

/*This function computes the topological number associated with NBH and a certain connectivity.
This function is a "hack": a more elegant/efficient implementation
should be done...
NBH is a binary object (0-1)
connectivity should be equal to 1,2,3, or 4 according to the
convention above*/
int checkTn(NBH *nbh_src,NBH *nbh_dst,int connectivity);

/*This function checks if a point is simple or not
NBH is a binary object (0-1)
connectivity should be equal to 1,2,3, or 4 according to the
convention above*/
int checkSimple(NBH *nbh,int connectivity);

#endif
