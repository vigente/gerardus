/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_connected.h
 * @brief  connected compontent detector
 * @author Martin Bergner
 *
 * The detector will detect block diagonal matrix structures as wells as generalized
 * set partitioning or covering master problems.
 *
 * It works as follows:
 * - It implicitly builds a graph with one vertex for every constraint and edges between constraints that
 *   share a node
 * - All vertices belonging to constraints of the form \f$\sum x_i = a \f$ for
 *   \f$x_i \in \mathbb Z, a\in \mathbb Z\f$ or of the form \f$\sum x_i \geq 1 \f$
 *   for \f$x_i \in \{0,1\} \f$ are removed
 * - The pricing problems correspond to connected components in the remaining graph
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_DEC_CONNECTED_H__
#define GCG_DEC_CONNECTED_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for connected constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeDetectionConnected(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
