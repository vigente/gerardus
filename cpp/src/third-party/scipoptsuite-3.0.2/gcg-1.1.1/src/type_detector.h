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

/**@file   type_detector.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for detectors in GCG projects
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_DETECTOR_H__
#define GCG_TYPE_DETECTOR_H__

#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "type_decomp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DEC_Detector DEC_DETECTOR;
typedef struct DEC_DetectorData DEC_DETECTORDATA;

/**
 * detector initialization method. This method is called when detection is
 * about to begin. It can be used to fill the detector data with needed
 * information. The implementation is optional.
 *
 * input:
 *  - scip            : SCIP data structure
 *  - detector        : detector data structure
 */
#define DEC_DECL_INITDETECTOR(x) SCIP_RETCODE x (SCIP* scip, DEC_DETECTOR* detector)

/**
 * detector deinitialization method. This method is called when the detection
 * is finished. It can be used to clean up the data created in
 * DEC_DECL_INITDETECTOR. The implementation is optional.
 *
 * input:
 *  - scip            : SCIP data structure
 *  - detector        : detector data structure
 */
#define DEC_DECL_EXITDETECTOR(x) SCIP_RETCODE x (SCIP* scip, DEC_DETECTOR* detector)

/**
 * detects the structure of a the problem. This mandatory method is called
 * when the detector should detect the structure.
 *
 * input:
 *  - scip            : SCIP data structure
 *  - detectordata    : detector data data structure
 *  - decdecomps      : a pointer to an array where detected decompositions
 *                      should be saved. The array needs to be created in this
 *                      method.
 *  - ndecdecomps     : pointer where the number of detected decompositions is
 *                      stored
 *  - result          : pointer where to store the result
 *
 * possible return values for result:
 *  - SCIP_SUCCESS    : the method completed and found decompositions
 *  - SCIP_DIDNOTFIND : the method completed without finding a decomposition
 *  - SCIP_DIDNOTRUN  : the method did not run
 */
#define DEC_DECL_DETECTSTRUCTURE(x) SCIP_RETCODE x (SCIP* scip, DEC_DETECTORDATA* detectordata, DEC_DECOMP*** decdecomps, int* ndecdecomps, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
