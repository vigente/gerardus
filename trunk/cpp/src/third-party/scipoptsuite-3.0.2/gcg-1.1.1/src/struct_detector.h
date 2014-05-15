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

/**@file   struct_detector.h
 * @brief  data structures for detectors
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_DETECTOR_H__
#define GCG_STRUCT_DETECTOR_H__

#include "type_detector.h"

#ifdef __cplusplus
extern "C" {
#endif

/** detector data structure */
struct DEC_Detector {
   const char*           name;               /**< name of the detector */
   DEC_DETECTORDATA*     decdata;            /**< custom data structure of the detectors */
   char                  decchar;            /**< display character of detector */
   const char*           description;        /**< description of the detector */
   int                   priority;           /**< detector priority */
   SCIP_Bool             enabled;            /**< flag to indicate whether detector is enabled */

   DEC_DECL_INITDETECTOR((*initDetection));  /**< initialization method of detector */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)); /**< structure detection method of detector */
   DEC_DECL_EXITDETECTOR((*exitDetection));  /**< deinitialization method of detector */
};

#ifdef __cplusplus
}
#endif

#endif
