/*
 * MatlabImageHeader.h
 *
 * Class to provide an interface to read the metadata from a an input
 * image. This allows to create a header class that can be used to
 * deal with the image more easily.
 *
 */

 /*
  * Author: Ramon Casero <rcasero@gmail.com>
  * Copyright © 2012 University of Oxford
  * Version: 0.1.0
  * $Rev$
  * $Date$
  *
  * University of Oxford means the Chancellor, Masters and Scholars of
  * the University of Oxford, having an administrative office at
  * Wellington Square, Oxford OX1 2JD, UK. 
  *
  * This file is part of Gerardus.
  *
  * This program is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * (at your option) any later version.
  *
  * This program is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details. The offer of this
  * program under the terms of the License is subject to the License
  * being interpreted in accordance with English Law and subject to any
  * action against the University of Oxford being under the jurisdiction
  * of the English Courts.
  *
  * You should have received a copy of the GNU General Public License
  * along with this program.  If not, see
  * <http://www.gnu.org/licenses/>.
  */

#ifndef MATLABIMAGEHEADER_H
#define MATLABIMAGEHEADER_H

class MatlabImageHeader {

 public:

  mxArray *data; // pointer to the image voxels
  mxClassID type; // pixel type
  std::vector<unsigned int> size;
  std::vector<double> spacing, origin;

  MatlabImageHeader(const mxArray *arg, std::string paramName);

  // get number of dimensions of the image
  size_t GetNumberOfDimensions() {
    return this->size.size();
  }
  
};

#endif /* MATLABIMAGEHEADER_H */
