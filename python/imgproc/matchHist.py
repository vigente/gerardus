#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
imgproc/matchHist.py
@package: gerardus
@subpackage: imgproc
@author: Ramon Casero <rcasero@gmail.com>
@copyright: © 2017 University of Oxford
@license: GPL v3

Version: 1.0.0

University of Oxford means the Chancellor, Masters and Scholars of
the University of Oxford, having an administrative office at
Wellington Square, Oxford OX1 2JD, UK. 

This file is part of Gerardus.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. The offer of this
program under the terms of the License is subject to the License
being interpreted in accordance with English Law and subject to any
action against the University of Oxford being under the jurisdiction
of the English Courts.

You should have received a copy of the GNU General Public License
along with this program.  If not, see
<http://www.gnu.org/licenses/>.
"""

import numpy as np

def matchHist(imref, im, maskref=np.ones(0, dtype=bool), mask=np.ones(0, dtype=bool), nbr_bins=256):
    """Modify image intensities to match the histogram of a reference image.
    
    imout = matchHist(imref, im)
    
    Args:
        imref, im: Grayscale or colour 2D images [row, col, channel]. im is the
                   image we want to modify so that it matches the histogram 
                   (channel by channel) of imref.

    Returns:                   
        imout: The modified version of im.
                   
    imout = matchHist(imref, im, maskref=[], mask=[], nbr_bins=256)
    
    Optional:
        maskref, mask: Bool masks for imref, im, respectively. The masks must 
                       have the same number of [rows,cols] as their respective
                       images, but only 1 channel. Pixels set to False in the
                       mask are completely ignored (default: no mask)
                       
        nbr_bins: Number of bins used to compute histograms (default: 256)
    """
    
    # duplicate inputs, to avoid modifying the objects they point to outside
    # this function
    imout = im.copy()
    imrefaux = imref.copy()
    
    # mask must be boolean
    if maskref.dtype != "bool":
        raise TypeError("maskref must be of type bool")
    if mask.dtype != "bool":
        raise TypeError("mask must be of type bool")
        
    # mask must have only one channel. The same is applied to each image 
    # channel
    if maskref.ndim > 2:
        raise ValueError("maskref can have at most one channel")
    if mask.ndim > 2:
        raise ValueError("mask can have at most one channel")
        
    # grayscale images will be treated as multi-channel images with 1 channel
    if len(imrefaux.shape) < 3:
        imrefaux = imrefaux[:,:,np.newaxis]
    if len(imout.shape) < 3:
        imout = imout[:,:,np.newaxis]
    
    # compute histogram of each channel of the reference and converted images
    for i in range(0, imrefaux.shape[2]):
        
        # extract channel from the image
        chanref = imrefaux[:, :, i]
        chan = imout[:, :, i]
        
        # extract masked pixels, if masks are provided. Otherwise, use all 
        # pixels flattening the channel
        if len(maskref) > 0:
            chanref_flat = chanref[maskref]
        else:
            chanref_flat = chanref.flatten()
            
        if len(mask) > 0:
            chan_flat = chan[mask]
        else:
            chan_flat = chan.flatten()

        # compute histograms
        imhistref, binsref = np.histogram(chanref_flat, nbr_bins, normed=True)
        imhist, bins = np.histogram(chan_flat, nbr_bins, normed=True)
            
        # cumulative distribution function
        cdfhistref = imhistref.cumsum()
        cdfhist = imhist.cumsum()
        
        # bin centers
        cbinsref = (binsref[:-1] + binsref[1:]) / 2.0
        cbins = (bins[:-1] + bins[1:]) / 2.0
        
        # map intensity values in current channel so that they match the 
        # reference histogram
        chan_flat_to_cdf = np.interp(chan_flat, cbins, cdfhist)
        chan_flat_mapped = np.interp(chan_flat_to_cdf, cdfhistref, cbinsref)
        
        # tranfer corrected pixels to image
        if len(mask) > 0:
            chan[mask] = chan_flat_mapped
        else:
            chan = np.reshape(chan_flat_mapped, chan.shape)
            
        imout[:, :, i] = chan
        
    # return corrected image
    return imout
    
      
###############################################################################
## TEST
###############################################################################

# if module is executed as scrip instead of imported, run test    
if __name__ == "__main__":

    import os 
    from scipy import misc

    import matplotlib.pyplot as plt

    # directory where this module lives    
    module_path = os.path.dirname(os.path.realpath(__file__))
    
    # path and name of test files
    imref_name = os.path.join(module_path, ".." + os.sep + "testdata", "right.png")
    im_name = os.path.join(module_path, ".." + os.sep + "testdata", "left.png")
    maskref_name = os.path.join(module_path, ".." + os.sep + "testdata", "right_mask.png")
    mask_name = os.path.join(module_path, ".." + os.sep + "testdata", "left_mask.png")

    # read test images and their masks
    imref = misc.imread(imref_name)
    im = misc.imread(im_name)
    maskref = misc.imread(maskref_name)
    mask = misc.imread(mask_name)
    maskref = maskref[:, :, 1]==255
    mask = mask[:, :, 1]==255

    # plot images
    plt.close('all')
    fig, ax = plt.subplots(2, 2)
    ax[0, 0].imshow(imref)
    ax[0, 0].set_title("Ref image")
    ax[1, 0].imshow(im)
    ax[1, 0].set_title("Image to be corrected")
    ax[0, 1].imshow(maskref)
    ax[0, 1].set_title("Ref mask")
    ax[1, 1].imshow(mask)
    ax[1, 1].set_title("Mask of image to be corrected")
 
    # match histogram of im to imref, only taking into account the pixels==1 in
    # the masks
    im_matched = matchHist(imref, im, maskref=maskref, mask=mask)
    
    # plot images
    fig, ax = plt.subplots(2, 2)
    ax[0, 0].imshow(imref)
    ax[0, 0].set_title("Ref image")
    ax[1, 0].imshow(im_matched)
    ax[1, 0].set_title("Corrected image")
    ax[0, 1].imshow(maskref)
    ax[0, 1].set_title("Ref mask")
    ax[1, 1].imshow(mask)
    ax[1, 1].set_title("Mask of corrected image")

    plt.show()
    
