This folder contains the source code from:

“Prospective acceleration of diffusion tensor imaging with compressed sensing using adaptive dictionaries”

The paper is available here:

http://onlinelibrary.wiley.com/doi/10.1002/mrm.25876/abstract

The 2D phantom experiment was used to generate figure 3 - it might be slightly simpler to understand than the 3D source code, but the algorithms are the same.

In the 3D code, CS_dictionary_3D is the one to look at first - it won’t run because all of the file paths have been redacted.

You might need to add the gerardus repository to your matlab path to get functions like “compute_helix_angle_on_image.m”, but most others should be already there.

You will need to get genPDF.m from Michael Lustig’s compressed sensing toolbox.

Feel free to contact me with any questions or comments at darryl.mcclymont@gmail.com

18 Jan 2016