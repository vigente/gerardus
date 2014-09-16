function blockface_correct_frame_shifts(indir, files, t, idxprop, idxnoprop, outdir)
% blockface_correct_frame_shifts  Correct shift between blockface image
% frames.
%
% When the wax block is sliced with a microtome, and photographs taken of
% the blockface, artifacts invariable happen. Artifacts are typically
% caused by the microtome tray not being pushed all the way back
% (translation artifacts), or by the camera, mirrors or any other elements
% in the light path being moved a little (similarity transform artifacts).
%
% This function assumes that the corrective transforms between consecutive
% frames have been already computed with elastix. The user provides a list
% of frames that need to be corrected, and blockface_correct_frame_shifts
% computes the accumulated transform for each frame, and applies it with
% transformix. Corrected images are saved with the same name as input
% images, but to a different directory.
%
% blockface_correct_frame_shifts(INDIR, FILES, T, IDXPROP, IDXNOPROP, OUTDIR)
%
%   INDIR is a string with the directory where the input files are kept.
%
%   FILES is the result of a dir() command with a list of the files to
%   correct.
%
%   T is a struct array where T(I) is the correction that registers image I
%   to image I-1. T can be obtained with blockface_intraframe_reg.
%
%   IDXPROP is a vector with a list of indices of frames that need to be
%   corrected, and whose correction has to propagate to subsequent frames.
%   For example, if the camera is moved a bit at frame I=30, that
%   displacement is going to show up in all frames from I=30 onwards.
%
%   IDXNOPROP is a vector with a list of indices of frames that need to be
%   corrected, but whose correction does not propagate. For example, if the
%   microtome tray is not pushed all the way in in frame I=30, that doesn't
%   mean that that problem will appear in frame I = 31. This is useful if
%   e.g. the camera has been very stable, but now and then the human
%   operator doesn't push the tray all the way back, as it avoids
%   unncessary small transforms and resampling of images.
%
%   OUTDIR is a string with the directory for output images. INDIR and
%   OUTDIR must be different to avoid overwriting the input files.
%
% See also: blockface_intraframe_reg.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.1.1
% $Rev$
% $Date$
% 
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(6, 6);
nargoutchk(0, 0);

% safety to prevent us from overwriting the source files
if (strcmp(indir, outdir))
    error('Input and output directory are the same, and files would be overwritten')
end

% if destination directory doesn't exist, create it
if (isempty(dir(outdir)))
    ok = mkdir(outdir);
    if (~ok)
        error(['Cannot create output directory: ' outdir])
    end
end

% first frame will be used as a reference for the rest
% copytogreyfile([pathstr filesep files(1).name], ...
%     [outdir filesep files(1).name]);

% initialize array for the accumulated transform for each frame
tAcc = t(1);

% clear up transform for first frame
tAcc(1).TransformParameters = [1 0 0 0];

% accumulate transforms that propagate
for I = 2:length(t)
    
    % if this frame needs to be corrected
    if (~isempty(idxprop) && I == idxprop(1))
        
        % update the accumulated transform with the transform for this
        % frame
        tAcc(I) = elastix_compose_afftransf(tAcc(I-1), t(I));
        
        % remove the first element of the index queue
        idxprop = idxprop(2:end);
        
    else 
        
        tAcc(I) = tAcc(I-1);
        
    end
    
end

% add transforms that do not propagate
for I = idxnoprop
   
    tAcc(I) = elastix_compose_afftransf(tAcc(I), t(I));
    
end

% correct frames
for I = 1:length(files)
    
    % is the transform the identity?
    if (all(tAcc(I).TransformParameters == [1 0 0 0]))
        
        % simply convert to grayscale and copy 
        copytogreyfile([indir filesep files(I).name], ...
            [outdir filesep files(I).name]);
        
    else
        
        % affine transformation of image to correct frame shifts
        fileout = transformix(tAcc(I), [indir filesep files(I).name]);
        
        % copy result to file
        copytogreyfile(fileout, ...
            [outdir filesep files(I).name]);
        
    end
    
end

end

function copytogreyfile(filein, fileout)

% check that the image is grayscale. Otherwise, we convert it,
% because elastix uses the first channel and ignores the rest. As
% all transformed images will be converted to grayscale, we also
% convert images that are not shifted, so that they are all
% grayscale
info = imfinfo(filein);

switch (info.ColorType)
    case 'truecolor'
        
        % load image
        im = imread(filein);
        
        % convert to grayscale
        im = rgb2gray(im);
        
        % write image to destination
        imwrite(im, fileout);
        
    case 'grayscale'
        
        % if the image is already grayscale, we can directly copy
        % it
        ok = copyfile(filein, fileout);
        if (~ok)
            error(['Cannot copy image ' filein ' to file ' fileout])
        end
       
    otherwise
        
        error('Image in input file is neither RGB nor grayscale')
end

end
