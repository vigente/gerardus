function out = scimat_resample(inp,out,kind)
% SCIMAT_RESAMPLE resamples a scimat file given an input and output scimat.
%
% Intended for mapping 2D slices in 3D volume or extract 2D
% slices from 3D volume. Also works nicely between 1D and 2D to extract
% profiles, for example.
%
% OUT = SCIMAT_RESAMPLE(INP, OUT, KIND)
%
%   INP is the input scimat file
%
%   OUT is the target scimat file
%
%   KIND is a string, deciding on the kind of interpolation/averaging

% Author: Nicolas Basty <nicolas.basty@eng.ox.ac.uk>
% Copyright © 2016 University of Oxford
% Version: 0.1.1
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
% _________________________________________________________________

[x,y,z] = inp.axis.spacing; %get in and output resolutions
[x_o,y_o,z_o] = out.axis.spacing;

if x_o>x %taking lower resolution
    x = x_o;
    y = y_o;
    z = z_o;
end

sigmaa = [x,y,inf]/2;
ws =  sigmaa*3; %window size
ws(3) = z/2;
rrrr=r;
if strcmp(kind,'gaussian') % get the weighting function
    W = @(xl,yl,zl) exp(-(xl.^2)/(2*sigmaa(1)^2)...
        -(yl.^2)/(2*sigmaa(2)^2)...
        -(zl.^2)/(2*sigmaa(3)^2));
elseif strcmp(kind,'sinc')
    W = @(xl,yl,zl) (sin(sigmaa(1)*xl)./(sigmaa(1)*xl))...
        .* (sin(sigmaa(2)*yl)./(sigmaa(2)*yl))...
        .* (sin(sigmaa(3)*zl)./(sigmaa(3)*zl));
elseif strcmp(kind,'nearest')
    %does something later, l84
else
    fprintf('Weighting method expected to be gaussian, sinc or nearest\n')
    return
end

To = scimat2transform(inp); %gets transform from input image mapping IDX to RW
Ti = scimat2transform(out);  %gets transform from target image mapping IDX to RW
T = pinv(To) * Ti; %combining gives IDXin to IDXout

% create grid of output�s indices
if numel(out.data)>2
    [X, Y, Z] = ndgrid(1:size(out.data,1), 1:size(out.data,2), 1:size(out.data,3));
elseif  numel(out.data)<3
    [X, Y, Z] = ndgrid(1:size(out.data,1), 1:size(out.data,2), 1);
end
% for each of the output index, find its corresponding input index
XYZin = T * [X(:), Y(:), Z(:), ones(size(X(:)))]';

% round them to do nearest neightbour assignment
Xin = round(XYZin(1,:));
Yin = round(XYZin(2,:));
Zin = round(XYZin(3,:));

% find which ones are valid (input space)
valid_indices = (Xin > 0) & (Xin <= size(inp.data,1)) & (Yin > 0) & (Yin <= size(inp.data,2)) & (Zin > 0) & (Zin <= size(inp.data,3));

%%
inpt = double(inp.data); %get data out of the structures
outp = double(out.data);

allind = 1:numel(valid_indices);
allind = allind(valid_indices);

[xi,yi,zi] = scimat_ndgrid(inp); % grid of input rw
[xo,yo,zo] = scimat_ndgrid(out); % grid of output rw

vox_out = [xo(:),yo(:),zo(:)];

[Xi, Yi, Zi] = ndgrid(1:size(inp.data,1), 1:size(inp.data,2), 1:size(inp.data,3));
xyzin = round(pinv(T)*[Xi(:) Yi(:) Zi(:) ones(size(Xi(:)))]');
xin = xyzin(1,:);
yin = xyzin(2,:);
zin = xyzin(3,:);
%using indices to get valid voxels, but the other part works better in
%rw... window size gets messed up when converting to index
valid_indices_out = (xin > 0) & (xin <= size(out.data,1)) & (yin > 0) & (yin <= size(out.data,2)) & (zin > 0) & (zin <= size(out.data,3));

inpt_val = inpt(valid_indices_out);
%keep voxels relevant to the averaging
xi = xi(valid_indices_out);
yi = yi(valid_indices_out);
zi = zi(valid_indices_out);

if strcmp(kind,'nearest')
    source_indices = sub2ind(size(inp.data), Xin(valid_indices), Yin(valid_indices), Zin(valid_indices));
    valid_intensities_input = inp.data(source_indices);
    out.data(valid_indices) = valid_intensities_input;
    return
end
outpdummy = zeros(numel(allind),1);
% xyz = [xi(:) yi(:) zi(:)];
% par parfor is actually slower for the data I am using so I commented it
% out
for r = 1:numel(allind)
    p =  vox_out(allind(r),:)';
    % get a neighbourhood of vox_out around this point (rw)
    mask = (abs(xi-p(1)) <= ws(1)) & (abs(yi - p(2)) <= ws(2)) & (abs(zi - p(3)) <= ws(3));
    %     mask = bsxfun(@le, (abs(bsxfun(@minus, xyz, p'))), ws);
    %     mask = sum(mask,2)==3;
    % get a shortlist of voxels in the window & distances for weighting
    xlist = xi(mask) - p(1);
    ylist = yi(mask) - p(2);
    zlist = zi(mask) - p(3);
    Ilist = inpt_val(mask);
    
    % get a weighted average
    outpdummy(r) = sum(W(xlist,ylist,zlist).*Ilist) / (sum(W(xlist,ylist,zlist)+eps));
end
outp(allind) = outpdummy;
outp(isnan(outp)) = 0;
out.data = outp;
end