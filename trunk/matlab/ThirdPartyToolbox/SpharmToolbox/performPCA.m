%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Harmonic Modeling and Analysis Toolkit (SPHARM-MAT) is a 3D 
% shape modeling and analysis toolkit. 
% It is a software package developed at Shenlab in Center for Neuroimaging, 
% Indiana University (SpharmMat@gmail.com, http://www.iupui.edu/~shenlab/)
% It is available to the scientific community as copyright freeware 
% under the terms of the GNU General Public Licence.
% 
% Copyright 2009, 2010, ShenLab, Center for Neuroimaging, Indiana University
% 
% This file is part of SPHARM-MAT.
% 
% SPHARM-MAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SPHARM-MAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SPHARM-MAT. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function performPCA(confs, objs)

numSbj = length(objs);
groupID = confs.GroupID;

% set parameters
fvecs = [];	% fvecs is accumulator of coefficients from each species.
keep = numSbj-1;	% number of eigenvectors to keep

% do processing of each file
for i = 1:numSbj
	% read in file.
	file = objs{i};
	load(file);
	[path, name, ext] = fileparts(file);
	clear('faces', 'vertices', 'sph_verts');

	[nrows ncols] = size(fvec);
	temp = fvec.';
	temp = reshape(temp, (nrows*3), 1);
	fvecs = [fvecs temp];
end
clear('fvec', 'temp', 'name');

fvecs = fvecs.';
[eigenvecs, scores, eigenvals] = princomp(fvecs);
% trace_eig = trace(eigenvals);
perc_variance_explained = eigenvals/sum(eigenvals);

% drop cells that are to not be kept
eigenvecs = eigenvecs(:,1:keep);
%scores = scores(:,1:keep);
eigenvals = eigenvals(1:keep,1);
perc_variance_explained = perc_variance_explained(1:keep,:);
cum_percent_explained = cumsum(perc_variance_explained);

new_name = [confs.OutDirectory '/' confs.OutputName];
if exist(new_name,'file')
    prompt = {'Enter new filename:'};
    dlg_title = 'New File Name';
    num_lines = 1;
    def = {new_name};
    answer = inputdlg(prompt,dlg_title,num_lines,def);    
    new_name = answer{1};
end

save(new_name, 'fvecs', 'eigenvecs', 'eigenvals', 'perc_variance_explained', ... 
    'cum_percent_explained', 'scores','groupID');

return;