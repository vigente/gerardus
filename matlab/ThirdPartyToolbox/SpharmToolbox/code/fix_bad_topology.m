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

% Goal: Fix a binary image by removing bad connectivities
% inObject: filename string of a binary object in .MAT format
% info = varargin{1}: connectivity 
% info can be one of four values between 1 and 4 - '1=(6+,18), 2=(18,6+),
% 3=(6,26), 4=(26,6)'
% epsilon = varargin{2}: hole size, holes bigger than epsilon won't be
% filled
%
% 04/15/2002 - created by Li Shen
% Modified by Sungeun Kim (10-09-08)

function [bim, origin, vxsize, new_name] = fix_bad_topology(inObject, confs)

epsilon = confs.Epsilon;
switch deblank(char(confs.Connectivity))
    case '(6+,18)'
        conn = 1;
    case '(18,6+)'
        conn = 2;        
    case '(6,26)'
        conn = 3;        
    case '(26,6)'
        conn = 4;        
end

% Read input object dataset
if ischar(inObject)
    [path,name,ext,ver] = fileparts(inObject);

    if ~strcmp(ext, '.mat')
        disp('Input object data should be Matlab data file (*.mat)');
        return;
    end

    load(inObject);
    roi = bim; 
else
    disp('Input Object should be the filename of a binary object (_bim.mat).');
    return;
end

tic;

infostr{1} = 'outlier or hole';
infostr{2} = 'vertex conn.';
infostr{3} = 'edge conn.';
infostr{4} = 'ring or hole';
infostr{5} = 'single ring';

% display information
disp('Make a simply connected volume by removing bad voxels');
disp(sprintf('(1) %s, (2) %s, (3) %s, (4) %s, (5) s', ...
              infostr{1}, infostr{2}, infostr{3}, infostr{4}, infostr{5}));

% Make binary image, in case reslice generates 2's
ind = find(roi>0); roi(ind) = 1;
vol = length(ind);
ind = find(roi<1); roi(ind) = 0;

% routine fix
roi = routine_fix(roi);
[vertices, faces] = gen_surf_data(roi,origin,vxsize);

% fill 3d holes
for k=1:3
    if (size(vertices,1)-size(faces,1)==2)
        disp('===> No 3D holes, skip');
        break;
    end
    disp('===> fixing 3D holes ...');
    roi = fill3dholes(roi,conn,epsilon,name);
    disp('===> routine fix again ...');
    roi = routine_fix(roi);
    [vertices, faces] = gen_surf_data(roi,origin,vxsize);
end

fvn = length(find(roi~=bim));
bim = roi;

disp(sprintf('fix_vox_num (%d) / vol (%d) = %f%%, save ...', ...
                      fvn, vol, fvn*100/vol));
                  
% Save surface object to a new file
postfix = name(end-2:end);
if strcmp(postfix, 'bim') | strcmp(postfix,'fix')
    new_name = [confs.OutDirectory '/' name(1:end-3) 'fix'];
else
    new_name = [confs.OutDirectory '/' name '_fix'];
end   

if exist(new_name,'file')
    prompt = {'Enter new filename:'};
    dlg_title = 'New File Name';
    num_lines = 1;
    def = {new_name};
    answer = inputdlg(prompt,dlg_title,num_lines,def);    
    new_name = answer{1};
end

save(new_name, 'bim', 'origin', 'vxsize');

seconds = toc;
hours = seconds/3600;
disp(['seconds=', num2str(seconds), ' hours=', num2str(hours)]);

return;


%
% routine fix
%

function roi = routine_fix(roi)

DIM = size(roi);

fixset = [];

for j = 1:10
    [roi, fs1] = outl_hole(roi, DIM,1);
    [roi, fs2] = vertex_conn(roi, DIM,2);
    [roi, fs3] = edge_conn(roi, DIM,3);
    if (isempty(fs1) & isempty(fs2) & isempty(fs3))
        break;
    else
        if (~isempty(fs1))
            fixset(end+1:end+size(fs1,1),:) = fs1;
        end
        if (~isempty(fs2))
            fixset(end+1:end+size(fs2,1),:) = fs2;
        end
        if (~isempty(fs3))
            fixset(end+1:end+size(fs3,1),:) = fs3;
        end
    end
end

[roi, fs4] = ring_hole(roi, DIM, 4);
if (~isempty(fs4))
    fixset(end+1:end+size(fs4,1),:) = fs4;        
    for j = 1:10
        [roi, fs1] = outl_hole(roi, DIM,1);
        [roi, fs2] = vertex_conn(roi, DIM,2);
        [roi, fs3] = edge_conn(roi, DIM,3);
        [roi, fs4] = ring_hole(roi, DIM,4);
        if (isempty(fs1) & isempty(fs2) & isempty(fs3) & isempty(fs4))
            break;
        else
            if (~isempty(fs1))
                fixset(end+1:end+size(fs1,1),:) = fs1;
            end
            if (~isempty(fs2))
                fixset(end+1:end+size(fs2,1),:) = fs2;
            end
            if (~isempty(fs3))
                fixset(end+1:end+size(fs3,1),:) = fs3;
            end
            if (~isempty(fs4))
                fixset(end+1:end+size(fs4,1),:) = fs4;
            end
        end
    end
end

fix_vox_num = size(fixset,1);

return;


%
% remove disconnected small components and holes
%

function [roi, fixset] = outl_hole(roi, d, infoid)

roi = reshape(roi,d);

fixset = [];

% check seperated small components (k=1)
% and then check holes
val = [0 1];
for k = 1:2
	[L, num] = bwlabeln(roi,6);
	
	if num>1
		% find the background
		for i = 1:num
            cnt(i) = length(find(L == i));
		end
		[bkgd_size, bkgd_ind] = max(cnt);
        sigind = find(L~=0 & L~=bkgd_ind);
        roi(sigind) = 0;
        [xs,ys,zs] = ind2sub(d, sigind);        
        fixset(end+1:end+length(xs),:) = [[xs,ys,zs], xs*0+val(k), xs*0+infoid];
	end
	
	roi = 1-roi;
end

return;


%
% fix bad vertex connectivities
%

function [roi, fixset] = vertex_conn(roi, d, infoid)

roi = reshape(roi,d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% work on 0 vertex connectivies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a work area so that all border voxels belong to the background
w = ones(d+2);
w(2:d(1)+1,2:d(2)+1,2:d(3)+1) = roi;
wd = d+2;
ind = find(w == 0);
w = zeros(d+2);
w(2:d(1)+1,2:d(2)+1,2:d(3)+1) = roi;

[xs,ys,zs] = ind2sub(wd,ind);

h = sum([w(sub2ind(wd,xs-1,ys,zs)), ...
         w(sub2ind(wd,xs+1,ys,zs)), ...
         w(sub2ind(wd,xs,ys-1,zs)), ...
         w(sub2ind(wd,xs,ys+1,zs)), ...
         w(sub2ind(wd,xs,ys,zs-1)), ...
         w(sub2ind(wd,xs,ys,zs+1))]')';

fixset = [];

% fix ring 0 voxels (single ring)
ind1 = ind(find(h==4 & ...
                w(sub2ind(wd,xs-1,ys,zs))==0 & ...
                w(sub2ind(wd,xs+1,ys,zs))==0 ...
           ));
ind2 = ind(find(h==4 & ...
                w(sub2ind(wd,xs,ys-1,zs))==0 & ...
                w(sub2ind(wd,xs,ys+1,zs))==0 ...
           ));
ind3 = ind(find(h==4 & ...
                w(sub2ind(wd,xs,ys,zs-1))==0 & ...
                w(sub2ind(wd,xs,ys,zs+1))==0 ...
           ));
ind1 = [ind1;ind2;ind3];
n4 = length(ind1);
if (n4~=0)	
	w(ind1) = 1;
    [a, b, c] = ind2sub(wd,ind1);
	fixset(end+1:end+n4,:) = [[a, b, c]-1, ind1*0+1, ind1*0+5];
end

% right up forward, right up backward, right down forward, right down backward
ruf = find(w(sub2ind(wd,xs-1,ys+1,zs+1)) == 0 & ...
           sum([w(sub2ind(wd,xs-1,ys,zs)), ...
                w(sub2ind(wd,xs,ys+1,zs)), ...
                w(sub2ind(wd,xs-1,ys+1,zs)), ...
                w(sub2ind(wd,xs,ys,zs+1)), ...
                w(sub2ind(wd,xs-1,ys,zs+1)), ...
                w(sub2ind(wd,xs,ys+1,zs+1))]')' == 6);

rub = find(w(sub2ind(wd,xs-1,ys+1,zs-1)) == 0 & ...
           sum([w(sub2ind(wd,xs-1,ys,zs)), ...
                w(sub2ind(wd,xs,ys+1,zs)), ...
                w(sub2ind(wd,xs-1,ys+1,zs)), ...
                w(sub2ind(wd,xs,ys,zs-1)), ...
                w(sub2ind(wd,xs-1,ys,zs-1)), ...
                w(sub2ind(wd,xs,ys+1,zs-1))]')' == 6);

rdf = find(w(sub2ind(wd,xs+1,ys+1,zs+1)) == 0 & ...
           sum([w(sub2ind(wd,xs+1,ys,zs)), ...
                w(sub2ind(wd,xs,ys+1,zs)), ...
                w(sub2ind(wd,xs+1,ys+1,zs)), ...
                w(sub2ind(wd,xs,ys,zs+1)), ...
                w(sub2ind(wd,xs+1,ys,zs+1)), ...
                w(sub2ind(wd,xs,ys+1,zs+1))]')' == 6);

rdb = find(w(sub2ind(wd,xs+1,ys+1,zs-1)) == 0 & ...
           sum([w(sub2ind(wd,xs+1,ys,zs)), ...
                w(sub2ind(wd,xs,ys+1,zs)), ...
                w(sub2ind(wd,xs+1,ys+1,zs)), ...
                w(sub2ind(wd,xs,ys,zs-1)), ...
                w(sub2ind(wd,xs+1,ys,zs-1)), ...
                w(sub2ind(wd,xs,ys+1,zs-1))]')' == 6);

cube = [];

if (~isempty(ruf))
    cube(end+1:end+length(ruf),:) = ...
        [xs(ruf),ys(ruf),zs(ruf), ...
         xs(ruf)-1,ys(ruf)+1,zs(ruf)+1];
end
    
if (~isempty(rub))
    cube(end+1:end+length(rub),:) = ...
        [xs(rub),ys(rub),zs(rub), ...
         xs(rub)-1,ys(rub)+1,zs(rub)-1];
end

if (~isempty(rdf))
    cube(end+1:end+length(rdf),:) = ...
        [xs(rdf),ys(rdf),zs(rdf), ...
         xs(rdf)+1,ys(rdf)+1,zs(rdf)+1];
end

if (~isempty(rdb))
    cube(end+1:end+length(rdb),:) = ...
        [xs(rdb),ys(rdb),zs(rdb), ...
         xs(rdb)+1,ys(rdb)+1,zs(rdb)-1];
end

% adjust according to cube (3: fix bad 0 vertex connectivities)
if (~isempty(cube))
    len = size(cube,1);
    for i = 1:len
        k = cube(i,:);
        if (w(k(1),k(2),k(3))==0 & w(k(4),k(5),k(6))==0)
            if (nb_sum(w,k(1),k(2),k(3)) >= nb_sum(w,k(4),k(5),k(6)))
                w(k(1),k(2),k(3)) = 1;                
                fixset(end+1,:) = [[k(1), k(2), k(3)]-1, 1, infoid];
            else
                w(k(4),k(5),k(6)) = 1;                
                fixset(end+1,:) = [[k(4), k(5), k(6)]-1, 1, infoid];
            end
        end
    end
end

cnts = [length(ruf), length(rub), length(rdf), length(rdb)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% work on 1 vertex connectivies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = find(w);
[xs,ys,zs] = ind2sub(wd,ind);

% right up forward, right up backward, right down forward, right down backward
ruf = find(w(sub2ind(wd,xs-1,ys+1,zs+1)) == 1 & ...
           sum([w(sub2ind(wd,xs-1,ys,zs)), ...
                w(sub2ind(wd,xs,ys+1,zs)), ...
                w(sub2ind(wd,xs-1,ys+1,zs)), ...
                w(sub2ind(wd,xs,ys,zs+1)), ...
                w(sub2ind(wd,xs-1,ys,zs+1)), ...
                w(sub2ind(wd,xs,ys+1,zs+1))]')' == 0);

rub = find(w(sub2ind(wd,xs-1,ys+1,zs-1)) == 1 & ...
           sum([w(sub2ind(wd,xs-1,ys,zs)), ...
                w(sub2ind(wd,xs,ys+1,zs)), ...
                w(sub2ind(wd,xs-1,ys+1,zs)), ...
                w(sub2ind(wd,xs,ys,zs-1)), ...
                w(sub2ind(wd,xs-1,ys,zs-1)), ...
                w(sub2ind(wd,xs,ys+1,zs-1))]')' == 0);
        
rdf = find(w(sub2ind(wd,xs+1,ys+1,zs+1)) == 1 & ...
           sum([w(sub2ind(wd,xs+1,ys,zs)), ...
                w(sub2ind(wd,xs,ys+1,zs)), ...
                w(sub2ind(wd,xs+1,ys+1,zs)), ...
                w(sub2ind(wd,xs,ys,zs+1)), ...
                w(sub2ind(wd,xs+1,ys,zs+1)), ...
                w(sub2ind(wd,xs,ys+1,zs+1))]')' == 0);

rdb = find(w(sub2ind(wd,xs+1,ys+1,zs-1)) == 1 & ...
           sum([w(sub2ind(wd,xs+1,ys,zs)), ...
                w(sub2ind(wd,xs,ys+1,zs)), ...
                w(sub2ind(wd,xs+1,ys+1,zs)), ...
                w(sub2ind(wd,xs,ys,zs-1)), ...
                w(sub2ind(wd,xs+1,ys,zs-1)), ...
                w(sub2ind(wd,xs,ys+1,zs-1))]')' == 0);

cube = [];

if (~isempty(ruf))
    cube(end+1:end+length(ruf),:) = ...
        [xs(ruf),ys(ruf),zs(ruf), ...
         xs(ruf)-1,ys(ruf)+1,zs(ruf)+1];
end
    
if (~isempty(rub))
    cube(end+1:end+length(rub),:) = ...
        [xs(rub),ys(rub),zs(rub), ...
         xs(rub)-1,ys(rub)+1,zs(rub)-1];
end

if (~isempty(rdf))
    cube(end+1:end+length(rdf),:) = ...
        [xs(rdf),ys(rdf),zs(rdf), ...
         xs(rdf)+1,ys(rdf)+1,zs(rdf)+1];
end

if (~isempty(rdb))
    cube(end+1:end+length(rdb),:) = ...
        [xs(rdb),ys(rdb),zs(rdb), ...
         xs(rdb)+1,ys(rdb)+1,zs(rdb)-1];
end

% adjust according to cube (4: remove bad 1 voxel connectivities)
if (~isempty(cube))
    len = size(cube,1);
    for i = 1:len
        k = cube(i,:);
        if (w(k(1),k(2),k(3))==1 & w(k(4),k(5),k(6))==1)
            if (nb_sum(w,k(1),k(2),k(3)) <= nb_sum(w,k(4),k(5),k(6)))
                w(k(1),k(2),k(3)) = 0;                
                fixset(end+1,:) = [[k(1), k(2), k(3)]-1, 0, infoid];
            else
                w(k(4),k(5),k(6)) = 0;                
                fixset(end+1,:) = [[k(4), k(5), k(6)]-1, 0, infoid];
            end
        end
    end
end

cnts = [length(ruf), length(rub), length(rdf), length(rdb)];

% make roi consistent with w
roi = w(2:d(1)+1,2:d(2)+1,2:d(3)+1);

return;


%
% fix bad edge connectivities
%

function [roi, fixset] = edge_conn(roi, d, infoid)

roi = reshape(roi,d);

% make a work area so that all border voxels belong to the background
w = zeros(d+2);
w(2:d(1)+1,2:d(2)+1,2:d(3)+1) = roi;
wd = d+2;

ind = find(w);
[xs,ys,zs] = ind2sub(wd,ind);

squa = [];

% right up and right down
ru = find(w(sub2ind(wd,xs-1,ys,zs))   == 0 & ...
          w(sub2ind(wd,xs,ys+1,zs))   == 0 & ...
          w(sub2ind(wd,xs-1,ys+1,zs)) == 1);
      
if (~isempty(ru))
    squa(end+1:end+length(ru),:) = ... 
                [xs(ru),ys(ru),zs(ru),...
                 xs(ru)-1,ys(ru)+1,zs(ru),...
                 xs(ru)-1,ys(ru),zs(ru),...
                 xs(ru),ys(ru)+1,zs(ru)];
end

rd = find(w(sub2ind(wd,xs+1,ys,zs))   == 0 & ...
          w(sub2ind(wd,xs,ys+1,zs))   == 0 & ...
          w(sub2ind(wd,xs+1,ys+1,zs)) == 1);

if (~isempty(rd))
    squa(end+1:end+length(rd),:) = ... 
                [xs(rd),ys(rd),zs(rd),...
                 xs(rd)+1,ys(rd)+1,zs(rd),...
                 xs(rd)+1,ys(rd),zs(rd),...
                 xs(rd),ys(rd)+1,zs(rd)];
end

% forward up and forward down
fu = find(w(sub2ind(wd,xs-1,ys,zs))   == 0 & ...
          w(sub2ind(wd,xs,ys,zs+1))   == 0 & ...
          w(sub2ind(wd,xs-1,ys,zs+1)) == 1);

if (~isempty(fu))
    squa(end+1:end+length(fu),:) = ... 
                [xs(fu),ys(fu),zs(fu),...
                 xs(fu)-1,ys(fu),zs(fu)+1,...
                 xs(fu)-1,ys(fu),zs(fu),...
                 xs(fu),ys(fu),zs(fu)+1];
end
     
fd = find(w(sub2ind(wd,xs+1,ys,zs))   == 0 & ...
          w(sub2ind(wd,xs,ys,zs+1))   == 0 & ...
          w(sub2ind(wd,xs+1,ys,zs+1)) == 1);
 
if (~isempty(fd))
    squa(end+1:end+length(fd),:) = ... 
                [xs(fd),ys(fd),zs(fd),...
                 xs(fd)+1,ys(fd),zs(fd)+1,...
                 xs(fd)+1,ys(fd),zs(fd),...
                 xs(fd),ys(fd),zs(fd)+1];
end
     
% left forward and right forward
lf = find(w(sub2ind(wd,xs,ys-1,zs))   == 0 & ...
          w(sub2ind(wd,xs,ys,zs+1))   == 0 & ...
          w(sub2ind(wd,xs,ys-1,zs+1)) == 1);

if (~isempty(lf))
    squa(end+1:end+length(lf),:) = ... 
                [xs(lf),ys(lf),zs(lf),...
                 xs(lf),ys(lf)-1,zs(lf)+1,...
                 xs(lf),ys(lf)-1,zs(lf),...
                 xs(lf),ys(lf),zs(lf)+1];
end
     
rf = find(w(sub2ind(wd,xs,ys+1,zs))   == 0 & ...
          w(sub2ind(wd,xs,ys,zs+1))   == 0 & ...
          w(sub2ind(wd,xs,ys+1,zs+1)) == 1);

if (~isempty(rf))
    squa(end+1:end+length(rf),:) = ... 
                [xs(rf),ys(rf),zs(rf),...
                 xs(rf),ys(rf)+1,zs(rf)+1,...
                 xs(rf),ys(rf)+1,zs(rf),...
                 xs(rf),ys(rf),zs(rf)+1];
end

% adjust according to squa (5: remove bad edge connectivities)
fixset = [];
if (~isempty(squa))
    len = size(squa,1);
    for i = 1:len
        k = squa(i,:);
        if (sum([w(k(1),k(2),k(3)) ...
                 w(k(4),k(5),k(6)) ...
                 w(k(7),k(8),k(9)) ...
                 w(k(10),k(11),k(12)) ...
                 ]==[1 1 0 0])==4)
            [start, val] = decide_squa(w,k);
            w(k(start),k(start+1),k(start+2)) = val;
            fixset(end+1,:) = [[k(start), k(start+1), k(start+2)]-1, val, infoid];
        end
    end
end

cnts = [length(ru),length(rd),length(fu),length(fd),length(lf),length(rf)];

% make roi consistent with w
roi = w(2:d(1)+1,2:d(2)+1,2:d(3)+1);

return;


%
% check ring on one slice, and hole between two slices
%

function [roi, fixset] = ring_hole(roi, d, infoid)

roi = reshape(roi,d);

% make a work area so that all border voxels belong to the background
w = zeros(d+2);
w(2:d(1)+1,2:d(2)+1,2:d(3)+1) = roi;
wd = d+2;

% exchange background and object
w = 1-w;

fixset = [];

for k = 1:3
    % label each slice
	for i = 1:wd(k)
        switch k
        case 1
            im = reshape(w(i,:,:),wd(2),wd(3));        
        case 2
            im = reshape(w(:,i,:),wd(1),wd(3));        
        case 3
            im = reshape(w(:,:,i),wd(1),wd(2));        
        end
        [lab, num(i)] = bwlabel(im,4);
        L{i} = lab;        
	end
    % original roi background should be labeled as '1'
    ind = find(num>1);
    for i = 1:length(ind);
        lab = L{ind(i)};
        for j = 1:num(ind(i))
            cnt(j) = length(find(lab == j));
        end
        [mval, mind] = max(cnt);
        if (mind > 1)
            disp('Make background labeled as 1');
            ind1 = find(lab == 1);
            ind2 = find(lab == mind);
            lab(ind2) = 1; % background
            lab(ind1) = mind;
            L{ind(i)} = lab;
        end
    end
    % check multiple component slices
    for i = 1:length(ind);
        lab = L{ind(i)};
        for j = 2:num(ind(i))
            % connect to background in the previous slice
            con_prev = find(lab == j & L{ind(i)-1} == 1);
            % connect to background in the next slice
            con_next = find(lab == j & L{ind(i)+1} == 1);
            if (~isempty(con_prev) & ~isempty(con_next))
                [xs, ys] = find(lab == j); 
                switch k
                case 1
                    im = reshape(w(ind(i),:,:),wd(2),wd(3));
                    sigind = sub2ind(wd,xs*0+ind(i),xs,ys);
                case 2
                    im = reshape(w(:,ind(i),:),wd(1),wd(3));        
                    sigind = sub2ind(wd,xs,xs*0+ind(i),ys);
                case 3
                    im = reshape(w(:,:,ind(i)),wd(1),wd(2));        
                    sigind = sub2ind(wd,xs,ys,xs*0+ind(i));
                end
                w(sigind) = 0; % note that really it is assigned 1 to roi
                [xs,ys,zs] = ind2sub(wd,sigind);
                fixset(end+1:end+length(xs),:) = [[xs,ys,zs]-1, xs*0+1, xs*0+infoid];
            end
        end        
    end
    % check hole between slices
    for i = 2:(wd(k)-2)
        overlap = L{i}*0;
        overlap_bkgd = find(L{i} == 1 & L{i+1} == 1);
        overlap(overlap_bkgd) = 1;
        [lab, m] = bwlabel(overlap,4);
        if (m>1)
            % deal with only the minimum one;
            for j = 1:m
                cnt(j) = length(find(lab == j));
            end
            [mval, mind] = min(cnt);
            % connect to background in the previous slice
            con_prev = find(lab == mind & L{i-1} == 1);
            % connect to background in the current slice
            con_curr = find(lab == mind & L{i} == 1);
            % connect to background in the next slice
            con_next = find(lab == mind & L{i+1} == 1);
            % connect to background in the next next slice
            con_next_next = find(lab == mind & L{i+2} == 1);
            if ((~isempty(con_curr) | ~isempty(con_prev)) & ...
                (~isempty(con_next) | ~isempty(con_next_next)))
                [xs, ys] = find(lab == mind); 
                switch k
                case 1
                    sigind = [sub2ind(wd,xs*0+i,xs,ys);sub2ind(wd,xs*0+i+1,xs,ys)];
                case 2
                    sigind = [sub2ind(wd,xs,xs*0+i,ys);sub2ind(wd,xs,xs*0+i+1,ys)];
                case 3
                    sigind = [sub2ind(wd,xs,ys,xs*0+i);sub2ind(wd,xs,ys,xs*0+i+1)];
                end
                w(sigind) = 0; % note that really it is assigned 1 to roi
                [xs,ys,zs] = ind2sub(wd,sigind);
                fixset(end+1:end+length(xs),:) = [[xs,ys,zs]-1, xs*0+1, xs*0+infoid];
            end
        end    
    end
    L = [];
    num = [];
end

% make roi consistent with w
roi = 1 - w(2:d(1)+1,2:d(2)+1,2:d(3)+1);

return;


%
% select a point to change its value 
% according to maximum number of different neighbours
%

function [start, val] = decide_squa(w,k)

for i = 1:3:10
    cnt((i-1)/3+1) = nb_sum(w, k(i), k(i+1), k(i+2));
end

cnt(1:2) = 6 - cnt(1:2);

[v, in] = max(cnt);

start = (in-1)*3+1;
if (in < 3)
    val = 0;
else
    val = 1;
end

return;


%
% calculate the sum of the direct neighbour
%
function b = nb_sum(w, x, y, z)

b = sum([w(x-1,y,z), ...
         w(x+1,y,z), ...
         w(x,y-1,z), ...
         w(x,y+1,z), ...
         w(x,y,z-1), ...
         w(x,y,z+1)]);

return;





