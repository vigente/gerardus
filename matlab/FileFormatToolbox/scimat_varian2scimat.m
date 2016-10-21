function scimat = scimat_varian2scimat(path)

% SCIMAT_VARIAN2SCIMAT Create SCIMAT struct from data obtained from Varian MR system.
%
% This function creates a struct with the scimat format (see "help scimat"
% for details) for data from Varian MR scanners. The output is a vector of scimat structures,
% one per 2D slice or a single scimat for 3D images. Depends on the size of the input
% image data. It follows some odd conventions for the rotations that only 
% became clear with the manual at hand which is not available online.
%
% SCIMAT = SCIMAT_VARIAN2SCIMAT(PATH)
%
%   PATH is a string of the path pointing the folder containing the image
%   data and the metainformation which is in the "procpar" file. The path
%   should end in .fid and the folder should contain a file called procpar
%   and another called image_mag.

% Author(s): Nicolas Basty <nicolas.basty@eng.ox.ac.uk>
% Darryl McClymont <darryl.mcclymont@cardiov.ox.ac.uk>
% Irvin Teh <irvin.teh@cardiov.ox.ac.uk>
% Copyright © 2016 University of Oxford
% Version: 0.1
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

% Read in the properties from the procpar file
props = init_params_nested(path); % Get properties from procpar
n = load_nii([path filesep 'image_mag.nii']); % Load nifti, function in third party toolbox
img = n.img; % getting image data from the nifti

% Determine (fairly arbitrarily) if the data is a 2D stack or 3D image
if size(img,1)/size(img,3) <5 && size(img,2)/size(img,3) < 5
    dim = 3;
else
    dim = 2;
end

% Euler angles & generate rotation matrix
psi = props.psi;
phi = props.phi;
theta = props.theta;
rotm = rotatematrix_varian(psi,theta,phi);

% Resolution
di = props.x_res;
dj = props.y_res;
dk = props.z_res;

% Field of view size in mm
FOVi = props.x_dim * di;
FOVj = props.y_dim * dj;
FOVk = props.z_dim * dk;

% rotate bottom corner of image to real world coordinates
min_xyz = rotm * (-[FOVi, FOVj, FOVk] / 2)';

% indices
[I, J, K] = ndgrid(1:size(img,1), 1:size(img,2), 1:size(img,3));
IJK = [I(:), J(:), K(:)];

% every pixel in the image gets its coordinates converted to real world
% coordinates. The first of every slice then gets passed on to the scimat
% min + rot * (index * resolution)
XYZ = bsxfun(@plus, min_xyz, rotm * (bsxfun(@times, IJK - 0.5, [di dj dk]))')';

X = reshape(XYZ(:,1), size(img));
Y = reshape(XYZ(:,2), size(img));
Z = reshape(XYZ(:,3), size(img));

% image data and rotation matrix are transposed here because of the varian
% convention

if dim == 2
    
    for i = 1:size(img,3)
        offset = [Y(1,1,i) X(1,1,i) Z(1,1,i)];
        scimat(i) = scimat_im2scimat(((img(:, :, i).')), [di dj dk], offset, rotm');
    end
    
elseif dim == 3
    
    % in case we do want a 3D image. (my HR images are 3D and I don't want a stack)
    offset = [Y(1,1,1) X(1,1,1) Z(1,1,1)];
    scimat = scimat_im2scimat(permute(img,[2 1 3]), [di dj dk], offset, rotm');
    
end

end

function rot = rotatematrix_varian(psi,theta,phi)
% Get rotation matrix for Varian convention
% Created N. Basty September 2016
% Define basic rotation matrices
rz = @(ang) [cosd(ang) sind(ang) 0; -sind(ang) cosd(ang) 0; 0 0 1];
% rx = @(ang) [1 0 0; 0 cosd(ang) sind(ang); 0 -sind(ang) cosd(ang)];
ry = @(ang) [cosd(ang) 0 sind(ang); 0 1 0; -sind(ang) 0 cosd(ang)];

rot = rz(psi)*ry(theta)*rz(-phi);
% according to the scanner's manual it should be X and not Y but this works
% in Matlab

end

function params = init_params_nested(path_full)

% Read procpar to setup reconstruction
% Add new variables to read in Items A & B
% Created I.Teh 21 Nov 2012

params = search_procpar_nested([path_full '/procpar'],...
    'seqfil','petable','pelist','pe_order','sampling','navigator','array','rcvrs','gain','ident',...
    'field_strength','orient','psf','v','v_str','ssc','simulation','seedarr','randseed',...
    'nseg','ns','np','nv','nv2','nv3','lro','lpe','lpe2','thk','te','etl','esp','pss','effTE','nsa','ne','nt',...
    'dwi','b0n','b1n','diffscheme','dro','dpe','dsl','diffpol','gdiff','tdelta','tDELTA',...
    'propeller','petable_split','philips_dwrot','nblades','undersample_pe3','undersample_blades',...
    'epiref_type','nread','nphase','nphase2','nnav','ss','zpe_flag',...
    'bidirectional','diff','psi','phi','theta','pescheme',...
    'rf1off','rf2off','rf3off','flip1','flip2','flip3','nshells',...
    'cs_flag','petable_array','petable_array2','nvnom','nvnom2','idx_array',...
    'bvalrr','bvalpp','bvalss','bvalrp','bvalrs','bvalsp','bvalue','cs_speedup','max_bval','mean_bval','image',...
    'pro','ppe','ppe2','nnav','dro2','dpe2','dsl2','anglelist','numblades','bladewiden'...
    );

params.nro=params.np/2;
if isempty(params.etl)
    params.etl=1;
end
params.npe=params.nv;
params.shots=params.nv/params.etl;                % # shots
params.rcvrs=length(strfind(params.rcvrs,'y'));   % # receivers
if isempty(params.nnav)
    params.nnav=0;
end

% Find length of arrays
params.image_n_all=1;
try
    if ~isempty(params.array)
        params.array(~isletter(params.array))=' ';
        params.array=regexp(params.array, '([^ ]*)', 'match');
        tmp=search_procpar_nested([path_full '/procpar'],params.array{1});
        tmp_name=fieldnames(tmp);
        params.image_n_all=length(getfield(tmp,tmp_name{1}));
    end
end

% Initialise base pulse sequence
k=strfind(params.seqfil,'_');
tmp=params.seqfil;
tmp(k)=[];
seqfil_base_array={'sems','semsdw','se3d','fsems','fse3d','epi','epip','me3d',...
    'radial','gems','ge2d','ge2dT1','ge3d','ge3dseg','mems','mgems','spuls','simulation','prescanpower'};
for idx=length(seqfil_base_array):-1:1
    k=strfind(tmp,seqfil_base_array{idx});
    if ~isempty(k)
        params.seqfil_base=seqfil_base_array{idx};
        break
    end
end

% Initialise sequence type based on seqfil (ie. 2D or 3D)
k=strfind(params.seqfil,'3d');
if ~isempty(k)
    params.seqfil_type='3d';
else
    params.seqfil_type='2d';
end

% EPI only
if ~isempty(params.nread)
    params.np=params.nread;
    params.nv=params.nphase;
    params.nv2=params.nphase2;
end

% 1D
if isempty(params.nv)
    params.nv=1;
end

% Get dimensions
if strcmp(params.cs_flag,'y')
    params.x_res = params.lro*10/(params.np/2);
    params.y_res = params.lpe*10/params.nvnom;
    params.x_dim = params.lro ./ params.x_res * 10;
    params.y_dim = params.lpe ./ params.y_res * 10;
    if params.nv2>0 %3D
        params.z_res = params.lpe2*10/params.nvnom2;
        params.z_dim = params.lpe2 ./ params.z_res * 10;
    else %2D
        params.z_res = params.thk;
        params.z_dim = params.ns;
    end
else
    params.x_res = params.lro*10/(params.np/2);
    params.y_res = params.lpe*10/params.nv;
    params.x_dim = params.lro ./ params.x_res * 10;
    params.y_dim = params.lpe ./ params.y_res * 10;
    if params.nv2>0 %3D
        params.z_res = params.lpe2*10/params.nv2;
        params.z_dim = params.lpe2 ./ params.z_res * 10;
    else %2D
        params.z_res = params.thk;
        params.z_dim = params.ns;
    end
end

% 1D
if isempty(params.bladewiden)
    params.bladewiden=1;
end

end

function params = search_procpar_nested(procpar_file,varargin)

% Search Varian procpar and extract as many variables as specified.
%
% Example:
% params=search_procpar(procpar_file,'np','nv','nv2');
% cellfun(@(n,v) assignin('caller',n,v),fieldnames(params),struct2cell(params));
% NB: The second line expands the 'params' structure into separate variables (Optional)
%
% Created, I.Teh, 23 Feb 2007
% Convert output variables into structure, I.Teh, 27 Nov 2012

% Load procpar
fid = fopen(procpar_file);
procpar = textscan(fid, '%s');
fclose(fid);
procpar = procpar{1};

% Go to correct idx depending whether data is from Varian or not
tmp = strcmp(procpar,'decpat2');
[C I] = find(tmp);
if ~isempty(C)  %if Varian data
    idx_add = 10;
else            % if not Varian data
    idx_add = 0;
end

% Find parameters
for idx = 1:length(varargin)
    variable = varargin{idx};
    tmp = strcmp(procpar,variable);
    [C I]=find(tmp);
    if ~isempty(C)
        C1=C(1)+1+idx_add;
        C2=C(1)+2+idx_add;
        num_elem = str2double(cell2mat(procpar(C1)));
        start = C2;
        for idx2 = start:start+num_elem-1
            tmp = cell2mat(procpar(idx2));
            if strcmp(tmp(1),'"')
                value = tmp(2:length(tmp)-1);
            else
                value(idx2-start+1) = str2double(tmp);
            end
        end
        varargout{idx} = value;
        clear value
    else
        varargout{idx} = [];
    end
end

% Convert strings to variables and assign values
for idx = 1:length(varargin)
    eval(sprintf('params.%s=varargout{idx};',varargin{idx}));
end

end

