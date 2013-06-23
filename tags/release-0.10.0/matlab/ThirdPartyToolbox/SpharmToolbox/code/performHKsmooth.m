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

% ============================================
% performHKsmooth.m
%
% Goal: Heat kernal smoothing, preprocessing for t-test
% input_signal: signal to smooth.
% tri         : triangle indices
% coord       : node coordinates
% nbr         : 1st neighbor nodes list
% sigma       : bandwidth
% n_smooth    : number of iterations
%
% EXAMPLE: bandwidth sigma=1, number of iteration=50
% output_signal=hk_smooth(input_signal,tri,coord,nbr,1,50);
%
% Mesh topology:
% There are total V + F rows.
% V = number of vertices
% F = number of faces
% It is assumed that the input mesh is topologically equvalent to a sphere, i.e.
% 2V - F = 4: Euler characteristic equation.
%
% Li Shen 
% 01/10/2005 - create

function output_signal = performHKsmooth(input_signal,vertices,faces,nbr,FWHM)

if isempty(FWHM) | FWHM<=0
    output_signal = input_signal;
    return;
end

sigma = FWHM/(2*sqrt(log(4)));
disp(sprintf('FWHM %f => sigma %f',FWHM,sigma));

if ismember(size(vertices,1),[42,162,642,2562,10242])
    output_signal = smooth_one_side(input_signal,vertices,faces,nbr,sigma);
else
    vn = size(vertices,1)/2; fn = size(faces,1)/2;
    output_signal = [smooth_one_side(input_signal(1:vn),     vertices(1:vn,:),     faces(1:fn,:),nbr,sigma); ...
                     smooth_one_side(input_signal((1:vn)+vn),vertices((1:vn)+vn,:),faces(1:fn,:),nbr,sigma)];
end

return;

%
% process one side
%

function output_signal = smooth_one_side(input_signal,coord,tri,nbr,sigma)

n_points = size(coord,1);
n_tri = size(tri,1);

n_smooth = 50; % number of iterations
sigma = sigma/sqrt(n_smooth);
output_signal = hk_smooth(input_signal,tri,coord',nbr,sigma,n_smooth);

return;

%
% Moo K. Chung's code
%

function output_signal=hk_smooth(input_signal,tri,coord,nbr,sigma,n_smooth)
%output_signal=hk_smooth(input_signal,tri,coord,nbr,sigma,n_smooth)
%
% input_signal: signal to smooth.
% tri         : triangle indices
% coord       : node coordinates
% nbr         : 1st neighbor nodes list
% sigma       : bandwidth
% n_smooth    : number of iterations
%
% EXAMPLE: bandwidth sigma=1, number of iteration=50
% output_signal=hk_smooth(input_signal,tri,coord,nbr,1,50);
%
% (C) Moo K. Chung
% Update history Feb 5. 2004. July 20, 2004. Oct 28, 2004.
% mchung@stat.wisc.edu
% http://www.stat.wisc.edu/softwares/hk/hk.html
%
% If you use this code, please reference one of the following papers. The details on
% the mathematical basis of of the algorithm can be found in these papers.
%
% [1] Chung, M.K., Robbins,S., Dalton, K.M., Davidson, R.J., Evans, A.C. (2004) 
%     Cortical thickness analysis in autism via heat kernel smoothing. NeuroImage, submitted. 
%     http://www.stat.wisc.edu/~mchung/papers/ni_heatkernel.pdf
%
% [2] Chung, M.K. (2004) Heat kernel smoothing and its application to cortical manifolds. 
%     Technical Report 1090. Department of Statististics, Universisty of Wisconsin-Madison. 
%     http://www.stat.wisc.edu/~mchung/papers/heatkernel_tech.pdf

%heat kernel
K=inline('exp(-x/(2*sigma^2))/sum(exp(-x/(2*sigma^2)))');

%smoothing weight computation
n_vertex=size(nbr,1);
weight=zeros(n_vertex,7);
for i=1:n_vertex
    if (i<=12)
        current_nbr=nbr(i,1:5);
        pipo=coord(:,current_nbr(1:5))-kron(ones(1,5),coord(:,i));
        distance=sum(pipo.*pipo);
        weight(i,:)=[K(sigma,[0 distance]) 0];
    else
        current_nbr=nbr(i,1:6);
        pipo=coord(:,current_nbr(1:6))-kron(ones(1,6),coord(:,i));
        distance=sum(pipo.*pipo);
        weight(i,:)=K(sigma,[0 distance]);
    end;
end;

%iterated smoothing
signal=input_signal;
for j=1:n_smooth
    Y=[[signal(nbr(1:12,1:5)) zeros(12,1)]; signal(nbr(13:n_vertex,:))];
    signal=sum(weight.*[signal Y],2);
end;

output_signal=signal;
return;

%
% get nbr (not done yet) 
% should be done before calling this function
%

function nbr = get_nbr(tri,coord)

n_points = size(coord,1);
n_tri = size(tri,1);

% compute the maximum degree of node
degree=zeros(n_points,1);
for j=1:n_tri
degree(tri(j,:))=degree(tri(j,:))+1;
end
max_degree=max(degree);

% find out the 1st neighbor nodes
nbr=zeros(n_points,max_degree);
for i_tri=1:n_tri
    for j=1:3
        cur_point = tri(i_tri,j);
        for k=1:3
            if (j ~= k)
                nbr_point= tri(i_tri,k);
                if find(nbr(cur_point,:)==nbr_point)
                    ;
                else
                    n_nbr = min(find(nbr(cur_point,:) == 0));
                    nbr(cur_point,n_nbr) = nbr_point;
                end;
            end;
        end;
    end;
end;

return;

%
% get sigma (not done yet) 
% should be done before calling this function
%

function sigma = get_sigma(vs,fs,sigfactor)

volume = hbm01_surfacevol(vs,fs); 
if size(vs,1)>642 volume = volume/2; end;
% sigma = volume^(1/3)/15;
% sigma = volume^(1/3)/30;
sigma = volume^(1/3)/sigfactor;
disp(sprintf('volume: %f; sigma: %f',volume,sigma));

return;