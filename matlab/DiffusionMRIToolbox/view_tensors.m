function view_tensors(D, HA, delta, dim_order)
% VIEW_TENSORS    Renders diffusion tensors in 3D, colour coded by HA
%
%
% Inputs:
%
%   D is the diffusion tensor from fit_DT in 2 dimensions only (i.e. a 3D
%   volume with [x, y, DT])
%
%   HA is the helix angle (in 2D). HA values of NaN are coloured grey. 
%
%   delta is the spacing between tensors (default 1.5E-3)
%
%   dim_order can be used (for example) when the tensors are referenced to
%   sagittal slices, whereas you want to display an axial slice. (default
%   [1 2 3])


% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright © 2014-2015 University of Oxford
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

% This function is an extension of plotDTI in the fanDTasia toolbox

% check arguments
narginchk(2, 4);
nargoutchk(0, 1);


if nargin < 4
    dim_order = [1 2 3]; % dimension
end
if nargin < 3
    delta = 1.5E-3; % spacing between tensors
end



% convert from 6x1 to 3x3 tensor representation
DT_33 = zeros([3 3 size(D,1) size(D,2)]);

DT_33(dim_order(1),dim_order(1),:,:) = D(:,:,1);
DT_33(dim_order(1),dim_order(2),:,:) = D(:,:,2);
DT_33(dim_order(1),dim_order(3),:,:) = D(:,:,3);
DT_33(dim_order(2),dim_order(1),:,:) = D(:,:,2);
DT_33(dim_order(2),dim_order(2),:,:) = D(:,:,4);
DT_33(dim_order(2),dim_order(3),:,:) = D(:,:,5);
DT_33(dim_order(3),dim_order(1),:,:) = D(:,:,3);
DT_33(dim_order(3),dim_order(2),:,:) = D(:,:,5);
DT_33(dim_order(3),dim_order(3),:,:) = D(:,:,6);


sz=size(DT_33);
if length(sz)==2
    nx=1;ny=1;
elseif length(sz)==3
    nx=sz(3);ny=1;
elseif length(sz)==4
    nx=sz(3);ny=sz(4);
end

ha = HA / 90;

hold on;
for i=1:nx
    for j=1:ny
        
        % Jet colourmap
        if isnan(ha(i,j))
            C = [0.8 0.8 0.8];
        else
            C = [red(ha(i,j)), green(ha(i,j)), blue(ha(i,j))];
        end
        
        
        % get the eigenvectors of the tensor
        [v,l]=eig(round(DT_33(:,:,i,j)*1E7)/1E7);
        
        % generate the ellipsoid
        [X,Y,Z]=ellipsoid(0,0,0,l(1,1),l(2,2),l(3,3),10);
        sz=size(X);
        for x=1:sz(1)
            for y=1:sz(2)
                A=[X(x,y) Y(x,y) Z(x,y)]';
                A=v*A;
                X(x,y)=A(1);
                Y(x,y)=A(2);
                Z(x,y)=A(3);
            end
        end
        X=X+(i-1)*delta*2;
        Y=Y+(j-1)*delta*2;
        
        % display the ellipsoid
        surf(real(X),real(Y),real(Z), 'faceColor', C, 'EdgeAlpha', 0, 'EdgeColor', 'none');
            
    end
end

axis equal
view([0 90]);
set(gca,'GridLineStyle','none')
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])

view([-1.2 0.5 0.8])
lighting phong
light('Position',[0 0 1],'Style','infinite','Color',[ 0.8 0.8 0.8]);

end


% these functions are to get the 'jet' colourmap (thanks stackexchange)
function val_out = interpolate( val,  y0,  x0,  y1,  x1 ) 
    val_out = (val-x0).*(y1-y0)./(x1-x0) + y0;
end

function val_out = base( val ) 
    if ( val <= -0.75 ) 
        val_out = 0;
    elseif ( val <= -0.25 ) 
         val_out = interpolate( val, 0.0, -0.75, 1.0, -0.25 );
    elseif ( val <= 0.25 ) 
        val_out =  1.0;
    elseif ( val <= 0.75 ) 
        val_out = interpolate( val, 1.0, 0.25, 0.0, 0.75 );
    else
        val_out = 0.0;
    end
end


function v = red( gray ) 
    v = base( gray - 0.5 );
end

function v = green( gray )
    v = base( gray );
end

function v = blue( gray )
    v = base( gray + 0.5 );
end