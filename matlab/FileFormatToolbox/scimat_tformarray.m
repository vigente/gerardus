function [ scimat_out ] = scimat_tformarray( Sci_HR, new_rotmat, new_size, new_spacing, new_min)
%UNTITLED4 transforms a 3D scimat image using tformarray


% new_rotmat = old_rotmat;
% new_size = old_size;
% new_spacing = old_spacing;
% new_min = old_min;

old_rotmat = Sci_HR.rotmat;

old_spacing_x = Sci_HR.axis(1).spacing;
old_spacing_y = Sci_HR.axis(2).spacing;
old_spacing_z = Sci_HR.axis(3).spacing;
old_spacing = [old_spacing_x, old_spacing_y, old_spacing_z];

old_min_x = Sci_HR.axis(1).min;
old_min_y = Sci_HR.axis(2).min;
old_min_z = Sci_HR.axis(3).min;
old_min = [old_min_x, old_min_y, old_min_z];

%%

% tformarray indexes from 0, not 1
T0 = eye(4);
T0(1:3,4) = -1;

% convert to mm:
S1 = diag([old_spacing([2 1 3]) 1]);

% convert the min field, which is the top left coordinate, into voxel
% centers (offset)
% (y, x, z) -> (x, y, z)
old_offset = old_min([2 1 3]) + old_spacing([2 1 3])/2;

% translate to origin:
T1 = eye(4);
T1(1:3, 4) = old_offset;

% rotate input:
R1 = pinv(old_rotmat);
R1(4,4) = 1;

T_old = T1 * R1 * S1 * T0;

%%
% scale of new image:
S2 = diag([new_spacing([2 1 3]) 1]);

% rotate to new plane:
R2 = pinv(new_rotmat);
R2(4,4) = 1;

% new transalation
new_offset = new_min([2 1 3]) + new_spacing([2 1 3])/2;
T2 = eye(4);
T2(1:3, 4) = new_offset;

T_new = T2 * R2 * S2 * T0;

%%
% this assumes multiplying by the left - TF * [i j k 1]'
TF = pinv(T_new) * T_old;

if max(abs(TF(4, 1:3))) > 1E-10
    disp('Transformation function might be broken, T = ')
    disp(TF)
end
TF(end, :) = [0 0 0 1];

% maketform assumes multiplying from the right, so we need to transpose the
% transformation
T = maketform('affine', TF');
Res = makeresampler('nearest','bound');

% the images are stored permuted as YXZ, reverse this
data = permute(Sci_HR.data, [2 1 3]);

rotated_image = tformarray(data, T, Res, [1 2 3], [1 2 3], new_size([2 1 3]), [], []);

% and re-permute when saving again
scimat_out.data = permute(rotated_image, [2 1 3]);
scimat_out.rotmat = new_rotmat;

for i = 1:3
    scimat_out.axis(i).size = new_size(i);
    scimat_out.axis(i).spacing = new_spacing(i);
    scimat_out.axis(i).min = new_min(i);
end


