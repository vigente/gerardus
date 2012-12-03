function [eig1, eig2, eig3] = structure_tensor(im)

% init memory for output
eig1 = zeros(size(im));
eig2 = zeros(size(im));
eig3 = zeros(size(im));

% compute image gradients in x, y, and z directions
[dx, dy, dz] = gradient(im);

% debug: loop to compute eigenvalues using eig() on the structure tensor
for I = 1:numel(im)
    % structure tensor
    s = [...
        dx(I)*dx(I) dx(I)*dy(I) dx(I)*dz(I); ...
        dy(I)*dx(I) dy(I)*dy(I) dy(I)*dz(I); ...
        dz(I)*dx(I) dz(I)*dy(I) dz(I)*dz(I) ...
        ];
    
    % compute eigenvalues
    [~, d] = eig(s);
    
    % sort eigenvalues
    d = sort(diag(d), 'descend');
    
    % assign eigenvalues
    eig1(I) = d(1);
    eig2(I) = d(2);
    eig3(I) = d(3);
end

% gradient crossed terms
p = (dx .* dy .* dx .* dy) ...
    + (dx .* dz .* dx .* dz) ...
    + (dy .* dz .* dy .* dz);

% special case when p==0 has to be computed differently
pidx = (p == 0);

% compute eigenvalues for special case p==0
eig1(pidx) = dx(pidx) .* dx(pidx);
eig2(pidx) = dy(pidx) .* dy(pidx);
eig3(pidx) = dz(pidx) .* dz(pidx);

% save memory and speed computations up by keeping only the elements that
% are not special case p==0
pidx = ~pidx;
p = p(pidx);
dx = dx(pidx);
dy = dy(pidx);
dz = dz(pidx);

% trace of the gradient matrix for each voxel / 3
q = (dx + dy + dz)/3;

% auxiliary variable
p = ((dx-q) .* (dx-q) ...
    + (dy-q) .* (dy-q) ...
    + (dz-q) .* (dz-q)) + 2 * p;
p = sqrt(p / 6);

% B = (1 / p) * (A - q * I)
% det(B)
r = (dx.*dx + dy.*dy + dz.*dz - q) .* q .* q ./ (2 * p);

% angle
phi = acos(r) / 3;

% fix numerical errors so that  -1 <= r <= 1
phi(r <= -1) = pi / 3;
phi(r >= 1) = 0;

% compute eigenvalues
eig1(pidx) = q + 2 .* p .* cos(phi);
eig3(pidx) = q + 2 .* p .* cos(phi + 2*pi/3);
eig2(pidx) = 3 * q - eig1(pidx) - eig3(pidx);
