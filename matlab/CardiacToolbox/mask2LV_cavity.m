function [ LV, LV_cavity ] = mask2LV_cavity( mask )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

mask_orig = mask;

D = bwdist(mask);

sz = size(mask);

xmin = floor(sz(1)*3/8);
xmax = floor(sz(1)*5/8);
ymin = floor(sz(2)*3/8);
ymax = floor(sz(2)*5/8);
zmin = floor(sz(3)*3/8);
zmax = floor(sz(3)*5/8);

D_small = D(xmin:xmax, ymin:ymax, zmin:zmax);

[~, idx] = max(D_small(:));

[r, c, s] = ind2sub(size(D_small), idx);

r = r + xmin - 1;
c = c + ymin - 1;
s = s + zmin - 1; % the center of the blood pool

L = bwlabeln(1-mask);

ll = 0;
while L(r, c, s) == L(1,1,1)
    mask = imdilate(mask, ones([3 3 3]));
    
    L = bwlabeln(1-mask);
    
    ll = ll + 1;
    
end

LV_cavity = L == L(r, c, s);

air = L == L(1,1,1);

for i = 1:ll
    air = imdilate(air, ones([3 3 3]));
end

for l = 1:3
    LV_cavity = imdilate(LV_cavity, ones([3 3 3]));
    LV_cavity(mask_orig) = 0;
    LV_cavity(air) = 0;
end

LV = mask_orig .* (1-LV_cavity);


end

