function X=RandSampleSphere(N,spl)
% Generate a uniform or stratified sampling of a unit sphere.
%
% INPUT ARGUMENTS:
%   - N   : desired number of point samples. N=200 is default.
%   - spl : can be 'uniform' or 'stratified'. The former setting is the
%           default.
%
% OUTPUT:
%   - X  : N-by-3 array of sample point coordinates.
%
% REFERENCE:
%   - Shao & Badler, 1996, Spherical Sampling by Archimedes' Theorem
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: May.2013
%

% Default arguments
if nargin<1 || isempty(N), N=200; end
if nargin<2 || isempty(spl), spl='uniform'; end

% Basic error checking
chk=strcmp(spl,{'uniform','stratified'});
if sum(chk)==0
    error('Invalid sampling option')
end

N=round(N);
if numel(N)~=1 || ~isnumeric(N) || N<1
    error('Invalid entry for 1st input argument')
end
if N<3, spl='uniform'; end

% Sample the unfolded right cylinder
if strcmp(spl,'stratified')
    
    % Partition the [-1,1]x[0,2*pi] domain into ceil(sqrt(N))^2 subdomains
    % and then draw a random sample for each
    n=ceil(sqrt(N));
    ds=2/n;
    [Xc,Yc]=meshgrid((-1+ds/2):ds:(1-ds/2));
    
    x=ds*(rand(n^2,1)-0.5);
    y=ds*(rand(n^2,1)-0.5);
    
    x=x+Xc(:);
    y=y+Yc(:);
    clear Xc Yc
    
    % Remove excess samples
    R=n^2-N;
    if R>0
        idx=randperm(n^2,R);
        x(idx)=[];
        y(idx)=[];
    end
    
    lon=(x+1)*pi;
    z=y;
    
else
    z=2*rand(N,1)-1;
    lon=2*pi*rand(N,1);
end

% Convert z to latitude
lat=acos(z);

% Convert spherical to rectangular co-ords
x=cos(lon).*sin(lat);
y=sin(lon).*sin(lat);
        
X=[x,y,z];
