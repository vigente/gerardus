%*****************************************************************************
% DSDP5:  Dual-Scaling Algorithm for Positive Semidefinite Programming
% Copyright (c) 2004 by
% S. J. Benson, Y. Ye
% Last modified: 20 October 2004
%*****************************************************************************
% check;
%
% This script tests the DSDP solver using the examples. Compare the output
% to the file check.out
% 
% diary on;
ntrials=3;
BB=cell(1);
for trials = [1:ntrials];
  N = 10;
  B = zeros(N); 
  for i=1:N,
    for j=1:N, B(i,j)=mod(i+j,trials+1); end;
  end;
  fprintf('MAXCUT \n');
  [y,X,obj] = maxcut(B); 
  fprintf('GPP \n');
  [y,X,obj] = gpp(B,1); 
  fprintf('ETP \n');
  [obj,d]=etp(B*B'+speye(N)');
  fprintf('THETA\n');
  [obj,X]=thetaproblem(B);
end;

%[X,y] = randdinf([10 4 3], 20);

OPTIONS=doptions;
OPTIONS.print=10;

fprintf('CONTROL1\n');
[AC,b] = readsdpa('control1.dat-s');
[STAT,y,X]=dsdp(b,AC,OPTIONS);

fprintf('ARCH0\n');
[AC,b] = readsdpa('arch0.dat-s');
[STAT,y,X]=dsdp(b,AC,OPTIONS);

% diary off;


%%
clc
b = [ 1 1 ]';
AAC = [ [ 1.0 0 0 ]' [ 0 0 1.0 ]' [ 4.0 -1.0 5.0 ]'];
AC = [];
AC{1,1} = 'SDP';        %type of cone
AC{1,2} = [2];          %dimension of cone
AC{1,3} = sparse(AAC);  %cone data
[STAT,y,X]=dsdp(b,AC);
y
XX=dmat(X{1})

%%
clc
b = [ 1 1 ]';
AAC = [ [ 1.0 0 0 ]' [ 0 0 1.0 ]' [ 4.0 -1.0 5.0 ]'];
AAC1 = [ [ 2.0 0 0 ]' [ 0 0 2.0 ]' [ 4.0 -1.0 6.0 ]'];
AC = [];
AC{1,1} = 'SDP';        %type of cone
AC{1,2} = [2 2];          %dimension of cone
AC{1,3} = [sparse(AAC); sparse(AAC1)]  %cone data
% AC{2,1} = 'SDP';
% AC{2,2} = 2;
% AC{2,3} = sparse(AAC);
[STAT,y,X]=dsdp_old(b,AC);
y
XX=dmat(X{1})

%% min x, where [x sqrt2; sqrt2 x] >= 0 (pos sd)
clc
b = -[ 1 ]';
AAC = [ -[ 1 0 1 ]' [ 0 sqrt(2) 0 ]'];
AC = [];
AC{1,1} = 'SDP';
AC{1,2} = [2];
AC{1,3} = sparse(AAC);
[STAT,y,X]=dsdp(b,AC);
y
XX=dmat(X{1})
[pobj,dobj,err]=derror(STAT,y,X,b,AC)

S = [0 sqrt(2); sqrt(2) 0] - full(dmat([-1 0 -1]'))*y

%%
clc
f = -1;
% sdp = {triu(-1*speye(2)); triu(rot90(sqrt(2)*speye(2)))};

sdp = [];
sdp.dim = 2;
sdp.blk = sparse([[0 sqrt(2) 0]' -[1 0 1]'])

[x,~,~,e] = dsdp(f,[],[],[],[],[],[],sdp)

%% repeat as nl problem
f = @(x) x;
c = @(x) eig([x sqrt(2); sqrt(2) x]);

opt = opti('fun',f,'nl',c,[0;0],[Inf;Inf],'ndec',1,'x0',1);

[x,f,e,i] = solve(opt)

%%
clc
C{1} = [1 0; 
           0 1];
At{1,1} = [0 1; 
           1 0]; 
At{1,2} = [1 1; 
           1 1];
b = [1; 2]; 
blk{1,1} = 's'; blk{1,2} = 2; 

AC{1,1} = 'SDP';
AC{1,2} = [2];
AC{1,3} = [dvec(At{1}) dvec(At{2}) dvec(C{1})];

[STAT,y,X]=dsdp(b,AC); 
y
XX=dmat(X{1})
[pobj,dobj,err]=derror(STAT,y,X,b,AC)

%%
clc
clear
DELTA = 1e-4;
C{1} = [ 0 1/2 0;
        1/2 DELTA 0;
        0 0 DELTA ];
A{1,1} = [ 0 -1/2 0;
        -1/2 0 0;
        0 0 0];
 A{2,1} = [ 1 0 0 ;
        0 0 0;
        0 0 0];
A{3,1} = [ 0 0 1;
        0 0 0;
        1 0 0];
A{4,1} = [ 0 0 0;
            0 0 1;
            0 1 0];
b = [1; 2*DELTA; 0; 0];

AC{1,1} = 'SDP'
AC{1,2} = [3];
AC{1,3} = [dvec(A{1}) dvec(A{2}) dvec(A{3}) dvec(A{4}) dvec(C{1})]

blk{1,1} = 's'; blk{1,2} = 3;

[STAT,y,X]=dsdp(b,AC);
y
XX=dmat(X{1})
[pobj,dobj,err]=derror(STAT,y,X,b,AC)


%% Design of a better format
clc
clear

% Standard form SDP
% SDP [Primal]
% min c'x
% st. X = F1*x1 + F2*x2 + ... + Fn*xn - F0
%     X >= 0 [positive semidefinite]

% SDP [Dual]    
% max <F0,Y>                                <> = trace of product of args [inner product]
% st. <F1,Y> = b1
%     <F2,Y> = b2
%       ...
%     Y >= 0 [positive semidefinite]   


% Alternative Primal-Dual Pair of SDP [when F1..Fn is not linearly independent, use dual form below]
% -A = F
% -b = c
%  X = Y, y = x, Z = X

% SDP [Primal]
% min <C,X>
% st. <A1,X> = b1
%     <A2,X> = b2
%       ...
%     X >= 0 [positive semidefinite]

% SDP [Dual] - form used by DSDP
% max b'y
% st. A1*y1 + A2*y2 + .. + An*yn + Z = A0
%     Z >= 0 [positive semidefinite]     


% p number of cones
% m number of variables

% b = m x 1
% X = m x m x p 

% New optiprob arg: sdcone
% sdcone = A cell array {F0;F1;F2;...;Fn}, each column is a new cone
% Alternatively written as {C;A1;A2;...;An} (automatically converted)

% New Calling Format [preprocessed sdcone into sparse, triu args]
% [y,X,info] = dsdp(f,A,b,Aeq,beq,lb,ub,sdcone,x0,opts)

%% YALMIP Example
x = sdpvar(1,1);
z = sdpvar(1,1);

X = [x 2;2 z];

F = set('X>=0');
% F = F+set('x==0');
F = F+set('z>=0');
F = F+set('x<=10');
F = F+set('z<=10');
sol = solvesdp(F,x+z)

%%
clc
yalmip('clear')
x = sdpvar(1,1);

X = [x sqrt(2); sqrt(2) x];
F = set('X>=0');

sol = solvesdp(F,x)

double(x)
double(X)

%%
clc
yalmip('clear')
x = sdpvar(1,1);

X = [x sqrt(2); sqrt(2) x];
F = set('X>=0');
F = F+set('1.5<=x<=1.5');

sol = solvesdp(F,x)

double(x)
double(X)

%% Just SDP
clc
f = -1;
% sdp.dims = 2'
% sdp.F = {triu(1*speye(2)); triu(rot90(sqrt(2)*speye(2)))};

ind = tril(ones(2))==1;

A = speye(2);
C = rot90(sqrt(2)*speye(2));

sdp = sparse([C(ind) -A(ind)])
opts.display=2;
[x,ff,e,i] = dsdp(f,[],[],[],[],sdp,[],opts)


%% Add Bounds
clc
f = -1;
% sdp.dims = 2'
% sdp.F = {triu(1*speye(2)); triu(rot90(sqrt(2)*speye(2)))};

ind = tril(ones(2))==1;

A = speye(2);
C = rot90(sqrt(2)*speye(2));

lb = 0;
ub = 1.42;

sdp = sparse([C(ind) -A(ind)])
[x,ff,~,e] = dsdp(f,[],[],[],[],lb,ub,sdp)

%% Add Linear Constraints
clc
f = -1;
% sdp.dims = 2'
% sdp.F = {triu(1*speye(2)); triu(rot90(sqrt(2)*speye(2)))};

ind = tril(ones(2))==1;

A = speye(2);
C = rot90(sqrt(2)*speye(2));

Ain = sparse(-1);
bin = -1.415;

lb = 0;
ub = 1.42;

sdp = sparse([C(ind) -A(ind)])
[x,ff,~,e] = dsdp(f,Ain,bin,[],[],lb,ub,sdp)

%% Add Linear Eq Constraints
clc
f = -1;
% sdp.dims = 2'
% sdp.F = {triu(1*speye(2)); triu(rot90(sqrt(2)*speye(2)))};

ind = tril(ones(2))==1;

A = speye(2);
C = rot90(sqrt(2)*speye(2));

Ain = sparse([-1]);
bin = -1.415;

Aeq = sparse(1);
beq = 1.5;

lb = 0;
ub = 2.42;

sdp = sparse([C(ind) -A(ind)])
% for i = 1:1000
[x,p,d,e] = dsdp(f,Ain,bin,Aeq,beq,lb,ub,sdp)
% end

