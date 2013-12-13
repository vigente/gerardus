%% SOS Testing
%% LP Solve Example SOS1
clc
f = [-1 -1 -3 -2 -2]';
A = [-1 -1 1 1 0;
     1 0 1 -3 0];
b = [30;30];
ub = [40;1;inf;inf;1];

sostype = '12';
sosind = {[1:2]' [3:5]'};
soswt = {[1:2]' [3:5]'};

% [x,fval,ef,stat] = cplexmilp(f,A,b,[],[],sos,sosind,soswt,[],ub)

% SOS 1 OPTI Version
opts = optiset('solver','scip');

prob = optiprob('f',f,'ineq',A,b,'ub',ub,'sos',sostype,sosind,soswt,'options',opts)

sos = [];
sos.type = sostype;
sos.index = sosind;
sos.weight = soswt;
prob2 = optiprob('f',f,'ineq',A,b,'ub',ub,'sos',sos,'options',opts)

Opt = opti(prob,opts)

[x,fval,ef,stat] = solve(Opt)

%% LP Solve Example SOS2
clc
f = [-1 -1 -3 -2 -2]';
A = sparse([-1 -1 1 1 0;
             1 0 1 -3 0]);
b = [30;30];

lb = zeros(5,1);
ub = [40;1;inf;inf;1];

sos = '1';
sosind = (1:5)';
soswt = (1:5)';

% [x,fval,ef,stat] = cplexmilp(f,A,b,[],[],sos,sosind,soswt,lb,ub,'IIIII');
% x
% stat

% SOS 2 OPTI Version
opts = optiset('solver','cbc');

prob = optiprob('f',f,'ineq',A,b,'bounds',lb,ub,'sos',sos,sosind,soswt,'options',opts)

Opt = opti(prob,opts)

[x,fval,ef,stat] = solve(Opt)


% [x,fval,ef,stat] = solve(Opt);
% x
% fval
% stat
