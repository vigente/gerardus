function [fun,lb,ub,nlcon,cl,cu,xtype] = GetBARON(B)
%Build BARON arguments from SymBuilder Object
%
%   Called By SYMBUILDER Class

%   Copyright (C) 2012 Jonathan Currie (I2C2)

%Check object is built
if(~IsBuilt(B))
    error('Please build the object using Draft() or Build()');
end

%Check we have an objective (no good without!)
if(B.noObjs == 0)
    error('You can only get BARON arguments when you have specified an objective!');
elseif(B.noObjs > 1)
    error('You can only get BARON arguments with a single objective');
end

%Clear existing functions from memory
clear symb_grad symb_hess symb_nlcon symb_nljac symb_obj

%Objective
fun = GetAllObj(B);
%Bounds
lb = B.lb;
ub = B.ub;
%Nonlinear Constraints
opts.preallocate = false;
[nlcon,cl,cu] = GetAllCon(B,opts); 
%Integer Constraints
xtype = B.xtype;