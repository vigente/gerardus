function O = GetOPTI(B,opts)
%Build OPTI Object from SymBuilder Object
%
%   Called By SYMBUILDER Class

%   Copyright (C) 2012-2013 Jonathan Currie (I2C2)

if(nargin < 2), opts = symbset; else opts = symbset(opts); end

%Check object is built
if(~IsBuilt(B))
    error('Please build the object using Draft() or Build()');
end

%Check we have an objective (no good without!)
if(B.noObjs == 0)
    error('You can only construct an OPTI object when you have specified an objective!');
elseif(B.noObjs > 1)
    error('You can only construct an OPTI object with a single objective');
end

%Set Option Defaults
use1stDerivs = true;
use2ndDerivs = true;
preallocate = true;

%Read Options
if(strcmpi(opts.use1stDerivs,'no')), use1stDerivs = false; end
if(strcmpi(opts.use2ndDerivs,'no')), use2ndDerivs = false; end
if(strcmpi(opts.preallocate,'no')), preallocate = false; end

%Clear existing functions from memory
clear symb_obj symb_grad symb_nlcon symb_nljac symb_hess

%Separate objective and constraint indices
olin = B.objlin(B.indobj);
clin = B.conlin; clin(B.indobj) = [];

%Get Bounds(common to all problems)
if(all(isinf(B.lb))), lb = []; else lb = B.lb; end
if(all(isinf(B.ub))), ub = []; else ub = B.ub; end

%Default options
sOpts = opts.solverOpts;
opts = optiset('warnings','off','solver',opts.solver,'display',opts.display);

%Check for Linear Problem
if(olin == 1 && all(clin == 1))
    %Get linear objective and constraints
    f = GetLinObj(B);
    [A,rl,ru] = GetLinCon(B);
    %Build OPTI object
    O = opti('f',f,'lin',A,rl,ru,'bounds',lb,ub,'xtype',B.xtype,'options',opts);
%Check for Quadratic Pprogram    
elseif(olin == 2 && all(clin == 1))
    %Get quadratic objective and linear constraints
    [H,f] = GetQuadObj(B);
    [A,rl,ru] = GetLinCon(B);
    %Build OPTI object
    O = opti('qp',H,f,'lin',A,rl,ru,'bounds',lb,ub,'xtype',B.xtype,'options',opts);
%Check for QCQP
elseif(olin == 2 && all(clin <= 2))
    warning('OPTI:SymbQCQP','This interface does not currently support QCQPs explicity, treating as NLP');
    %Get nonlinear objective and nonlinear constraints
    [obj,grad] = GetAllObj(B);
    [nlcon,cl,cu,nljac,nljacstr] = GetAllCon(B);
    %Check if we want derivatives
    if(~use1stDerivs)
        grad = []; nljac = []; nljacstr = []; lH = []; lHstr = [];
    elseif(~use2ndDerivs)
        lH = []; lHstr = [];
    else
        [lH,lHstr] = GetHessLag(B);
    end
    %Build OPTI object
    O = opti('obj',obj,'grad',grad,'nl',nlcon,cl,cu,'nljac',nljac,'nljacstr',nljacstr,'hess',lH,'Hstr',lHstr,'bounds',lb,ub,'xtype',B.xtype,'options',opts);
%Check for NLP
elseif(olin <= 3 && all(clin <= 3))
    %Get nonlinear objective & Hessian
    if(use1stDerivs)
        fprintf('Generating Objective & Gradient....');
        [obj,grad] = GetAllObj(B);
    else
        fprintf('Generating Objective....');
        obj = GetAllObj(B);
        grad = [];
    end
    fprintf('Done\n');
    %Check if we want derivatives
    if(~use1stDerivs)
        grad = []; lH = []; lHstr = [];
    elseif(~use2ndDerivs)
        lH = []; lHstr = [];
    else
        fprintf('Generating Hessian....');
        [lH,lHstr] = GetHessLag(B);
        fprintf('Done\n');
    end
    %Check for variations of constraints
    if(all(clin == 0)) %no constraints
        O = opti('obj',obj,'grad',grad,'hess',lH,'Hstr',lHstr,'bounds',lb,ub,'xtype',B.xtype,'options',opts);
    elseif(all(clin == 1)) %only linear
        [A,rl,ru] = GetLinCon(B); %OPTI takes care of this gradient
        O = opti('obj',obj,'grad',grad,'lin',A,rl,ru,'hess',lH,'Hstr',lHstr,'bounds',lb,ub,'xtype',B.xtype,'options',opts); 
    else %quadratic, nonlinear, or mixed constraints
        if(use1stDerivs)
            fprintf('Generating Constraints & Jacobian....');
        else
            fprintf('Generating Constraints....');
        end
        if(preallocate)
            conopt.preallocate = 1;
        else
            conopt.preallocate = 0;
        end
        [nlcon,cl,cu,nljac,nljacstr] = GetAllCon(B,conopt); 
        fprintf('Done\n');
        if(~use1stDerivs), nljac = []; nljacstr = []; end
        %If we are building an MINLP, add BONMIN linear information
        if(any(B.xtype ~= 'C'))
            %Constraint Linearity
            cons_lin = clin;
            cons_lin(clin > 1) = 0;
            %Variable linearity
            varlin = ones(size(B.jac,2),1);
            nlvars = symvar(B.jac); %Variables in the jacobian are nonlinear WHAT ABOUT GRADIENT?
            for i = 1:length(nlvars)
                varlin(logical(nlvars(i) == B.vars)) = 0;
            end
            opts = optiset(opts,'warnings','off','solverOpts',bonminset(sOpts,'cons_lin',cons_lin,'var_lin',varlin));
        else
            opts = optiset(opts,'warnings','off');
        end
        O = opti('obj',obj,'grad',grad,'nl',nlcon,cl,cu,'nljac',nljac,'nljacstr',nljacstr,'hess',lH,'Hstr',lHstr,'bounds',lb,ub,'xtype',B.xtype,'options',opts);
    end
else
    error('problem type not yet supported');
end