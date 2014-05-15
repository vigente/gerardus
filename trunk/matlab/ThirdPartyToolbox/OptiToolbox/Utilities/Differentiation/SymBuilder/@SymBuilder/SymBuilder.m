classdef SymBuilder < handle
%SYMBUILDER  Create a SYMBUILDER object for creating symbolic optimization problems
%
%   B = SymBuilder() creates a blank optimization problem
%
%   See also SymBuilder.AddObj SymBuilder.AddCon SymBuilder.Build
%
%   Copyright (C) 2012 Jonathan Currie (www.i2c2.aut.ac.nz)
    
    properties(SetAccess=private)
        vars        %Symbolic array of symbolic variables
        sobj        %Symbolic matrix of equations
        jac         %Symbolic jacobian of equations
        hess        %Structure of hessians of each equation
        hesslag     %Symbolic hessian of the lagrangian  
        Opt         %OPTI object
    end
    
    properties(SetAccess=private)%,GetAccess=private)
        eqs         %Equations string
        cl          %Constraint Lower Bounds
        cu          %Constraint Upper Bounds
        constnt     %Cell array of constants
        exprsn      %Cell array of expressions
        bnds        %Cell array of bounds
        lb          %Double array of lower bounds
        ub          %Double array of upper bounds
        vartypes    %Cell array of variable types
        xtype       %Char array of variable types
        indobj      %Index of objectives in equations
        objlin      %Index of linear, quadratic and nonlinear objectives
        conlin      %Index of linear, quadratic and nonlinear constraints
        noObjs      %Counter of number objectives
        noCons      %Counter of number of constraints
        bldstat     %Build Status
        bldtime     %Build Time
        resgrp      %Cell array of result groups
        resexp      %Cell array of result expressions
    end
    
    methods
        function B = SymBuilder()
            %Initialization
            B.noObjs = 0;
            B.noCons = 0;
            Unbuild(B);            
        end
                
        function B = AddObj(B,str)
        %Add a General Objective (Linear or Nonlinear)
        
            %Concatenate with existing equations
            if(isempty(B.eqs))
                B.eqs = str;
                B.cl = NaN; %place holders
                B.cu = NaN;
            else
                B.eqs = sprintf('%s;%s',B.eqs,str);
                B.cl = [B.cl;NaN];
                B.cu = [B.cu;NaN];
            end
            if(isempty(B.indobj))
                B.indobj = length(B.cl);
            else
                B.indobj = [B.indobj;length(B.cl)];
            end
            %Add no objs
            B.noObjs = B.noObjs + 1;  
            %Indicate Rebuild Required
            Unbuild(B);
        end
                
        function B = AddCon(B,str_in)
        %Add a General Constraint (Linear or Nonlinear)
        
            %Parse constraint string
            [str,l,u] = SymBuilder.parseConstraint(str_in);            
            %Concatenate with existing equations if not empty
            if(~isempty(str))
                if(isempty(B.eqs))
                    B.eqs = str;
                    B.cl = l;
                    B.cu = u;
                else
                    B.eqs = sprintf('%s;%s',B.eqs,str);
                    B.cl = [B.cl;l];
                    B.cu = [B.cu;u];
                end            
                %Add no cons
                B.noCons = B.noCons + 1;  
                %Indicate Rebuild Required
                Unbuild(B);
            else
                optiwarn('symb:emptycon','The following constraint will be ignored by SymBuilder as both sides appear to be numeric.\n  %s\n',str_in);
            end
        end
        
        function B = AddConstraint(B,str)
        %Same as AddCon()
            B = AddCon(B,str);
        end
                
        function B = AddConstant(B,str,val,varargin)
        %Identify variable as a constant 
        
            %Concatenate with existing constants
            B.constnt = [B.constnt; {str val}];   
            %If user has passed multiple constants, process all
            if(nargin > 3 && ~isempty(varargin))
                if(mod(length(varargin),2))
                    error('You must supply a name and value for each constant');
                else
                    for i = 1:2:length(varargin)
                        B.constnt = [B.constnt; varargin(i:i+1)];
                    end
                end
            end
            %Indicate Rebuild Required
            Unbuild(B);
        end  
                
        function B = AddExpression(B,str,grp,ex)
        %Identify variable as an expression
        
            %Parse expression string
            [str,exp] = SymBuilder.parseExpression(str);
            %Concatenate with existing expressions
            B.exprsn = [B.exprsn; {str exp}]; 
            %If group passed, add to result group
            if(nargin > 3), B = AddResultExp(B,grp,exp,ex);
            elseif(nargin > 2), B = AddResultExp(B,grp,exp); end
            %Indicate Rebuild Required
            Unbuild(B);
        end
        
        function B = AddBound(B,str)
        %Add variable bound
        
            %Parse expression string
            [var,llb,lub] = SymBuilder.parseBound(str);
            %Concatenate with existing bounds
            B.bnds = [B.bnds; {var llb lub}];             
            %Indicate Rebuild Required
            Unbuild(B);        
        end
        
        function B = AddInteger(B,str)
        %Add integer constraint
        
            %Parse expression string
            [var,xx] = SymBuilder.parseInteger(str);
            %Concatenate with existing constants
            B.vartypes = [B.vartypes; {var xx}];             
            %Indicate Rebuild Required
            Unbuild(B);        
        end
        
        function B = AddResultGroup(B,name,str)
        %Add Result Group
            %Concatenate with existing groups
            B.resgrp = [B.resgrp; {name str}];        
        end
        
        function B = AddResultExp(B,name,str,bin)
        %Add Result Expression            
            %Parse expression string
            [group,name] = SymBuilder.parseResExp(name);
            %Concatenate with existing groups
            if(nargin > 3)
                B.resexp = [B.resexp; {group name str bin}];
            else
                B.resexp = [B.resexp; {group name str []}];        
            end
        end
                
        function B = Draft(B)
        %Draft the object by solving for the system Jacobian
            t = tic;
            %Build Symbolic representation of the equations
            fprintf('\nGenerating Symbolic Representation of Equations...');
            buildSymRep(B);
            fprintf('Done\n');
            %Generate Equations Jacobian
            fprintf('Generating Symbolic Jacobian...');
            buildJac(B);
            fprintf('Done\n');
            %Determine Linear & Nonlinear Equations
            fprintf('Generating Equation Linearity...');
            detLinearity(B);
            fprintf('Done\n');
            %Save build time
            B.bldtime = toc(t);
            %Assign build Status
            B.bldstat = 'draft';
        end
                
        function B = Build(B)
        %Build the object by solving for the system Hessian
            t = tic;
            %Build Symbolic representation of the equations
            fprintf('\nGenerating Symbolic Representation of Equations...');
            buildSymRep(B);
            fprintf('Done\n');
            %Generate Equations Jacobian
            fprintf('Generating Symbolic Jacobian...');
            buildJac(B);
            fprintf('Done\n');
            %Generate Equations Hessians (also determines linearity)
            fprintf('Generating Symbolic Hessian...');
            buildHess(B);
            fprintf('Done\n');
            %Save build time
            B.bldtime = toc(t);
            %Assign build Status
            B.bldstat = 'built';
        end
        
        function B = Unbuild(B)
        %Unbuild object (free memory)
        
            B.jac = [];
            B.hess = [];
            B.hesslag = [];
            B.bldstat = 'unbuilt';
            B.bldtime = 0;
            B.Opt = [];
            
            clear symb_grad symb_hess symb_nlcon symb_nljac symb_obj
        end
                
        function f = GetLinObj(B)
        %Get Linear Objective function [f'*x]
        
            %Check is built
            if(~IsBuilt(B))
                error('Please build the object using Draft() or Build()');
            end
            
            ind = B.objlin == 1;
            if(sum(ind) == 0)
                error('This object does not contain a linear objective');
            end
            f = double(B.jac(ind,:))';
        end
                
        function [obj,grad] = GetNLObj(B)
        %Get Nonlinear Objective Function [f(x)]
        
            %Check is built
            if(~IsBuilt(B))
                error('Please build the object using Draft() or Build()');
            end
        
            %Get Nonlinear Jacobian & Indices
            ind = B.objlin >= 2; 
            if(sum(ind) == 0)
                error('This object does not contain a nonlinear objective');
            end
            nljac = B.jac(ind,:);
            %Build NL Con Callback + row constraints
            B.buildMFun('symb_obj','obj',B.sobj(ind),B.vars);
            obj = @symb_obj;
            %Build NL Gradient Callback
            B.buildMFun('symb_grad','grad',nljac,B.vars);
            grad = @symb_grad;
        end
        
        function [obj,grad] = GetAllObj(B)
        %Get Objective Function [f(x)]
        
            %Check is built
            if(~IsBuilt(B))
                error('Please build the object using Draft() or Build()');
            end
        
            %Get All Objective Jacobian & Indices
            ind = B.objlin >= 1; 
            if(sum(ind) == 0)
                error('This object does not contain an objective');
            end
            ajac = B.jac(ind,:);
            %Build Con Callback + row constraints
            B.buildMFun('symb_obj','obj',B.sobj(ind),B.vars);
            obj = @symb_obj;

            %Optional Return Args
            if(nargout > 1)
                %Build Gradient Callback
                B.buildMFun('symb_grad','grad',ajac,B.vars);
                grad = @symb_grad;
            end
        end
                
        function [A,rl,ru] = GetLinCon(B)
        %Get Linear Constraint Function [rl <= A*x <= ru]
        
            %Check is built
            if(~IsBuilt(B))
                error('Please build the object using Draft() or Build()');
            end
            %Get Linear Indices
            ind = B.conlin == 1; 
            %Check we have some constraints
            if(B.noCons && sum(ind) > 0)
                %Get Linear Jacobian
                A = sparse(double(B.jac(ind,:)));
                rl = B.cl(ind);
                ru = B.cu(ind);
            else
                A = []; rl = []; ru = [];
            end
        end
        
        function [nlcon,cl,cu,nljac,nljacstr] = GetNLCon(B)   
        %Get Nonlinear Constraint Function [cl <= c(x) <= cu]
        
            %Check is built
            if(~IsBuilt(B))
                error('Please build the object using Draft() or Build()');
            end
            %Get Nonlinear Indices
            ind = B.conlin >= 2;
            %Check we have some constraints
            if(B.noCons && sum(ind) > 0)      
                %Get NL Jacobian
                nljac = B.jac(ind,:);
                %Build NL Con Callback + row constraints
                B.buildMFun('symb_nlcon','con',B.sobj(ind),B.vars);
                nlcon = @symb_nlcon;
                cl = B.cl(ind);
                cu = B.cu(ind);

                %Build NL Jacobian Structure Callback
                jacstr = zeros(size(nljac));
                jacstr(nljac ~= 0) = 1;
                nljacstr = @() sparse(jacstr);

                %Build NL Jacobian Callback
                B.buildMFun('symb_nljac','jac',nljac,B.vars);
                nljac = @symb_nljac;
            else
                nlcon = []; cl = []; cu = [];
                nljac = []; nljacstr = [];
            end            
        end
        
        function [nlcon,cl,cu,nljac,nljacstr] = GetAllCon(B,opts)
        %Get all constraints as nonlinear callback function
        
            %Check is built
            if(~IsBuilt(B))
                error('Please build the object using Draft() or Build()');
            end
            %Check we have some constraints
            if(B.noCons)
                %Get Full Jacobian & Indices
                ind = B.conlin > 0;
                ajac = B.jac(ind,:);

                %Check for options
                if(nargin < 2), opts = []; end

                %Build Con Callback + row constraints
                B.buildMFun('symb_nlcon','con',B.sobj(ind),B.vars,opts);
                nlcon = @symb_nlcon;
                cl = B.cl(ind);
                cu = B.cu(ind);

                %Optional return args
                if(nargout > 3)           
                    %Build Jacobian Structure Callback
                    jacstr = zeros(size(ajac));
                    jacstr(logical(ajac ~= 0)) = 1;
                    nljacstr = @() sparse(jacstr);

                    %Build Jacobian Callback
                    B.buildMFun('symb_nljac','jac',ajac,B.vars);
                    nljac = @symb_nljac;
                end
            else
                nlcon = []; cl = []; cu = [];
                nljac = []; nljacstr = [];
            end
        end
                
        function [H,f] = GetQuadObj(B)
        %Get Quadratic Objective Function [0.5*x'*H*x + f'*x]   
            
            %Check is built
            if(~IsBuilt(B))
                error('Please build the object using Draft() or Build()');
            end
        
            %Get quadratic objective indices
            ind = B.objlin == 2;
            if(isempty(ind))
                error('This object does not contain a quadratic objective');
            elseif(sum(ind) > 1)
                error('This interface only supports a single quadratic objective');
            end
            %Get Hessian
            H = 0.5*(sparse(double(B.hess(ind).H)));
            %Get Gradient, remove 2nd derivative parts
            f = B.jac(ind,:); v = symvar(f);
            f = double(subs(f,v,{zeros(size(v))}).');
        end
                
        function [lH,lHstr] = GetHessLag(B)
        %Get Hessian of the Lagrangian Function [@(x,sigma,lambda)]
        
            %Check is built
            if(~IsBuilt(B))
                error('Please build the object using Draft() or Build()');
            end
            %Check if we have done a full build
            if(isempty(B.hess))
                error('You can only obtain the Hessian of the Lagrangian if you called Build()\n\nDraft() only computes 1st Derivatives%s','.');
            end
            %Build HessLag
            buildHessLag(B);
            %Build NL Hessian Callback
            B.buildMFun('symb_hess','hess',B.hesslag,B.vars);
            lH = @symb_hess;
            
            lHstr = zeros(size(B.hesslag));
            lHstr(logical(B.hesslag ~= 0)) = 1;
            lHstr = @() sparse(lHstr);
        end       
        
        %Get Build Status
        function tf = IsBuilt(B)
            if(strcmp(B.bldstat,'unbuilt'))
                tf = false;
            else
                tf = true;
            end
        end
        
        function [x,fval,ef,info] = solve(B,x0,opts)
        %Lowercase version
            if(nargin < 3), opts = []; end
            if(nargin < 2), x0 = []; end
            [x,fval,ef,info] = Solve(B,x0,opts);
        end
        
        function [x,fval,ef,info] = Solve(B,x0,opts)
        %Solve Object using OPTI

            %Check is built
            if(~IsBuilt(B))
                error('Please build the object first using Draft() or Build()');
            end
            if(nargin < 3 || isempty(opts)), opts = symbset; end
            if(nargin < 2), x0 = []; end
            
            fprintf('\nGenerating OPTI Object....\n');
            %Cusomtize settings based on solver
            switch(lower(opts.solver))
                %White box solvers
                case {'scip','baron'}
                     opts = symbset(opts,'use1stDerivs','no','use2ndDerivs','no','preallocate','no');
                %Derivative free solvers
                case {'nomad','pswarm','gmatlab'}                    
                    opts = symbset(opts,'use1stDerivs','no','use2ndDerivs','no');
                %No Hessian Support                   
                case {'filtersd','lbfgsb','nlopt'}
                    opts = symbset(opts,'use2ndDerivs','no');
            end
            B.Opt = GetOPTI(B,opts);
            fprintf('Done\n\n');
            %Solve
            [x,fval,ef,info] = solve(B.Opt,x0);
        end
    end
    

    methods(Access=private)
                    
        function buildSymRep(B)
        %Build symbolic representation of the equations
        
            %Check if something new has been added
            if(strcmp(B.bldstat,'unbuilt'))
                %Build Symbolic Object
                digits(16);
                symobj = sym(sprintf('[%s]',B.eqs)); 
%                 %Substitute constants (1st round)
%                 if(~isempty(B.constnt))
%                     symobj = subs(symobj,B.constnt(:,1),B.constnt(:,2));
%                 end
                %Substitute expressions
                if(~isempty(B.exprsn))                                           
                    %Now subs into full equation system
                    symobj = subs(symobj,B.exprsn(:,1),B.exprsn(:,2));
                    %Have to repeat until all nested expressions are sub'd
                    n = 10;  %max depth
                    no = size(B.exprsn,1);
                    se = sym(B.exprsn(:,1));
                    while n > 0
                        v = symvar(symobj); alldone = 1;
                        for i = 1:no
                            if(any(se(i) == v))
                                symobj = subs(symobj,B.exprsn(i,1),B.exprsn(i,2));
                                alldone = 0;
                            else
                                if(i == no && alldone) %ensure we have checked them all
                                    n = -10;
                                    break;
                                end
                            end
                        end  
                        n = n - 1;
                    end
                    if(n == 0)
                        error('Maximum expression recursion depth reached!');
                    else
%                         fprintf('%d exp recursions required\n',n+10);
                    end
                end
                %Convert constants to decimal form
%                 digits(16); %set vpa digits
%                 BV = cell(size(B.constnt,1),1);
%                 for i = 1:size(B.constnt,1)
%                     BV{i} = sym(B.constnt{i,2},'d');
%                 end
                %Substitute constants
                if(~isempty(B.constnt))
                    symobj = subs(symobj,B.constnt(:,1),B.constnt(:,2));
                end
                %Save symbolic vector
                B.sobj = symobj;
                %Save Variables
                B.vars = symvar(B.sobj);
                
                %Build default bounds and xtype
                n = length(B.vars);
                B.lb = -1*Inf(n,1);
                B.ub = Inf(n,1);
                B.xtype = repmat('C',1,n);   
                %Get names of variables and their indices
                names = SymBuilder.detVarNames(B.vars);
                if(isempty(names) && size(B.bnds,1) > 0)
                    error('Cannot process bounds without declared variables!');
                end
                                
                %Process bound declarations
                for i = 1:size(B.bnds,1)
                    %Check if bounds are numerical
                    llb = str2double(B.bnds{i,2});
                    if(isnan(llb)) %not a number, check constants
                        ind = strcmp(B.constnt(:,1),B.bnds{i,2});
                        if(any(ind))
                            llb = B.constnt{ind,2};
                        else
                            error('Unknown lower bound: %s',B.bnds{i,2});
                        end
                    end
                    lub = str2double(B.bnds{i,3});
                    if(isnan(lub)) %not a number, check constants
                        ind = strcmp(B.constnt(:,1),B.bnds{i,3});
                        if(any(ind))
                            lub = B.constnt{ind,2};
                        else
                            error('Unknown upper bound: %s',B.bnds{i,3});
                        end
                    end
                    %Check if variable exists
                    ind = logical(B.vars == sym(B.bnds{i,1}));
                    if(any(ind))
                        if(~isinf(llb)), B.lb(ind) = llb; end
                        if(~isinf(lub)), B.ub(ind) = lub; end                    
                    else %could be applied to multiple variables
                        ind = strcmp(names(:,1),B.bnds{i,1});
                        if(isempty(ind) || all(ind  == 0))
                            error('Unknown bound: %s',B.bnds{i,1});
                        elseif(sum(ind) > 1)
                            error('Bound name is ambiguous and matches more than one variable name');
                        else                            
                            %Extract indices
                            start = names{ind,2}(1);
                            send = names{ind,2}(end);
                            if(~isinf(llb)), B.lb(start:send) = llb; end
                            if(~isinf(lub)), B.ub(start:send) = lub; end 
                        end
                    end
                end
                
                %Process integer declarations
                for i = 1:size(B.vartypes,1)
                    %Check if variable exists
                    ind = logical(B.vars == sym(B.vartypes{i,1}));
                    if(any(ind))
                         B.xtype(ind) = B.vartypes{i,2};                
                    else %could be applied to multiple variables
                        ind = strcmp(names(:,1),B.vartypes{i,1});
                        if(isempty(ind))
                            error('Unknown variable: %s',B.vartypes{i,1});
                        elseif(sum(ind) > 1)
                            error('Variable name is ambiguous and matches more than one variable name');
                        else                            
                            %Extract indices
                            start = names{ind,2}(1);
                            send = names{ind,2}(end);
                            B.xtype(start:send) = B.vartypes{i,2};
                        end
                    end
                end
            end
        end
                
        function buildJac(B)
        %Evaluate Jacobian of entire equation system
        
            if(isempty(B.sobj))
                error('You must build a symbolic representation of the system first!');
            end
            %Call SymToolbox jacobian function
            B.jac = jacobian(B.sobj);
        end  
                
        function buildHess(B)
        %Evaluate Hessian(s) of entire equation system
        
            if(isempty(B.jac))
                error('You must evaluate the Jacobian of the system first!');
            end
            noeq = B.noObjs + B.noCons;
            B.hess = struct('H',cell(noeq,1),'lin',cell(noeq,1));
            %Build Index Vectors
            B.conlin = zeros(noeq,1);
            B.objlin = zeros(noeq,1);
            %Go through each equation
            for i = 1:noeq
                J = B.jac(i,:);
                %Determine if linear
                if(isempty(symvar(J)))
                    B.hess(i).lin = 1;
                    B.hess(i).H = 0;
                else
                    %Hessian is the Jacobian of the Jacobian (in this case)
                    H = jacobian(J,B.vars);
                    s = symvar(H);
                    %Determine if quadratic
                    if(isempty(s))
                        B.hess(i).lin = 2;
                        B.hess(i).H = H;
                    else
                        %Otherwise must be nonlinear, save entire symobj (for now)
                        B.hess(i).lin = 3;
                        B.hess(i).H = H;
                    end
                end
                if(any(B.indobj == i))
                    %If the objective is linear or quadratic, it must not contain a constant (strictly), so check here
                    if(B.hess(i).lin <= 2)
                        eq = B.sobj(i);
                        v = symvar(eq);
                        if(double(subs(eq,v,{zeros(size(v))})) ~= 0) %not a LP/QP
                            B.hess(i).lin = 3;
                        end
                    end                        
                    B.objlin(i) = B.hess(i).lin;
                else
                    B.conlin(i) = B.hess(i).lin;
                end
            end
        end
                
        function buildHessLag(B)
        %Build Lagrangian Hessian assuming IPOPT (tril) format 
        
            if(isempty(B.hess))
                error('You must evaluate the Hessian of the system first!');
            end
            %Start with objective
            nvars = length(B.vars);
            ind = B.objlin >= 2;
            if(any(ind)) %check we have a nonlinear (or quadratic) objective  
                if(sum(ind) > 1)
                    error('This interface only support single objective problems');
                end
                s = sym('sigma');
                objH = B.hess(ind).H;
                H = s*objH;
            else
                H = zeros(nvars);
            end
            %Now each nonlinear (or quadratic) constraint
            ind = B.conlin > 1;
            %Two sets of indexes, one for lambda, one for B.hess
            index_BH = find(ind);
            ind(B.indobj) = [];
            index_L = find(ind);
            for i = 1:length(index_BH)
                l = sym(sprintf('lambda(%d)',index_L(i)));
                H = H + l*B.hess(index_BH(i)).H;
            end  
                
            %Save resulting Hessian
            B.hesslag = tril(H);
        end
                
        function detLinearity(B)
        %Determine linearity of the supplied equations (without Hessian)
        
            if(isempty(B.jac))
                error('You must evaluate the Jacobian of the system first!');
            end
            %Build Index Vectors
            noeq = B.noCons + B.noObjs;
            B.conlin = zeros(noeq,1);
            B.objlin = zeros(noeq,1);
            %Determine Linearity (based on existence of symvars)
            for i = 1:noeq %must search all equations
                %Gather first derivative vars
                s = symvar(B.jac(i,:));
                %Objective
                if(any(B.indobj == i))
                    if(isempty(s))
                        B.objlin(i) = 1; %lin
                    else
                        B.objlin(i) = 3; %nonlin
                    end
                else                   
                    if(isempty(s))
                        B.conlin(i) = 1;
                    else
                        B.conlin(i) = 3;
                    end
                end
            end            
        end 
    end
    
    
    methods(Static)        
        %Parse Constraint String
        [str,l,u] = parseConstraint(str);       
        %Extract RHS of expression string, remove from eq, and return LHS
        [str,exp] = parseExpression(str);
        %Extract Group and Name from Result Expression Name
        [group,name] = parseResExp(name);
        %Parse bound string
        [var,lb,ub] = parseBound(str);
        %Parse integer string
        [var,xtype] = parseInteger(str);
        %Determine Variable names
        varn = detVarNames(svar);
        %Build MATLAB function file
        buildMFun(name,mode,sobj,svar,opts);
        %Convert symbolic expression into matlab function
        fun = sym2fun(sobj,svar);        
    end
    
end

