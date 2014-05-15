function buildMFun(name,mode,sobj,svar,opts)
% BUILDMFUN Build a MATLAB function for nonlinear callbacks

%Determine callback type
switch(mode)
    case 'obj'
        title = 'Objective Function Callback';
        fcall = sprintf('function j = %s(x)\n',name);
        var = 'j';
    case 'grad'
        title = 'Objective Gradient Callback';
        fcall = sprintf('function g = %s(x)\n',name);
        var = 'g';
    case 'con'
        title = 'Constraint Function Callback';
        fcall = sprintf('function c = %s(x)\n',name);
        var = 'c';
    case 'jac'
        title = 'Constraint Jacobian Callback';
        fcall = sprintf('function J = %s(x)\n',name);
        var = 'J';
    case 'hess'
        title = 'Hessian of the Lagrangian Callback';
        fcall = sprintf('function H = %s(x,sigma,lambda)\n',name);
        var = 'H';
    otherwise
        error('Unknown MFun mode: %s',mode);
end

%Open file
fp = fopen([name '.m'],'w');
%Write Header            
fprintf(fp,fcall);
fprintf(fp,'%% %s\n',upper(name));
fprintf(fp,'%%\n%% %s\n\n',title);
fprintf(fp,'%% Symbolic Builder Auto-Generated Callback Function\n');
fprintf(fp,'%% Generated %s\n\n',datestr(now));

%Distribute Vars            
v = SymBuilder.detVarNames(svar);
if(~strcmp(v{1,1},'x') || size(v,1) > 1)
    fprintf(fp,'%% Slice Input Array\n');
    for i = 1:size(v,1)
        if(v{i,2}(end)-v{i,2}(1) == 0)
            fprintf(fp,'%s = x(%d);\n',v{i,1},v{i,2}(1));
        else
            fprintf(fp,'%s = x(%d:%d);\n',v{i,1},v{i,2}(1),v{i,2}(end));
        end
    end
end

%Check for options
if(nargin < 5 || ~isstruct(opts))
    preallocate = true;
else
    if(isfield(opts,'preallocate'))
        preallocate = opts.preallocate;
    end
end

%Build cell array of indexed variable strings
ivar = buildVarIndex(v,length(svar));      
%Ensure both columns
if(size(svar,1) > 1), svar = svar.'; end
if(size(ivar,1) > 1), ivar = ivar'; end
%Subs out individual symbolic variables into our indexed list
eq = subs(sobj,svar,ivar);

%Get equation size (matrices treated differently)
[r,c] = size(eq);

%Enter Equations (Var Type Dictates Entry Type)
switch(var)
    case 'j' %scalar
        fprintf(fp,'\n%% Equation:\n');
        fprintf(fp,'%s(1,1) = %s;\n',var,char(eq));
    case {'g','c'} %vector, dense    
        %HACK
        if(~isa(eq,'sym')), eq = sym(eq); end
        %Check for row
        if(c > 1), eq = eq.'; tr = 1; else tr = 0; end
        %Convert equations to string
        str = char(eq);
        %Remove 'matrix' if is found
        ind = strfind(str,'matrix');
        if(ind)
            str = str(ind+7:end-1);
        end
        if(preallocate)
            fprintf(fp,'\n%% Preallocate\n');
            if(tr)
                fprintf(fp,'%s = zeros(1,%d);\n',var,length(eq));
            else
                fprintf(fp,'%s = zeros(%d,1);\n',var,length(eq));
            end
        end
        if(strcmp(str(1),'[')), str = str(2:end); end %remove extra square brackets
        if(strcmp(str(end-1:end),']]')), str = str(1:end-2); end
        ss = regexp(str,'], ','split');
        fprintf(fp,'\n%% Equations:\n');
        if(tr)
            for i = 1:length(ss)
                if(~strcmp(ss{i}(2:end),'0'))
                    if(strcmp(ss{i}(1),'[')), ss{i} = ss{i}(2:end); end %remove extra square brackets
                    fprintf(fp,'%s(1,%d) = %s;\n',var,i,ss{i});
                end
            end
        else
            for i = 1:length(ss)
                if(~strcmp(ss{i}(2:end),'0'))
                    if(strcmp(ss{i}(1),'[')), ss{i} = ss{i}(2:end); end %remove extra square brackets
                    fprintf(fp,'%s(%d,1) = %s;\n',var,i,ss{i});
                end
            end   
        end
    case {'J','H'} %matrix, sparse
        %Get NonZero elements from equation
        nzel = logical(eq ~= 0);               
        [rows,cols] = find(sparse(nzel));
        nz = nnz(nzel);                
        nzsym = eq(nzel);

        %OLD WAY (Inefficient?)
        if(preallocate)
            fprintf(fp,'\n%% Preallocate\n');
            fprintf(fp,'%s = spalloc(%d,%d,%d);\n',var,r,c,nz);
        end
        fprintf(fp,'\n%% Sparse Matrix:\n');
        for i = 1:nz
            fprintf(fp,'%s(%d,%d) = %s;\n',var,rows(i),cols(i),char(nzsym(i)));
        end

        %FASTER BUT LESS CLEAR
%         %Row & Col Indices
%         fprintf(fp,'\n%% Row & Col Indices\n');
%         fprintf(fp,'rows = [');
%         for i = 1:length(rows)
%             if(i == length(rows))
%                 fprintf(fp,'%d];\n',rows(i));
%             else
%                 fprintf(fp,'%d, ',rows(i));
%             end
%         end
%         fprintf(fp,'cols = [');
%         for i = 1:length(cols)
%             if(i == length(cols))
%                 fprintf(fp,'%d];\n',cols(i));
%             else
%                 fprintf(fp,'%d, ',cols(i));
%             end
%         end
%         fprintf(fp,'\n%% Preallocate\n');
%         fprintf(fp,'vals = zeros(%d,1);\n',nz);
%         
%         fprintf(fp,'\n%% Nonzero Entries\n');
%         for i = 1:nz
%             fprintf(fp,'vals(%d,1) = %s;\n',i,char(nzsym(i)));
%         end
%         
%         fprintf(fp,'\n%% Create Sparse Return Matrix\n');
%         fprintf(fp,'%s = sparse(rows,cols,vals);\n',var);
end

%Close file
fclose(fp);


%Convert variables from names to indexed strings [nominally m1 to m(1)]
function ivar = buildVarIndex(v,no)   

%Create cell to store strings in
ivar = cell(no,1); n = 1;

%Process v cell matrix, building strings as we go
for i = 1:size(v,1)
    name = v{i,1};
    indices = v{i,2};
    %If single variable, easy!
    if(length(indices) == 1)
        ivar{n} = name; n = n + 1;
    %Otherwise have to run through indexing list
    else
        for j = 1:length(indices)
            ivar{n} = sprintf('%s(%d)',name,j); n = n + 1;
        end
    end
end      