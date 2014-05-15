function res = Results(B,names)
%RESULTS  Pretty Print Results of a Solved SymBuilder Object

if(isempty(B.Opt) || isempty(B.Opt.sol))
    error('You cannot print results from an object which has not been solved!');
end

fprintf('\n');
disp('------------------------------------------------------');
fprintf('SymBuilder Optimization Results\n');

%Overall solve status
switch(B.Opt.ef)
    case 1
        fprintf(' SOLVED in %1.5gs\n',B.Opt.info.Time);
    case 0
        fprintf(' SOLVE LIMIT REACHED in %1.5gs\n',B.Opt.info.Time);
    otherwise
        fprintf(' NOT SOLVED! STATUS: %s\n\n',B.Opt.info.Status);
        disp('------------------------------------------------------');
        return;
end
%Print Status
fprintf(' STATUS: %s\n',B.Opt.info.Status);
%Print Iterations
if(isfield(B.Opt.info,'Iterations'))
    fprintf(' ITERATIONS: %d\n',B.Opt.info.Iterations);
end
if(isfield(B.Opt.info,'Nodes'))
    fprintf(' NODES: %d\n',B.Opt.info.Nodes);
end
if(isfield(B.Opt.info,'AbsGap'))
    fprintf(' ABS GAP: %1.6g\n',B.Opt.info.AbsGap);
end
%Print Cost
fprintf('\n COST: %1.6g\n',B.Opt.obj);
disp('------------------------------------------------------');


%For each variable in the model, assign its value from the solution
for i = 1:length(B.vars)
    eval(sprintf('%s = %1.15g;',char(B.vars(i)),B.Opt.sol(i)));
end
%Eval constants
for i = 1:length(B.constnt)
    eval(sprintf('%s = %1.15g;',B.constnt{i,1},B.constnt{i,2}));
end
%Eval expressions until no errors (missing vars occur) (poor form Jonny)
done = false;
while ~done
    done = true;
    %Eval constants
    for i = 1:length(B.exprsn)
        try
            eval(sprintf('%s = %1.15g;',B.exprsn{i,1},eval(B.exprsn{i,2})));
        catch ME
            done = false;
        end
    end        
end

%Iterate through each solution group, printing results as we go
if(isempty(B.resgrp))
    fprintf('NO result groups exist, please add result groups and variables in order to display information here.\n\n');
elseif(isempty(B.resexp))
    fprintf('NO result expressions exist, please add result groups and variables in order to display information here.\n\n');
else
    %For each group
    for i = 1:size(B.resgrp,1)
        fprintf('%s: %s\n',B.resgrp{i,1},B.resgrp{i,2}); %display title
        %Find expressions / resexp to display
        group = B.resgrp{i,1};
        ind = ismember(B.resexp(:,1),group);
        %For each match, eval and print
        if(any(ind))
            ind = find(ind);
            [~,sind] = sort(B.resexp(ind,2));
            ind = ind(sind);
            for j = 1:length(ind)
                k = ind(j);
                if(~isempty(B.resexp{k,4}))
                    ex = eval(B.resexp{k,4});
                    if(ex == 1 || ex == 0)
                        if(ex == 1), str = 'On'; else str = 'Off'; end
                        fprintf('  - %-12s = %-13s  [%s]\n',B.resexp{k,2},getNum(eval(B.resexp{k,3})),str);
                    else
                        fprintf('  - %-12s = %-13s  %-13s\n',B.resexp{k,2},getNum(eval(B.resexp{k,3})),getNum(ex));
                    end
                else
                    fprintf('  - %-12s = %-13s\n',B.resexp{k,2},getNum(eval(B.resexp{k,3})));
                end
            end
        end
        fprintf('\n');
    end    
end
disp('------------------------------------------------------');

%If names provided for each variable display
if(nargin > 1)
    res = cell2struct(num2cell(B.Opt.sol),names);
else
    res = B.Opt.info;
end

function str = getNum(num)
if(abs(num) < 1e-8)
    str = '0';
else
    str = sprintf('%-13.5g',num);
end