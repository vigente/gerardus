function opts = optiRemoveDefaults(opts,defs)
%OPTIREMOVEDEFAULTS  Remove Default Options from an Options Structure
oFn = fieldnames(opts);
for i = 1:length(oFn)
    label = oFn{i};
    if(isfield(defs,label))
        if(ischar(opts.(label)))
            if(strcmpi(defs.(label),opts.(label)))
                opts = rmfield(opts,label);
            end
        else
            if(defs.(label) == opts.(label))
                opts = rmfield(opts,label);
            end
        end
    end
end