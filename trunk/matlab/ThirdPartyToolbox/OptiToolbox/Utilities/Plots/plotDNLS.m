function plotDNLS(prob,opts,xb)
%plotDNLS Plot Parameter Estimation Problem

%   Copyright (C) 2013 Jonathan Currie (I2C2)

%Measurement Plot Color
measC = [0.4 0.4 0.4];
%Measurement Plot Style
measS = 'o';
%Initial Condition plot style
icS = 'sq';

% Insert initial conditions if also solved
ind = isnan(prob.odez0);
if(any(ind))
    len = sum(ind);
    estZ0 = xb(end-len+1:end);
    prob.odez0(ind) = estZ0;
    xb = xb(1:end-len);
end

% Generate smooth plot 
dopts = optidynset(opts.dynamicOpts);
tspan = [prob.xdata(1) prob.xdata(end)];
ode = @(t,z) prob.ode(t,z,xb);
switch(dopts.integrator)
    case 'ode45'
        [t,z] = ode45(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode15s'
        [t,z] = ode15s(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode23'
        [t,z] = ode23(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode113'
        [t,z] = ode113(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode23t'
        [t,z] = ode23t(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode23tb'
        [t,z] = ode23tb(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode23s'
        [t,z] = ode23s(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode15i'
        [t,z] = ode15i(ode,tspan,prob.odez0,dopts.odeOpts);
    otherwise
        [t,z] = ode45(ode,tspan,prob.odez0,dopts.odeOpts);
end

%Plot smooth sim
plot(t,z);
title('Dynamic Parameter Estimation Solution');
xlabel('time'); ylabel('z');
hold on;

%Plot Solved Initial Conditions
if(any(ind))
    for i = 1:sum(ind)
        plot(tspan(1),estZ0(i),icS,'Color',measC);
    end
end

%If we have xdata_old, then should be in cell format, plot each cell
if(isfield(prob,'xdata_old') && ~isempty(prob.xdata_old))
    if(iscell(prob.xdata_old) && iscell(prob.ydata_old))
        for i = 1:length(prob.xdata_old)
            plot(prob.xdata_old{i},prob.ydata_old{i},measS,'Color',measC);
        end 
    elseif(~iscell(prob.xdata_old) && ~iscell(prob.ydata_old))
        for i = 1:size(prob.ydata_old,2)
            plot(prob.xdata_old,prob.ydata_old(:,i),measS,'Color',measC);  
        end
    else
        error('Expected xdata_old and ydata_old to be cell arrays');
    end    
else
    %Otherwise reshape ydata (always a column in OPTI) and plot
    if(~isempty(dopts.stateIndex))
        nstates = length(prob.odez0(dopts.stateIndex));
    else
        nstates = length(prob.odez0);
    end
    ydata = reshape(prob.ydata,length(prob.xdata),nstates);
    %Plot
    plot(prob.xdata,ydata,measS,'Color',measC);
end

hold off;