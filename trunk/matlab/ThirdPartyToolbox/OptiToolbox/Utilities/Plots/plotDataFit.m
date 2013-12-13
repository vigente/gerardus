function plotDataFit(prob,xb,doLog)
%PLOTDATAFIT Add new Plot with Data + Fit

%   Copyright (C) 2011 Jonathan Currie (I2C2)
    
if(nargin(prob.fun) == 2)
    if(~isempty(prob.xdata))
        %If xdata is a matrix, then just plot 1:no points
        if(size(prob.xdata,1) > 1 && size(prob.xdata,2) > 1)
            if(size(prob.xdata,1) ~= length(prob.ydata))
                optiwarn('opti:dfit','Cannot plot data fit as size(xdata,1) and length(ydata) are not the same length!');
            end    
            t = 1:length(prob.ydata);
            y = prob.fun(xb,prob.xdata);
            plot(t,prob.ydata,'o',t,y);            
        else
            if(length(prob.xdata) ~= length(prob.ydata))
                optiwarn('opti:dfit','Cannot plot data fit as xdata and ydata are not the same length!');
            end    
            %Generate Smooth Data
            x = linspace(prob.xdata(1),prob.xdata(end),1e3);
            y = prob.fun(xb,x);
            plot(prob.xdata,prob.ydata,'o',x,y);
        end
        title(['NLS Curve Fit - SSE: ' num2str(sum((prob.fun(xb,prob.xdata)-prob.ydata).^2))]);
        xlabel('x'); ylabel('y');
        legend('Original Data','NLS Fit');
    end
else
    t = 1:length(prob.ydata);
    plot(t,prob.ydata,'o',t,prob.fun(xb)); ylabel('y');
    title(['NLS Curve Fit - SSE: ' num2str(sum((prob.fun(xb)-prob.ydata).^2))]);
    xlabel('x'); 
    legend('Original Data','NLS Fit');
end




