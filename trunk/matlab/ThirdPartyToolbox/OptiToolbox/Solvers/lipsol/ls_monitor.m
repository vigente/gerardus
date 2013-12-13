function ls_monitor(iter,rb,rc,ru,dgap)
% MONITOR     - Graphic monitor for iteration progress.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County
% Modified J.Currie AUT May 2013

global rb_prev rc_prev ru_prev dg_prev

if iter == 0
   clf;
   drawnow;
   ymax = min(16,ceil(1+log10(max([rb,rc,ru,dgap]))));
   axis([0 50 -16 ymax]);
   title('LIPSOL Monitor'); 
   xlabel('Iteration'); 
   ylabel('Log10(Residuals)');
   hold on;
   rb   = max(1.e-32,rb);   rb_prev = rb;
   rc   = max(1.e-32,rc);   rc_prev = rc;
   ru   = max(1.e-32,ru);   ru_prev = ru;
   dgap = max(1.e-32,dgap); dg_prev = dgap;
else
    rb   = max(1.e-32,rb);
    rc   = max(1.e-32,rc);
    ru   = max(1.e-32,ru);
    dgap = max(1.e-32,dgap);
    plot([iter-1 iter], log10([rb_prev rb]),   ':', ...
        [iter-1 iter], log10([rc_prev rc]),   '-.', ...
        [iter-1 iter], log10([ru_prev ru]),   '--', ...
        [iter-1 iter], log10([dg_prev dgap]), '-' );
    if(iter==1)
        legend('Duality Gap','P-Infeasibility','D-Infeasibility','Upper Bounds');
    end
    drawnow;
    rb_prev = rb;
    rc_prev = rc;
    ru_prev = ru;
    dg_prev = dgap;
end;
