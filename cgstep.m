
% reference: Clarbout 2004 p139-p143
% more details can be found in github.com/jbrussell

function [m,r,Dm,ds] = cgstep(niter,m,dm,Dm,r,dr,ds)

% -------------------------------------------------------------------------------------------------------------
% output£º
% m  [x] - model solution
% r  [R] - residual
% Dm [s] - previous descent step
% ds [S] - descent vector

% input£º
% niter
% m  [x] - model solution
% dm [g] - gradient vector
% Dm [s] - previous descent step
% r  [R] - residual
% dr [G] - conjugate gradient vector
% ds [S] - descent vector
% -------------------------------------------------------------------------------------------------------------

if niter == 0
    Dm = zeros(size(m));
    ds = zeros(size(r));
    if r'*r == 0
        error('r = 0');
    end
    alpha = (dr'*r)/(dr'*dr);
    beta = 0;
else
    dr_dr = dr'*dr;
    ds_ds = ds'*ds;
    dr_ds = dr'*ds;
    determ = dr_dr*ds_ds - dr_ds*dr_ds + (0.00001*(dr_dr*ds_ds)+1e-15);
    dr_r = dr'*r;
    ds_r = ds'*r;
    alpha = (ds_ds*dr_r - dr_ds*ds_r)/determ;
    beta = (-dr_ds*dr_r + dr_dr*ds_r)/determ;
end

Dm = alpha*dm + beta*Dm;

ds = alpha*dr + beta*ds;

m = m + Dm;

r = r - ds;

end
