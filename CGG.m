
% conjugate gradient method
% reference: Conjugate guided gradient (CGG) method for robust inversion and its application to velocity-stack inversion.  junji 2006
% more details can be found in github.com/jbrussell

function [m,r,niter,fiter] = CGG(m0,LL,d,maxiter,rthresh)

% -------------------------------------------------------------------------------------------------------------
% output£º
% m             Radon solution in frequency domain
% r             residual vector
% niter         number of iterations
% fiter         the final state of the iteration: 1, the maximum number of iterations is reached 0, the residual threshold is met -1, abnormal

% input£º
% m0            Radon solution in frequency domain(initial)
% LL            operator matrix L
% d             left half Fourier transform of seismic gather matrix
% maxiter       the maximum number of iterations
% rthresh       residual threshold
% -------------------------------------------------------------------------------------------------------------

m = zeros(size(m0));
r = d;
Dm = [];
ds = [];

epsilon = max(abs(d))/100;                               % junji 2006 p5

niter = 0;
while niter < maxiter && norm(r)/norm(d) >= rthresh
    
    % Weighting Matrices
    if niter == 0
        Wm = eye(length(m));
        Wr = eye(length(r));
    else
        % junji 2006 Formula 9
        Wm = diag(abs(m).^(1/2));           % L1 norm
        
        % junji 2006 Formula 7
        Wr = diag(abs(r).^(-1/2));          % L1 norm
        
        Wr(diag(abs(r))<=epsilon & diag(abs(r))~=0) = epsilon;
    end
    
    dm = Wm'*LL'*Wr'*r;
    dr = LL*dm;
    [m,r,Dm,ds] = cgstep(niter,m,dm,Dm,r,dr,ds);
    niter = niter + 1;
end

if niter == maxiter
    fiter = 1;
elseif niter < maxiter
    fiter = 0;
elseif niter > maxiter
    fiter = -1;
end

end
