
% forward Radon transform

function [ Rfft,f,iF,rfft,allniter,allfiter] = Radon_forward(p,t,M,delta,maxiter,rthresh,method)

% -------------------------------------------------------------------------------------------------------------
% Algorithms from Ji 2006 & Claerbout 1992
% org ed: github.com/jbrussell
% -------------------------------------------------------------------------------------------------------------

% Define some array/matrices lengths.
it = length(t);
iF = pow2(nextpow2(it)+1);
iDelta = length(delta);
ip = length(p);
Mfft = fft(M,iF,2);
dF = 1/(t(1)-t(2));
f  = (((1:floor((iF+1)/2))-1)/iF)*dF*-1;

% Define blocks
delta_block = repmat(delta,ip,1)';
p_block = repmat(p,iDelta,1);
m0 = zeros(length(p),1);
Rfft = zeros(ip,iF);
rfft = zeros(size(Mfft));

parfor j = 1:length(f)
    
    exp_arg = -1i*2*pi*f(j).*delta_block.*p_block;
    L = exp(exp_arg);
    d = Mfft(:,j);
    
    switch lower(method)
        case 'cgg'
            [m,r,niter,fiter] = CGG(m0,L,d,maxiter,rthresh);
        otherwise
            error('Incorrect method name.');
    end
    
    Rfft(:,j) = m;
    rfft(:,j) = r;
    allniter(:,j) = niter;
    allfiter(:,j) = fiter;
    
    if (mod(j,100)==0)
        fprintf("%d/%d freq has been finished...\n",j,length(f));
    end
    
end

