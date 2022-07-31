
% inverse Radon transform

function M = Radon_inverse(t,p,R,delta,ref_dist,line_model)

% -------------------------------------------------------------------------------------------------------------
% Algorithms from Ji 2006 & Claerbout 1992
% org ed: github.com/jbrussell
% -------------------------------------------------------------------------------------------------------------

% Define some array/matrices lengths.
it = length(t);
iF = pow2(nextpow2(it)+1);
iDelta = length(delta);
ip = length(p);

% Exit if inconsistent data is input.
if(min([ip,it] ~= size(R)))
    fprintf('Dimensions inconsistent!\nsize(R)~=[length(p),length(t)]\n');
    M = 0;
    return;
end

% Preallocate space in memory.
Mfft = zeros(iDelta, iF);
A = zeros(iDelta,ip);
Tshift = A;

% Define some values.
Dist_array = delta - ref_dist;
dF = 1/(t(1)-t(2));
Rfft = fft(R,iF,2);

% Populate ray parameter then distance data in time shift matrix.
for j = 1:iDelta
    Tshift(j,:) = p;
end

for k = 1:ip
    Tshift(:,k) = Tshift(:,k).*Dist_array';
end

% Loop through each frequency.
for i=1:floor((iF+1)/2)
    
    % Make time-shift matrix, A.
    f = ((i-1)/iF)*dF;
    A = exp(  (2i*pi*f).*Tshift  );
    
    % Apply Radon operator.
    Mfft(:,i) = A*Rfft(:,i);
    
    % Assuming Hermitian symmetry of the fft make negative frequencies the complex conjugate of current solution.
    if(i ~= 1)
        Mfft(:,iF-i+2) = conj(Mfft(:,i));
    end
end

M = ifft(Mfft,iF,2, 'symmetric');
M = M(:,1:it);
return;
