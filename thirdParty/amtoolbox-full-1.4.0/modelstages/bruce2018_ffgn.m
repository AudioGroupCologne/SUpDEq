function y = bruce2018_ffgn(N, tdres, Hinput, noiseType, mu, sigma)
%BRUCE2018_FFGN  Fast (exact) fractional Gaussian noise and Brownian motion generator
%   Usage: y = bruce2018_ffgn(N, tdres, Hinput, noiseType, mu, sigma)
%     
%   Input parameters:
%     N             : is the length of the output sequence
%     tdres         : is the time resolution (1/sampling rate)
%     Hinput        : is the Hurst index of the resultant noise (0 < H <= 2). For 0 < H <= 1,the output will be fractional Gaussian noise with Hurst index H.  For 
%                     1 < H <= 2, the output will be fractional Brownian motion with Hurst
%                     index H-1.  Either way, the power spectral density of the output will
%                     be nominally proportional to 1/f^(2H-1)
%     noiseType     : is 0 for fixed fGn noise and 1 for variable fGn 
%     mu            : is the mean of the noise. [default = 0]
%     sigma         : is the standard deviation of the noise [default = 1]
%
%   Output parameters:
%     y         : a sequence of fractional Gaussian noise with a mean of zero and a 
%                 standard deviation of one or fractional Brownian motion derived from such
%                 fractional Gaussian noise.
%
%   BRUCE2018_FFGN returns a vector containing a sequence of fractional Gaussian 
%   noise or fractional Brownian motion.  The generation process uses an FFT which makes it 
%   very fast.
%
%   References:
%     R. Davies and D. Harte. Tests for hurst effect. Biometrika, 74(1):95 --
%     101, 1987.
%     
%     J. Beran. Statistics for long-memory processes, volume 61. CRC Press,
%     1994.
%     
%     J. Bardet. Statistical study of the wavelet analysis of fractional
%     brownian motion. Information Theory, IEEE Transactions on,
%     48(4):991--999, 2002.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/bruce2018_ffgn.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB MEX M-Signal
%   #Author: B. Scott Jackson (2005)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% Check to see if running under Matlab or Octave
if exist ('OCTAVE_VERSION', 'builtin') ~= 0
    pkg load signal;
end

%---- Check input arguments ---------- %

if ( (nargin < 5) || (nargin > 6) )
	error('Requires Five to Six input arguments.')
end

if (prod(size(N)) ~= 1) || (prod(size(Hinput)) ~= 1) || ~isnumeric(N) || ~isnumeric(Hinput) ...
        || ~isreal(N) || ~isreal(Hinput) || ~isfinite(N) || ~isfinite(Hinput)
	error('All input arguments must be finite real scalars.')
end

if (N <= 0)
	error('Length of the return vector must be positive.')
end

if (tdres > 1)
	error('Original sampling rate should be checked.')
end 

if (Hinput < 0) || (Hinput > 2)
	error('The Hurst parameter must be in the interval (0,2].')
end

if (nargin > 4)
	if (prod(size(mu)) ~= 1) || ~isnumeric(mu) || ~isreal(mu) || ~isfinite(mu)
		error('All input arguments must be finite real scalars.')
	end
end
	
if (nargin > 5)
	if (prod(size(sigma)) ~= 1) || ~isnumeric(sigma) || ~isreal(sigma) || ~isfinite(sigma)
		error('All input arguments must be finite real scalars.')
	end
	if (sigma <= 0)
		error('Standard deviation must be greater than zero.')
	end
end

% Downsampling No. of points to match with those of Scott jackson (tau 1e-1)
resamp = ceil(1e-1/tdres);
nop = N; N = ceil(N/resamp)+1; 
if (N<10)
    N = 10;
end

% Determine whether fGn or fBn should be produced.
if ( Hinput <= 1 )
	H = Hinput;
	fBn = 0;
else
	H = Hinput - 1;
	fBn = 1;
end

% Calculate the fGn.
if (H == 0.5)
	y = randn(1, N);  % If H=0.5, then fGn is equivalent to white Gaussian noise.
else
    % If this function was already in memory before being called this time,
    % AND the values for N and H are the same as the last time it was
    % called, then the following (persistent) variables do not need to be
    % recalculated.  This was done to improve the speed of this function,
    % especially when many samples of a single fGn (or fBn) process are
    % needed by the calling function.
    persistent Zmag Nfft Nlast Hlast
    if isempty(Zmag) || isempty(Nfft) || isempty(Nlast) ||isempty(Hlast) || N ~= Nlast || H ~= Hlast
		% The persistent variables must be (re-)calculated.
        Nfft = 2^ceil(log2(2*(N-1)));
		NfftHalf = round(Nfft/2);
		
		k = [0:NfftHalf, (NfftHalf-1):-1:1];
		Zmag = 0.5 .* ( (k+1).^(2.*H) - 2.*k.^(2.*H) + (abs(k-1)).^(2.*H) );
		clear k
		
		Zmag = real(fft(Zmag));
		if ( any(Zmag < 0) )
			error('The fast Fourier transform of the circulant covariance had negative values.');
		end
        Zmag = sqrt(Zmag);
        
        % Store N and H values in persistent variables for use during subsequent calls to this function.
        Nlast = N;
        Hlast = H;
    end
    if noiseType == 0 % for fixed fGn
%         rng(16); % fixed seed from MATLAB
        randn('seed',37) % fixed seed from MATLAB
    end
    
	Z = Zmag.*(randn(1,Nfft) + i.*randn(1,Nfft));
	
	y = real(ifft(Z)) .* sqrt(Nfft);
	clear Z
	
	y((N+1):end) = [];
end

% Convert the fGn to fBn, if necessary.
if (fBn)
	y = cumsum(y);
end

% Resampling back to original (1/tdres): match with the AN model
y = resample(y,resamp,1);  % Resampling to match with the AN model

% define standard deviation
if (nargin < 6)

    if mu<0.2
        sigma = 1;%5  
    else

        if mu<20

            sigma = 10;
        else
             sigma = mu/2;
        end
    end
end
y = y*sigma;

y = y(1:nop);











