function y = zilany2014_ffgn(N, tdres, Hinput, noiseType, mu, sigma)
%ZILANY2014_FFGN  Fast (exact) fractional Gaussian noise and Brownian motion generator
%   Usage: Y = zilany2014_ffGn(N, tdres, Hinput)
%
%   Input parameters:
%     N       : is the length of the output sequence.
%     tdres   : is the time resolution (1/sampling rate)
%     Hinput  : is the "Hurst" index of the resultant noise (0 < H <= 2).  For 0 < H <= 1, 
%               the output will be fractional Gaussian noise with Hurst index H.  For 
%               1 < H <= 2, the output will be fractional Brownian motion with Hurst
%               index H-1.  Either way, the power spectral density of the output will
%               be nominally proportional to 1/f^(2H-1).
%
%   zilany2014_ffGn(...) returns a vector containing a sequence of fractional Gaussian 
%   noise or fractional Brownian motion.  The generation process uses an FFT 
%   which makes it very fast. This method is based on an embedding of the 
%   covariance matrix in a circulant matrix.
%
%
%   ZILANY2014_FFGN accepts the following optional parameters:
%
%     'noiseType'    is 0 for fixed fGn noise and 1 for variable fGn. [default = 1]
%
%     'mu'           is the mean of the noise. [default = 0]
%
%     'sigma'        is the standard deviation of the noise. [default = 1]
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/zilany2014_ffgn.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB MEX M-Signal
%   #Author: Muhammad Zilany 
%   #Author: B. Scott Jackson (2005)
%   #Author: Robert Baumgartner: adapted to the AMT
%   #Author: Clara Hollomey (2020): adapted to AMT 1.0
%   #Author: Piotr Majdak (2021): C1 and C2 outputs

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

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
        % rng(16); % fixed seed from MATLAB
        randn('seed',37) % fixed seed from MATLAB
    end
    
	Z = Zmag.*(randn(1,Nfft) + 1i.*randn(1,Nfft));
	
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
    if mu<0.5
        sigma = 3;%5  
    else
        if mu<18
            sigma = 30;%50   % 7 when added after powerlaw
        else
            sigma = 200;  % 40 when added after powerlaw        
        end
    end
end
y = y*sigma;

y = y(1:nop);


