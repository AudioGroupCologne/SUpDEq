function [Mixsigmin, Mixsigmax, ECparams4Opt] = hauth2020_fftcon(MixsigL,MixsigR,fc,fs,sigmadelta0,Delta0,bin_inaccuracy)
%HAUTH2020_FFTCON Cancellation process
%
%   Usage: 
%     [Mixsigmin, Mixsigmax, ECparams4Opt] = hauth2020_fftcon(MixsigL,MixsigR,fc,fs,sigmadelta0,Delta0,bin_inaccuracy)
%
%   Input parameters: 
%     MixsigL        : Left ear signal
%     MixsigR        : Right ear signal
%     fc             : center frequency of frequency channel
%     fs             : sampling frequency
%     sigmadelta0    : term to calculate binaural processing inaccuracy
%     Delta0         : term to calculate binaural processing inaccuracy
%     bin_inaccuracy : flag indicating the use of (1) binaural processing
%                      inaccuracies or (0) assuming binaural processing to be deterministic
%
%   Output parameters:
%     Mixsigmin      : EC processed signal using the minimization strategy
%     Mixsigmax      : EC processed signal using the maximization strategy
%     ECparams4Opt   : Structure containing EC parameters (Delay and Uncertainties)
%
%   First, left and right ear signal are transformed to
%   frequency domain in order to calculate the cross power spectral density.
%   By calculating the angle of the CPSD and, the frequency dependent ITD can
%   be obtained by dividing the angle by the angular frequency. In order to
%   obtain the final ITD value for the frequency band, an integration window
%   derived from the absolute value of the CPSD is applied.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/hauth2020_fftcon.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Author: Christopher F. Hauth (2020)
%   #Author: Dr. Thomas Brand (2020)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% caluclate number of fft points for computation of correlation
lenSig = length(MixsigL);
fftpoints = (length(MixsigL)+length(MixsigR));
% compute fft representation
MixsigLFFT = fft(MixsigL,fftpoints); 
MixsigRFFT = fft(MixsigR,fftpoints);

%% Frequency domain fractionaly delay estimation
  % Apply ERB window in the frequency domain
ERB       = 24.7.*(4.37.*fc./1000 + 1);
 % compute CrossPowerSpectralDensity
PHY_LR    = MixsigLFFT.*conj(MixsigRFFT);

Frequency = linspace(0, 44100, fftpoints);
PHY_phase = angle(PHY_LR);

% Make sure phase is consistent in case of pi (+-pi)
PHY_phase_unwrapped = PHY_phase;
PHY_phase_unwrapped(PHY_phase<0) = PHY_phase(PHY_phase<0)+2*pi;
% Convert IPD to ITD
PHY_time      = PHY_phase./(2.*pi.*Frequency');
PHY_time_unwrapped      = PHY_phase_unwrapped./(2.*pi.*Frequency');

% Normalize
PHY_abs_norm = abs(PHY_LR)./max(abs(PHY_LR));
window = find(Frequency>=fc-ERB/2&Frequency<=fc+ERB/2);

% compute standard deviation of tau 
EC_std = std(PHY_time(window));
% do the same for unwrapped phase:
EC_std_unwrapped = std(PHY_time_unwrapped(window));
% compute standard deviation of the delay processign inaccuracy
if  EC_std <= EC_std_unwrapped
    EC_tau = sum(PHY_time(window).*(PHY_abs_norm(window)/sum(PHY_abs_norm(window))));
else
    EC_tau = sum(PHY_time_unwrapped(window).*(PHY_abs_norm(window)/sum(PHY_abs_norm(window)))); %EC_tau_unwrapped 
end
std_min_frac = sigmadelta0.*(1+(abs(EC_tau)./(Delta0)));
% compute a set of 2 Gaussian variables for each ear
errorLR = 1.*std_min_frac.*randn(2,1);
% Save EC params in a struct. It will also be used to EC process optional
% signals
ECparams4Opt.errorLR = errorLR;
ECparams4Opt.fftpoints = fftpoints; 
ECparams4Opt.EC_tau = EC_tau; 
% apply EC mechanism in the frequency domain:
[Mixsigmin, Mixsigmax] =local_ecfftprocess(MixsigLFFT,MixsigRFFT,EC_tau,errorLR,fs,fftpoints,lenSig,bin_inaccuracy);
end

function [Mixsigmin, Mixsigmax] = local_ecfftprocess(MixsigLFFT,MixsigRFFT,EC_tau,errorLR,fs,fftpoints,orig_len,bin_inaccuracy)
% Usage: [Mixsigmin Mixsigmax] = EC_FFTprocess(MixsigLFFT,MixsigRFFT,EC_tau,errorLR,fs,fftpoints,orig_len,bin_inaccuracy)
% This function applies the Equalization-Cancellation process in the frequency
% domain.
% The Equalization is applied two both ears symmetrically.
% Input: 
% MixsigLFFT - FFT representation of left ear signal 
% MixsigRFFT - FFT representation of right ear signal
% EC_tau     - Delay for the EC process
% errorLR    - Pair of binaural processing inaccuracies
% fs         - sampling frequency
% fftpoints  - fftpoints in the fft
% orig_len   - original length of the signal
% bin_inaccuracy - flag indicating the use of (1) binaural processing
% inaccuracies or (0) assuming binaural processing to be deterministic

% Output: 
% Mixsigmin - EC processed signal applying minimization strategy
% Mixsigmax - EC processed signal applying maximization strategy
% Authors: Christopher Hauth <christopher.hauth@uni-oldenburg.de>
%          Dr. Thomas Brand  <thomas.brand@uni-oldenburg.de>
% Date: 22.10 2020
% Version: 1.0.


 f = linspace(0,fs/2,floor(fftpoints./2+1))';
 % invert frequency vector for complex conjugate
 f_inv = f(end-1:-1:2); 
 % frequency vector in rad
 omega = 2.*pi*[f;f_inv];
 % copy signals to generate buffer for processed signals
 FFTMixsigLeq = MixsigLFFT;
 FFTMixsigReq = MixsigRFFT;
%--------------------------------------------------------------------------
                %% Apply Equalization Mechanism %%
% calculate phase factor for left and right ear with a set of gaussian
% distributed RVs as binaural processing inaccuracies
if bin_inaccuracy
    for kk = 1:size(errorLR,2)
        phaseL(:,kk) = exp(-1j.*omega(1:(end/2)+1).*(EC_tau./2 + errorLR(1,kk)));
        phaseR(:,kk) = exp(+1j.*omega(1:(end/2)+1).*(EC_tau./2 + errorLR(2,kk)));
        % complex conjugate
        phaseLcc(:,kk) = exp(+1j.*omega((end/2)+2:end).*(EC_tau./2 + errorLR(1,kk)));
        phaseRcc(:,kk) = exp(-1j.*omega((end/2)+2:end).*(EC_tau./2 + errorLR(2,kk)));
    end
    % take the mean over all phase terms to obtain averaged processing
    % binaural processing inaccuracy
    phaseL   = mean(phaseL,2);
    phaseLcc = mean(phaseLcc,2);
    phaseR   = mean(phaseR,2);
    phaseRcc = mean(phaseRcc,2);
else
        phaseL = exp(-1j.*omega(1:(end/2)+1).*(EC_tau./2));
        phaseR = exp(+1j.*omega(1:(end/2)+1).*(EC_tau./2));
        % complex conjugate
        phaseLcc = exp(+1j.*omega((end/2)+2:end).*(EC_tau./2));
        phaseRcc = exp(-1j.*omega((end/2)+2:end).*(EC_tau./2));
end
% equalize phases of left and right ear signals
% Left ear
FFTMixsigLeq(1:(end/2)+1)  = FFTMixsigLeq(1:(end/2)+1).*phaseL;
FFTMixsigLeq((end/2)+2:end)= FFTMixsigLeq((end/2)+2:end).*phaseLcc;
   
% Right ear 
FFTMixsigReq(1:(end/2)+1)  = FFTMixsigReq(1:(end/2)+1).*phaseR;
FFTMixsigReq((end/2)+2:end)= FFTMixsigReq((end/2)+2:end).*phaseRcc;
%--------------------------------------------------------------------------
                    %% Apply Cancellation %%
FFTMixsigMin = FFTMixsigLeq - FFTMixsigReq; % destructive interference
FFTMixsigMax = FFTMixsigLeq + FFTMixsigReq; % constructive interference
%--------------------------------------------------------------------------
                    %% Apply IFFT and cut signals
 % Level minimization
 Mixsigmin = real(ifft(FFTMixsigMin));
 Mixsigmin = Mixsigmin(1:orig_len); 
 % Level maximization
 Mixsigmax = real(ifft(FFTMixsigMax));
 Mixsigmax = Mixsigmax(1:orig_len);    
end
%--------------------Licence ---------------------------------------------
% Copyright (c) <2020> Christopher F. Hauth
% Dept. Medical Physics and Acoustics
% Carl von Ossietzky University Oldenburg 
% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to 
% permit persons to whom the Software is furnished to do so, subject 
% to the following conditions:
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
% END OF FILE


