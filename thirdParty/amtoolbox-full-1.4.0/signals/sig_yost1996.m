function outsig = sig_yost1996(d,iterations,gn,siglen,fs)
%SIG_YOST1996	Generate iterated rippled noise from Yost (1996)
%   Usage: outsig=sig_yost1996(delay,iterations,gn,siglen,fs)
%
%   Input parameters:
%      d          : delay in ms of the time-shifted noise adding process
%      iterations : number of iterations of the adding process
%      gn         : relative gain of irn
%      siglen	  : signal length of irn in samples
%      fs         : sampling rate in Hz
%
%   SIG_YOST1996(d,iterations,gn,siglen,fs) generates a signal consisting of
%   white noise, with iterations added to itself with a delay of d (in
%   ms).
%
%   An example:
%
%     fs = 44100;
%     signal = sig_yost1996(4,6,1,fs,fs);
%     sound(signal,fs)
%
%   References:
%     W. A. Yost. Pitch strength of iterated rippled noise. The Journal of
%     the Acoustical Society of America, 100(5):3329--3335, Nov. 1996.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_yost1996.php


%   #Author: Hagen Wierstorf
%   #Author: Daniel Pressnitzer
%   #Author: Stefan Uppenkamp
%   #Author: Piotr Majdak (2017)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% ------ Checking of input parameters ---------

narginchk(5,5);

% ------ Computation --------------------------

% Frequency to which the delay corresponds
freq = 1000/d;
% Number of samples for the noise (slightly longer to avoid circular iterations 
% (S. Uppenkamp)
noiselen = siglen + round(iterations*fs/freq);
% Number of samples for the delay
delaylen = round(fs/freq);
% Create white noise
noisesig = randn(1,noiselen);

% Iterate delay and add n times
for ii = 1:iterations
  dnoise(delaylen+1:noiselen) = noisesig(1:noiselen-delaylen);
  dnoise(1:delaylen) = noisesig(noiselen-delaylen+1:noiselen);
  noisesig = noisesig + gn*dnoise;
end

% Take first bit of result as IRN
outsig = noisesig(1:siglen);
% Scale to RMS of 1
outsig = outsig./rms(outsig);



