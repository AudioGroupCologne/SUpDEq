function outsig = sig_yost1996(d,iterations,gn,siglen,fs);
%sig_yost1996	Generate iterated rippled noise from Yost (1996)
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
%     the Acoustical Society of America, 100(5):3329-3335, Nov. 1996.
%     
%     
%
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/signals/sig_yost1996.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% AUTHOR: Hagen Wierstorf, Daniel Pressnitzer, Stefan Uppenkamp
%         09.07.2017 Piotr Majdak

% ------ Checking of input parameters ---------

error(nargchk(5,5,nargin));

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
outsig = outsig/rms(outsig);


