function signalout = sig_breebaart2001(inttype,fc,sl,sdur,sphase,nbw,nl,ndur,nphase,hannramp,fs)
%sig_breebaart2001 computes the signals for the breebaart2001 experiment
%
%   Usage:   signalout = sig_breebaart2001(inttype,fc,sl,sdur,sphase,nbw,nl,ndur,nphase,hannramp,fs)
%
%   signalout = SIG_BREEBAART2001(inttype,fc,sl,sdur,sphase,nbw,nl,ndur,nphase,hannramp,fs)
%   generates the signals for the breebaart2001 experiment  where
%   sinusoidal signals are masked by noise. 'inttyp' determines if it is a 
%   target interval (signal + noise) or just a reference interval (noise 
%   only). The sinusoidal has a center frequency of 'fc' (in Hz), an 
%   overall level of 'sl' (in dB SPL), a duration of 'sdur' (in seconds) 
%   and an interaural phase difference of 'sphase' (in rad).  The masker is 
%   generated at a center frequency 'fc' in Hz,  with a  bandwidth of nbw
%   (in Hz), duration ndur (in sec), overall level nl (in dB SPL) and 
%   an interaural phase difference of 'nphase'. Both signal and masker are 
%   gated with hanning ramps with a duration of 'hannramp' (in s). The 
%   sampling rate is defined with fs (in Hz).
%
%   Input parameters:
%       inttyp     : 'target' or 'reference' defines the interval type
%       fc         : center frequency
%       sl         : overall signal level
%       sl         : duration of the signal
%       sphase     : interaural phase difference of signal
%       nbw        : bandwidth of the noise
%       nl         : overall noise level
%       ndur       : duration of the noise
%       nphase     : interaural phase difference of noise
%       hannramp   : duration of the hanning ramps
%       fs         : sampling frequency
%
%   The following parameters for 'nphase' are possible:
%       'nphase' = 0      no phase difference, correlation = 1; 
%       'nphase' = pi     masker interaurally phase reversed, correlation = -1;
%       0 < 'nphase' < 1  defines interaural correlation 
%
%   Output parameters:
%       signalout:  the generated experimental signal which can be feeded 
%                   directly into the model
%
%   References:
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. II. Dependence on spectral
%     parameters. J. Acoust. Soc. Am., 110:1089-1104, August 2001.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/signals/sig_breebaart2001.php

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

%   Author:     Martina Kreuzbichler

% generate noise
noise(:,1) = sig_bandpassnoise(fc,fs,ndur,nl,nbw);

% set diotic phase difference for noise masker
if nphase == 0
    noise(:,2) = noise(:,1);
elseif nphase == pi;
    noise(:,2) = -noise(:,1);
else % set correlation
    noise(:,2) = sig_bandpassnoise(fc,fs,ndur,nl,nbw);
    noiseL = 0.5*sqrt(2)*sqrt(1+nphase)*noise(:,1) + 0.5*sqrt(2)*...
        sqrt(1-nphase)*noise(:,2);
    noiseR = 0.5*sqrt(2)*sqrt(1+nphase)*noise(:,1) - 0.5*sqrt(2)*...
        sqrt(1-nphase)*noise(:,2);
    noise = [noiseL, noiseR];
end

% apply hanning ramps to noise
n_ramp = round(hannramp*fs); 
noise= rampsignal(noise,n_ramp);

% set noiselevel
noise = setdbspl(noise,nl);

if strcmp(inttype,'target')

    signal(:,1) = sin(2*pi*(0:(sdur*fs)-1)'*fc/fs);

    % set diotic phase difference for signal
    if sphase == 0
        signal(:,2) = signal(:,1);
    else
        signal(:,2) = sin(2*pi*(0:(sdur*fs)-1)'*fc/fs-sphase);
    end

    % apply hanning ramps to signal  
    signal = rampsignal(signal,n_ramp);

    % set signal level
    signal = setdbspl(signal,sl);

    if sdur < ndur
        zerolength = round((ndur-sdur)*fs/2);
        signal = [zeros(zerolength,2); signal; zeros(zerolength,2)];
    end
    
    signalout = signal + noise;

elseif strcmp(inttype,'reference')
    signalout = noise;
end
    


