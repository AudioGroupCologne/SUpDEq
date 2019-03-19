function outsig = sig_itdildsin(fc,itd,ild,fs)
%sig_itdildsin Generate a sinusoid with a interaural time difference
%   Usage: outsig = sig_itdildsin(fc,itd,ild,fs)
%
%   Input parameters:
%       fc      : carrier frequency of the sinusoid (Hz)
%       itd     : ITD of the left signal, this can be positive or negative (ms)
%       ild     : ILD of the right signal, this can be positive or negative (dB)
%       fs      : sampling rate (Hz)
%
%   Output parameters:
%       outsig  : two channel 1 s long sinusoid
%
%   SIG_ITDILDSIN(fc,itd,ild,fs) generates a sinusoid with a interaural time 
%   difference of itd, a interaural level difference of ild and a frequency of 
%   fc.
%
%   The output is scaled to have a maximum value of 1-eps.  
%
%   References:
%     B. C. J. Moore. An Introduction to the Psychology of Hearing. Academic
%     Press, 5th edition, 2003.
%     
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/signals/sig_itdildsin.php

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

% AUTHOR: Hagen Wierstorf


% ------ Checking of input parameters ---------

error(nargchk(4,4,nargin));

if ~isnumeric(fc) || ~isscalar(fc) || fc<0
    error('%s: f must be a positive scalar.',upper(mfilename));
end

if ~isnumeric(itd) || ~isscalar(itd)
    error('%s: itd must be a scalar.',upper(mfilename));
end

if ~isnumeric(ild) || ~isscalar(ild)
    error('%s: ild must be a scalar.',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs must be a positive scalar!',upper(mfilename));
end


% ------ Computation --------------------------

% Create a one second time 
t = (1:fs)/fs;
% Right signal
sigr = sin(2*pi*fc.*t);
% Time shift in samples
itdsamples = ceil(fs * abs(itd)/1000);
% Left signal with ITD shift
sigl = [zeros(1,itdsamples) sin(2*pi*fc.*t(1:end-itdsamples))];
% Combine left and right signal to outsig
% Check if we have a positive or negative ITD and switch left and right signal
% for negative ITD
if itd<0
    % Apply ILD
    sigl = gaindb(sigl,ild);
    outsig = [sigr' sigl'];
else
    % Apply ILD
    sigr = gaindb(sigr,ild);
    outsig = [sigl' sigr'];
end
% Scale outsig
outsig = outsig / (max(abs(outsig(:)))+eps);

