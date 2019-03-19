function outsig = sig_ildsin(fc,ild,fs)
%sig_ildsin Sinusoid with a interaural level difference (ILD)
%   Usage: outsig = sig_itdsin(fc,ild,fs)
%
%   Input parameters:
%       fc      : carrier frequency of the sinusoid (Hz)
%       ild     : ILD of the right signal, positive or negative (dB)
%       fs      : sampling rate (Hz)
%
%   Output parameters:
%       outsig  : two channel 1 s long sinusoid 
%
%   SIG_ILDSIN(fc,ild,fs) generates a sinusoid with a interaural level difference
%   of ild and a frequency of fc.
%
%   The output is scaled to have a maximum value of 1-eps.
%
%   References:
%     B. C. J. Moore. An Introduction to the Psychology of Hearing. Academic
%     Press, 5th edition, 2003.
%     
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/signals/sig_ildsin.php

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

error(nargchk(3,3,nargin));

if ~isnumeric(fc) || ~isscalar(fc) || fc<0
    error('%s: fc must be a positive scalar.',upper(mfilename));
end

if ~isnumeric(ild) || ~isscalar(ild)
    error('%s: itd must be a scalar.',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs must be a positive scalar!',upper(mfilename));
end


% ------ Computation --------------------------

% Create a one second time 
t = (1:fs)/fs;
% Left signal
sigl = sin(2*pi*fc.*t);
% Right signal with level difference of ILD
sigr = gaindb(sigl,ild);
% Combine left and right channel
outsig = [sigl' sigr'];
% Scale outsig
outsig = outsig / (max(abs(outsig(:)))+eps);


