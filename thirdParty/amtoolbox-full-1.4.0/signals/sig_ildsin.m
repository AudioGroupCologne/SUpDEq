function outsig = sig_ildsin(fc,ild,fs)
%SIG_ILDSIN Sinusoid with a interaural level difference (ILD)
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_ildsin.php


%   #Author: Hagen Wierstorf

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% ------ Checking of input parameters ---------

narginchk(3,3);

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



