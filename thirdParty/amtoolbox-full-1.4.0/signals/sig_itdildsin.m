function outsig = sig_itdildsin(fc,itd,ild,fs)
%SIG_ITDILDSIN Generate a sinusoid with a interaural time difference
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_itdildsin.php


%   #Author: Hagen Wierstorf

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% ------ Checking of input parameters ---------

narginchk(4,4);

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


