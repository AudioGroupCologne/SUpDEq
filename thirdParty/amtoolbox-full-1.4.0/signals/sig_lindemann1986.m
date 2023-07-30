function outsig = sig_lindemann1986(fc,mf,fs)
%SIG_lINDEMANN1986 Generate a binaural modulated sinus for the Lindemann (1986) model
%   Usage: outsig = sig_lindemann1986(fc,mf,fs)
%
%   Input parameters:
%       fc  : carrier frequency of the sinus (Hz)
%       mf  : binaural modulation frequency (Hz)
%       fs  : sampling rate (Hz)
%
%   Output parameters:
%       outsig  : fs x 2 sinusoid signal
%
%   SIG_LINDEMANN1986(fc,mf,fs) generates an binaural modulated sinusoid with a
%   carrier frequency of f and a frequency moving around the two ears of
%   mf.
%
%   See also: lindemann1986 demo_lindemann1986
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_lindemann1986.php


%   #Author: Hagen Wierstorf
%   #Author: Peter Soendergaard
%   #Author: Piotr Majdak

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
    error('%s: f has to be a positive scalar.',upper(mfilename));
end

if ~isnumeric(mf) || ~isscalar(mf) || mf<=0
    error('%s: mf has to be a positive scalar.',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end


% ------ Computation --------------------------
% Create a one second time 
t = (1:fs)/fs;
% Left signal
sigl = sin(2*pi*fc.*t);
% Right signal with amplitude modulation
sigr = sin(2*pi*fc.*t + sin(2*pi*mf.*t));
outsig = [sigl' sigr'];
% Scale outsig
outsig = outsig / (max(abs(outsig(:)))+eps);



