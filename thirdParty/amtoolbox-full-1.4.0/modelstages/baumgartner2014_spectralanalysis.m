function varargout = baumgartner2014_spectralanalysis(sig,varargin)
%BAUMGARTNER2014_SPECTRALANALYSIS Approximation of spectral analysis by auditory periphery
%   Usage:  [mp,fc] = baumgartner2014_spectralanalysis(sig)
%
%   Input parameters:
%     sig     : incoming time-domain signal
%
%   Output parameters:
%     mp      : spectral magintude profile
%     fc      : center frequencies of auditory filters
%
%   BAUMGARTNER2014_SPECTRALANALYSIS(...) computes temporally integrated
%   spectral magnitude profiles.
%
%   BAUMGARTNER2014_SPECTRALANALYSIS accepts the following optional parameters:
%
%     'flow',flow    Set the lowest frequency in the filterbank to
%                    flow. Default value is 700 Hz.
%
%     'fhigh',fhigh  Set the highest frequency in the filterbank to
%                    fhigh. Default value is 18000 Hz.
%
%     'fs',fs        Define the sampling rate of the impulse responses. 
%                    Default value is 48000 Hz.
%
%     'space',sp     Set spacing of auditory filter bands (i.e., distance 
%                    between neighbouring bands) to sp in number of
%                    equivalent rectangular bandwidths (ERBs). 
%                    Default value is 1 ERB.
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2014_spectralanalysis.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA CircStat M-SIGNAL M-Stats O-Statistics
%   #Author: Robert Baumgartner (2014), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


definput.import={'baumgartner2014'};
[flags,kv]=ltfatarghelper({},definput,varargin);

%% Spectral Analysis, Eq.(2)

if kv.space == 1 % Standard spacing of 1 ERB
  [mp,fc] = auditoryfilterbank(sig(:,:),kv.fs,'flow',kv.flow,'fhigh',kv.fhigh);
else
  fc = audspacebw(kv.flow,kv.fhigh,kv.space,'erb');
  [bgt,agt] = gammatone(fc,kv.fs,'complex');
  mp = 2*real(ufilterbankz(bgt,agt,sig(:,:)));  % channel (3rd) dimension resolved
end
Nfc = length(fc);   % # of bands

% Set back the channel dimension
mp = reshape(mp,[size(sig,1),Nfc,size(sig,2),size(sig,3)]);

% Averaging over time (RMS)
mp = 20*log10(squeeze(rms(mp)));      % in dB

if size(mp,2) ~= size(sig,2) % retreive polar dimension if squeezed out
    mp = reshape(mp,[size(mp,1),size(sig,2),size(sig,3)]);
end

varargout{1} = mp;
if nargout == 2
  varargout{2} = fc;
end

end


