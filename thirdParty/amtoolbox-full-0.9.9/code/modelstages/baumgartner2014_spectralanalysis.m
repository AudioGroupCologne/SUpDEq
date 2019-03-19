function varargout = baumgartner2014_spectralanalysis(sig,varargin)
%baumgartner2014_spectralanalysis - Approximation of spectral analysis by auditory periphery
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
%     Acoustical Society of America, 136(2):791-802, 2014.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2014_spectralanalysis.php

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

% AUTHOR: Robert Baumgartner

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
