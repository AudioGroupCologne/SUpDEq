function difSTD=georganti2013(signal,P)
%GEORGANTI2013 Binaural spectral-magnitude difference standard deviation according to Georganti et al., (2013)
%   Usage: difSTD = georganti2013(signal, P)
%
%   Input parameters:
%       signal : binaural input signal
%
%       P.fs: sampling rate in Hz
%
%       P.timeFr: Frame size in seconds
%
%       P.fmin: lower frequency (Hz) for the BSDM STD calculation
%
%       P.fmax: upper frequency (Hz) for the BSDM STD calculation
%
%   Output parameters:
%       difSTD: Binaural spectral-magnitude difference standard deviation (dB)
%
%   See also: exp_georganti2013
%
%   References:
%     E. Georganti, T. May, S. van de Par, and J. Mourjopoulos. Extracting
%     sound-source-distance information from binaural signals. In J. Blauert,
%     editor, The Technology of Binaural Listening, Modern Acoustics and
%     Signal Processing, pages 171-199. Springer Berlin Heidelberg, 2013.
%     [1]http ]
%     
%     E. Georganti, T. May, S. van de Par, and J. Mourjopoulos. Sound source
%     distance estimation in rooms based on statistical properties of
%     binaural signals. Audio, Speech, and Language Processing, IEEE
%     Transactions on, 21(8):1727-1741, Aug 2013.
%     
%     References
%     
%     1. http://dx.doi.org/10.1007/978-3-642-37762-4_7
%     
%
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/georganti2013.php

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

%
%   Developed with Matlab 7.1.0.584 (R2010b) v.0.1
%   Updated with Matlab R2011b (7.13.0.564) v.1.0
%
%   Please send bug reports to:
%     Eleftheria Georganti
%     Postdoctoral Researcher
%     Experimental Audiology, ENT
%     University Hospital of Zurich/University of Zurich
%     Zurich, Switzerland
%     eleftheria.georganti@uzh.ch
%
%



if ~exist('P','var'), P=[]; end
  
if ~isfield(P,'fmin'), P.fmin=20; end % lower frequency in Hz - default value
if ~isfield(P,'fmax'), P.fmax=23000; end  % upper frequency in Hz - default value
if ~isfield(P,'fs'), P.fs=44100; end
if ~isfield(P,'timeFr'), P.timeFr=1; end

P.sampleFr = P.fs * P.timeFr; % Frame size in samples
if ~isfield(P,'hop'), P.hop = P.sampleFr/2; end % Overlap
P.nFFT = P.sampleFr; % FFT points

P.freq = (P.fs/P.nFFT)*(0:(P.nFFT/2-1)); % Frequency index


fmin_id = min(find((P.freq>P.fmin)));
fmax_id = min(find((P.freq>P.fmax)));

difSTD=zeros(1,length(1:P.hop:length(signal)-P.hop));
idx = 1;

for kk = 1:P.hop:length(signal)-P.hop

    % Calculate magnitude spectrums in dB of the left & right signals
    leftFFT  = 20*log10(abs(fft(signal(kk:kk+P.hop-1,1))));
    rightFFT = 20*log10(abs(fft(signal(kk:kk+P.hop-1,2))));

    % Subtract the magnitude spectrums
    specDIF  = leftFFT(1:end/2)-rightFFT(1:end/2);

    % Calculate the differential standard deviation for the
    % frequency range of interest
    difSTD(1,idx) = std(specDIF(fmin_id:fmax_id));

    idx = idx+1;

end
    


