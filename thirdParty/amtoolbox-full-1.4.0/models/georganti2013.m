function difSTD=georganti2013(signal,P)
%GEORGANTI2013 Distance estimation
%   Usage: difSTD = georganti2013(signal, P)
%
%   Input parameters:
%       signal   : binaural input signal
%
%       P        : input struct with the following fields
%
%                  .fs     : sampling rate in Hz
%
%                  .timeFr : Frame size in seconds
%
%                  .fmin   : lower frequency (Hz) for the BSDM STD calculation
%
%                  .fmax   : upper frequency (Hz) for the BSDM STD calculation
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
%     Signal Processing, pages 171--199. Springer Berlin Heidelberg, 2013.
%     [1]http ]
%     
%     E. Georganti, T. May, S. van de Par, and J. Mourjopoulos. Sound source
%     distance estimation in rooms based on statistical properties of
%     binaural signals. Audio, Speech, and Language Processing, IEEE
%     Transactions on, 21(8):1727--1741, Aug 2013.
%     
%     References
%     
%     1. http://dx.doi.org/10.1007/978-3-642-37762-4_7
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/georganti2013.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Author: Eleftheria Georganti (2013)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

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
fmax_id = min(find((P.freq>=P.fmax)));

difSTD=zeros(1,length(1:P.hop:length(signal)-P.hop));
idx = 1;

for kk = 1:P.hop:length(signal)-P.hop

    % Calculate magnitude spectrums in dB of the left & right signals
    leftFFT  = 20*log10(abs(fft(signal(kk:kk+P.hop-1,1),P.nFFT)));
    rightFFT = 20*log10(abs(fft(signal(kk:kk+P.hop-1,2),P.nFFT)));

    % Subtract the magnitude spectrums
    specDIF  = leftFFT(1:end/2)-rightFFT(1:end/2); 

    % Calculate the differential standard deviation for the
    % frequency range of interest
    difSTD(1,idx) = std(specDIF(fmin_id:fmax_id));

    idx = idx+1;

end


