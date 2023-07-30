function [cqmag,fc,cqmaghr,fvec] = langendijk2002_spectralanalysis( insig,varargin )
%LANGENDIJK2002_SPECTRALANALYSIS FFT-based filter bank with constant relative bandwidth according
%   Usage:  [cqmag] = langendijk2002_spectralanalysis( insig )
%           [cqmag,fc,cqmaghr,fvec] = langendijk2002_spectralanalysis( insig,fs,flow,fhigh,bw )
%
%   Input parameters:
%      insig   : Impulse response or complex spectrum
%      fs      : Sampling rate, default is 48kHz.
%      flow    : Lowest frequency, minimum: 0.5kHz, default is 2kHz
%      fhigh   : Highest frequency, default is, default is 16kHz  
%      bw      : bandwidth, possible values 3,6,9,12, default is 6.
%
%   Output parameters:
%      cqmag   : mean magnitudes of CQ-bands in dB
%      fc      : center frequencies of bands (geo. mean of corners)
%      cqmaghr : same as cqmag but for all freq. bins (high resolution)
%      fvec    : freq. vector according to FFT-resolution
%
%   LANGENDIJK2002_SPECTRALANALYSIS(insig) approximates a constant-Q filter bank by averaging the
%   magnitude bins of a DFT. LANGENDIJK2002_SPECTRALANALYSIS results in 'bw' dB-magnitudes per octave.
%
%   References:
%     E. Langendijk and A. Bronkhorst. Contribution of spectral cues to human
%     sound localization. J. Acoust. Soc. Am., 112:1583--1596, 2002.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/langendijk2002_spectralanalysis.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: none
%   #Author : Robert Baumgartner (2013), OEAW Acoustical Research Institute

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.keyvals.fs=48000;
definput.keyvals.flow=2000;
definput.keyvals.fhigh=16000;
definput.keyvals.bw=6;

[flags,kv]  = ltfatarghelper({'fs','flow','fhigh','bw'},definput,varargin);


% input signal given in time or frequency domain?
if isreal(insig)    % -> TD
    nfft = 2^12;%max(2^12,size(insig,1));
    y = abs(fft(insig,nfft));
else                % -> FD
    y = abs(insig);
    nfft = size(insig,1);
end
fvec = 0:kv.fs/nfft:kv.fs-kv.fs/nfft;

octs = log2(kv.fhigh/kv.flow);    % # of octaves
jj = 0:octs*kv.bw; 
n = round(2.^((jj)/kv.bw)*kv.flow/kv.fs*nfft);       % startbins
fc = zeros(length(jj)-1,1);                       % center frequencies
cqmag = zeros(length(jj)-1,size(y,2),size(y,3));  % mean magnitudes of CQ-bands
cqmaghr = zeros(size(y));                         % same but for all freq. bins (high resolution)
for ind = jj(1)+1:jj(end)
    nj = n(ind+1)-n(ind);
    idn = n(ind):n(ind+1)-1;
    fc(ind) = sqrt(fvec(n(ind))*fvec(n(ind+1)));  % geometric mean
    cqmag(ind,:,:) = sqrt(1/(nj)*sum(y(idn,:,:).^2,1));
    cqmaghr(idn,:,:) = repmat(cqmag(ind,:,:),[length(idn),1,1]);
end

cqmag = 20*log10(cqmag);
cqmaghr(nfft/2+2:end,:,:) = cqmaghr(nfft/2:-1:2,:,:);
cqmaghr = 20*log10(cqmaghr);

end



