function [sig_out,gt_out,fc] = bischof2023_filterbank(sig_in,N,len,fs,f_range)
%BISCHOF2023_FILTERBANK filterbank used in model bischof2023
%   Usage: [sig_out,gt_out,fc] = bischof2023_filterbank(sig_in,N,len,fs,f_range)
%
%
%   Input parameters:
%        sig_in      : Input time signal to be filtered
%        N           : Filter order of the gammatone filters (default: 4)
%        len         : Filter length in samples (default: 512)
%        fs          : Sampling frequency in Hz(default: 44100)
%        f_range     : (optional) Defines the upper and lower frequencies for
%                      the bark filterbank and chooses the right channel
%
%   Output parameters:
%       sig_out     : Matrix with filtered input signal. The matrix dimensions are:
%
%                     size(sig_out,1): length(signal)
%
%                     size(sig_out,2): numer of used bark bands(either defined by fs of f_range)
%
%                     size(sig_out,3): number of channels of sig_in
%
%       gt_out      : Matrix with all gammatone filters
%       fc          : Center frequencies of used bark filters
%
%
%   BISCHOF2023_FILTERBANK
%   creates BARK-scaled gammatone filters and filters the input signal.
%   Corner frequencies of Bark Scale according to Zwicker (1961).
%
%   See also: bischof2023 exp_bischof2023
%
%
%   References:
%     E. Zwicker. Subdivision of the audible frequency range into critical
%     bands (frequenzgruppen). J. Acoust. Soc. Am., 33(2):248--248, 1961.
%     
%     T. Gunawan and E. Ambikairajah. Speech enhancement using temporal
%     masking and fractional bark gammatone filters. In 10th International
%     Conference on Speech Science & Technology, pages 420--425, Sydney,
%     Australia, 2004.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/bischof2023_filterbank.php


%   #Author: Norbert F. Bischof (2023)
%   #Author: Pierre G. Aublin
%   #Author: Bernhard Seeber (2023)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.

%%
% Corner Frequencies of Bark Scale according to
% E. Zwicker. Subdivision of the Audible Frequency Range into Critical 
% Bands (Frequenzgruppen). J. Acoust. Soc. Am., 33(2):248, 1961.
BarkCornerFreq = [10,100,200,300,400,510,630,770,920,1080,1270,1480,1720,...
                  2000,2320,2700,3150,3700,4400,5300,6400,7770,9500,12000,15500];

% define default variables
if nargin < 5; f_range = [10 15500]; end % default frequency range from 10 Hz to 15.5 kHz
if nargin < 4; fs = 44100; end % default sampling frequency of 44.1 kHz
if nargin < 3; len = 512; end % default length of gammatone filter in samples
if nargin < 2; N = 4; end % default order for the exponential filter

% limits given f_range to Bark corner frequencies from Zwicker
if f_range(1) > f_range(2)
    error('Lower frequency bound must be smaller than upper frequency bound!');
end
if f_range(2) > BarkCornerFreq(end)
    f_range(2) = BarkCornerFreq(end);
end
if f_range(1) < BarkCornerFreq(1)
    f_range(1) = BarkCornerFreq(1);
end

idx_l = find((BarkCornerFreq-f_range(1)) >= 0);
idx_u = find((BarkCornerFreq-f_range(2)) >= 0);

if (BarkCornerFreq(idx_l(1))-f_range(1)) == 0
    f_range(1) = BarkCornerFreq(idx_l(1));
else
    f_range(1) = BarkCornerFreq(idx_l(1)-1);
end
f_range(2) = BarkCornerFreq(idx_u(1));

% time vector for filter in seconds
t = (1/fs:1/fs:len/fs)';

% determine used Bark Bands
tmp = find(BarkCornerFreq <= f_range(1));
idx_start = tmp(end);
if ~ismember(f_range(1),BarkCornerFreq)
    idx_start = idx_start +1;
end

tmp = find(BarkCornerFreq >= f_range(2));
idx_end = tmp(1);
if ~ismember(f_range(2),BarkCornerFreq)
    idx_end = idx_end-1;
end

% derives band width and center frequency of the selected auditory filters
% defined with upper and lower corner frequencies
CorFreq = BarkCornerFreq(idx_start:idx_end);
CenterFreq = zeros(1,length(CorFreq)-1);
BandWidth = zeros(1,length(CorFreq)-1);
for ii = 1:length(CenterFreq)
    BandWidth(1,ii) = CorFreq(1,ii+1) - CorFreq(1,ii);
    CenterFreq(1,ii) = CorFreq(1,ii) + BandWidth(1,ii)/2;
end
fc = CenterFreq; % for output variable

% constanc b taken from
% Gunawan,T.S. & Ambikairajah,E.(2004).Speech Enhancement Using Temporal
% Masking and fractional bark gammatone filters. In proceedings of the 10th
% International Conference on Speech Science & Technology, pp.420-425, Sydney.
b = 1.65;
a0 = 1e5;
phi = 0;

gt1 = zeros(len,length(CenterFreq));
g = zeros(1,length(CenterFreq));
gt_out = zeros(size(gt1));
for ii = 1:length(CenterFreq)
    % generate gammatone IR for each bark band
    gt1(:,ii) = a0*t.^(N-1).*exp(-2*pi*b*BandWidth(ii).*t).*cos(2*pi*CenterFreq(ii).*t + phi);
    
    % normalize each gammatone filter at center frequency to 0 dB gain
    g(1,ii) = sqrt((sum(gt1(:,ii).*cos((2*pi*CenterFreq(1,ii))/fs*(0:len-1)'),1)).^2 +...
        (sum(gt1(:,ii).*sin((2*pi*CenterFreq(1,ii))/fs*(0:len-1)'),1)).^2);
    
    gt_out(:,ii) = gt1(:,ii)./g(1,ii);
end

% fiter input signal with gammatone filter
sig_out = zeros(length(sig_in),length(CenterFreq),size(sig_in,2));
for ii = 1:size(sig_in,2)
    for jj = 1:length(CenterFreq)
        sig_out(:,jj,ii) = filter(gt_out(:,jj),1,sig_in(:,ii));
    end
end

