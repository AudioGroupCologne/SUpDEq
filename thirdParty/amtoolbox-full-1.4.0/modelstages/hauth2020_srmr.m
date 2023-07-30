function [ratio, energy] = hauth2020_srmr(s,varargin)
%HAUTH2020_SRMR Computes the speech-to-reverberation modulation energy ratio of the given signal
%
%   Usage:
%     [ratio, energy] = hauth2020_srmr(s, fs, 'fast', 0, 'norm', 0, 'minCF', 4, 'maxCF', 128)
%
%   Input parameters: 
%     s:            either the path to a WAV file or an array containing a single-channel speech sentence.
%     fs:           sampling rate of the data in s. If s is the path to a WAV file,  
%                   this parameter has to be omitted.
%
%   Output parameters:
%     ratio:        the SRMR score.
%     energy:       a 3D matrix with the per-frame modulation spectrum extracted from the input.
%
%
%
%   HAUTH2020_SRMR calculates the speech-to-reverberation modulation energy ratio using
%   the modulation filterbank described by Ewert and Dau in "Characterizing frequency 
%   selectivity for envelope fluctuations" (2000).
%
%   The code has been derived from the SRMR toolbox where it has been published under the MIT license.
%
%
%   Optional parameters:
%
%     'fast',F    flag to activate (F = 1)/deactivate (F = 0) the fast implementation. 
%                 The default is 'fast', 0 (this can be omitted).
%
%     'norm',N    flag to activate (N = 1)/deactivate (N = 0) the normalization step in the 
%                 modulation spectrum representation, used for variability reduction. The default is 'norm', 0.
%
%     'minCF',cf1   value of the center frequency of the first filter in the modulation filterbank.
%                   The default value is 4 Hz.
%
%     'maxCF',cf8   value of the center frequency of the first filter in the modulation filterbank. 
%                   The default value is 128 Hz if the normalization is off and 30 Hz if normalization is on.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/hauth2020_srmr.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Author: Christopher F. Hauth (2020)
%   #Author: Dr. Thomas Brand (2020)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% Load (if needed) and preprocess file
if ischar(s)
    fs = [];
    args = {varargin{:}};
else
    if isempty(varargin)
        error('Second argument must be the sampling rate if input is a vector');
    else
        if ~isnumeric(varargin{1})
            error('Second argument must be the sampling rate if input is a vector');
        else
            fs = varargin{1};
        end
    end
    args = {varargin{2:end}};
end

fast = 0;
norm = 0;

% Parameter parsing
for i = 1 : 2 : length(args)
    name = args{i};
    value = args{i+1};
    switch name
        case 'fast'
            fast = value;
        case 'norm'
            norm = value;
        case 'minCF'
            minCF = value;
        case 'maxCF'
            maxCF = value;
        case 'single'
            single = value;
        case 'window'
            window = value;
        otherwise
            error('Wrong parameter in parameter list');
    end
end

%% Fixed parameters:
% Modulation filterbank
nModFilters = 8;
if single
    nModFilters = 1;
end

if ~exist('minCF','var')
    minCF = 4;
end
if ~exist('maxCF', 'var')
    if norm == 1
        maxCF = 30;
    else
        maxCF = 128;
    end
end

%wLengthS = 0.256; % Window length in seconds. 
wLengthS = 0.250;
wLengthS = window;

%wIncS = 0.064; % Window increment in seconds;
wIncS = window;
wIncS = 0;
% Sampling rate for the modulation spectrum representation
if fast
    mfs = 400;
else
    mfs = fs;
end

wLength = ceil(wLengthS*mfs); % Window length in samples
wInc = ceil(wIncS*mfs); % window increment in samples

%% Cochlear filterbank/envelope computation

if fast
    % Compute acoustic band envelopes using gammatonegram
    %[tempEnv, ~] = gammatonegram(s, fs, 0.010, 0.0025, nCochlearFilters, lowFreq, fs/2);
   % tempEnv = flipud(tempEnv); % make the frequencies have the same order as the time-domain implementation
else
    % Pass the signal through the cochlear filter bank to produce the cochlear
    % outputs at each critical band.
    % The dimensions will be: cochlearOutputs(nCochlearFilters, nSamples), with 
    % the last being lowest frequency, and the first being highest frequency. 
%    cochlearOutputs = cochlearFilterBank(fs, nCochlearFilters, lowFreq, s);

    % Compute the temporal envelope for each critical band.  The dimensions will
    % be: temporalEnvelopes(nCochlearFilters, nSamples)
   % tempEnv = abs(hilbert(cochlearOutputs'))';    
end
tempEnv = abs(hilbert(s));    
tempEn = rms(tempEnv);
%tempEnv = tempEnv/tempEn;
%% Modulation spectrum
modFilterCFs = local_computemodulationcfs(minCF, maxCF, nModFilters);  
if single
    modFilterCFs = local_computemodulationcfs(minCF, maxCF, 1);  
end
%w = hamming(wLength);

    modulationOutput = local_modulationfilterbank(tempEnv(:,:), modFilterCFs, mfs, 2);%(tempenv(k,:))
    for m=1:nModFilters
        % Window frames with Hamming window
       % modOutFrame = buffer(modulationOutput(m,:), wLength,% (wLength-wInc));
        modOutFrame = modulationOutput(m,:);
        energy(1,m,:) = sum(modOutFrame.^ 2);
    end


%% Modulation energy thresholding
if norm
    peak_energy = max(max(mean(energy)));
    min_energy = peak_energy*0.00001;
    energy(energy < min_energy) = min_energy;
    energy(energy > peak_energy) = peak_energy;
end

%% Computation of K*

avg_energy = mean(energy,3);

total_energy = sum(sum(avg_energy));

AC_energy = sum(avg_energy,2);
%AC_perc = AC_energy*100./total_energy;

%AC_perc_cumsum=cumsum(flipud(AC_perc));
%K90perc = find(AC_perc_cumsum>90);

%BW = cochFilt_BW(K90perc(1));

%cutoffs = calc_cutoffs(modFilterCFs, fs, 2);

%if (BW > cutoffs(5)) && (BW < cutoffs(6))
 %   Kstar=5;
%elseif (BW > cutoffs(6)) && (BW < cutoffs(7))
 %   Kstar=6;
%elseif (BW > cutoffs(7)) && (BW < cutoffs(8))
%   Kstar=7;
%elseif (BW > cutoffs(8))
 %   Kstar=8;
%end

%% Modulation energy ratio
if single
    ratio = sum(avg_energy);
else
    ratio = sum(sum(avg_energy(:,1:4)))/sum(sum(avg_energy(:,5:8)));
    %ratio = var(sum(modulationOutput));
    %ratio = ratio./tempEn;
end
end



function out = local_modulationfilterbank(x, mcf, fs, q)
% y = modulationFilterBank(x, mcf, fs, q)
% -------------------------------------------------------------------------
% in: The input signal to be filtered, expected to be in row-vector form.
% mcf: A vector containing the center frequencies, in Hz, of each filter in
%   the filterbank. (Example: [2, 4, 8, 16, 32, 64, 128, 256])
% fs: The sampling rate, in Hz.
% q: The desired Q-value of the filters.
%
% out: The outputs of the modulation filterbank, organized as a matrix of
%   dimensions (length(mcf), length(x)).
% -----------------------------------------------
% Filters the input signal through the modulation filterbank described by
% Ewert and Dau in "Characterizing frequency selectivity for envelope
% fluctuations" (2000).  The original comment blocks are included. 

out = zeros(length(mcf),length(x));

B = zeros(length(mcf),3);
A = zeros(length(mcf),3);

for i = 1:length(mcf)

    w0 = 2*pi*mcf(i)/fs;
    [b3,a3] = local_makemodulationfilter(w0,q);
    B(i,:) = b3; A(i,:) = a3;
    out(i,:)= filter(b3, a3, x);   
end

end

function cfs = local_computemodulationcfs(minCF, maxCF, nModFilters)
% cfs = computeModulationCFs(minCF, maxCF, nModFilters)
% Computes the center frequencies of the filters needed for the modulation
% filterbank used on the temporal envelope (or modulation spectrum) of the
% cochlear channels.
% -----------------------------------------------
% minCF: Center frequency of the first modulation filter
% maxCF: Center frequency of the last modulation filter
% nModFilters: Number of modulation filters between minCF and maxCF
%
% cfs: The center frequencies of the filters needed for the modulation
% filterbank.

% Spacing factor between filters.  Assumes constant (logarithmic) spacing.
spacingFactor = (maxCF/minCF)^(1/(nModFilters-1));

% Computes the center frequencies
cfs = zeros(nModFilters, 1);
cfs(1) = minCF;
for i = 2:nModFilters
  cfs(i) = cfs(i - 1)*spacingFactor;
end

end

% [b, a] = makeModulationFilter(w0, q)
% -------------------------------------------------------------------------
% w0: normalized center frequency of the 2nd order bandpass filter
% q: The desired Q-value
%
% b, a: filter coefficients

function [b,a] = local_makemodulationfilter(w0,Q)
% w0 is cf of 2nd-order bandpass filter
% Q is the Q of the filter
W0 = tan(w0/2);
B0 = W0/Q;
 
b = [B0; 0; -B0];
a = [1 + B0 + W0^2; 2*W0^2 - 2; 1 - B0 + W0^2];

b = b/a(1);
a = a/a(1);

end


