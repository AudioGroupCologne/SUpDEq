% [HpIRs_target, HpIRs_icd, f_icd] = AKpHeadphone(HpIRs, target, targetName, fLim, fs)
% 
% plots headphone impulse responses (HpIRs)
% a) averaged left and right ear HpIRs
% b) HpIRs relative to a target response (e.g. diffuse field HRIR)
% c) Inter channel differences in auditory filters
%
% See AKplotDemo.m for examples
%
% I N P U T:
% HpIRs    - SOFA file ('SimpleHeadphoneIR' convention), or matrix of
%            size [N x M x R], where N: number of samples per HpIR.
%            M: number of HpIRs per ear, R: number of ears, must be 2
% target   - target response of size [N x 1]. N can be different from the N
%            in the HpIRs. If target is not passed, the mean HpIR is used
%            for this.
% fLim, fs - see AKerbError.m for documentation and default values
%
% O U T P U T:
% HpIRs_target - averaged HpIRs relative to the target
% HpIRs_icd    - Inter channel differences in auditory filters at
%                frequencies f_icd
% f_icd        - center frequencies of auditory filters [Hz]
%
%
% 9/2016 v.1 fabian.brinkmann@tu-berlin.de

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License.

function [HpIRs_target, HpIRs_icd, f_icd] = AKpHeadphone(HpIRs, target, targetName, fLim, fs)

if nargin < 3
    targetName = 'target';
end
if nargin < 4
    fLim = [50 20000];
end
if nargin < 5
    fs = 44100;
end

% check if we have SOFA input
if isstruct(HpIRs)
    fs    = HpIRs.Data.SamplingRate;
    HpIRs = shiftdim(HpIRs.Data.IR, 2);
end

% check size of input IRs and zero pad if neccessary
if size(HpIRs,1)>size(target,1)
    target(end+1:size(HpIRs,1),:) = 0;
else
    HpIRs(end+1:size(target,1),:) = 0;
end

% average HpIRs
if size(HpIRs, 2) > 1
    HpIRs_avg = ifft(mean(fft(HpIRs), 2), 'symmetric');
else
    HpIRs_avg = HpIRs;
end

HpIRs_avg = squeeze(HpIRs_avg);

% check if we have a target function
if ~exist('target', 'var')
    target = ifft(mean(fft(HpIRs)), 'symmetric');
end

% apply target
HpIRs_target = AKm(HpIRs_avg, target, '/', 'ms');

% average deviation should be 0 dB
tmp = mean(mean(AKerbError(HpIRs_target,AKdirac(size(HpIRs_target,1)), [max(fLim(1), 50) min(fLim(2), 16000)], fs)));
HpIRs_target = HpIRs_target * 10^(-tmp/20);

% get inter-channel differences in auditory filters
[HpIRs_icd,f_icd] = AKerbError(HpIRs_avg(:,2),HpIRs_avg(:,1), fLim, fs);

% plot results
% make figure if none exists
if isempty(findall(0,'Type','Figure'))
    AKf(20,30)
end
subplot(3,1,1)
    AKp([target HpIRs_avg], 'm2d', 'c', [0 0 0; .2 .2 .8; .8 .2 .2], 'N', fs/20, 'dr', 60, 'lw', 1.2, 'x', [20 20000])
    legend('target', 'left', 'right', 'location', 'SouthWest')
    title 'Averaged HpIRs'
    AKgrid
subplot(3,1,2)
    AKp(HpIRs_target, 'm2d', 'c', [.2 .2 .8; .8 .2 .2], 'N', fs/20, 'lw', 1.2, 'x', [20 20000])
    title(['HpIRs relative to ' targetName])
    
    l = db(fft(HpIRs_target));
    l = max(max( l( round(fLim(1)/fs*size(HpIRs,1))+1:round(fLim(2)/fs*size(HpIRs,1)),:)));
    l = l + 10-mod(l,10);
    
    ylim([l-40 l]);
    
    AKgrid
subplot(3,1,3)
    
    set(gca, 'xscale', 'log', 'xlim', [20 20000], 'xTick', [10:10:100 200:100:1000 2000:1000:10000 20000], 'xTickLabel', {.01 '' '' '' '' '' '' '' '' .1 '' '' '' '' '' '' '' '' 1 '' '' '' '' '' '' '' '' 10 20})
    if max(abs(HpIRs_icd)) <= 5
        set(gca, 'yTick', -5:5, 'ylim', [-5 5])
    elseif max(abs(HpIRs_icd)) <= 10
        set(gca, 'yTick', -10:2:10, 'ylim', [-10 10])
    else
        set(gca, 'yTick', -20:5:20, 'ylim', [-20 20])
    end
    
    title 'Inter-channel differences in auditory filters'
    xlabel 'f in kHz'
    ylabel 'magnitude in dB'
    
    AKgrid
    
    plot(f_icd, HpIRs_icd, '.-k', 'MarkerSize', 10)
    
if ~nargout
    clear HpIRs_target HpIRs_icd f_icd
end
