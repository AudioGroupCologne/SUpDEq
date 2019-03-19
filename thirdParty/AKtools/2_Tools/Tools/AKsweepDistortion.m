% [groupDelay_k, f] = AKsweepDistortion(groupDelay, N, fs, K)
%
% numerically calculates, and plots the frequency dependend group delay
% of non-harmonic distortion products of impulse responses measured with
% swept sines. This is usefull, when non-harmonic distortion products
% should be analyzed.
%
% See AKsweepFD.m for example (parameter do_advanced_plot)
%
% I N P U T
% groupDelay   - group delay as passed from AKsweepSynth.m, i.e. single
%                sided spectrum.
% N            - length of original sweep in samples (If NFFT_double=1
%                during the call of AKsweepSynth, pass 2*N to this
%                function)
% fs           - sampling frequency in Hz
% K            - highest order of non-harmonic distortion products that is
%                analyzed, and plottet (default = 5)
%
% O U T P U T
% gourpDelay_k - frequency dependent gourp delay of non-harmonic distortion
%                products in seconds
% f            - frequencies of groupDelay_k
%
% v1 2013/02 fabian.brinkmann@tu-berlin.de, Audio Communication Group,
%            TU Berlin
%
% [1] Weinzierl et al. (2009): 'Generalized multiple sweep measurement', In
% 126th Audio. Eng. Soc Convention (Paper 7767). Munich, Germany.

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
function [groupDelay_k, f] = AKsweepDistortion(groupDelay, N, fs, K)

if ~exist('K', 'var')
    K = 5;
end

% discard 0 Hz
groupDelay = groupDelay(2:end);

% frequency vector
f  = (fs/N:fs/N:fs/2)';

% get group delay of impulse responses of non-harmonic products [1, eq. 2]
% (numerical calculation)
groupDelay_k = zeros(size(f,1), K-1);

for k = 2:K
    groupDelay_k(:,k-1) = interp1(f(:), groupDelay(:), f(:)/k, 'spline');
    groupDelay_k(:,k-1) = groupDelay_k(:,k-1) - groupDelay;
end

% make figure of none exists
if isempty(findall(0,'Type','Figure'))
    % set size and ranges also for plottinf
    figure('units','normalized','outerposition',[0 0 1 1], ...
           'PaperUnits', 'normalized', 'PaperSize', [1 1], 'PaperPosition', [0 0 1 1])
    % set background color to white
    set(gcf, 'color', [1 1 1])
end

% plot
semilogy(groupDelay_k, f, 'LineWidth', 2)
ylim([10 20000])
grid on

% legend, and label
leg = '';
for k = 2:K
    leg = [leg '''k' num2str(k) ''', ']; %#ok<AGROW>
end
leg = [leg ' ''location'', ''best'''];

eval(['legend(' leg ')'])
ylabel('f in Hz'); xlabel('grp. delay in s');
title('Group delay of non-harmonic distortion products vs frequency')

if nargout == 0
    clear groupDelay_k f
end
