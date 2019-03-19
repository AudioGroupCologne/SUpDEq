% both_sided = AKsingle2bothSidedSpectrum(single_sided, is_even)
%
% can be used to switch back and forth between single and both sided
% spectra, e.g.
% y = AKboth2singleSidedSpectrum(x);
% x = AKsingle2bothSidedSpectrum(y);
%
% Note that only the real part of the frequency bins at 0 Hz and Nyquist
% (i.e. half the sampling rate) are considered before generating the both
% sided spectrum.
%
% I N P U T
% single-sided   - single sided spectrum , of size [N M C], where N is the
%                  number of frequency bins, M the number of measurements
%                  and C the number of channels.
%                  N must correspond to frequencies of 
%                  0 <= f <= fs/2, where f is the sampling frequency
% is_even        - true if both sided spectrum had even number of
%                  taps (default).
%                  if is_even>1, it denotes the number of samples of the
%                  both sided spectrum (default = 1)
%
% O  U T P U T
% both-sided spectrum.
%
% F. Brinkmann, Audio Communication Group, TU Berlin, 04/2013

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
function both_sided = AKsingle2bothSidedSpectrum(single_sided, is_even)

if ~exist('is_even', 'var')
    is_even = 1;
end

if is_even>1
    is_even = 1-mod(is_even,2);
end

N = size(single_sided,1);

if is_even
    % make sure that the bin at nyquist frequency is real
    % (there might be rounding errors that produce small immaginary parts)
    single_sided(end,:,:) = real(single_sided(end,:,:));
    % mirror the spectrum
    both_sided = [single_sided; flipud(conj(single_sided(2:N-1,:, :)))];    
else
    % mirror the spectrum
    both_sided = [single_sided; flipud(conj(single_sided(2:N,:, :)))];
end

% make sure that the bin at 0 Hz is real
% (there might be rounding errors that produce small immaginary parts)
both_sided(1,:,:) = real(both_sided(1,:,:));