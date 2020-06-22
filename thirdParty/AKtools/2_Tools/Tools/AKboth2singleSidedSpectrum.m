% [single_sided, is_even] = AKboth2singleSidedSpectrum(both_sided)
% cuts frequencies larger than half the sampling rate from a both-sided
% spectrum and returns only the single-sided spectrum
%
% can be used to switch back and forth between single and both sided
% spectra, e.g.
% y = AKboth2singleSidedSpectrum(x);
% x = AKsingle2bothSidedSpectrum(y);
%
% I N P U T:
% both-sided spectrum, of size [N M C], where N is the number of frequency
% bins, M the number of measurements and C the number of channels. N must
% correspond to frequencies 0 <= f <= fs, where fs is the sampling
% frequency.
%
% O U T P U T:
% single-sided spectrum.
% is_even is a boolean that specifies wether or not the both-sided input
% spectra had an even number of frequency bins. This is needed for
% recovering the both-sided spectrum -> see AKsingle2bothSidedSpectrum.m
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
function [single_sided, is_even] = AKboth2singleSidedSpectrum(both_sided)

NFFT    = size(both_sided,1);
is_even = ~mod(NFFT,2);

if is_even
    single_sided = both_sided(1:NFFT/2+1,:,:);
else
    single_sided = both_sided(1:ceil(NFFT/2),:,:);
end
