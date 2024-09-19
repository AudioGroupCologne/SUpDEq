function y = AKbinauralCoherence(x, fs, N)
% Applys binaural coherence to a two channel input signal based on [1]
%
% For example usage see
% AKroomSimulationDemo.m
% AKdiffuseReverbTailDemo.m
%
% I N P U T
% x  - two channel input time signal of size [M x 2] where M is the number
%      of samples
% fs - sampling frequency in Hz (default = 44100)
% N  - filter length in samples (default = 2^12)
%
% O U T P U T
% y  - filtered input signal, also in the time domain
%
%
% [1] Christian Borss and Rainer Martin (2009): "An Improved Parametric
%     Model for Perception-Based Design of Virtual Acoustics." In: Proc.
%     of the 35th International AES Conference: Audio for Games. London
%
% 2014/11 fabian.brinkmann@tu-berlin.de

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 

if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('N', 'var')
    N = 2^12;
end

% in case of non-stationary filtering for generation of the reverb, we want
% to avoid having a bin at fs/2 in our fft. This bin had to be real and
% would commonly be set to 0 (-inf dB) and thus we would have a dip in the
% spectrum...
if ~mod(N, 2)
    N = N+1;
end

% frequency vector
f = (0:fs/N:fs/2)';

% interaural coherence from Borss2009 eq.(8)
gamma    = sin(pi*f/550)./(pi*f/550) .* max([zeros(size(f)) 1-f/2700], [], 2);
gamma(1) = 1;

% filter for coherence adjustment from Borss2009 eq.(17)(18)
H_b          = sqrt(1/2 * (1-sqrt(1-gamma.^2)));
H_b(gamma<0) = -H_b(gamma<0);

H_a = sqrt(1-H_b.^2);

% get double sided spectra and impulse responses
H_b = AKsingle2bothSidedSpectrum(H_b, 1-mod(N,2));
H_a = AKsingle2bothSidedSpectrum(H_a, 1-mod(N,2));

h_b = ifft(H_b, 'symmetric');
h_a = ifft(H_a, 'symmetric');

h_b = circshift(h_b, [round(N/2), 0]);
h_a = circshift(h_a, [round(N/2), 0]);

% plot
% subplot(2,1,1)
% hp([h_b h_a], 'ir2d', 'c', 'cyc')
% subplot(2,1,2)
% hp([h_b h_a], 's2d', 'c', 'cyc')

% filter
y(:,1) = fftfilt(h_a, [x(:,1); zeros(N-1,1)]) + fftfilt(h_b, [x(:,2); zeros(N-1,1)]);
y(:,2) = fftfilt(h_a, [x(:,2); zeros(N-1,1)]) + fftfilt(h_b, [x(:,1); zeros(N-1,1)]);

% circshift to account for the delay of the filter
y = circshift(y, [-round(N/2) 0]);
y = y(1:end-N+1,:);

