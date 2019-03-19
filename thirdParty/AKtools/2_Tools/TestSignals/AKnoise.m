% noise = AKnoise(samples, channels, noise_color)
% generates a noise signal
%
% See AKtestSignalsDemo.m for examples
%
% I N P U T:
% samples     - duration of noise signal in samlpes (default = 44100)
% channels    - number of channels (default = 1)
% noise_color - noise color 'white', or 'pink' (default = pink)
%
% O U T P U T:
% noise       - noise signal of size [samples x channels]
%
% v1 12/2015 fabian.brinkmann@tu-berlin.de, Audio Communication Group,
%            TU Berlin

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
function noise = AKnoise(samples, channels, noise_color)

% set default values
if ~exist('samples', 'var')
    samples = 44100;
end
if ~exist('channels', 'var')
    channels = 1;
end
if ~exist('noise_color', 'var')
    noise_color = 'pink';
end

% generate noise
switch lower(noise_color)
    case 'white'
        noise = randn(samples,channels) - .5;
        noise = noise / max(abs(noise(:)));
    case 'pink'
        noise = pinknoise(samples*channels)';
        noise = reshape(noise, [samples channels]);
end
    
   

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Pink Noise Generation with MATLAB Implementation   %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov       07/30/13  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function: y = pinknoise(N) 
% N - number of samples to be returned in row vector
% y - row vector of pink (flicker) noise samples
%
% The function generates a sequence of pink (flicker) noise samples. 
% Pink noise has equal energy in all octaves (or similar log bundles) of frequency.
% In terms of power at a constant bandwidth, pink noise falls off at 3 dB per octave.

% Copyright (c) 2013, Hristo Zhivomirov
% Copyright (c) 2012, Przemyslaw Baranski
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
function y = pinknoise(N)

% difine the length of the vector
% ensure that the M is even
if rem(N,2)
    M = N+1;
else
    M = N;
end

% generate white noise with sigma = 1, mu = 0
x = randn(1, M);

% FFT
X = fft(x);

% prepare a vector for 1/f multiplication
NumUniquePts = M/2 + 1;
n = 1:NumUniquePts;
n = sqrt(n);

% multiplicate the left half of the spectrum so the power spectral density
% is inversely proportional to the frequency by factor 1/f, i.e. the
% amplitudes are inversely proportional to 1/sqrt(f)
X(1:NumUniquePts) = X(1:NumUniquePts)./n;

% prepare a right half of the spectrum - a copy of the left one,
% except the DC component and Nyquist frequency - they are unique
X(NumUniquePts+1:M) = real(X(M/2:-1:2)) -1i*imag(X(M/2:-1:2));

% IFFT
y = ifft(X);

% prepare output vector y
y = real(y(1, 1:N));

% normalise
y = y./max(abs(y));

end
