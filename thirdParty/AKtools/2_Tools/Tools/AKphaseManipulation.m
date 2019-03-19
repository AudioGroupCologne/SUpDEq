% [y, dev_max] = AKphaseManipulation(x, fs, type, type_param, verbose)
%
% manipulates the phase of input signal x, to be zero-, linear- or
% minimum-phase.
% 
% See AKphaseManipulationDemo.m for examples
%
% INPUT (default):
% x           - input time signal [rows=samples x columns=channels]
% fs          - sampling frequency in Hz
% type        - desired phase behaviour 'lin', 'min', 'zero' ('min_phase etc. does work as well')
% type_param  - type='min': artifacts from the hilbert transform produce
%                           errors that usually get smaller if the input is
%                           zero padded. Thus, the length of x is changed
%                           by zero padding according to type_param.
%                           1=doubled, 2=quadroupled, etc.
%                           (Default = 1)
%               type='lin': sets the group dealy samples. This affects the
%                           point where the maximum energy in y occurs.
%                           A 0 sample delay will result in a zero phase
%                           signal that has the maximum at the irst sample.
%                           The default of (size(x,1)-1)/2 ensures that the
%                           max. is in the center of the output signal y.
%                           NOTE that parts of the energy will show up at
%                           the end of y if the group delay is too small.
%               type='zero': in this case type_param has no effect on y.
% verbose (1) - Differences between original and min phase magnitude
%               spectra is printed to screen. It is given separately for
%               frequencies that are 0-60dB, 60-100dB, and 100-inf dB below
%               the maximum. Set 1 or 0.
%
%
% OUTPUT:
% y           - impulse responses showing the desired phase behaviour
% dev_max     - maximum spectral deviations between original and min-phase
%               magnitude response in dB. Values are given for three ranges
%               of the spectrum. First value is the maximal deviation in
%               the upper 60 dB of the magnitude response, the second value
%               in the -60 to -100 dB range, and the third value covers
%               everything below -100 dB
%
% fabian.brinmann@tu-berlin.de, alexander.lindau@tu-berlin.de
% Audio Communication Group, TU Berlin, 09/2013

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
function [y, dev_max] = AKphaseManipulation(x, fs, type, type_param, verbose)

% auxiallary variables
type = lower(type);
N    = size(x,1);

% set the default parameters
id = find(type == '_');
if id
    type = type(1:id-1);
end
if ~exist('type_param', 'var')
    if strcmpi(type, 'min')
        type_param = 1;
    else
        type_param = (N-1)/2;
    end
end
if ~exist('verbose', 'var')
    verbose = 1;
end
if ~isreal(x)
    error('AKphaseManipulation:input', 'Input signal x must be real')
elseif isnan(x)
    error('AKphaseManipulation:input', 'Input signal x must not contain any NaN values')
end

% phase manipulation
switch type
    case 'zero'
        % results in zerophase impulse response
        % [Moeser (2008):"Theroetische Akustik." Skript zur Vorlesung S.25]
        y = ifft(abs(fft(x)), 'symmetric');
        
        dev_max = zeros(3,1);
    case 'lin'
        % desired group delay in seconds
        group_delay = type_param / fs;
        % spacing between frequency bins
        delta_omega = 2*pi* fs/N;
        % gradient of phase response
        delta_phi   = -group_delay * delta_omega;
        
        % get magnitude spectrum
        Y_abs = abs(fft(x));
        % generate phase spectrum
        Y_ang = linspace(0, (N-1)*delta_phi, N)';
        Y_ang = repmat(Y_ang, 1, size(Y_abs,2));
        
        % construct complex spectrum
        Y = Y_abs .* exp(1j*Y_ang);
        
        % obtain impulse response (matlab takes care of symmetry)
        y = ifft(Y, 'symmetric');
        
        dev_max = zeros(3,1);
    case 'min'
        % set the FFT order
        Nfft = nextpow2(N)+type_param;
        
        % pre-processing for removal of hilbert-artifacts: windowing, zeropadding
        y = ifft(abs(fft(x)));
        y = circshift(y, floor(N/2));
        
        % symmetrical zeropadding to avoid artifact from hilbert transform
        y = [zeros(floor(2^Nfft-N/2),size(y,2)); y; zeros(ceil(2^Nfft-N/2),size(y,2))];
        
        % create Cepstral/Hilbert minmumphase
        % use the AKrcpes function that can handle zeros in the spectrum
        % and multi-channel data
        y = AKrceps(y);
        % previous version:
        %         for n = 1:size(x, 2)
        %             [~, y(:, n)] = rceps(y(:, n));
        %         end
        
        % cut out "Hilbert-echo" .. dabei geht Sperrdaempfung in Arsch!
        y = y(1:N, :);
        
        
        % get difference between minimun-phase and original-phase
        % magnitude responses
        X = 20*log10(abs(fft(x)));
        deviation = 20*log10(abs(fft(y))) - X;
        % set deviation = 0, for bins where the magnitude response is more
        % than x dB below the maximum
        dev_max = zeros(2, size(x,2));
        for m = 1:2
            dB_val = [60 100];
            for n = 1:size(x,2)
                tmp    = abs(deviation(:,n));
                dB_max = max(X(:,n));
                dB_id  = X(:,n) > dB_max-dB_val(m);
                dev_max(m,n) = max(tmp(dB_id));
            end
        end
        dev_max(m+1,:) = max(abs(deviation));
        dev_max        = max(dev_max,[],2);
        
        if verbose
            disp('Difference between original and minimum-phase magnitude response:')
            disp([ num2str(round(dev_max(1)*100)/100) ' dB (0-60 dB dynamic)'])
            disp([ num2str(round(dev_max(2)*100)/100) ' dB (60-100 dB dynamic)'])
            disp([ num2str(round(dev_max(3)*100)/100) ' dB (100-inf dB dynamic)'])
        end
        
        if nargout == 1
            clear dev_max
        end
        
        
    otherwise
        error(['phase_type ' type ' not known. Use ''zero'', ''min'' or ''lin''']);
end

end


% ----------------------------------------------------------------------- %
% implementation of rceps that copes with zeros in the spectrum that might
% occur at Nyquist, and multi-channel data
function xMin = AKrceps(x)

N             = size(x,1);
Xabs          = abs(fft(x));
Xabs(Xabs==0) = eps;
xLog          = real(ifft(log(Xabs)));
isOdd         = mod(N,2);
xWindow       = [1; 2*ones((N+isOdd)/2-1,1) ; ones(1-rem(N,2),1); zeros((N+isOdd)/2-1,1)];
xWindow       = AKm(xWindow, xLog, '*');
xMin          = real(ifft(exp(fft(xWindow))));

end
