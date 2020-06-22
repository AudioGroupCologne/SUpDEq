% [H, h] = AKinterpolateSpectrum(H, f, N, extrapMode, fs, phaseType, logInterp, doPlot)
% generates impulse response from a sparse magnitude spectrum by
% interpolation and inverse Fourier transform
%
% For example usage see AKism.m
%
% I N P U T
% H           - magnitude spectrum of size [F x W] where F is the
%               number of frequencies and W the number of walls
% f           - frequencies in Hz corresponding to values in H
% N           - desired impulse response length (default = 128)
% extrapMode  - Method for inter- and extrapolating the absorption
%               coefficients below, inside and above tha range given by f.
%               E.g. {'linear' 'linear', 'linear'} which is the default
%               performs linear interpolation in all ranges.
% fs          - sampling frequency in Hz (default = 44100)
% phaseType   - 'min' or 'lin' to obtain minimum or linear phase impuse
%               repsonses (default = 'min')
% logInterp   - interpolate on db(H) and log(f), default = false
% doPlot      - plot result (default = false)
%
% O U T P U T
% H           - single sided spectrum of size [N/2+1 x w]
% h           - impulse response(s) of size [N x W]. If h is returned, H
%               has a minimum length of 513 frequency bins.
%
% 2017/05 - fabian.brinkmann@tu-berlin.de

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
function [H, h] = AKinterpolateSpectrum(H, f, N, extrapMode, fs, phaseType, logInterp, doPlot)

% default values
if ~exist('N', 'var')
    N = 128;
end
if ~exist('extrap_mode', 'var')
    extrapMode = {'lin' 'lin' 'lin'};
end
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('phaseType', 'var')
    phaseType = 'min';
end
if ~exist('logInterp', 'var')
    logInterp = false;
end
if ~exist('doPlot', 'var')
    doPlot = false;
end

f = reshape(f, [numel(f) 1]);

% discard valus above fs/2
H = H(f<=fs/2,:);
f = f(f<=fs/2,:);

H_cp = H;

% number of channels
C = size(H, 2);

% design impulse responses with a minimum length of 1024 samples
if nargout == 2
    if mod(N, 2)
        M = max(N, 1025);
    else
        M = max(N, 1024);
    end
else
    M = N;
end

% get bin of lowest and highest specified T values
f_lim = [ceil(f(1)/fs*M+1) floor(f(end)/fs*M+1)];

% target frequencies
f_interp = (0:fs/M:fs/2)';

% allocate space for output
h = zeros(numel(f_interp), C);

% use the log spectrum
if logInterp
    H        = db(H);
    f        = log(f);
    f_interp = log(f_interp);
end

for cc = 1:C
    % interpolate H values
    h(f_lim(1):f_lim(2), cc) = interp1(log(f), H(:,cc), log(f_interp(f_lim(1):f_lim(2))), extrapMode{2});
    
    % extrapolate H values
    h(2:f_lim(1)-1, cc) = interp1(log(f), H(:,cc), log(f_interp(2:f_lim(1)-1)), extrapMode{1}, 'extrap');
    h(f_lim(2)+1:end, cc) = interp1(log(f), H(:,cc), log(f_interp(f_lim(2)+1:end)), extrapMode{3}, 'extrap');  
end


if logInterp
    % we can not interpolate to 0 Hz on a log scale
    h(1,:) = h(2,:);
    % de-log
    h = 10.^(h/20);
end

H = h;

if nargout == 2
    % get single sided spectra
    h = AKsingle2bothSidedSpectrum(H, 1-mod(M,2));
    
    % get zero phase IRs
    h = ifft(h, 'symmetric');
    
    % get minimum/linear phase IRs
    NFFTdouble = 1;
    [h, dev] = AKphaseManipulation(h, fs, phaseType, NFFTdouble, false);
    while dev(1)>1 && 2^NFFTdouble * M < 2^18
        NFFTdouble = NFFTdouble + 1;
        [h, dev] = AKphaseManipulation(h, fs, phaseType, NFFTdouble, false);
    end
    
    if N < M
        if strcmpi(phaseType, 'min')
            h = AKfade(h, N, 0, 10);
        else
            if mod(N, 2)
                h = h(ceil(M/2)-floor(N/2):ceil(M/2)+floor(N/2),:);
            else
                h = h(M/2-N/2:M/2+N/2-1,:);
            end
            
            h = AKfade(h, [], 5, 5);
        end
    end
    
end

% plot
if doPlot
    AKf
    if nargout == 2
        for cc = 1:C
            subplot(2, C, cc)
                AKp(h(:,cc), 'etc2d', 'xu', 'n', 'x', [-10 N])
                title(['Channel ' num2str(cc) ': IR'])
            subplot(2, C, cc+C)
                AKp(h(:,cc), 'm2d', 'N', fs/5)
                ylabel Magnitude
                hold on
                plot(f, db(H_cp(:,cc)), 'xr')
                title(['Channel ' num2str(cc) ': Spectrum'])
        end
    else
        semilogx(f_interp, db(H(:,cc)), 'k')
        hold on
        semilogx(f, db(H_cp(:,cc)), 'or')
        xlabel frequency
        ylabel Magnitude
        AKfrequencyTicks
        title(['Channel ' num2str(cc) ': Spectrum'])
    end
end