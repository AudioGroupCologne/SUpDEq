% [f_nm, f, is_even, YnmInv] = AKsht(x, doFFT, multi, N, SHTmode, fs, compact, SHmode, tikhEps)
%
% does a discrtete spherical harmonics/fourier transform (SHT) with
% sampling weights according to eq. (3.2) in [1] or without sampling
% weights according to (3.34) in [1]. The real valued transfrom is done
% according to eq. (233) in [2].
%
% See AKsphericalHarmonicsTransformDemo.m for examples
%
% I N P U T:
% x       - input signal. One column holds one channel (e.g. impulse 
%           response). Any kind of data sampled on a sphere of uniform
%           radius can be passed.
% doFFT   - if true (default), an FFT is applied to x before spherical
%           harmonics transform. Use this if you pass impulse responses in
%           x. After the FFT only the first half of the spectrum is used.
% multi   - a) the spatial sampling grid in form of a Qx2 matrix where the
%           first column holds the azimuth and the second the elevation. In
%           this case the SHT is done without sampling weights.
%           b) a Qx3 matrix where the third column holds the sampling
%           weights
%           c) the inverted Ynm matrix for calculating the SHT without
%           sampling weights
% N       - Order of the spherical harmonics transform. Only necessary of
%           passing a spatial sampling grid in multi.
% SHTmode - 'db_unwrap'  separate transfrom on logarithmic absolute and
%                        unwrapped phase data (default)
%           'abs_unwrap' separate transform on magnitude and unwrapped
%                        phase data
%           'complex'    transform on complex data
%           'real_imag'  separate transform on real and imaginary part
%           'db'         transform on logarithmic absolut input data
%           'abs'        transform on absolute input data
% fs      - sampling frequncy in Hz (default = 44100)
% compact - false: Two SH transforms are carried out if SHTmode is
%                  'db_unwrap', 'abs_unwrap', or 'real_imag'.
%           true: A single SH transform is carried out, e.g. on
%                 magnitude + 1j*unwrapPhase.
%           This leads to identical results, but the compact mode will only
%           return a single set of f_nm coefficients (default = false)
% SHmode  - 'complex' to use complex valued spherical harmonics (default), 
%           or 'real' to use real valued spherical harmonics which results
%           in less data. If compact = false it can be used with SHTmodes
%           db_unwrap, abs_unwrap, and real_imag. It can be used with
%           SHTmodes db , and abs regardless of the compact parameter.
% tikhEps - Define epsilon of Tikhonov regularization if regularization 
%           should be applied [3]. Only appicable to least-squares solution
%           (without weights)
%           Default: 0 (no Tikhonov regularization)
%
% O U T P U T:
% f_nm    - spherical harmonics coefficients. One column holds coefficients
%           for a single frequency specified in f, or for one row in x if
%           doFFT=false. The third dimension holds f_nm coefficients
%           according to SHTmode, e.g. if SHTmode='db_unwrap', f_nm(:,:,1)
%           will hold the coefficients for the log. absolute values and
%           f_nm(:,:,2) will hold the coefficients for the unwrapped phase.
%           If compact=true f_nm the third dimension will be singular.
% f       - frequencies corresponding to f_nm if doFFT=true.
% is_even - specifies if the input signal was of even length, if doFFT=true
% YnmInv  - matrix with inverted spherical harmonics coefficients. Can be
%           passed to AKsht to save proecessing time in case multiple
%           transformations are done 
%
% Note: If you apply the transformation on the phase values as well, the
%       impulse reponses should be free of a common delay.
%
% [1] Boaz Rafaely: Fundamentals of spherical array processing. In.
%     Springer topics in signal processing. Benesty, J.; Kellermann, W.
%     (Eds.), Springer, Heidelberg et al. (2015).
% [2] Franz Zotter: Analysis and synthesis of sound-radiation with 
%     spherical arrays. Ph.D. dissertation, University of Music and
%     Performing arts (2009).
% [3] Ramani Duraiswami, Dimitry N. Zotkin, and Nail A. Gumerov: "Inter-
%     polation and range extrapolation of HRTFs." IEEE Int. Conf. Acoustics
%     , Speech, and Signal Processing (ICASSP), Montreal, Canada, May 2004,
%     p. 45-48, doi: 10.1109/ICASSP.2004.1326759.

%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group, TU Berlin
% 09/2015 - inital dev. (fabian.brinkmann@tu-berlin.de)
% 07/2020 - added Tikhonov regularization (Johannes.Arend@th-koeln.de)

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
function [f_nm, f, is_even, YnmInv] = AKsht(x, doFFT, multi, N, SHTmode, fs, compact, SHmode, tikhEps)

% -------------------------------------------------- set default parameters
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('SHTmode', 'var')
    SHTmode = 'db_unwrap';
end
if ~exist('compact', 'var')
    compact = false;
end
if ~exist('SHmode', 'var')
    SHmode = 'complex';
end
if ~exist('tikhEps', 'var')
    tikhEps = 0;
end

% ------------------------------- transfer into frequency domain if desired
if doFFT
    f            = (0:fs/size(x,1):fs/2)';
    [x, is_even] = AKboth2singleSidedSpectrum(fft(x));
else
    is_even = [];
    f       = [];
end

M = size(x, 1);

% ------------------- get inverted Ynm matrix if not passed to the function
if size(multi, 2) <= 3
    if ~exist('N', 'var')
        error('You have to specify N if you pass a spatial sampling grid in ''multi''')
    end
    % get SH matrix
    [Ynm, n, m] = AKsh(N, [], multi(:,1), multi(:,2), SHmode);
    
    if size(multi, 2) == 2
        %Check if Tikhonov regularization should be applied and perform
        %least-squares solution with Tikhonov
        if tikhEps == 0
            % discrete SHT without sampling weights and without regularization: calculate inverse matrix
            YnmInv = pinv(Ynm);
        else
            %Create diagonal matrix according to [3]
            I = eye(size(n,2));
            D = diag(1 + n.*(n+1)) .* I;
            
            % Inverse SH matrix for Least-Squares SH transform with Tikhonov regularization
            YnmInv = (Ynm' * Ynm + tikhEps*D)^-1 * Ynm';
        end
    else
        % discrete SHT using sampling weights: calculate conjugate matrix
        YnmConj = zeros(size(Ynm));
        
        for nn = 0:N
            id = find(n == nn);
            % eq. (1.10) in [1]
            YnmConj(:, id) = AKm( (-1).^(m(id)), Ynm(:, fliplr(id)), '*');
        end
        
        YnmInv = 'YnmInv is not calculated, if sampling weights are passed';
    end
    
else
    % discrete SHT using matrix inversion - inverse matrix was passed
    YnmInv = multi;
end

% --------------------------------- pre process x to match desired SHT mode
switch SHTmode
    case 'db_unwrap'
        X        = zeros(size(x,1), size(x,2), 2);
        X(:,:,1) = db(abs(max(x, eps)));
        X(:,:,2) = unwrap(angle(x));
    case 'abs_unwrap'
        X        = zeros(size(x,1), size(x,2), 2);
        X(:,:,1) = abs(x);
        X(:,:,2) = unwrap(angle(x));
    case 'complex'
        X        = x;
    case 'real_imag'
        X        = zeros(size(x,1), size(x,2), 2);
        X(:,:,1) = real(x);
        X(:,:,2) = imag(x);
    case 'db'
        X        = db(abs(x));
    case 'abs'
        X        = abs(x);
    otherwise
        error(['SHTmode ''' SHTmode ''' unknown'])
end

if compact && size(X,3)==2
    X = X(:,:,1) + 1j*X(:,:,2);
end

if strcmpi(SHmode, 'real') &&  ~all( isreal( X(:) ) )
    error('AKsht:SHmode', 'Set SHmode must be ''complex'' if transform data is complex')
end

% ----------------------------------------------------- do the SH transform

if isnumeric(YnmInv)
    f_nm = zeros(size(YnmInv,1), size(X,1), size(X,3));
    
    % SHT using inverse SH matrix -> eq. (3.34) in [1]
    for m = 1:M
        for c = 1:size(X,3)
            f_nm(:,m,c) = YnmInv * X(m,:,c).';
        end
    end
else
    f_nm = zeros(size(Ynm,2), size(X,1), size(X,3));
    
    if strcmpi(SHmode, 'complex')
        % SHT using conjugate SH matrix -> eq. (3.2) in [1]
        X = AKm(multi(:,3)', X, '*');
        
        for c = 1:size(X, 3)
            for m = 1:size(X, 1)
                f_nm(:,m,c) = sum( repmat( X(m,:,c).', 1, (N+1)^2) .* YnmConj ).';
            end
        end
    else
        % SHT using SH matrix -> eq. (233) in [2]
        Ynm_tw = ( Ynm.' * diag(multi(:,3)) );
        
        for c = 1:size(X, 3)
            for m = 1:size(X, 1)
                f_nm(:,m,c) = Ynm_tw * X(m,:,c)';
            end
        end
    end

    f_nm = 4*pi*f_nm;
end
