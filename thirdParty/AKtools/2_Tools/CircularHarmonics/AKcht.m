% [d_n, f, is_even] = AKcht(x, doFFT, phi, N, CHTmode, fs, compact, CHmode)
%
% does a discrtete circular harmonics/fourier transform (CHT) with
% sampling weights similar to eq. (3.2) in [2] if the angular values passed
% by phi are equally spaced or without sampling weights according to (14)
% in [1], which is similar to (3.34) in [1].
%
% See AKcircularHarmonicsDemo.m for examples
%
% I N P U T:
% x       - input signal. One column holds one channel (e.g. impulse 
%           response). Any kind of data sampled on a sphere of uniform
%           radius can be passed.
% doFFT   - if true (default), an FFT is applied to x before spherical
%           harmonics transform. Use this if you pass impulse responses in
%           x. After the FFT only the first half of the spectrum is used.
% phi     - the spatial sampling grid in form of a Qx1 vector holding
%           the angular values. Coordinate convention according to AKch.m
% N       - Order of the spherical harmonics transform.
% CHTmode - 'db_unwrap'  separate transfrom on logarithmic absolute and
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
% CHmode  - 'complex' to use complex valued circular harmonics (default), 
%           or 'real' to use real valued circular harmonics which results
%           in less data. If compact = false it can be used with CHTmodes
%           db_unwrap, abs_unwrap, and real_imag. It can be used with
%           SHTmodes db , and abs regardless of the compact parameter.
%
% O U T P U T:
% d_n     - circular harmonics coefficients. One column holds coefficients
%           for a single frequency specified in f, or for one row in x if
%           doFFT=false. The third dimension holds f_nm coefficients
%           according to SHTmode, e.g. if SHTmode='db_unwrap', f_nm(:,:,1)
%           will hold the coefficients for the log. absolute values and
%           f_nm(:,:,2) will hold the coefficients for the unwrapped phase.
%           If compact=true f_nm the third dimension will be singular.
% f       - frequencies corresponding to f_nm if doFFT=true.
% is_even - specifies if the input signal was of even length, if doFFT=true
%
% Note: If you apply the transformation on the phase values as well, the
%       impulse reponses should be free of a common delay.
%
% [1] Noam R Shabtai, Gottfried Behler, and Michael Vorlaender: "Generation
%     of a reference radiation pattern of string instruments using
%     automatic excitation and acoustic centering." J. Acoust. Soc. Am.,
%     138(5): EL479-EL456, (2015).
% [2] Boaz Rafaely: Fundamentals of spherical array processing. In.
%     Springer topics in signal processing. Benesty, J.; Kellermann, W.
%     (Eds.), Springer, Heidelberg et al. (2015).
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group, TU Berlin
% 02/2018

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
function [d_n, f, is_even] = AKcht(x, doFFT, phi, N, CHTmode, fs, compact, CHmode)

% -------------------------------------------------- set default parameters
narginchk(4,8)

if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('CHTmode', 'var')
    CHTmode = 'db_unwrap';
end
if ~exist('compact', 'var')
    compact = false;
end
if ~exist('CHmode', 'var')
    CHmode = 'complex';
end

% check format of angular values
phi = reshape(phi, [numel(phi), 1]);
% discard doubel entries
phi = unique( mod(phi, 360) );

% ------------------------------- transfer into frequency domain if desired
if doFFT
    f            = (0:fs/size(x,1):fs/2)';
    [x, is_even] = AKboth2singleSidedSpectrum(fft(x));
else
    is_even = [];
    f       = [];
end

M = size(x, 1);

% --------------------- check if we have evenly distributed  angular values
ang_diff = diff([phi; phi(1)+360]);
ang_diff = ang_diff - ang_diff(1);

if all(abs(ang_diff) < 1e-3)
    regular = true;
    weight  = 1 / numel(phi);
else
    regular = false;
end

clear ang_diff

% --------------------------------- get e_n matrix (inverted or conjugated)
% get CH matrix
e_n = AKch(N, phi, CHmode);
if regular
    % discrete CHT using sampling weights: calculate conjugate matrix
    e_nConj = conj(e_n);
else
    % discrete CHT without sampling weights: calculate inverse matrix
    e_nInv = pinv(e_n);
end

% --------------------------------- pre process x to match desired CHT mode
switch CHTmode
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
        error(['CHTmode ''' CHTmode ''' unknown'])
end

if compact && size(X,3)==2
    X = X(:,:,1) + 1j*X(:,:,2);
end

if strcmpi(CHmode, 'real') &&  ~all( isreal( X(:) ) )
    error('AKcht:CHmode', 'Set CHmode must be ''complex'' if transform data is complex')
end
% ----------------------------------------------------- do the CH transform
if regular
    d_n = zeros(size(e_nConj,2), size(X,1), size(X,3));
    
    % CHT using conjugate CH matrix -> similar to eq. (3.2) in [1]
    X =  weight * X;
    
    for c = 1:size(X, 3)
        for m = 1:size(X, 1)
            d_n(:,m,c) = X(m,:,c) * e_nConj;
        end
    end
    
    % normalize coefficients - this is a bit weird and I did not find the
    % reason why it has to be done in the time that I took.
    d_n = d_n * 8 * pi^2; 
    if strcmpi(CHmode, 'complex')
        d_n = d_n / 2;
    else
        d_n( ceil(numel(d_n)/2) ) = d_n( ceil(numel(d_n)/2) ) / 2;
    end

else
    d_n = zeros(size(e_nInv,1), size(X,1), size(X,3));
    
    % CHT using inverse CH matrix -> similar to eq. (3.34) in [1]
    for m = 1:M
        for c = 1:size(X,3)
            d_n(:,m,c) = e_nInv * X(m,:,c).';
        end
    end
end
