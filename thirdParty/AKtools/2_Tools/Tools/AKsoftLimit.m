% H = AKsoftLimit(H, limit_dB, knee, fRange, fs, isEven)
%
% applies soft limiting/compression to the complex, single sided spectra H.
% Soft limiting can either be applied using an arcus tanges function
% according to [1]. equation (3.25), or using classic limiting/compression
% logic according to [2], equation (4). Both algorithms give similar
% results if passing knee = 30 to AKsoftLimit.
%
%
% Example:
% Hin  = logspace(0, 3, 100)';          -> generate a real input spectrum
% Hout = AKsoftLimit(Hin, 20, 'atan');  -> arcus tanges limiting at 20 dB
% Hout = AKsoftLimit(Hin, 20, 10);      -> soft knee limiting with 10 dB
%                                          knee width
% Note:
% See AKboth2singleSidedSpectrum.m, and AKsingle2bothSidedSpectrum.m for
% converting spectra.
%
%
% I N P U T
% H        - single sided real or complex spectra of size [N C M], where N
%            is the number of frequency bins, and C and M a positive integer
% limit_dB - the maximum gain in dB that is allowed in the spectrum.
% knee     - 'atan': to apply arcus tangens soft limiting according to [1],
%                    equation (3.25) - this is the functions default.
%             W    : Applies soft knee limiting with knee width W in dB
%                    according to [2], equation (4). To disable the soft
%                    knee, pass W=0.
%            [W R] : Applies soft knee compression with ratio R in dB
% fRange   - two element vector [f_low f_high ]that specifies a frequency
%            range in Hz where the limiting is applied. If fRange is not
%            passed limiting is applied to the entire spectrum by default
% fs       - sampling frequency in Hz. Only needed if fRange is passed.
%            Default = 44100.
% isEven   - true or false to denote whether or not the both sided spectrum
%            has an even number of bins. Only needed if fRange is passed.
%
%
% O U T P U T
% H        - complex, single sided, and limited spectrum
%
%
% R E F E R E C E S
% [1] B. Bernschuetz (2016). Microphone arrays and sound field decomposi-
%     tion for dynamic binaural synthesis. Ph.D. dissertation, TU Berlin,
%     Berlin, Germany.
% [2] D. Giannoulis, M. Massberg, and J.D. Reiss (2012). "Digital Dynamic
%     Range Compressor Design?A Tutorial and Analysis." J. Audio Eng. Soc,
%     60(6), 399--408. 
%
%
% 2017/06 - fabian.brinkmann@tu-berlin.de

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
function H = AKsoftLimit(H, limit_dB, knee, fRange, fs, isEven)

% ------------------------------------------------------ set default values
if ~exist('knee', 'var')
    knee = 'atan';
end

if ~exist('isEven', 'var')
    isEven = true;
end

if ~exist('fs', 'var')
    fs = 44100;
end

if ~exist('fRange', 'var')
    fRange = false;
end

% -------------------------------------------------- select frequency range
if any(fRange)
    % number of frequency bins from the complete spectrum
    if isEven
        N = 2* ( size(H, 1) - 1 );
    else
        N = 2* size(H, 1) - 1;
    end
    
    % frequency resolution
    df = fs/N;
    
    % frequency bins where soft limiting is applied
    fID    = [ round( min(fRange)/df + 1 ) round( max(fRange)/df + 1 )];
    fID(2) = min( fID(2),  size(H, 1) );
    
else
    fID = [1 size(H, 1)];
end

% corresponding frequency bins
Hf = H( fID(1):fID(2), : , :);

% -------------------------------------- pre compute frequently used values
% set the maximum gain
limit_lin = 10^(limit_dB/20);
% get absolute values within frequency range of interest
Hf_abs = abs(Hf);

% ----------------------------------------------------- apply soft limiting
if ischar(knee)
    % arcus tanges soft limitng according to [1]
    
    alpha = Hf_abs / limit_lin;
    H( fID(1):fID(2), : , :) = 2/pi * Hf./alpha .* atan(pi/2 * alpha);
else
    % soft knee limiting according to [2]
    if numel(knee) == 1
        knee = [knee inf];
    end
    
    % precompute values
    Hf_dB = db(Hf_abs);
    
    % select different ranges for soft knee limiting - [2], eq. (4)
    idSoft = 2 * abs(Hf_dB - limit_dB) <= knee(1);
    idHard = 2 *    (Hf_dB - limit_dB) >  knee(1);
    
    % calculate and apply the soft limit gains
    gainSoft = (1/knee(2)-1) * ( Hf_dB(idSoft)-limit_dB + knee(1)/2 ).^2 / (2*knee(1));
    gainSoft = 10.^(gainSoft/20);
    
    Hf(idSoft) = Hf(idSoft) .* gainSoft;
    
    % calculate and apply the compression/limiting
    if knee(2) ~= inf
        gainHard = limit_dB + (Hf_dB(idHard)-limit_dB)/knee(2);
        gainHard = 10.^(gainHard/20);
        
        Hf(idHard) = Hf(idHard)./Hf_abs(idHard) .* gainHard;
    else
        Hf(idHard) = Hf(idHard)./Hf_abs(idHard) * limit_lin;
    end
    
    % change output
    H( fID(1):fID(2), : , :) = Hf;
end
