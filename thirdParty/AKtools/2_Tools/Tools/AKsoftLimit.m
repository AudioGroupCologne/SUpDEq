% H = AKsoftLimit(H, limit_dB, knee, fRange, fs, isEven, target, lim_dir)
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
%            is the number of frequency bins, and C and M a positive
%            integers
% limit_dB - the maximum gain in dB that is allowed in the spectrum of size
%            [N C M] where, N, C and M are 1 or equal the size of H.
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
% target   - reference of size [N C M] that is applied to H before soft
%            limiting (H./target). N, C, and M must be 1 or match the size
%            of H. The default target=false leafes H as it is.
% lim_dir  - 'upper': limits H in case db(H./target)>limit_dB
%            'lower': limits H in case db(1./(H./target))>limit_dB
%            'both' : applies both.
%            The default is 'upper'.
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
% 2017/06 - initial development (fabian.brinkmann@tu-berlin.de)
% 2020/07 - added 'target' and 'lim_dir' and possibility to pass limit_dB
%           as an array (fabian.brinkmann@tu-berlin.de)

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
function H = AKsoftLimit(H, limit_dB, knee, fRange, fs, isEven, target, lim_dir)

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

if ~exist('target', 'var')
    target = false;
end

if ~exist('lim_dir', 'var')
    lim_dir = 'upper';
end

% prepare limit_direction for looping
if strcmpi(lim_dir, 'both')
    lim_dir = {'upper' 'lower'};
else
    lim_dir = {lim_dir};
end

% ----------------------------------------------------------- prepare input
% make sure target is real
if ~islogical(target)
    target = abs(target);
end

%ceck dimensions of limit_dB and target
if ndims(H) < ndims(limit_dB)
    error('limit_dB can not have more dimensions than H.')
end
if ndims(H) < ndims(target)
    error('target can not have more dimensions than H.')
end

% match dimensions of H and limit_dB
reps = nan(1,ndims(H));
for dd = 1:ndims(H)
    if size(limit_dB,dd)==1 && size(H,dd)>1
        reps(dd) = size(H,dd);
    elseif size(limit_dB,dd) == size(H,dd)
        reps(dd) = 1;
    else
        error('The number of elements limit_dB has in each dimension must be 1 or equal H.')
    end
end
limit_dB = repmat(limit_dB, reps);


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
    
    % turn into mask
    f = false(size(H,1),1);
    f(fID(1):fID(2)) = true;
    % match size of mask
    Hsize = size(H);
    f = repmat(f, [1 Hsize(2:end)]);
    
else
    f = true(size(H));
end

% ----------------------------------------------------------- soft limiting
for ll = 1:numel(lim_dir)
    
    % invert data if lower soft limit is applied ------
    if strcmpi(lim_dir{ll}, 'lower')
        H(H==0)           = eps;
        H                 = 1./H;
        target(target==0) = eps;
        target            = 1./target;
    end
    
    % apply the target ------
    if ~islogical(target)
        target(target==0) = eps;
        H = H ./ target;
    end
       
    % pre compute frequently used values ------
    % set the maximum gain
    limit_lin = 10.^(limit_dB/20);
    % get absolute values within frequency range of interest
    H_abs = abs(H);
    
    % apply soft limiting ------
    if ischar(knee)
        % arcus tanges soft limiting according to [1]
        alpha = H_abs ./ limit_lin;
        H(f)  = 2/pi * H(f)./alpha(f) .* atan(pi/2 * alpha(f));
    else
        % soft knee limiting according to [2]
        if numel(knee) == 1
            knee(2) = inf;
        end
        
        % precompute values
        H_dB = 20*log10(H_abs);
        
        % select different ranges for soft knee limiting - [2], eq. (4)
        idSoft = 2 * abs(H_dB - limit_dB) <= knee(1) & f;
        idHard = 2 *    (H_dB - limit_dB) >  knee(1) & f;
        
        % calculate and apply the soft limit gains
        gainSoft = (1/knee(2)-1) * ( H_dB(idSoft)-limit_dB(idSoft) + ...
            knee(1)/2 ).^2 / (2*knee(1));
        gainSoft = 10.^(gainSoft/20);
        
        H(idSoft) = H(idSoft) .* gainSoft;
        
        % calculate and apply the compression/limiting
        if knee(2) ~= inf
            gainHard = limit_dB(idHard) + (H_dB(idHard)-limit_dB(idHard))/knee(2);
            gainHard = 10.^(gainHard/20);
            
            H(idHard) = H(idHard)./H_abs(idHard) .* gainHard(idHard);
        else
            H(idHard) = H(idHard)./H_abs(idHard) .* limit_lin(idHard);
        end
    end
    
    %remove the target ------
    if ~islogical(target)
        H = H .* target;
    end
    
    %revert data if lower soft limiting was applied ------
    if strcmpi(lim_dir{ll}, 'lower')
        H      = 1./H;
        target = 1./target;
    end
end