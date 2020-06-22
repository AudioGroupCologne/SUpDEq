% [l, r, az, el, HATO] = AKhrirInterpolation(az, el, HATO, type)
% interpolates head-related impulse responses (HRIRs) for desired source
% position and head-above-torso orientation (HATO). HATO is interpolated
% using inverse distance weighting on magnitude and unwrapped phase
% separately (c.f. [1]). Source positions can be interpolated in the
% spherical harmonics (SH) domain, or using HRIRs directly.
%
% See AKhrirInterpolationDemo.m for use cases
%
% I N P U T:
% az     - Azimuth angle(s) of desired source position(s) in degree
%          (0=front, 90=left, 180/-180=back, 270/-90=right). az can be a
%          scalar or vector. The size of az will be matched, in case az is
%          scalar, and el or HATO are vectors.
% el     - Elevation angle(s) of desired source position(s) in degree
%          (90=North Pole, 0=front, -90=South Pole). el can be a
%          scalar or vector. The size of el will be matched, in case el is
%          scalar, and az or HATO are vectors.
% HATO   - Head-above-torso orientation(s) for source position(s) in 
%          degree, between 50 and -50 (Same coordinate convention as az). 
%          HATO can be a scalar or vector. The size of HATO will be
%          matched, in case HATO is scalar, and el or HATO are vectors.
% type   - a string that specifies the data that is used for interpolation.
%          It consitst of up to three keywords that are spearated by
%          underscores:
%          'measured', or 'modeled': to use acuoustically measured, or
%                                    numerically simulted impulse responses
%          'sh', or 'ir'           : to use HRIRs stored in the spherical
%                                    harmonics (sh) domain or time domain
%                                    (ir)
%          'hrir', or 'dir'        : to use HRIRs or directional IRs where
%                                    the common transfer function was
%                                    removed.
%          e.g. 'measured_sh_hrir' which is the default, interpolates the
%          measured HRIRs in the SH domain.
%
% O U T P U T:
% l      - Left ear HRIRs at deisred source position and HATO
% r      - Right ear HRIRs at deisred source position and HATO
% az     - Source azimuth
% el     - Source elevation
% HATO   - HATOs clipped to the range between 50 and -50 degree.
%
%
% [1] F Brinkmann, R Roden, A Lindau, S Weinzierl: Audibility and
%     interpolation of head-above-torso orientation in binaural technology.
%     IEEE J. Sel. Topics Signal Process., 9(5):931-942 (2015)
%
%
%
% v1 2016/06 fabian.brinkmann@tu-berlin.de, Audio Communication Group, TU
%            Berlin

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
function [l, r, az, el, HATO] = AKhrirInterpolation(az, el, HATO, type)

%% ------------------------------------------------------ 1. HRIR parameter

% resolution of available HATOs
dHATO = 10;
% range of available HATOs
rHATO = [310 50];
% number of frequencies coded by SH coefficients
Nfreq = 129;
% Number of SH coefficients per frequency
Nnm = (35+1)^2;
% number of HRIR samples
Nsamp = 256;
% number of grid points in impulse response data
Ngrid = 11950;


%% --------------------------------------------------------- 2. input check
if ~exist('type', 'var')
    type = 'measured_sh_hrir';
end

% separate data type and interpolation domain
type = strsplit(type, '_');
for nn = 1:numel(type)
    if any( strcmpi( type{nn}, {'measured' 'modeled'} ) )
        type_acquisition = type{nn};
    elseif any( strcmpi( type{nn}, {'sh' 'ir'} ) )
        type_domain = type{nn};
    elseif any( strcmpi( type{nn}, {'hrir' 'dir'} ) )
        type_data = type{nn};
    else
        error('AKhrirInterpolation:type', [type{nn} 'is not a valid flag for ''type'''])
    end
end

if ~exist('type_acquisition', 'var')
    type_acquisition = 'measured';
end
if ~exist('type_domain', 'var')
    type_domain = 'sh';
end
if ~exist('type_data', 'var')
    type_data = 'hrir';
end

clear type

% number of outpt HRIRs
Nhrir = max([numel(el) numel(az) numel(HATO)]);

clear id


%% ---------------------------------------------- 3. Format input arguments
az   = reshape(mod(az,360), numel(az), 1);
el   = reshape(el, numel(el), 1);
HATO = reshape(mod(HATO,360), numel(HATO), 1);

if numel(az) ~= Nhrir
    az = repmat(az(1), Nhrir, 1);
end
if numel(el) ~= Nhrir
    el = repmat(el(1), Nhrir, 1);
end
if numel(HATO) ~= Nhrir
    HATO = repmat(HATO(1), Nhrir, 1);
end

%% -------------------- 4. List source positions and HATO for interpolation
% clip to range of available HATOs
HATO(HATO>rHATO(2) & HATO <=180) = 50;
HATO(HATO>180 & HATO <rHATO(1))  = 310;

% get HATOs for interpolation (eq. (4) in [1])
HATOint      = zeros(Nhrir, 2);
HATOint(:,1) = mod(floor(HATO/dHATO)*dHATO, 360);
HATOint(:,2) = mod((floor(HATO/dHATO)+1)*dHATO, 360);

% don't interpolate HATO if data exists
HATOint(~mod(HATO,dHATO),2) = NaN;

% get needed source positions (eq. (8) in [1])
azInt = mod(HATOint-repmat(HATO, [1 2])+repmat(az, [1 2]), 360);

% get HATOs needed for interpolation
HATOunique = unique(HATOint(:));
HATOunique = HATOunique(~isnan(HATOunique));


%% ----------------------------------------------------------- 5. Load data   
% needed HATOs
hato = HATOunique;

if strcmpi(type_domain, 'sh')
    % allocate memory
    data_l = zeros(Nnm, Nfreq, numel(hato));
    data_r = data_l;
    
    % load
    for nn = 1:numel(hato)
        try
            load(['FABIAN_HRIR_' type_acquisition '_HATO_' num2str(hato(nn))], 'SH');
            data_l(:,:,nn) = SH.coeffLeft;
            data_r(:,:,nn) = SH.coeffRight;
        catch
            AKdependencies('FABIAN')
        end
    end
else
    % allocate memory
    data_l = zeros(Nfreq, Ngrid, numel(hato));
    data_r = data_l;
    
    % load
    for nn = 1:numel(hato)
        try
            tmp = SOFAload(['FABIAN_HRIR_' type_acquisition '_HATO_' num2str(hato(nn)) '.sofa']);
            [tmp, ~, gridAll] = AKsofa2ak(tmp);
            gridAz            = gridAll(:,1);
            gridEl            = gridAll(:,2);
            data_l(:,:,nn) = AKboth2singleSidedSpectrum(fft(tmp(:,:,1)));
            data_r(:,:,nn) = AKboth2singleSidedSpectrum(fft(tmp(:,:,2)));
            
        catch
            AKdependencies('FABIAN')
        end
    end
end

%% ----------------------------------------- 6. interpolate source position
% allocate space for HRTFs needed for interpolation
lInt = zeros(Nfreq, Nhrir, 2);
rInt = lInt;

for nn = 1:numel(hato)
    for ii = 1:2
        id = HATOint(:,ii)==hato(nn);
        if any(id)
            if strcmpi(type_domain, 'sh') % -> get HRIRs from SH coefficients
                lInt(:,id,ii) = AKisht(data_l(:,:,nn), false, [azInt(id,ii) 90-el(id)], 'complex');
                rInt(:,id,ii) = AKisht(data_r(:,:,nn), false, [azInt(id,ii) 90-el(id)], 'complex');
            else % -> get HRIRs from IR data
                tmp = irInterpolate([azInt(id,ii) el(id)], data_l(:,:,nn), data_r(:,:,nn), gridAz, gridEl);
                lInt(:,id,ii) = tmp.left;
                rInt(:,id,ii) = tmp.right;
            end
        end
    end
end

clear tmp
%% ---------------------------------------------------- 7. interpolate HATO
% only interpolate when necessary
idInt = ~isnan(HATOint(:,2));

% allocate space for final HRIRs
l = zeros(Nfreq, numel(HATO));
r = l;

if any(idInt)
    % get weights for interpolation (eq. (5), and (4) in [1])
    HATOweight = acosd(cosd(HATOint(idInt,:)-repmat(HATO(idInt), [1 2])));
    HATOweight = HATOweight ./ repmat(sum(HATOweight,2), [1 2]);

    HATOweightA = repmat(HATOweight(:,2)', [Nfreq 1]);
    HATOweightB = repmat(HATOweight(:,1)', [Nfreq 1]);


    % interpolate HATO separately for magnitude and unwrapped phase
    l(:,idInt) = ( abs(lInt(:,idInt,1)).*HATOweightA + abs(lInt(:,idInt,2)).*HATOweightB ) .* ...
                   exp(1j * ( unwrap(angle(lInt(:,idInt,1))).*HATOweightA + unwrap(angle(lInt(:,idInt,2))).*HATOweightB ) );

    r(:,idInt) = ( abs(rInt(:,idInt,1)).*HATOweightA + abs(rInt(:,idInt,2)).*HATOweightB ) .* ...
                   exp(1j * ( unwrap(angle(rInt(:,idInt,1))).*HATOweightA + unwrap(angle(rInt(:,idInt,2))).*HATOweightB ) );
end

% don't interpolare HATO if not neccessary           
l(:,~idInt) = lInt(:,~idInt,1);
r(:,~idInt) = rInt(:,~idInt,1);

%% ----------------------------------------------------------- 8. get HRIRs
% set bin at fs/2 to 0
l(end,:) = 0;
r(end,:) = 0;

% get HRIRs
l = ifft(AKsingle2bothSidedSpectrum(l), 'symmetric');
r = ifft(AKsingle2bothSidedSpectrum(r), 'symmetric');

% get DIRs
if strcmpi(type_data, 'dir')
    ctf = SOFAload(['FABIAN_CTF_' type_acquisition '_inverted_smoothed.sofa']);
    ctf = AKsofa2ak(ctf);
    
    l = fftfilt(ctf, l);
    r = fftfilt(ctf, r);
end

% WARNING: an earlier version of AKtools reversed the phase of the HRIRs
% because AKsht and AKist used the ' operator for transposing a complex
% vector. However, Matlab calculates the complex conjugate transpose in
% this case. The current version now uses the .' operator.
if strcmpi(type_domain, 'sh')
    if sum(abs(l(1:Nsamp/2,1))) < sum(abs(l(Nsamp/2:end,1)))
        l = flipud(l);
        r = flipud(r);
    end
end

end


function data = irInterpolate(pos, left, right, az, el)
% function to obtain HRTF(s). Uses Matlabs 'matfile' function to load
% single HRTFs without having to load the whole data set. 
% However, it only works with Matlab 2012a and later.
%
% INPUT
% pos   - positions to be obtained
%         [source_azimuth1 source_elevation1;
%          source_azimuth2 source_elevation2; ...
% left  - left ear single sided complex spectra
% right - right ear single sided complex spectra
% az    - source azimuth
% el    = source elevation
%
% OUTPUT
% data  - HRIRs
%
%
% v1.  2012/07 fabian.brinkmann.tu-berlin.de, Audio Communication Group,
%              TU Berlin
% v1.1 2016/06 adapted to new data format and interpolation in frequency
%              domain

% length of HRIRs and sampling frequency
N_ir = length(left(:,1));

% allocate memory for output
data.left  = zeros(N_ir, size(pos,1));
data.right = zeros(N_ir, size(pos,1));

for n = 1:size(pos,1)
    % --------------------------- get indices of closest available point(s)
    % check if source position is included in the grid
    if any(az == pos(n,1) & el == pos(n,2))
        % -- no interpolation between sources --
        pos_i.ir_id = find(az == pos(n,1) & el == pos(n,2));
        pos_i.ir_d  = 1;
        pos_match = true;
    else
        pos_match = false;
    end
    
    if pos_match % -> source position included in grid
        data.left(:,n)  = left(:, pos_i.ir_id);
        data.right(:,n) = right(:, pos_i.ir_id);
        
    else % -> source position not included in grid
        % -- interpolation between sources --
        % calculate great circle distance from desired position to every point
        % in the grid
        cur_pos = repmat(pos(n,1:2), size(az, 1), 1);
        d = acosd(sind(cur_pos(:,2)) .* sind(el) +...
                  cosd(cur_pos(:,2)) .* cosd(el) .*...
                  cosd(cur_pos(:,1) - az));
        % sort
        [d, d_id] = sort(d); % ...dann suche k?rzeste Orthodome zu dem gew?nschten Punkt
        
        
        % get points for interpolation
        if ~pos_match
            % interpolate between neighboring points of same elevation
            if ~mod(pos(n,2), 2) % wenn die gew?nschte Quellposition auf einem grid-Kreis liegt (also die Elevation im Grid liegt)...
                N = 2;
                pos_i.ir_id = d_id(1:N); % ...dann werden die zwei n?chsten Punkte auf diesem grid-Kreis zur Interpolation verwendet
                pos_i.ir_d  = d(1:N);
                % make sure that the two point are on the same
                % elevation
                while el(pos_i.ir_id(1)) ~= el(pos_i.ir_id(2))
                    N = N+1;
                    pos_i.ir_id(2) = d_id(N);
                    pos_i.ir_d(2)  = d(N);
                end
                % interpolate within triangle of neighboring points
            else % wenn die gew?nschte Quellposition auf NICHT auf einem grid-Kreis liegt (also die Elevation NICHT im Grid liegt)...
                N = 3;
                pos_i.ir_id = d_id(1:3);% ...dann werden die drei n?chsten Punkte im Grid zur Interpolation verwendet
                pos_i.ir_d  = d(1:3);
                % wobei folgendes sichergestellt werden muss:
                % make sure that
                % - one of the three points is from a differing
                %   elevation
                % - there are two points below and one above the given
                %   point (or the other way around)
                while ( ...
                        el(pos_i.ir_id(1)) == el(pos_i.ir_id(2)) && ...
                        el(pos_i.ir_id(1)) == el(pos_i.ir_id(3)) ...
                        ) ...
                        || ...
                        ( ...
                        ( el(pos_i.ir_id(1))<=pos(n,2) && el(pos_i.ir_id(2))<=pos(n,2) && el(pos_i.ir_id(3))<=pos(n,2) ) ...
                        || ...
                        ( el(pos_i.ir_id(1))>=pos(n,2) && el(pos_i.ir_id(2))>=pos(n,2) && el(pos_i.ir_id(3))>=pos(n,2) ) ...
                        )
                    
                    N = N+1;
                    pos_i.ir_id(3) = d_id(N);
                    pos_i.ir_d(3)  = d(N);
                end
            end
        end
        
        % -------------------------- inverse distance interpolation
        % allocate memory for interpolated HRIRs
        ir_l_m = zeros(N_ir, 1);
        ir_l_a = ir_l_m;
        ir_r_m = ir_l_m;
        ir_r_a = ir_l_m;
        
        
        for ir_id = 1:length(pos_i.ir_id)
            % get HRIR
            ir_tmp_l   = left(:, pos_i.ir_id(ir_id));
            ir_tmp_r   = right(:, pos_i.ir_id(ir_id));
            
            % interpolate magnitude and unwrapped phase
            ir_l_m = ir_l_m + abs(ir_tmp_l)           / pos_i.ir_d(ir_id);
            ir_l_a = ir_l_a + unwrap(angle(ir_tmp_l)) / pos_i.ir_d(ir_id);
            
            ir_r_m = ir_r_m + abs(ir_tmp_r)           / pos_i.ir_d(ir_id);
            ir_r_a = ir_r_a + unwrap(angle(ir_tmp_r)) / pos_i.ir_d(ir_id);
        end
        % normalize HRIRs by inverse distance weights
        data.left(:,n)  = ir_l_m./sum(1./pos_i.ir_d) .* exp(1j*ir_l_a./sum(1./pos_i.ir_d));
        data.right(:,n) = ir_r_m./sum(1./pos_i.ir_d) .* exp(1j*ir_r_a./sum(1./pos_i.ir_d));
    end
end

end
