function [hISM, hSR, hHybrid, ISM, SR] = AKroomSimulation(rs)
% [hISM, hSR, hHybrid, ISM, SR] = AKroomSimulation(rs)
% room simulation of a recangular room
%
% for example usage see
% AKroomSimulationDemo.m
%
% I N P U T
% rs      - input data, see AKroomSimulationDemo.m for explanation
%
%
% O U T P U T
% hISM    - impulse response holding the early reflection calculated with
%           the image source model for a rectangular room.
%           Size [N x RO x RC x S]
%           N : length in samples
%           RO: number of receiver orientations
%           RC: number of reciever channels
%           S : number of sources
% hSR     - impulse response holding the late/stochastic reverberation
%           calulated based on colored and decaying noise.
%           Size [M x RC x S]
%           M : lengths in samples
%           RC: nuber of receiver channels
%           S : number of sources
% hHybrid - combination of hISM and hSR
% ISM     - Struct holding the output of AKism for each source position
%           d       : distances between receiver and sources in ascending
%                     order [M x 1]
%           A       : amplitudes at the F frequencies of rs.f [M x F]
%           Rec_az  : source exit angle (azimuth) [M x 1]
%           Rec_el  : source exit angle (elevation) [M x 1]
%           Rec_AZ  : receiver incident angle (azimuth) [M x 1]
%           Rec_EL  : receiver incident angle (elevation) [M x 1]
%           Src_az  : receiver incident angle after rotation the receiver
%                     according to rs.recRot (azimuth) [M x R]
%           Src_el  : receiver incident angle after rotation the receiver
%                     according to rs.recRot (elevation) [M x R]
%           X       : absolute image source positions [M x 3] x/y/z
%           X_rel   : image source positions relative to the receiver
%                     [M x 3] x/y/z
%           N       : image source order [M x 1]
%           wallLog : struct that contains information of the reflection
%                     paths of all reflections up to the 3rd order
%                     ID       : gives the walls that a ray was reflected
%                                at in reverse order
%                     position : intersection points of rays and walls
%                                x/y/z
%                     N        : image source order
% SR      - Struct holding additional parameters for the calculation of the
%           stochastic reverberation.
%           V              : room volume [m^3]
%           S              : surface of each wall [m^2]
%           alpha_mean     : mean absorption for each frequency in rs.f
%                            divided by the surface of the room
%           T              : reverberation time calculated from the
%                            parameters above
%           dtf_src        : diffuse field transfer function of the source
%           dtf_rec        : diffuse field transfer function of the
%                            receiver
%           dtf            : combination of the above
%           T_interpolated : reverberation time interpolated to the
%                            frequencies that are needed to for calculating
%                            the stochastic reverberation
%           f_interpolated : frequencies for the above
%           mixingTime     : perceptual mixing time estimated from the room
%                            volume [ms]
%           gain           : gain correction for the stochastic
%                            reverberation to match the level of the ISM at
%                            the mixing time
%
%
% [1] Allen, J. B. & Berkley, D. A.: "Image method for efficiently
%     simulating small-room acoustics." J. Acoust. Soc. Am., 65(4),
%     943-950 (1979).
% [2] Lehmann, E. A. & Johansson, A. M.: "Prediction of energy decay in
%     room impulse responses simulated with an image-source model."
%     J. Acoust. Soc. Am., 124(1), 269-277 (2008).
% [3] Borss, C. & Martin, R.: "An improved parametric Model for perception-
%     based design of virtual acoustics." AES 35th International Conference
%     London, UK, 2009.
% [4] Brinkmann, F. & Erbes, V. & Weinzierl, S.: "Extending the closed form
%     image source model for source directivty." Fortschritte der Akustik -
%     DAGA 2018, Munich, Germany, March 2018.
%
%
% P R O C E S S I N G
% 0. pre-processing
% 1. load directivity data
% 2. image source model (ISM)
% 3. stochastic reverberation (SR)
% 4. apply inverse diffuse field transfer function of the receiver
% 5. combine ISM and SR impulse responses
%
% 2018/02 fabian.brinkmann@tu-berlin.de

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

%% ------------------------------------------------------ 0. pre-processing

if rs.verbose
    disp('AKroomSimulation: Checking input data')
end

% --- check if the source view was passed --- %
if isempty(rs.srcView)
    
    % get source position in spherical coordinates
    if strcmpi(rs.srcCoordinates, 'spherical')
        tmp = rs.srcPos;
    else
        [az, el] = cart2sph(rs.srcPos(:,1)-rs.recPos(1), rs.srcPos(:,2)-rs.recPos(2), rs.srcPos(:,3)-rs.recPos(3));
        tmp      = [az el];
        tmp      = tmp*180/pi;
    end
    
    % wrap to 360 deg
    tmp(:,1) = mod(tmp(:,1), 360);
    
    % mirror to get source rotation
    rs.srcView = [mod(180+tmp(:,1), 360) -tmp(:,2)];
    clear az el tmp
    
end

% --- get the source position(s) in carthesian coordinates --- %
if strcmpi(rs.srcCoordinates, 'spherical')
    [x, y, z] = sph2cart(rs.srcPos(:,1)/180*pi, rs.srcPos(:,2)/180*pi, rs.srcPos(:,3));
    rs.srcPos = [x + rs.recPos(1) y + rs.recPos(2) z + rs.recPos(3)];
    
    rs.srcCoordinates = 'carthesian';
    
    clear x y z
end

% --- check if receiver, and source(s) are inside the room --- %
if any(rs.recPos < 0)    || rs.recPos(1) > rs.L(1)        || rs.recPos(2) > rs.L(2)        || rs.recPos(3) > rs.L(3)      || ...
        any(rs.srcPos(:) < 0) || any(rs.srcPos(:,1) > rs.L(1)) || any(rs.srcPos(:,1) > rs.L(1)) || any(rs.srcPos(:,1) > rs.L(1)) 
    error('AKroomSimulation:input', 'At least one receiver or source position is outside the room.')
end

%% ----------------------------------------------- 1. load direcitvity data

if rs.verbose
    disp('AKroomSimulation: Loading directivity data')
end

if strcmpi(rs.src, 'QSC-K8')
    H.src  = load('QSC_K8_directivity'); % spherical harmonics coefficients
    H.srcN = 1024;                                          % length of impulse responses
    
    % give the mean time of arrival, i.e. number of leading zeros in
    % the impulse responses (source directivity is minimum phase)
    H.srcTOA = 0;
else
    H.srcN   = 0;
    H.srcTOA = 0;
end

if strcmpi(rs.rec, 'FABIAN')
    H.rec  = load('FABIAN_HRIR_measured_HATO_0');           % spherical harmonics coefficients
    H.recN = 256;                                           % length of impulse response
    
    % estimate the mean time of arrival, i.e. number of leading zeros in
    % the impulse responses (use left ear data only)
    hTmp     = AKisht(H.rec.SH.coeffLeft, true, [[0 0 90 180 270 0]' [0 90 90 90  90 180]'], 'complex');
    H.recTOA = AKonsetDetect(hTmp, 0, 6);
    H.recTOA = round(mean(H.recTOA));
    clear hTmp
    
    rs.recHRTF = true;
else
    H.recN   = 0;
    H.recTOA = 0;
    
    rs.recHRTF = false;
end


%% -------------------------------------------------- 2. image source model

if rs.ISM
    
    % allocate output struct
    ISM = struct('d',             [], ...
                 'A',             [], ...
                 'Rec_az',        [], ...
                 'Rec_el',        [], ...
                 'Rec_AZ',        [], ...
                 'Rec_EL',        [], ...
                 'Src_az',        [], ...
                 'Src_el',        [], ...
                 'X',             [], ...
                 'X_rel',         [], ...
                 'N',             [], ...
                 'wallLog',       []);
    ISM = repmat(ISM, [size(rs.srcPos, 1) 1]);
    
    % track the maximum distance between receiver and image sources
    dMax = 0;
    dMin = NaN;
    
    % --- a: ISM positions --- %
    for nn = 1:size(rs.srcPos, 1)
        
        if rs.verbose
            disp(['AKroomSimulation - ISM: Calculating image source model (source ' num2str(nn) '/' num2str(size(rs.srcPos, 1)) ')'])
        end
        
        % get the image sources: distance, amplitude, and exit/incident angles
        [ISM(nn,1).d, ISM(nn,1).A, ISM(nn,1).Rec_az, ISM(nn,1).Rec_el, ISM(nn,1).Src_az, ISM(nn,1).Src_el, ISM(nn,1).X, ISM(nn,1).X_rel, ISM(nn,1).N, ISM(nn,1).wallLog] ...
            = AKism(rs.L, rs.srcPos(nn,:), rs.srcView(nn,:), rs.recPos, rs.recView, rs.alpha, rs.ISMtruncation, rs.c);
        
        % update maximum distance
        dMax = max( dMax, max(ISM(nn,1).d) );
        dMin = min( dMin, min(ISM(nn,1).d) );
        
    end
    
    % --- b: apply receiver rotation to incident angles --- %
    % (only needed if the receiver directivity is not OMNI)
    if ~strcmpi(rs.rec, 'OMNI') && any(rs.recRot(:))
        
        if rs.verbose
            disp('AKroomSimulation - ISM: Applying head rotation')
        end
        
        % loop source positions
        for nn = 1:size(rs.srcPos, 1)
            % apply the rotation
            [ISM(nn).Rec_AZ, ISM(nn).Rec_EL] = AKroomSimulationRotation(ISM(nn).Rec_az, ISM(nn).Rec_el, rs.recRot(:,1), rs.recRot(:,2));
        end
    else
        for nn = 1:size(rs.srcPos, 1)
            ISM(nn).Rec_AZ = ISM(nn).Rec_az;
            ISM(nn).Rec_EL = ISM(nn).Rec_el;
        end
    end
    
    % --- c: calculate impulse responses --- %
    
    % disable air absorption, if distances becoe to large. The current
    % implementation gets inefficient in this case
    if dMax > 500 && rs.airAbsorption
        disp('AKroomSimulation - ISM: Air absorption disabeled because image-source-to-receiver distances of > 500 m were found')
        rs.airAbsorption = false;
    end
    
    % get length of each reflection in samples
    if rs.airAbsorption || size(rs.alpha,1) > 1
        H.N = max([H.srcN H.recN 512]);
    else
        H.N = max([H.srcN H.recN 1]);
    end
    % make sure its even
    if H.N > 1
        H.N = H.N + mod(H.N,2);
    end
    
    % get frequency vector and air absorption
    f         = (0:rs.fs/H.N:rs.fs/2)';
    alpha_air = AKairAttenuation(f, rs.h_r, rs.T, rs.p_a);
    
    % allocate space for ISM impulse responses
    if rs.recHRTF
        C = 2;
    else
        C = 1;
    end
    if islogical(rs.recRot)
        hISM = zeros(round(dMax/rs.c*rs.fs)+H.N, 1,                  C, size(rs.srcPos, 1));
    else
        hISM = zeros(round(dMax/rs.c*rs.fs)+H.N, numel(rs.recRot)/2, C, size(rs.srcPos, 1));
    end
    clear C
    
    % loop across source positions
    for ss = 1:size(rs.srcPos, 1)
        
        if rs.verbose
            disp(['AKroomSimulation - ISM: Calculating impulse response(s) from image source model (source ' num2str(ss) '/' num2str(size(rs.srcPos, 1)) ')'])
        end
        
        % loop across image sources / reflections
        for ii = 1:numel(ISM(ss).d)
            
            % apply wall, and distance induced damping
            if size(rs.alpha,1) == 1
                h    = ISM(ss).A(ii) * ones(floor(H.N/2)+1,1);
            else
                h    = AKinterpolateSpectrum( ISM(ss).A(ii,:)', rs.f', H.N, {'nearest' 'lin' 'lin'}, rs.fs);
            end
            
            % apply air absorption
            if rs.airAbsorption
                h = h .* 10.^(-alpha_air*ISM(ss).d(ii)/20);
            end
            
            % filter with source directivity
            if H.srcN
                % get source directivity
                hTmp = AKisht(H.src.SH.coeff, H.src.SH.doFFT, [ISM(ss).Src_az(ii) 90-ISM(ss).Src_el(ii)], H.src.SH.SHTmode, H.src.SH.isEven, H.src.SH.compact, H.src.SH.SHmode);
                % discard imaginary small part, due to rounding erros
                hTmp = real(hTmp);
                % interpolate to desired frequencies
                hTmp = AKinterpolateSpectrum(hTmp, H.src.N, H.N, {'nearest' 'linear' 'nearest'}, rs.fs);
                
                % apply
                h = h .* hTmp;
            end
            
            % get min-phase IR
            if H.N > 1
                h(end,:,:) = abs(h(end,:,:));
                h = AKsingle2bothSidedSpectrum(h);
                h = ifft(h, 'symmetric');
                h = AKphaseManipulation(h, rs.fs, 'min', 4, false);
            end
            
            
            % filter with receiver directivity
            if H.recN
                % SH values
                Ynm = AKsh(H.rec.SH.order, [], ISM(ss).Rec_AZ(ii,:)',  90-ISM(ss).Rec_EL(ii,:)');
                
                hL = AKisht(H.rec.SH.coeffLeft,  true, Ynm, 'complex');
                hR = AKisht(H.rec.SH.coeffRight, true, Ynm, 'complex');
                                
                % filter
                if size(hL, 1) < size(h,1)
                    hL = fftfilt(hL, h);
                    hR = fftfilt(hR, h);
                else
                    hL = fftfilt(h, hL);
                    hR = fftfilt(h, hR);
                end
                
                % put together
                h = cat(3, hL, hR);
            end
            
            % add the final reflection to the impulse response
            n = round(ISM(ss).d(ii)/rs.c * rs.fs) + 1;
            hISM(n:n+H.N-1, :, :, ss) = hISM(n:n+H.N-1, :, :, ss) + h;
            
        end
        
    end
    
    % compensate the leading zeros from the directivity impulse responses
    if rs.verbose
        disp('AKroomSimulation - ISM: Removing leading zeros from directivity files')
    end
    
    H.totalTOA = H.srcTOA + H.recTOA;
    if round(dMin/rs.c * rs.fs) > H.totalTOA + 20  % add a safety margin of 20 samples
        hISM = hISM( max(H.totalTOA,1):end ,:,:,:);
        H.totalTOAcompensated = true;
    else
        H.totalTOAcompensated = false;
        disp('AKroomSimulation - ISM: Leading zeros could not be removed from directivity files, because the initial delay in the impulse response is too short') 
    end
    
    clear h hL hR hTmp n nn mm ss ii Ynm f alpha_air
    
else
    ISM  = [];
    hISM = [];
end


%% --------------------------------------------------- 3. stochastic reverb

if rs.SR
    
    % --- a: estimate the reverberation time from the absorption coefficients and room size --- %
    if rs.verbose
        disp('AKroomSimulation - SR: Estimating reverberation time from rs.alpha')
    end
    % room volume
    SR.V = prod(rs.L);
    % surface of each wall
    SR.S = [rs.L(2)*rs.L(3) rs.L(2)*rs.L(3) ...
            rs.L(1)*rs.L(3) rs.L(1)*rs.L(3) ...
            rs.L(1)*rs.L(2) rs.L(1)*rs.L(2)];
    % mean absorption
    SR.alpha_mean = AKm(rs.alpha, SR.S, '*');
    SR.alpha_mean = sum(SR.alpha_mean, 2) / sum(SR.S); 
    % take air absorption into account
    [~, m] = AKairAttenuation(rs.f);
    % estimate reverberation time using (Eyring's) Sabine's formula
    % including air damping
    if rs.airAbsorption
        % SR.T = .161 * SR.V ./ (-sum(SR.S) * log(1-SR.alpha_mean) + 4*m*SR.V);
        SR.T = .161 * SR.V ./ (sum(SR.S) * SR.alpha_mean  + 4*m*SR.V);
    else
        % SR.T = .161 * SR.V ./ (-sum(SR.S) * log(1-SR.alpha_mean) );
        SR.T = .161 * SR.V ./ (sum(SR.S) * SR.alpha_mean);
    end
    
    clear m
    
    
    % --- b: calcualte diffuse field transfer function --- %
    if rs.verbose
        disp('AKroomSimulation - SR: Calculate diffuse field transfer function')
    end
    if strcmpi(rs.src, 'QSC-K8')
        SR.dtf_src = sqrt( AKshEnergy(H.src.SH.coeff) )';
        SR.dtf_src = AKinterpolateSpectrum(SR.dtf_src, H.src.N, H.N, {'nearest' 'linear' 'nearest'}, rs.fs);
        SR.dtf_src = AKsingle2bothSidedSpectrum(SR.dtf_src, H.srcN);
        SR.dtf_src = ifft(SR.dtf_src, 'symmetric');
        SR.dtf_src = AKphaseManipulation(SR.dtf_src, rs.fs, 'min', 4, false);
    else
        SR.dtf_src = false;
    end
    
    if strcmpi(rs.rec, 'FABIAN')
        SR.dtf_rec = mean( [...
                         sqrt( AKshEnergy(H.rec.SH.coeffLeft) )'  ...
                         sqrt( AKshEnergy(H.rec.SH.coeffRight) )' ...
                         ], 2);
        SR.dtf_rec = AKfractOctSmooth(SR.dtf_rec, 'welti', rs.fs, 3);
        SR.dtf_rec = AKsingle2bothSidedSpectrum(SR.dtf_rec, H.srcN);
        SR.dtf_rec = ifft(SR.dtf_rec, 'symmetric');
        SR.dtf_rec = AKphaseManipulation(SR.dtf_rec, rs.fs, 'min', 4, false);
        if islogical(SR.dtf_src)
            SR.dtf = SR.dtf_rec;
        else
            if H.recN > H.srcN
                SR.dtf = fftfilt(SR.dtf_src, SR.dtf_rec);
            else
                SR.dtf = fftfilt(SR.dtf_rec, SR.dtf_src);
            end
        end
    else
        if islogical(SR.dtf_src)
            SR.dtf = false;
        else
            SR.dtf = SR.dtf_src;
        end
    end
    
    
    % --- c: calculate the stochastic reverb --- %
    if rs.recHRTF
        
        if rs.verbose
            disp('AKroomSimulation - SR: Calculate the stochastic reverberation for a binaural receiver')
        end
        
        for nn = 1:size(rs.srcPos, 1)
            % calculate stochastic reverberation
            [reverb_tail, SR.T_interpolated, SR.f_interpolated] = AKdiffuseReverbTail(SR.T', rs.f', 'V', SR.V, 'c', rs.c, 'fs', rs.fs, ...
                                                                                      'extrap_mode', {'nearest', 'spline', 'linear'}, 'decay_mode', rs.SRdecayMode, 'decay_range', rs.SRdynamic, ...
                                                                                      'bin_coherence', true, 'dtf', SR.dtf, 'do_plot', false);
            
            % allocate space for output
            if nn == 1
                hSR = zeros(size(reverb_tail, 1), 2, size(rs.srcPos, 1));
            end
            
            % save to output
            hSR(:,1,nn) = reverb_tail(:,1);
            hSR(:,2,nn) = reverb_tail(:,2);
            
        end
        
    else
        
        if rs.verbose
            disp('AKroomSimulation - SR: Calculate the stochastic reverberation for a non-binaural receiver')
        end
        
        % calculate stochastic reverberation
        [reverb_tail, SR.T_interpolated, SR.f_interpolated] = AKdiffuseReverbTail(SR.T', rs.f', 'V', SR.V, 'c', rs.c, 'fs', rs.fs, ...
                                                                          'extrap_mode', {'nearest', 'spline', 'linear'}, 'decay_mode', rs.SRdecayMode, 'decay_range', rs.SRdynamic, ...
                                                                          'bin_coherence', false, 'dtf', SR.dtf, 'rev_channel', size(rs.srcPos, 1), 'do_plot', false);
        
        % convert to desired format
        hSR = zeros(size(reverb_tail, 1), 1, size(rs.srcPos, 1));
        for nn = 1:size(rs.srcPos, 1)
            hSR(:,:,nn) = reverb_tail(:,nn);
        end
    end
    
    % --- d: post-process the stochastic reverberation --- %
    if rs.ISM
        if rs.verbose
            disp('AKroomSimulation - SR: Delay, level, and fade-in')
        end
        
        % ---delay the SR ---
        SR.mixingTime = dMax/rs.c;
        N_mix   = round(dMax/rs.c * rs.fs);
        N_start = round(dMin/rs.c * rs.fs);
        hSR     = [zeros( N_start, size(hSR, 2), size(hSR, 3) ); hSR];
        
        % --- level the SR ---
        % get mean the level of the SR around the truncation time
        if N_mix+256 < size(hSR,1)
            L_SR = max( hSR( N_mix-256:N_mix+256, :,:,: ) );
            L_SR = mean( L_SR(:) );
        else
            disp('AKroomSimulation - SR: hISM is longer than hSR. The level of the SR could not be set')
            L_SR = 1;
        end
        
        % estimate the level of the ISM at the truncation time between 500
        % and 1000 Hz
        if size(rs.alpha,1) == 1
            id = 1;
        else
            id = rs.f >=500 & rs.f <= 1e3;
            if ~any(id)
                id = rs.f <= 500;
                if ~any(id)
                    id = 1:size(rs.alpha,1);
                end
            end
        end
        
        L_ISM = zeros(size(rs.srcPos, 1), sum(id));
        for nn = 1:size(rs.srcPos, 1)
            L_ISM(nn,:) = max( ISM(nn).A( max(1, size(ISM(nn).A,1)-10):end, id) );
        end
        L_ISM = mean(L_ISM(:));
        
        % estimate the level change caused by the directivity files
        if any(SR.dtf)
            [~, L_DTF] = AKnormalize(SR.dtf, 'abs', 'mean', 'mean', 1, [500 1000], false, rs.fs);
        else
            L_DTF = 1;
        end
        
        % get the overall gain
        SR.gain = L_ISM*L_DTF / L_SR;
        % apply gain correction
        SR.gain = SR.gain * 10^(rs.SRgain/20);
        % level the stochastic reverberation
        hSR     = hSR * SR.gain;

        
        % --- fade in the SR ---
        
        % get the sample where the fade in starts
        if isnumeric(rs.SRfadeDuration)
            % taken from input parameters
            N_fade = max(1, N_mix - round(rs.SRfadeDuration*rs.fs));
        else
            % detect the time of arrival of the first reflection
            N_fade = zeros(size(rs.srcPos, 1), 1);
            for nn = 1:size(rs.srcPos, 1)
                N_fade(nn) = ISM(nn).d( min( numel(ISM(nn).d), 2 ) );
            end
            N_fade = mean(N_fade) / rs.c * rs.fs;
            N_fade = round(N_fade);
            
            % add time of arrival detected in the direcitvity impulse
            % responses
            if ~ H.totalTOAcompensated
                N_fade = N_fade + H.totalTOA;
            end
            
        end
        
        % fade duration
        N = N_mix - N_fade + 1;
        
        % generate the fade
        if N > 1
            if strcmpi(rs.SRfadeType, 'lin')
                % linear fade
                fade = linspace(0, 1, N)';
            else
                % sine fade
                id = strfind(rs.SRfadeType, '_');
                if ~isempty(id)
                    sine_power = str2double( rs.SRfadeType(id+1:end) );
                else
                    sine_power = 1;
                end
                
                fade = sin( linspace(0, pi/2, N)' ).^sine_power;
            end
            
            % remove unwanted part before N_fade
            if N_fade < size(hSR,1)
                hSR(1:N_fade, :,:,:) = 0;
            end
            % apply the fade
            if N_mix < size(hSR,1)
                hSR(N_fade:N_mix, :,:,:) = AKm( hSR(N_fade:N_mix, :,:,:), fade, '*' );
            end
        else
            % remove unwanted part before N_mix
            hSR(1:N_mix, :,:,:) = 0;
        end
        
    end
    
    clear dMax fade id L_DTF L_ISM L_SR N N_fade N_mix nn reverb_tail nn
    
else
    SR  = [];
    hSR = [];
end


%% ------- 4. apply inverse diffuse field transfer function of the receiver
if rs.recHRTF && rs.recDTFcompensation
    if rs.verbose
        disp('AKroomSimulation: Applying diffuse field compensation')
    end
    
    % get the inverse diffuse field transfer function
    if strcmpi(rs.rec, 'FABIAN')
        dtf = mean( [...
            sqrt( AKshEnergy(H.rec.SH.coeffLeft) )'  ...
            sqrt( AKshEnergy(H.rec.SH.coeffRight) )' ...
            ], 2);
        dtf = 1./dtf;
        dtf = AKsoftLimit(dtf, 0, 5, [0 rs.fs/2], rs.fs, true);
        dtf = AKfractOctSmooth(dtf, 'welti', rs.fs, 3);
        dtf = AKsingle2bothSidedSpectrum(dtf, H.srcN);
        dtf = ifft(dtf, 'symmetric');
        dtf = AKphaseManipulation(dtf, rs.fs, 'min', 1, false);
    else
        dtf = false;
    end
    
    % apply the inverse diffuse field transfer function
    if rs.ISM && any(dtf)
        for nn = 1:size(rs.srcPos, 1)
            for mm = 1:size(hISM, 3)
                hISM(:, :, mm, nn) = fftfilt(dtf, hISM(:, :, mm, nn));
            end
        end
    end
    
    if rs.SR && any(dtf)
        for nn = 1:size(rs.srcPos, 1)
            hSR(:, :, nn) = fftfilt(dtf, hSR(:, :, nn));
        end
    end
    
    if ~any(dtf)
        if rs.verbose
            disp('AKroomSimulation: Unknown receiver. Diffuse field compensation aborted!')
        end
    end
    
end


%% -------------------------------- 5. combine ISM and SR impulse responses

if rs.ISM && rs.SR
    
    if rs.verbose
        disp('AKroomSimulation: Combining ISM and SR')
    end
    
    % allocate space
    Nism = size(hISM,1);
    Nsr  = size(hSR, 1);
    hHybrid = zeros(max(Nism, Nsr), size(hISM,2), size(hISM,3), size(hISM,4));
    
    % move SR IRs into place
    for mm = 1:size(hSR,3)
        for nn = 1:size(hSR,2)
            hHybrid(1:Nsr,:,nn,mm) = repmat(hSR(1:Nsr,nn,mm), [1 size(hISM,2)]);
        end
    end
    
    % add ISM
    hHybrid(1:size(hISM,1),:,:,:) = hHybrid(1:size(hISM,1),:,:,:) + hISM;

else
    hHybrid = [];
end