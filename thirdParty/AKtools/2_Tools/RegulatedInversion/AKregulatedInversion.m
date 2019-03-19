function data = AKregulatedInversion(s, raw_data)
% frequency dependant regulated inversion
%
% See AKregulatedInversionDemo.m for examples
%
% I N P U T:
% s                       - struct containing the parameters for inversion
%                           (see AKregulatedInversionDemo.m)
% raw_data (optional)     - raw_data as impulse responses
%                           [rows=samples columns=measurements 3rd dimension=channel]
%                           If raw_data is not passed, wav-files from
%                           s.data_dir are used instead
% 
% O U T P U T:
% Final and premilinary results are stored in the 'data' struct wich is
% returned by the function, it includes:
% raw                     - data after processing step 2
%                           [samples x measurements x channels]
% prefilter               - data after processing step 3
% aligned                 - data after processing step 4
% truncated               - data after processing step 5
% averaged                - data after processing step 6
%                           [samples x channels]
% regularization          - regularization function as calculated in
%                           processing step 8
% target                  - target function as calculated in
%                           processing step 9
% average_sm              - data after processing step 10
% average_norm            - data after processing step 11
% compensation_filter_raw - compensation filter as calculated in proccesing
%                           step 12
% compensation_filter     - compensation filter as calculated in proccesing
%                           step 14 and saved in step 15
% compensation_result     - simulated compensation result from plot 7 in
%                           step 16
% ERB                     - energetic difference in auditory filters from
%                           plot 7 in step 16
%                           1. dimension = number of measurements
%                           2. dimension = number of filters
%                           3. dimension = number of channles
% f_c                     - center frequencies of auditory filters
%
%
% PROCESSING:
% 1.  directory and data handling:
%     create directories for saving (in-between) results of processing steps
%     and documentation
% 2.  read in raw transfer functions that should be inverted
% 3.  apply pre-filter to compensate for transducers in the measurement
%     chain (e.g. recording microphones)
% 4.  time align impulse responses based on an the time of arrival (estimated
%     by onset detection). This compensates for varying delays and time of
%     flight.
% 5.  Truncate impulse responses
% ------ plot1: time alignment and trunciation
% 6.  Average impulse responses
% ------ plot2: averaging and frequency range for normalization (step 11)
% 7.  decide if regularization or smoothing is will be applied (based on
%     s.inverion method)
% 8.  Caluclate regularization function (b,a coefficients and impulse response)
% ------ plot3: regularization function
% 9.  Calculate target function
% ------ plot4: target function
% 10. Fract. octave smoothing and min-phase of averaged transfer functions
%     (optional)
% 11. Normalize averaged transfer functions to match level of targer
%     function (see plot2 for normalization range)
% 12. Calculate compensation filter (AKffdInverseFilter.m)
% ------ plot5: compensation filtering result 1
%
%    AKffdInverseFilter.m
%       Returns the FIR-LMS-inverse using regularization and a target function with arbitrary phase
%       - implementation acc. to [1]
%       - length of filter N will match length of input x
%       - 3 phase methods:
%         a) Filter design with minimumphase target function (default)
%         b) Filter design with linearphase target (modelling delay = N/2)
%         c) Filter design with zerophase target
%
% 13. Apply post-filter to compensate transducers in playback-chain
% 14. Reconstruct interchannel gain differences (were discarded in step 11)
% ------ plot6: inverse filter
% 15. Save compensation filter
% 16. Plot simulated compensation
% ------ plot7: compensation filtering result 2
%
%
% (C) Fabian Brinkmann, Alexander Lindau, Zora Schaerer
% TU Berlin, Audio Communication Group, 2009-2013
%
%
% [1] Norcross(2006): "Inverse Filtering Design Using a Minimal-Phase Target
%     Function from Regularization", AES preprint no. 6929.
% [2] Schaerer, Lindau (2009): "Evaluation of Equalization Methods for Binaural Signals."
%     In Proc. 126th AES Convention, Munich, Germany.
% [3] Brinkmann (2011): "Individual headphone compensation for binaural
%     synthesis." MA Thesis, TU Berlin, Berlin, Germany.
%     http://www2.ak.tu-berlin.de/~akgroup/ak_pub/abschlussarbeiten/2011/Brinkmann_MagA.pdf

% issues:
%
% TBD: sinnlose Kombinationen von hp_len, fade_in, fade_out und margin
% werden nicht ueberall abgefangen
% TBD: avoid cases: hp_len > win_len & hp_fade_out/in, if APPLY_SMOOTHING frac > 0
% TBD: hp_filter_margin ... change to automatic set up?

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

%% ------------------------------------------------ 0 set default parameter
%  some parameter were added after the initial release of AKtools. The
%  default parameters are set for downwards compatability and ensure that
%  calls the the latest version produce identical results.

if ~isfield(s, 'g_norm')
    s.g_norm = 0;
end

s_copy = s;

%% ------------------------------------------ 1 directory and data handling

% directory for saving final compensation filters
if any(s.save_filter)
    s.dest_dir = fullfile(s.save_dir,'compensation_filters');
    if exist(s.dest_dir,'dir')~=7
        mkdir(s.dest_dir)
    end
end

% make filter design documentation directory
if exist(fullfile(s.save_dir,'filterdesign_doc'), 'dir') ~= 7
    mkdir(fullfile(s.save_dir,'filterdesign_doc'))
end

%% ----------------------------------------------------- 2 read in raw data
if exist('raw_data', 'var')
    data.raw              = raw_data;
    s.num_of_measurements = size(raw_data, 2);
    s.numChannel          = size(raw_data, 3);
else
    wav_list = dir(fullfile(s.data_dir, '*.wav'));
    s.num_of_measurements = numel(wav_list);
    
    if isempty(wav_list)
        error('No wav files found in s.data_dir')
    end
    
    for wav_idx = 1:length(wav_list)
        [tmp, s.fs] = audioread(fullfile(s.data_dir,wav_list(wav_idx).name));
        
        s.numChannel = size(tmp,2);
        
        for n = 1:s.numChannel
            data.raw(:,wav_idx,n) = tmp(:,n);
        end
    end
end

clear wav_list tmp n wav_idx raw_data


%% ----------------------------------------- 3 apply pre-filter to raw data
if s.prefilter
    % load prefilter
    load(s.prefilter, 'pre_post_filter', 'fs');
    
    if fs ~= s.fs
        error('sampling rates of pre-filter and raw data do not match!')
    end
    
    data.prefilter = zeros(size(data.raw,1)+size(pre_post_filter,1), size(data.raw,2), size(data.raw,3)); %#ok<NODEF>
    
    % check if number of channles match across prefilter and raw data
    if size(pre_post_filter,2) ~= s.numChannel
        pre_post_filter = repmat(pre_post_filter(:,1), [1 s.numChannel]);
        disp('Pre-filter: Number of channels does not match raw data - first channel of pre-filter used for all channels of raw data!');
    end
    
    % filter (linear convolution)
    for m = 1:s.num_of_measurements
        for n = 1:s.numChannel
            data.prefilter(:,m,n) = fftfilt(pre_post_filter(:,n), [data.raw(:,m,n); zeros(size(pre_post_filter(:,n)))]);
        end
    end
    
    data_current = data.prefilter;
else
    data_current = data.raw;
end

% temporary save for plots
data_pre_shortening = data_current;

clear pre_post_filter m n fs


%% ---------------------------------------- 4 time allign impulse responses
% (base on time of arrival of first channel of each impulse response)
if s.time_align
    data.aligned = zeros(size(data_current));
    
    % estimate time of arrival via onset detection
    s.toa = AKonsetDetect(data_current(:,:,1), 1, s.threshold);
    
    % circshift to desired position
    if min(s.toa) > s.peak_protect
        for n = 1:s.num_of_measurements
            data.aligned(:,n,:) = circshift(data_current(:,n,:), [-s.toa(n)+s.peak_protect 0 0]);
        end
    else
        for n = 1:s.num_of_measurements
            data.aligned(:,n,:) = [zeros(s.peak_protect-s.toa(n), 1, size(data.aligned, 3)); data_current(1:end-s.peak_protect+s.toa(n),n,:)];
        end
    end
    
    data_current = data.aligned;
end

clear n

%% ------------------------------------------- 5 truncate impulse responses
% (this and previous step)

s.Nsamples = min([s.Nsamples size(data_current,1)]);
data.truncated = data_current(1:s.Nsamples,:,:);

if s.win_len
    win = cos(linspace(0, pi/2, s.win_len)').^2;
    win = repmat(win, [1 s.num_of_measurements s.numChannel]);
    data.truncated(end-s.win_len+1:end,:,:) = data.truncated(end-s.win_len+1:end,:,:) .* win;
end

data_current = data.truncated;

if s.plot_shorten
    AKf(10*s.numChannel, 20)
    set(gcf, 'Name', 'Plot 1: TOA alignment and shortening', 'numberTitle', 'off')
    for n = 1:s.numChannel
        subplot(2,s.numChannel,n)
            AKp(data_pre_shortening(:,:,n), 'et2d', 'c', [.7 .7 .7], 'xu', 'n', 'fs', s.fs)
            AKp(data_current(:,:,n), 'et2d', 'xu', 'n', 'fs', s.fs)
            title({['ETC, channel ' num2str(n)] 'raw (grey), aligned & shortened (colored)'})
        subplot(2,s.numChannel,n+s.numChannel)
            AKp(data_pre_shortening(:,:,n), 'm2d', 'c', [.7 .7 .7], 'x', [20 s.fs/2], 'fs', s.fs)
            AKp(data_current(:,:,n), 'm2d', 'x', [20 s.fs/2], 'fs', s.fs)
            title({['Magnitude response, channel ' num2str(n)]  'raw (grey), aligned & shortened (colored)'})
    end
    
    if s.do_print
        print('-dpdf', '-r300', fullfile(s.save_dir,'filterdesign_doc','plot1_toaelim_timealign_shortening'));
    end    
end

clear win n data_pre_shortening

%% -------------------------------------------------------------- 6 average
% average complex spectra or time signals 
if strcmpi(s.avg_type, '_complex')
    data.averaged = ifft(mean(fft(data_current), 2), 'symmetric');
elseif strcmpi(s.avg_type, '_time')
    data.averaged = mean(data_current, 2);
else
    error(['''' s.avg_type ''' is not valid for s.avg_type. Try ''_complex'' or ''_time'''])
end

data.averaged = squeeze(data.averaged);

% get lower and upper frequency bounds for normalization
% (will be applied in step 11)
f_bounds = s.f_norm/2^(1/(2*s.f_norm_frac));
f_bounds = [f_bounds f_bounds*2^(1/s.f_norm_frac)];
% get corresponding indicees
f_bounds = round(f_bounds/(s.fs/s.Nsamples))+1;
tmp_bounds = (f_bounds-1) * (s.fs/s.Nsamples);

if s.plot_avg
    AKf(10*s.numChannel, 20)
    set(gcf, 'Name', 'Plot 2: Averaging', 'numberTitle', 'off')
    for n = 1:s.numChannel
        subplot(2,s.numChannel,n)
            AKp(data_current(:,:,n), 'et2d', 'c', [.7 .7 .7], 'xu', 'n', 'fs', s.fs)
            AKp(data.averaged(:,n), 'et2d', 'xu', 'n', 'fs', s.fs)
            title(['ETC, channel ' num2str(n) ' (Averaged=black)'])
        subplot(2,s.numChannel,n+s.numChannel)
            AKp(data_current(:,:,n), 'm2d', 'c', [.7 .7 .7], 'x', [20 s.fs/2], 'fs', s.fs)
            AKp(data.averaged(:,n), 'm2d', 'x', [20 s.fs/2], 'fs', s.fs)
            tmp = get(gca, 'ylim');
            hold on
            plot([tmp_bounds(1) tmp_bounds(1)], [tmp(1) tmp(2)], 'r--')
            plot([tmp_bounds(2) tmp_bounds(2)], [tmp(1) tmp(2)], 'r--')
            title({['Magnitude response, channel ' num2str(n)] 'freq. range for averaging (red)'})
    end
    
    if s.do_print
        print('-dpdf', '-r300', fullfile(s.save_dir,'filterdesign_doc','plot2_averaging'));
    end    
end

data_current = data.averaged;

clear tmp tmp_bounds


%% ------------- 7 decide if regularization or smoothing is will be applied
% (based on s.inversion method)
switch s.inversion_method
    case 1 
        s.APPLY_REG    = 1;
    case 2
        s.APPLY_REG    = 1;
    case 3 
        s.APPLY_REG    = 0;
        if s.frac
            s.APPLY_SMOOTHING = 1;
        end
    case 4
        s.APPLY_REG    = 1;
    case 5
        s.APPLY_REG       = 0;
        s.APPLY_SMOOTHING = 1;
        s.frac            = s.inv_data.frac;
    case 6
        s.APPLY_REG    = 1;
end

%% ------------------------------------ 8 calculate regularization function
% (filter coefficients, impulse response and plot)
switch s.inversion_method
    % LMS inversion of mean transfer function with shelve for amplitude
    % regularization
    case 1
        for idx_c = 1:s.numChannel
            [s.b_REG(idx_c, :), s.a_REG(idx_c, :)] = AKhighshelve2(s.reg_shelve_freq(idx_c), s.fs, s.reg_shelve_gain(idx_c), 1/sqrt(2), 1/sqrt(2), 'III');
            s.b_REG(idx_c, :) = 10^(-s.reg_shelve_gain(idx_c)/20)*s.b_REG(idx_c, :);   
        end
        
    % LMS inversion of mean transfer function with mean inverse of frac. oct.
    % smoothed (optional) transfer function for amplitude regularization    
    case 2
        AVG           = fft(data_current);
        [AVG, is_even] = AKboth2singleSidedSpectrum(AVG);
        AVG_SM = zeros(size(AVG));
        % fractional octave smoothing
        if s.reg_frac
            for idx_c = 1:s.numChannel
                AVG_SM(:, idx_c) = AKfractOctSmooth(AVG(:, idx_c), 'welti', s.fs, s.reg_frac);
            end
            AVG_SM = AKsingle2bothSidedSpectrum(AVG_SM, is_even);
            AVG = AVG_SM;
        else
            AVG = abs(AVG);
        end
        % set amplitude spectrum constant below eq_data.align(1)
        if s.reg_align(1)
            f_norm_reg2 = round(s.reg_align(1)/(s.fs/s.Nsamples))+1;
            for idx_a = 1:size(AVG, 2)
                AVG(1:f_norm_reg2, idx_a) = AVG(f_norm_reg2, idx_a);
                AVG(end-f_norm_reg2+2:end, idx_a) = AVG(f_norm_reg2, idx_a);
            end
        end
        % set amplitude spectrum constant above eq_data.align(2)
        if s.reg_align(2)
            f_norm_reg2 = round(s.reg_align(2)/(s.fs/s.Nsamples))+1;
            for idx_a = 1:size(AVG, 2)
                AVG(f_norm_reg2:end-f_norm_reg2+2, idx_a) = AVG(f_norm_reg2, idx_a);
            end
        end
        % limit dynamic of inversion
        AVG = 1./AVG;
        if s.reg_inv_dyn
            for idx_a = 1:size(AVG, 2)
                inv_dyn = min(AVG(:, idx_a)) * 10^(s.reg_inv_dyn/20);
                AVG(:, idx_a) = min(AVG(:, idx_a), inv_dyn);
            end
        end
        % normalize spectrum
        for idx_a = 1:size(AVG, 2)
            AVG(:, idx_a) = AVG(:, idx_a) / max(AVG(:, idx_a));
        end
        s.b_REG = ifft(AVG, 'symmetric')';
        s.a_REG = ones(size(s.b_REG, 1), 1);
    
    % LMS inversion of mean transfer function with oct. smoothed transfer
    % function - 1/x oct. smoothed mean transfer function is used for
    % regularization    
    case 4
        AVG           = abs(fft(data_current));
        [AVG, is_even] = AKboth2singleSidedSpectrum(AVG);
        AVG_rough = zeros(size(AVG));
        % rough fractional octave smoothing
        for idx_c = 1:s.numChannel
            AVG_rough(:, idx_c) = AKfractOctSmooth(AVG(:, idx_c), 'welti', s.fs, s.reg_frac_rough);
        end
        AVG_rough = AKsingle2bothSidedSpectrum(AVG_rough, is_even);
        % fine fractional octave smoothing
        if s.reg_frac_fine
            AVG_fine = zeros(size(AVG));
            for idx_c = 1:s.numChannel
                AVG_fine(:, idx_c) = AKfractOctSmooth(AVG(:, idx_c), 'welti', s.fs, s.reg_frac_fine);
            end
            AVG_fine = AKsingle2bothSidedSpectrum(AVG_fine, is_even);
        else
            AVG_fine = AVG;
        end
        B_REG = (AVG_rough ./ AVG_fine);
        % align and normalize
        for idx_c = 1:s.numChannel
            % set amplitude spectrum constant below eq_data.align(2)
            if s.reg_align(1)
                f_norm_reg5 = round(s.reg_align(1)/(s.fs/s.Nsamples))+1;
                B_REG(1:f_norm_reg5, idx_c) = B_REG(f_norm_reg5, idx_c);
                B_REG(end-f_norm_reg5+2:end, idx_c) = B_REG(f_norm_reg5, idx_c);
            end
            % set amplitude spectrum constant above eq_data.align(2)
            if s.reg_align(2)
                f_norm_reg5 = round(s.reg_align(2)/(s.fs/s.Nsamples))+1;
                B_REG(f_norm_reg5:end-f_norm_reg5+2, idx_c) = B_REG(f_norm_reg5, idx_c);
            end
            B_REG(:, idx_c) = B_REG(:, idx_c) / max(B_REG(:, idx_c));
        end
        s.b_REG = ifft(B_REG, 'symmetric')';
        s.a_REG = ones(size(s.b_REG, 2), 1);
    
    % LMS inversion with Param. EQ used for regularization    
    case 6
        for idx_c = 1:numel(s.reg_beta)
            % initialize matrix for impulse responses of single filters
            h_reg = [];
            
            % get parametrical EQs and shelves
            for idx = 1:size(s.reg_data, 1)
                if s.reg_data{idx, 1} == idx_c
                    % PEQs
                    if strcmpi(s.reg_data{idx, 2}, 'PEQ')
                    [b, a] = AKpeq2(s.reg_data{idx, 3}, s.fs, s.reg_data{idx, 4}, s.reg_data{idx, 5}, 'III');
                    elseif strcmpi(s.reg_data{idx, 2}, 'LS')
                        if isnan(s.reg_data{idx, 5})
                            % f specifies mid gain frequency
                            [b, a] = AKlowshelve2(s.reg_data{idx, 3}, s.fs, s.reg_data{idx, 4}, 1/sqrt(2), 1/sqrt(2), 'III');
                        else
                            % f specifies -3 db cut off frequency
                            [b, a] = AKlowshelve2(s.reg_data{idx, 3}, s.fs, s.reg_data{idx, 4}, 1/sqrt(2), 1/sqrt(2), 'I');
                        end
                    elseif strcmpi(s.reg_data{idx, 2}, 'HS')
                        if isnan(s.reg_data{idx, 5})
                            % f specifies mid gain frequency
                            [b, a] = AKhighshelve2(s.reg_data{idx, 3}, s.fs, s.reg_data{idx, 4}, 1/sqrt(2), 1/sqrt(2), 'III');
                        else
                            % f specifies -3 db cut off frequency
                            [b, a] = AKhighshelve2(s.reg_data{idx, 3}, s.fs, s.reg_data{idx, 4}, 1/sqrt(2), 1/sqrt(2), 'I');
                        end
                    end
                    h_reg(:, end+1) = impz(b, a, s.Nsamples);
                end
            end
            
            % get spectra and combine transfer functions
            B_REG = fft(h_reg);
            B_REG = prod(abs(B_REG), 2);
            
            % limit dynamic of inversion
            if s.reg_dyn
                hp_dyn = min(B_REG) * 10^(s.reg_dyn(idx_c)/20);
                B_REG  = min(B_REG, hp_dyn);
                B_REG  = B_REG / max(B_REG);
            end
            
            % get impulse Response
            s.b_REG(idx_c,:) = ifft(B_REG, 'symmetric')';
            s.a_REG(idx_c,1) = 1;
        end
    otherwise % cases 3 and 5 don't use amplitude regularization
        s.b_REG(1:s.numChannel,1) = 1;
        s.a_REG(1:s.numChannel,1) = 1;
end

% get impulse responses from filter coefficients
data.regularization = zeros(s.Nsamples, s.numChannel);
for n = 1:s.numChannel
    % check impulse response length
    if impzlength(s.b_REG(n,:),s.a_REG(n,:))>s.Nsamples
        warning('reg_inv:regularization_function', 'Filter length given by Nsmaples is to short to realize regularization function properly.');
    end
    % create impulse response
    data.regularization(:,n) = impz(s.b_REG(n,:),s.a_REG(n,:),s.Nsamples);
end
s.regularization = data.regularization;

%plot 3: reg-functions
if s.plot_regularization && s.inversion_method~=3 && s.inversion_method~=5
    AKf(20*s.numChannel, 12)
    set(gcf, 'Name', 'Plot 3: Regularization function', 'numberTitle', 'off')
    for n = 1:s.numChannel
        subplot(1,s.numChannel,n)
            AKp(data.regularization(:,n), 'm2d', 'fs', s.fs, 'dr', 40)
            title(['Regularization, channel ' num2str(n)])
    end
    
    if s.do_print
        print('-dpdf', '-r300', fullfile(s.save_dir,'filterdesign_doc','plot3_regularization'));
    end    
end

clear idx_c AVG AVG_SM n f_norm_reg2 idx_a inv_dyn AVG_fine AVG_rough B_REG f_norm_reg5 h_reg b a gain hp_dyn idx is_even

%% -------------------------------------------- 9 calculate target function
% (replicating BKamp-Preset-Targets, or arbitrary buttworth bandpass)
if any(s.target)
    
    if ischar(s.target)
        % design a suitable low pass characteristic for the target function
        % based on a BK211 measured on FABIAN HATS
        
        % or derivation of these target function parameters see
        % Erbes, Vera (2013): ?BKsystem Reference Manual?, audio communication
        % group, TU Berlin
        
        lp_order = 2; % was 6(?) in Erbes 2013, had to be adapted for more well-behaved inversion behavour
        lp_fcut = 16400;
        
        % design a suitable high pass characteristic for the target function
        % (replicating  the BK211's low-end roll-off as realized through
        % respective BKamp presets
        switch s.target
            case 'fr'
                hp_fcut        = 59;        % was 61.7Hz in Erbes (2013), had to be adapted for more well-behaved inversion behavour
                hp_order       = 24/6;      % was 18/6 in Erbes (2013), had to beadapted for more well-behaved inversion behavour
                s.target       = [4 59 2 16400];
            case 'sub50'
                hp_fcut        = 56.5;      % was 58.5 Hz in Erbes (2013), had to be adapted for more well-behaved inversion behavour
                hp_order       = 36/6;
                s.target       = [6 56.5 2 16400];
            case 'sub85'
                hp_fcut        = 86.1;
                hp_order       = 36/6;
                s.target       = [6 86.1 2 16400];
            case 'sub120'
                hp_fcut        = 106.3;
                hp_order       = 30/6;
                s.target       = [5 106.3 2 16400];
            case 'sub150'
                hp_fcut        = 180.3;
                hp_order       = 36/6;
                s.target       = [6 180.3 2 16400];
        end
    else
        hp_order = s.target(1);
        hp_fcut  = s.target(2);
        lp_order = s.target(3);
        lp_fcut  = s.target(4);
    end
    % assume butterworth characteristics for convenience
    if hp_order
        spec_hp = fdesign.highpass('N,F3dB', hp_order, hp_fcut, s.fs);
        h_hp    = design(spec_hp, 'butter');
    end
    
    if lp_order
        spec_lp = fdesign.lowpass('N,F3dB', lp_order, lp_fcut, s.fs);
        h_lp    = design(spec_lp, 'butter');
    end
        
    % design a target function as a cascade of hp and lp filter
    if hp_order && lp_order
        h_target = dfilt.cascade(h_hp,h_lp);
    elseif hp_order
        h_target = h_hp;
    elseif lp_order
        h_target = h_lp;
    else
        error('all filter orders were 0. At least one filter order must be 1 or higher.');
    end
    
    % test if filter is stable
    if ~isstable(h_target)
        error('Target function is not stable! Try again with lower orders or changed cut-off frequenies.')
    end
    
    if impzlength(h_target)>s.Nsamples
        warning('reg_inv:target_function', 'Filter length given by Nsmaples is to short to realize target function properly.');
    end
    
     % check minimum phase
    if isminphase(h_target)
        fprintf('\nRaw target function is minimum phase.\n')
    else
        fprintf('\nRaw target function is not minimum phase.\n')
    end
    
    % parameter are stored for results plot
    s.target_param = s.target;
    
    % create impulse response (AKffdInverseFilter expects it to be saved in s.target)
    data.target    = impz(h_target,s.Nsamples);
    s.target       = data.target;
else
    % target function will be dirac in this case
    s.target    = zeros(s.Nsamples,1);
    s.target(1) = 1;
    data.target = s.target;
    % cut-off frequencies are needed for plot6
    s.target_param = [0 20 0 s.fs/2];
end

% plot 4: show target function
if s.plot_target
    AKf(20,20)
    set(gcf, 'Name', 'Plot 4: Target function', 'numberTitle', 'off')
    subplot(2,1,1)
    AKp(data.target,'et2d', 'x', [-10 s.Nsamples/s.fs*1000], 'fs', s.fs);
    subplot(2,1,2)
    AKp(data.target,'m2d', 'fs', s.fs);
    
    if  s.do_print
        print('-dpdf', '-r300', fullfile(s.save_dir,'filterdesign_doc','plot4_target function'));
    end
end

clear hp_order hp_fcut lp_order lp_fcut spec_hp h_hp spec_lp h_lp h_target k po ze

%%  10 Fract. octave smoothing and min-phase of averaged transfer functions
% (optional)
if s.APPLY_SMOOTHING && s.frac
    AVG           = fft(data_current);
    [AVG, is_even] = AKboth2singleSidedSpectrum(AVG);
    AVG_SM        = zeros(size(AVG));
    for idx_c = 1:s.numChannel
        AVG_SM(:, idx_c) = AKfractOctSmooth(AVG(:, idx_c), 'welti', s.fs, s.frac);
    end
    AVG_SM = AKsingle2bothSidedSpectrum(AVG_SM, is_even);
    data.average_sm = ifft(AVG_SM, 'symmetric');
    disp('Applying min-phase after smoothing averaged transfer function')
    data.average_sm = AKphaseManipulation(data.average_sm, s.fs, 'min_phase', s.Nfft_double);
    
    data_current = data.average_sm;
end

clear AVG AVG_SM is_even idx_c

%% --------------------------------------------- 11  normalize avg @ f_norm
% and upper frequency bounds for normalization were calculated in step 6

mag_AVG = abs(fft(data_current));
s.norm = 1./mean(mag_AVG(f_bounds(1):f_bounds(2),:)) * 10^(s.g_norm/20);
for i=1:s.numChannel
    data.average_norm(:,i) = data_current(:,i)*s.norm(i);
end

data_current = data.average_norm;
s.avg        = data_current;

clear f_bounds i mag_AVG

%% --------------------------------------- 12 calculate compensation filter

if s.inversion_method ~= 5
    s = AKffdInverseFilter(s);
else
    AVG = abs(fft(data_current));
    for idx_c = 1:s.numChannel
        % inverse
        AVG(:, idx_c) = 1 ./ AVG(:, idx_c);
        % constrain inversion dynamic
        if s.inv_data.inv_dyn
            inv_dyn = min(AVG(:, idx_c)) * 10^(s.inv_data.inv_dyn/20);
            AVG(:, idx_c) = min(AVG(:, idx_c), inv_dyn);
        end
        % apply target function
        if s.target
            TARGET = abs(fft(data.target));
            AVG(:, idx_c) = AVG(:, idx_c) .* TARGET;
        end
        s.compensation_filter(:, idx_c) = ifft(AVG(:, idx_c), 'symmetric');
        % get desired phase (lin, min, zero)
        s.compensation_filter(:, idx_c) = AKphaseManipulation(s.compensation_filter(:, idx_c), s.fs, s.phase_type, 0);
    end
end

% save compensation filter to data struct
data.compensation_filter_raw = s.compensation_filter;

%-------------------------------------------------------------------------%
% Plot 5: inverse filtering result
if  s.plot_compensation1
    
    AKf(20*s.numChannel, 30);
    set(gcf, 'Name', 'Plot 5: Results (averaged)', 'numberTitle', 'off')
    
    for ch = 1:s.numChannel
        
        % estimate compensated system by convolution of inverse filter and
        % input data (linear convolution)
        comp = fftfilt(data.compensation_filter_raw(:,ch),[data_current(:,ch); zeros(size(data_current(:,ch)))]);
        
        % subplot 1: averaged data, compensation filter and regularization function
        subplot(4,s.numChannel,ch)
        if s.APPLY_REG && s.reg_beta(ch) ~= 0
            AKp(data.regularization(:,ch), 'm2d', 'fs', s.fs, 'c', 'b')
        end
            AKp(data.compensation_filter_raw(:,ch), 'm2d', 'fs', s.fs, 'c', 'r')
            AKp(data_current(:,ch), 'm2d', 'fs', s.fs, 'x', [20 s.fs/2])
        if s.APPLY_REG && s.reg_beta(ch) ~= 0
            title(['Ch ' num2str(ch) ': Input (black), filter (red), regularization (blue)'])
        else
            title(['Ch ' num2str(ch) ': Input (black), filter (red)'])
        end
        % subplot 2: target and compensated system
        subplot(4,s.numChannel,ch+s.numChannel)
            AKp(data.target, 'm2d', 'fs', s.fs, 'c', 'g')
            AKp(comp, 'm2d', 'fs', s.fs, 'c', 'r', 'dr', [-10 5], 'x', [20 s.fs/2])
            title(['Ch ' num2str(ch) ': target (green), compensation result (red)'])        
        % subplot 3: etc of compensated system
        subplot(4,s.numChannel,ch+2*s.numChannel)
            AKp(comp, 'et2d', 'fs', s.fs)
            title(['Ch ' num2str(ch) ': ETC of compensation result'])
                % subplot 3: etc of compensated system
        subplot(4,s.numChannel,ch+3*s.numChannel)
            AKp(data.compensation_filter_raw(:,ch), 'et2d', 'fs', s.fs)
            title(['Ch ' num2str(ch) ': ETC of compensation filter'])
        
        if ch == s.numChannel
            if s.do_print
                print('-dpdf', '-r300',fullfile(s.save_dir,'filterdesign_doc','plot5_inverse_filter_design_results'))
            end
        end
    end
end

clear ch comp

%% --------------------------------------------------- 13 apply post-filter
if s.postfilter
    % load postfilter
    load(s.postfilter, 'pre_post_filter', 'fs');
    
    if fs ~= s.fs
        error('sampling rates of post-filter and raw data do not match!')
    end
    
    data.compensation_filter = zeros(size(data.compensation_filter_raw,1)+size(pre_post_filter,1), size(data.compensation_filter_raw,2));
    
    % check if number of channles match across postfilter and raw data
    if size(pre_post_filter,2) ~= s.numChannel
        pre_post_filter = repmat(pre_post_filter(:,1), [1 s.numChannel]);
        disp('Pre-filter: Number of channels does not match raw data - first channel of pre-filter used for all channels of raw data!');
    end
    
    % filter (linear convolution)
    for n = 1:s.numChannel
        data.compensation_filter(:,n) = fftfilt(pre_post_filter(:,n), [data.compensation_filter_raw(:,n); zeros(size(pre_post_filter(:,n)))]);
    end
else
    data.compensation_filter = data.compensation_filter_raw;
end

clear pre_post_filter n fs

%% --------------------------- 14 reconstruct interchannel gain differences
% (They were eliminated in step 11)
% As the microphone sensitivities have been compensated before, and the
% ADDA-path characteristics should be near identical between channels
% (when  assuring identical mic-preamp gains!!!) the offset in s.norm
% should be small (~+-0.5 dB) and reflect differences in left and right
% headphone speaker.

for n = 2:s.numChannel
    interchannel_ratio = s.norm(n)/s.norm(1);
    data.compensation_filter(:,n) = data.compensation_filter(:,n) * interchannel_ratio;
end

clear n interchannel_ratio

%-------------------------------------------------------------------------%
% Plot 6: inverse filtering result
if  s.plot_inverseFilter
    
    AKf(20*s.numChannel, 12);
    set(gcf, 'Name', 'Plot 6: Inverse/compensation filter', 'numberTitle', 'off')
    
    for ch = 1:s.numChannel
        
        subplot(1,s.numChannel,ch)
        AKp(data.compensation_filter(:,ch), 'm2d', 'fs', s.fs, 'x', [20 s.fs/2], 'dr', [-20 20])
        title(['Ch ' num2str(ch) ': Inverse filter'])
        
        if ch == s.numChannel
            if s.do_print
                print('-dpdf', '-r300',fullfile(s.save_dir,'filterdesign_doc','plot6_InverseFilter'))
            end
        end
    end
end


%% -------------------------------------------- 15 save compensation filter
% (attenuate by -s.att dB before)
data.compensation_filter = data.compensation_filter * 10^(-abs(s.att)/20);

for nn = 1:numel(s.save_filter)
    switch s.save_filter(nn)
        case 1  %create fWonder format wav-files/folder structure
            if s.numChannel ~= 2
                error('Saving in fWonder format only works for data with 2 channels');
            end
            
            % check directories / create new
            if exist(fullfile(s.dest_dir,'source1'), 'dir')~=7
                mkdir(fullfile(s.dest_dir,'source1'))
            end
            if exist(fullfile(s.dest_dir,'source2'), 'dir')~=7
                mkdir(fullfile(s.dest_dir,'source2'))
            end
            
            % blank one side of binaural filter (fwonder requirement, as each binaural "source" has only a singular input)
            blank = zeros(size(data.compensation_filter(:, 1)));
            left_wav  = cat(2, data.compensation_filter(:, 1), blank);
            right_wav = cat(2, blank, data.compensation_filter(:, 2));
            
            %wavwrite(left_wav, s.fs, s.nbits, (fullfile(s.dest_dir,'source1',s.save_name)));
            %wavwrite(right_wav, s.fs, s.nbits, (fullfile(s.dest_dir,'source2',s.save_name)));
            
            audiowrite((fullfile(s.dest_dir,'source1',[s.save_name '.wav'])), left_wav, s.fs, 'BitsPerSample', s.nbits);
            audiowrite((fullfile(s.dest_dir,'source2',[s.save_name '.wav'])), right_wav, s.fs, 'BitsPerSample', s.nbits);
            
        case 2 % singular multichannel-wav file
            % wavwrite(data.compensation_filter, s.fs, s.nbits, (fullfile(s.dest_dir,s.save_name)));
            audiowrite((fullfile(s.dest_dir,[s.save_name '.wav'])), data.compensation_filter, s.fs, 'BitsPerSample', s.nbits);
        case 3
            compensation_filter = data.compensation_filter; %#ok<NASGU>
            save(fullfile(s.dest_dir, s.save_name), 'compensation_filter');
            clear compensation_filter
        case {4 5}
            
            % get convention and write convention specific meta data
            if s.save_filter(nn) == 4
                H                  = SOFAgetConventions('GeneralFIR');
                H.ReceiverPosition = [0 .0875 0; 0 -.0875 0];
            else
                H                           = SOFAgetConventions('SimpleFreeFieldHRIR');
                H.GLOBAL_DatabaseName       = s.save_name;
                H.GLOBAL_ListenerShortName  = s.save_name;
            end
            
            % common meta data
            H.GLOBAL_ApplicationName    = 'Matlab';
            H.GLOBAL_ApplicationVersion = version;
            H.GLOBAL_Comment            = 'The filter was designed with AKtools (www.ak.tu-berlin.de/aktools)';
            H.GLOBAL_References         = 'Fabian Brinkmann and Stefan Weinzierl: "AKtools---An Open Software Toolbox for Signal Acquisition, Processing, and Inspection in Acoustics." In: 142nd AES Convention, e-Brief 309, Berlin , Germany, 2017.';
            H.GLOBAL_Title              = 'Headphone filter';
            
            % it would be nice to save s_copy to the SOFA file, but User
            % defied structs are not supported so far
            if isfield(s_copy, 'save_filter')
                s_copy = rmfield(s_copy, {'save_filter' 'save_dir' 'script_dir' 'plot_shorten' 'plot_avg' 'plot_regularization' 'plot_target' 'plot_compensation1' 'plot_inverseFilter' 'plot_compensation2' 'do_print'});
            end
            
            % save
            if s.save_filter(nn) == 4
                % FIR data
                H.Data.SamplingRate = s.fs;
                H.Data.IR           = zeros(1, 2, s.Nsamples);
                H.Data.IR(:,1,:)    = data.compensation_filter(:,1);
                H.Data.IR(:,2,:)    = data.compensation_filter(:,2);
                H.Data.Delay        = zeros(1, 2);
                % save
                H = SOFAsave(fullfile(s.dest_dir, [s.save_name '.sofa']), H); %#ok<NASGU>
            else
                % left ear FIR data
                H.Data.SamplingRate = s.fs;
                H.Data.IR           = zeros(1, 2, s.Nsamples);
                H.Data.IR(:,1,:)    = data.compensation_filter(:,1);
                % save left ear
                H = SOFAsave(fullfile(s.dest_dir, [s.save_name '_left.sofa']), H);
                % right ear FIR data
                H.Data.IR           = zeros(1, 2, s.Nsamples);
                H.Data.IR(:,2,:)    = data.compensation_filter(:,2);
                % save right ear
                H = SOFAsave(fullfile(s.dest_dir, [s.save_name '_right.sofa']), H); %#ok<NASGU>
            end
    end
end

clear nn blank left_wav right_wav params

%% ----------------------------------------- 16 plot simulated compensation
% (compensation filter without post-filter applied to raw transfer functions)

if s.plot_compensation2
      
    % allocate memory
    comp_filter = data.compensation_filter_raw;
    comp_result = zeros(size(data.truncated,1)+size(comp_filter,1), size(data.truncated,2), size(data.truncated,3));
    
    for ch = 1:s.numChannel
        % apply interchannel gains and normalization values to compensation filter to match raw data
        % (they are not included in the raw compensation filter. they of course are included in the final compensation filter)
        if ch > 1
            interchannel_ratio = s.norm(ch)/s.norm(1);
            % remove interchannel gain diffs
            comp_filter(:,ch) = data.compensation_filter_raw(:,ch) * interchannel_ratio;
        end

        for n = 1:s.num_of_measurements
            % apply compensation filter and normalization gains to truncated data
            % (apply normalization gain of first channel to all filters, the
            % remaining was done by the interchannel ratio in the previous step)
            comp_result(:,n,ch) = fftfilt(comp_filter(:,ch), [data.truncated(:,n,ch)*s.norm(1); zeros(size(comp_filter(:,ch)))]);
        
        end
            
    end
    
    data.compensation_result = comp_result;
    
    % ERB-filter-bank error in dB
    [ERBerror, f_c] = AKerbError(reshape(comp_result(1:s.Nsamples,:,:), [s.Nsamples s.num_of_measurements*s.numChannel]), data.target, [s.target_param(2) s.target_param(4)], s.fs);
    ERBerror        = reshape(ERBerror, [numel(f_c) s.num_of_measurements s.numChannel]);
    
    clear n
    
    % discard error values above cut off-frequency
    n = find(f_c > s.target_param(4), 1, 'last');
    if n
        ERBerror = ERBerror(n:end,:,:);
        f_c      = f_c(n:end);
    end
    clear n
    
    % save to ouput
    data.ERB = ERBerror;
    data.f_c = f_c;
    
    % plot
    AKf(10*s.numChannel,20)
    set(gcf, 'Name', 'Plot 7: Compensation results', 'numberTitle', 'off')
    for ch = 1:s.numChannel
        % plot equalized raw data
        subplot(2,s.numChannel, ch)
            AKp(comp_result(1:s.Nsamples,:,ch), 'm2d', 'fs', s.fs, 'x', [20 s.fs/2], 'dr', [-10 10])
            title(['Ch ' num2str(ch) ': Compensated raw data'])
        % plot ERB error
        subplot(2,s.numChannel, ch+s.numChannel)
            hold on
            plot(f_c, ERBerror(:,:,ch), 'xk')
            set(gca, 'XTick', [100 1e3 10e3 20e3 40e3 80e3], 'XTickLabel', {'100' '1k' '10k' '20k' '40k' '80k'}, 'XScale', 'log', ...
                'YTick', -5:5)
            axis([20 s.fs/2 -5 5])
            grid on
            box on
            title(['Ch ' num2str(ch) ': Error in ERB-filters'])
            xlabel('frequency in Hz'); ylabel('amplitude in dB')
    end
    
    if s.do_print
        print('-dpdf', '-r300', fullfile(s.save_dir,'filterdesign_doc',['plot7_reliab_inv_method_' num2str(s.inversion_method) '_len_' num2str(s.Nsamples) '_marging_' num2str(s.margin) ]));
    end    
end

clear n ch err_tmp err comp_filter comp_result
