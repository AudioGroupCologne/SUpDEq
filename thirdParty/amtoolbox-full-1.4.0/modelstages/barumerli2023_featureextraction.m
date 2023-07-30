function [varargout] = barumerli2023_featureextraction(sofa_obj, varargin)
%BARUMERLI2023_FEATUREEXTRACTION extract binaural and monaural cues from SOFA object
%
%   Usage: [template, target] = barumerli2023_featureextraction(sofa_obj)
%
%   Input parameters:
%       sofa_obj: Struct in SOFA format with DTFs
%
%   Output parameters:
%     template : internal templates with specific feature points
%     target   : (optional) preprocessed target struct
%
%   BARUMERLI2023_FEATUREEXTRACTION(...) computes temporally integrated
%   spectral magnitude profiles, itd and ild.
%
%
%   Additional input parameters: 
%
%     'fs'          Sampling rate. Default: 48kHz. 
%
%     'flow'        Low frequency for auditory bands. Default: 700Hz
% 
%     'fhigh'       High frequency for auditory bands. Default: 18kHz
% 
%     'space'       Auditory bands spacing. Default: 1 ERB
%         
%     'targ_az'     column vector to select the binaural stimulus in the 
%                   templates with the specified azimuth in degree shall 
%                   be used as target. Use in combination with targ_el. 
%                   Default: [].
%         
%     'targ_el'     column vector to select the binaural stimulus in the 
%                   templates with the specified elevation in degree shall 
%                   be used as target. Use in combination with targ_az. 
%                   Default: [].
%                     
%     'source_ir'   Specify a custom sound source to be convolved with HRTFs. 
%                   The default consider broadband noise. Default: []
%         
%     'source_fs'   Sampling rate of the source. Default: 0 Hz
% 
%   
%   Further, cache flags (see amt_cache) and plot flags can be specified:
%
%     'template'          Compute only templates.
%
%     'target'            Compute only targets.
%
%     'pge'               Use spectral gradients as monaural cues. 
%
%     'dtf'               Use spectral amplitudes as monaural cues. 
%
%     'monaural_none'     Do not use monaural cues. 
%
%     'reijniers'         Compute feature space as in Reijniers et al. 2014
%       
%     'source'            Use sound source provided with the parameter source_ir 
%                         to compute targets.
%
%
%
%
%
%   See also: barumerli2023
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/barumerli2023_featureextraction.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: MATLAB SOFA M-STATISTICS M-Control M-Signal
%   #Author: Roberto Barumerli (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   References: barumerli2022 


    definput.import={'amt_cache', 'barumerli2023_featureextraction'};
    definput.keyvals.fs = sofa_obj.Data.SamplingRate;
    
    [flags, kv]  = ltfatarghelper({}, definput, varargin);
    % some checks... please if you need to change them be careful... 
    assert(~xor(flags.do_source, ~isempty(kv.source_ir)), ...
        'Please add the flag source to enable the source computation!')
    assert(xor(flags.do_all, flags.do_template | flags.do_target))
    if(flags.do_source)
        assert(logical(flags.do_target), ...
        'If you need to add a source, please specify the target flag!')
    end
    
    %% extract coordinates
    % Get directions from SOFA file
    coords = barumerli2023_coordinates(sofa_obj);
    % NOTE: assume position on a sphere with radius of 1 meter
    % if you change it there will be a problem with the cached points
    % for the sphere interpolation! 
    coords.normalize_distance();
    
    % do_all, do_template and do_target allow to reduce the number of times
    % that the features require to be computed
    
    if (flags.do_all || flags.do_template)
        features = local_computefeatures(sofa_obj, coords, kv, flags);
    end
    
    %% TEMPLATE
    if flags.do_template || flags.do_all
        %% SPHERICAL HARMONIC INTERPOLATION
        template.monaural = [];
        
        if ~strcmp(flags.feature_monaural, 'none')
            % split monaural
            pl = size(features.monaural, 2);
            monaural_left = features.monaural(:, 1:pl/2);
            monaural_left = local_resamplefeatures(monaural_left, features.coords);
            monaural_right = features.monaural(:, (pl/2+1):end);
            monaural_right = local_resamplefeatures(monaural_right, features.coords);

            template.monaural = [monaural_left, monaural_right];
        end
        
        [template.itd, template.coords] = ...
            local_resamplefeatures(features.itd, features.coords);
        template.ild = ...
            local_resamplefeatures(features.ild, features.coords);

        template.fc = features.fc;
    end

    %% TARGET
    if flags.do_target
        target = local_computefeatures(sofa_obj, coords, kv, flags);
    elseif flags.do_all 
        assert(logical(flags.do_source_broadband))
        assert(exist('template', 'var') == 1)
        
        target = template;
        if ~isempty(kv.targ_az)
            coords_search = barumerli2023_coordinates([kv.targ_az, kv.targ_el, ones(size(kv.targ_el))], 'spherical');
            [coords_new, idx] = extract_directions_from_coords(target.coords, coords_search);
            target.coords = coords_new;
            target.itd = target.itd(idx,:);
	    
            if ~isempty(target.ild)
                target.ild = target.ild(idx,:);
            end
            
            if ~isempty(target.monaural)
                target.monaural = target.monaural(idx,:);
            end
        end
    end

    % output parameters
    if flags.do_template
        varargout{1} = template;
    elseif flags.do_target
        varargout{1} = target;
    else
        varargout{1} = template;
        varargout{2} = target;
    end
end

function [feature] = local_computefeatures(sofa_obj, coords, kv, flags)  
    % normalize HRTF 
    sofa_frontal = local_extractdirections(sofa_obj, barumerli2023_coordinates([0,0,1], 'spherical'));
    sofaFData = sofa_frontal.Data.IR(1, :, :);
    sofa_obj.Data.IR = sofa_obj.Data.IR ./ (max(abs(sofaFData(:)))+eps);

    stimulus = sofa_obj.Data.IR;
    
    if flags.do_target
        if ~isempty(kv.targ_az)
            coords_search = barumerli2023_coordinates([kv.targ_az, kv.targ_el, ones(size(kv.targ_el))], 'spherical');
            [sofa_obj, coords] = local_extractdirections(sofa_obj, coords_search);
        end
        stimulus = sofa_obj.Data.IR;
        
        if flags.do_source
            % normalize sound source
            kv.source_ir = kv.source_ir ./ max(abs(kv.source_ir));
            stimulus = local_convolvesource(sofa_obj.Data.IR, sofa_obj.Data.SamplingRate, kv.source_ir, kv.source_fs);
        end
    end
        
    %% compute LATERAL
    % parameters to transform into the jnd scale
    % check Reijniers2014 for these magic numbers
    a = 32.5e-6;
    b = 0.095;
    
    % ITD
    itd = itdestimator(stimulus, 'MaxIACCe', 'fs', kv.fs);    
    % transform into the jnd scale - check Reijniers2014
    itd = sign(itd) .* ((log(a + b * abs(itd)) - log(a)) / b); 
    
    % ILD
    ild = (mag2db(squeeze(rms(stimulus(:,1,:), 'dim', 3))) - ...
                    mag2db(squeeze(rms(stimulus(:,2,:),'dim', 3))));
    
    feature.itd = itd;
    feature.ild = ild;
    
    %% compute POLAR
    % compute spectral analysis
    [dtf, fc] = local_spectralanalysis(stimulus, kv);
    for ch=1:size(dtf, 2)
        for side=1:2
            dtf(:,ch,side,:) = sqrt(max(squeeze(dtf(:,ch,side,:)),0));
        end
    end
    % Averaging over time (RMS)
    dtf = (rms(dtf, 'dim', 4)); 
    
    % convert in dB
    dtf = 20*log10(dtf);
    
    if strcmp(flags.feature_monaural, 'dtf')
        feature.monaural = dtf;
    elseif strcmp(flags.feature_monaural, 'pge')
        % compute Positive Gradient Extraction (PGE)
        pge = zeros(size(dtf, 1), size(dtf, 2)-1, size(dtf, 3));
        for i=1:size(dtf,1)
            [pge(i, :, :), gfc] = ...
                baumgartner2014_gradientextraction(squeeze(dtf(i,:,:)), fc);
        end
        fc = gfc;
        feature.monaural = pge;
    elseif strcmp(flags.feature_monaural, 'monaural_none')
        feature.monaural = [];
    else
        error(['The monaural feature ', ...
                flags.monaural_feature, ' was not recognized'])
    end
    
    feature.fc = fc;
    
    %% combine SPECTRAL
    if strcmp(flags.feature_monaural, 'monaural_none')
        feature.monaural = [];
    elseif flags.do_reijniers
        feature.ild = [squeeze(dtf(:,:,1)) - squeeze(dtf(:,:,2))];
        feature.monaural = [squeeze(dtf(:,:,1)) + squeeze(dtf(:,:,2))];
        warning(['Calibrations before 2021-12-31 are not valid', ...
            'since ild and monaural have been flipped'])
    elseif flags.do_dtf || flags.do_pge
        monaural_idx = find(fc>kv.monoaural_bw(1) & fc<kv.monoaural_bw(2));
        assert(~isempty(monaural_idx), 'monoaural frequency bands empty')
        feature.fc = fc(monaural_idx);
        feature.monaural = reshape(feature.monaural(:,monaural_idx,:), ...
                                size(feature.monaural(:,monaural_idx,:), 1),[]);
    else
        error(['The monaural feature was not recognized'])
    end
    
    feature.coords = coords;
end

function [sofa_obj, coords_new, idx] = local_extractdirections(sofa_obj, coords_search)
    coords = barumerli2023_coordinates(sofa_obj);
    coords.normalize_distance();

    [coords_new, idx] = extract_directions_from_coords(coords, coords_search);

    sofa_obj.API.('M') = length(idx);
    sofa_obj.Data.IR = sofa_obj.Data.IR(idx,:,:);
    sofa_obj.SourcePosition = sofa_obj.SourcePosition(idx,:);
end
    

function [coords_new, idx] = extract_directions_from_coords(coords, coords_search)
    if(coords_search.count_pos() ~= 0)
        % TODO: warning... some points are avoided because of numerical
        [idx, coords_new] = coords.find_positions(coords_search);
        
        if(numel(idx) ~= coords_search.count_pos())
            amt_disp(sprintf('Requested HRTF''s points: %i\nFound: %i', ...
                coords_search.count(), numel(idx)))
        end

    else
        coords_new = coords;
        idx = 1:coords_search.count_pos();
    end
end

function stimulus = local_convolvesource(hrir, hrir_fs, source, source_fs)
    if source_fs <= 0
        error('source_fs is zero')
    end
    if source_fs ~= hrir_fs
        fsgcd = gcd(hrir_fs, source_fs);
        source = resample(source, hrir_fs/fsgcd, source_fs/fsgcd);
    end
    
    stimulus = zeros(size(hrir,1), ...
                     size(hrir,2), ...
                     size(hrir,3) + length(source) - 1);
    for i = 1:size(hrir, 1)
        stimulus(i,:,:) = lconv(squeeze(hrir(i,:,:))',source)';
    end
end

function [dtf, fc] = local_spectralanalysis(stimulus, kv)
    % this function expect the stimulus organized as 
    % [directions, ear_channel, time_index]
    assert(size(stimulus, 2) == 2);
    [dir_len, ear_len, time_len] = size(stimulus);
    dir_idx = 1;
    ear_idx = 2;
    time_idx = 3;
    
    % permute in order to use ufilterbankz
    stimulus = permute(double(stimulus),[time_idx, dir_idx, ear_idx]);
    % pad to account for longer filters in the filterbank
    pad_len = 0.05; % secs
    pad_mat = zeros(pad_len*kv.fs - time_len, dir_len, ear_len);
    stimulus = cat(1, stimulus, pad_mat);
    
    % compute templates features
    if kv.space == 1 % Standard spacing of 1 ERB
      [dtf,fc] = auditoryfilterbank(stimulus(:,:), kv.fs, 'flow', ...
                                            kv.flow, 'fhigh', kv.fhigh);
    else
      fc = audspacebw(kv.flow, kv.fhigh, kv.space, 'erb');
      [bgt,agt] = gammatone(fc, kv.fs, 'complex');
      % channel (3rd) dimension resolved
      dtf = 2*real(ufilterbankz(bgt,agt, stimulus(:,:)));  
    end
    
    % restore 2 channels
    dtf_size = size(stimulus);
    dtf = reshape(dtf,[dtf_size(1),length(fc),dtf_size(2),dtf_size(3)]);
    
    dtf = permute(dtf, [3 2 4 1]);
end

function [feature_interp, coords_interp] = local_resamplefeatures(feature, coords)
    if isempty(feature)
        feature_interp = zeros(0);
        coords_interp = zeros(0);
        return
    end
    
    % sample uniformly over sphere with N is number of directions
    % NOTE: amt_cache('get', 'dirs') contains the sampled point on a unitary
    % sphere
    %dirs = amt_cache('get', 'dirs');
    dirs = amt_load('barumerli2023','dirs.mat');

    % remove the points from the unitary sphere below HRTF lowest elevation
    coords_init = coords.return_positions('cartesian');
    dirs = dirs.cache.value;
    dirs = dirs(dirs(:,3) > min(coords_init(:, 3)),:);
    coords_interp = barumerli2023_coordinates(dirs, 'cartesian');

    %% interpolate at uniformly distributed directions and update feature
    % calculate spherical harmonic coefficients of H and itd, using tikonov regularization
    sh_order = 15; % spherical harmonic order 
    Y_N_tik = local_SH(sh_order, coords); 
    
    % calculate SH coefficients of H and ITD, using tikonov regularization
    lambda = 4.;
    SIG = eye((sh_order+1)^2);
    SIG(1:(2+1)^2,1:(2+1)^2) = 0;
    
    % interpolate at uniformly distributed directions and update
    Y_N_interp = local_SH(sh_order, coords_interp); 
    
    for c = 1:size(feature, 3)
        c_feat = (Y_N_tik'*Y_N_tik+lambda*SIG)\Y_N_tik'*squeeze(feature(:,:,c));
        feature_interp(:,:,c) = Y_N_interp*c_feat;
    end
end

function Y_N = local_SH(N, coords)
% calculate spherical harmonics up to order N for directions dirs [azi ele;...] (in radiant)
% 
    dirs = coords.return_positions('spherical');
    dirs = [deg2rad(dirs(:,1)), deg2rad(dirs(:,2))];
    
    N_dirs = size(dirs, 1);
    N_SH = (N+1)^2;
	dirs(:,2) = pi/2 - dirs(:,2); % convert to inclinations

    assert(N_SH < N_dirs, ...
        ['Spherical harmonics: beware that the number of provided ',...
        'coordinates is too low to obtain a precise interpolation'])
    
    Y_N = zeros(N_SH, N_dirs);

	  % n = 0
	Lnm = legendre(0, cos(dirs(:,2)'));
	Nnm = sqrt(1./(4*pi)) * ones(1,N_dirs);
	CosSin = zeros(1,N_dirs);
	CosSin(1,:) = ones(1,size(dirs,1));
	Y_N(1, :) = Nnm .* Lnm .* CosSin;
	
	  % n > 0
	idx = 1;
    for n=1:N
        
        m = (0:n)';            

		Lnm = legendre(n, cos(dirs(:,2)'));
		condon = (-1).^[m(end:-1:2);m] * ones(1,N_dirs);
		Lnm = condon .* [Lnm(end:-1:2, :); Lnm];
		
		mag = sqrt( (2*n+1)*factorial(n-m) ./ (4*pi*factorial(n+m)) );
		Nnm = mag * ones(1,N_dirs);
		Nnm = [Nnm(end:-1:2, :); Nnm];
		
		CosSin = zeros(2*n+1,N_dirs);
			% m=0
		CosSin(n+1,:) = ones(1,size(dirs,1));
			% m>0
		CosSin(m(2:end)+n+1,:) = sqrt(2)*cos(m(2:end)*dirs(:,1)');
			% m<0
		CosSin(-m(end:-1:2)+n+1,:) = sqrt(2)*sin(m(end:-1:2)*dirs(:,1)');

		Ynm = Nnm .* Lnm .* CosSin;
        Y_N(idx+1:idx+(2*n+1), :) = Ynm;
        idx = idx + 2*n+1;
    end
    
    Y_N = Y_N.';
end


