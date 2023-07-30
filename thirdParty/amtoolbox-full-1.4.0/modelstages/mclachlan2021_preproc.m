function [template, target] = mclachlan2021_preproc(SOFAtemplate, varargin)
%MCLACHLAN2021_PREPROC - extract HRTF using gammatone frequency bands and ITDs from SOFA object
%   Usage: [template, target] = mclachlan2021_preproc(SOFAobj)
%
%   Input parameters:
%       SOFAtemplate: Struct in SOFA format with HRTFs
%
%   Output parameters:
%     template : template struct with spectral components
%     target   : template struct with spectral components
%
%   MCLACHLAN2021_PREPROC(...) computes temporally integrated
%   spectral magnitude profiles and itd.
%
%   MCLACHLAN2021_PREPROC accepts the following optional parameters:
% 
%     'source_ir',source_ir   Set the sound source's impulse reponse. 
%                             Default value a broadband sound source 
%                             with 0dB amplitude.
%
%     'fb_ch',fb_ch     Set the number of channels for the gammatone
%                       filterbank to fb_ch.
%                       Default value is 30.
%
%     'fb_low',fb_low   Set the lowest frequency in the filterbank to
%                       fb_low. Default value is 300 Hz.
%
%     'fb_high',fb_high   Set the highest frequency in the filterbank to
%                         fhigh. Default value is 15000 Hz.
%
%     'SNR',SNR         Set the signal to noise ratio corresponding to
%                       different sound source intensities.
%                       Default value is SNR = 75 [dB]
%
%     'targ_az',targ_az   Set the azimuth of a set of sound sources
%                         to targ_el. It can be a scalar or a column vector
%                         Default value is []: all target azimuths are
%                         used. Must have the same size of targ_el.
%
%     'targ_el',targ_el   Set the elevation of a set of sound sources
%                         to targ_el. It can be a scalar or a column vector
%                         Default value is []: all target elevations are
%                         used. Must have the same size of targ_az.
%
%   See also: exp_reijniers2014 reijniers2014
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/mclachlan2021_preproc.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: M-Signal M-Image
%   #Author: Glen McLachlan (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   References:
%     R. Barumerli, P. Majdak, R. Baumgartner, J. Reijniers, M. Geronazzo,
%     and F. Avanzini. Predicting directional sound-localization of human
%     listeners in both horizontal and vertical dimensions. In Audio
%     Engineering Society Convention 148. Audio Engineering Society, 2020.
%     
%     R. Barumerli, P. Majdak, R. Baumgartner, M. Geronazzo, and F. Avanzini.
%     Evaluation of a human sound localization model based on bayesian
%     inference. In Forum Acusticum, 2020.
%     
%     J. Reijniers, D. Vanderleist, C. Jin, C. S., and H. Peremans. An
%     ideal-observer model of human sound localization. Biological
%     Cybernetics, 108:169--181, 2014.    

%   AUTHOR: Michael Sattler and Roberto Barumerli (adapted from code provided by Jonas Reijniers)

definput.import={'amt_cache'};
definput.import={'mclachlan2021'};

[~, kv]  = ltfatarghelper({}, definput, varargin);

fs = SOFAtemplate.Data.SamplingRate; % sampling rate
ir_pad = 0.05*fs; % 50ms pad

% Get directions from SOFA file 
SOFAcoords_sph = SOFAcalculateAPV(SOFAtemplate);

% assume position on a sphere with radius of 1 meter
SOFAcoords_sph(:, 3) = 1;

% convert polar (radians) to cartesian
SOFAcoords_crt = zeros(size(SOFAtemplate.SourcePosition));
[SOFAcoords_crt(:,1), SOFAcoords_crt(:,2), SOFAcoords_crt(:,3)] = ...
    sph2cart(SOFAcoords_sph(:,1)*pi/180,SOFAcoords_sph(:,2)*pi/180,SOFAcoords_sph(:,3));

%% S computation
if isequal(kv.source_ir, 0)
    S = zeros(kv.fb_ch, 1);
else
    % if S is not the default value compute its spectrum
    % pad to account for longer filters in the filterbank
    kv.source_ir = padarray(kv.source_ir(:), ...
                                [abs(ir_pad - length(kv.source_ir)) 0],'post');
    S = 2*real(ufilterbankz(bgt,agt,kv.source_ir)); 
    % Averaging over time (RMS)
    S = 20*log10(squeeze(rms(S, 'dim', 1))+eps); 
end

%% sample uniformly over sphere with N is number of directions
% contains the sampled points on a unitary sphere in cartesian coords
%dirs=amt_load('reijniers2014','dirs.mat');
dangle = 1;

UNIcoords=amt_load('mclachlan2021', 'dirs_2000.mat');
UNIcoords=UNIcoords.dirs;
if(isempty(UNIcoords))
    error('New directions grid not available. Please check your internet connection!')
end

% remove the points from the unitary sphere below original HRTF lowest elevation
idx = find(UNIcoords(:,3) > min(SOFAcoords_crt(:, 3))); 
UNIcoords = UNIcoords(idx,:);
num_dirs = length(idx); 

% create 1 degree (yaw) rotated over  coordinates (for dITD)
UNIcoords_rot = mclachlan2021_rotatedirs(UNIcoords,dangle,kv.rot_type);

%% ITD computation 
% do alignment as was performed by Katz 2014
%template_itd = itdestimator(SOFAtemplate,'Threshold', 'lp', ...
%    'upper_cutfreq', 3000, 'butterpoly', 10, 'threshlvl', -10, 'silent');
template_itd = itdestimator(SOFAtemplate,'Threshold', 'lp', ...
    'upper_cutfreq', 3000, 'butterpoly', 10, 'threshlvl', -10);
%lat=pi/2-acos(SOFAcoords_crt(:,2));
%pol=sign(SOFAcoords_crt(:,3)).*acos(SOFAcoords_crt(:,1)./(cos(lat)));
%template_itd = woodworth(lat,pol);

%transformation to JND
a = 32.5e-6;
b = 0.095;
template_itd = sign(template_itd) .* ((log(a + b * abs(template_itd)) - log(a)) / b); 

%% Pad HRIR vector
time_idx = find(SOFAtemplate.API.Dimensions.Data.IR == 'N');
dir_idx = find(SOFAtemplate.API.Dimensions.Data.IR == 'M');
ear_idx = find(SOFAtemplate.API.Dimensions.Data.IR == 'R');

% permute in order to use ufilterbankz
hrir = permute(double(SOFAtemplate.Data.IR),[time_idx, dir_idx, ear_idx]);
% pad to account for longer filters in the filterbank
pad_mat = zeros(ir_pad - SOFAtemplate.API.('N'), SOFAtemplate.API.('M'), SOFAtemplate.API.('R'));
hrir = cat(1, hrir, pad_mat);

%% Gammatone filterbank
fc = fc_ERB(kv.fb_ch, kv.fb_low, kv.fb_high);
% if the number of channels exceed the fb_high 
% the vector will be shorter than kv.fb_ch
kv.fb_ch = length(fc);  

[bgt,agt] = gammatone(fc,fs,'complex');

%% H_L and H_R generation
template_hrtf = 2*real(ufilterbankz(bgt,agt,hrir(:,:))); 
hrtf_size = size(hrir);
template_hrtf = reshape(template_hrtf,[hrtf_size(1),kv.fb_ch,hrtf_size(2),hrtf_size(3)]);
clear hrir
% Averaging over time (RMS)
template_hrtf = 20*log10(squeeze(rms(template_hrtf, 'dim', 1))+eps);      % in dB

% normalize
template_hrtf = template_hrtf - max(template_hrtf(:));

% add source spectrum to target and to template
template_hrtf = template_hrtf + repmat(S(:), 1, size(template_hrtf, 2), 2);

% SHOULD TEMPLATE INCLUDE SOUND SOURCE SPECTRUM?

SNR = kv.SNR; % defined as maximal SNR (in interval 2kHz-7kHz)

% account for SNR and frequency-dependent hearing sensitivity (see section 2.1 in SI) 
template_hrtf = max(template_hrtf ,-SNR); 
template_hrtf(fc<=2000,:,:) = max(template_hrtf(fc<=2000,:,:),-SNR + 10);
template_hrtf(fc>=7000,:,:) = max(template_hrtf(fc>=7000,:,:),-SNR + 20);

%% interpolate to uniform distribution
%spherical harmonics basis functions
Y_N = local_SH(kv.SHorder, [SOFAcoords_sph(:,1)*pi/180, SOFAcoords_sph(:,2)*pi/180]); 

% tikonov regularisation
lambda = 4;
SIG = eye((kv.SHorder+1)^2);
SIG(1:(2+1)^2,1:(2+1)^2) = 0;

%expansion coefficients
expcoef_H(:,:,1) = transpose((Y_N'*Y_N+lambda*SIG)\Y_N'*squeeze(template_hrtf(:,:,1))');
expcoef_H(:,:,2) = transpose((Y_N'*Y_N+lambda*SIG)\Y_N'*squeeze(template_hrtf(:,:,2))');
expcoef_itd = (Y_N'*Y_N+lambda*SIG)\Y_N'*template_itd;

% interpolate itd values to original coords
[AZ,EL] = cart2sph(UNIcoords(:,1),UNIcoords(:,2),UNIcoords(:,3));
Y_N = local_SH(kv.SHorder, [AZ EL]); 
template_itd = (Y_N*expcoef_itd)';
template_hrtf = []; % clear and replace hrtfs to uniform grid
template_hrtf(:,:,1) = transpose(Y_N*squeeze(expcoef_H(:,:,1))');
template_hrtf(:,:,2) = transpose(Y_N*squeeze(expcoef_H(:,:,2))');  

%% dITD computation
% interpolate itd values to rotated coords
[AZ,EL] = cart2sph(UNIcoords_rot(:,1),UNIcoords_rot(:,2),UNIcoords_rot(:,3));
Y_N = local_SH(kv.SHorder, [AZ EL]); 
template_itdt = (Y_N*expcoef_itd)';

%compute dITD
template_ditd = (template_itdt-template_itd);

%dITD expansion coefficient
[AZ,EL] = cart2sph(UNIcoords(:,1),UNIcoords(:,2),UNIcoords(:,3));
Y_N = local_SH(kv.SHorder, [AZ EL]); %SH basis functions
expcoef_ditd = (Y_N'*Y_N+lambda*SIG)\Y_N'*template_ditd';

%plot the dITD values in JND per degree (take at SHorder 5)
%plot_reijniers2014(target.coords,(target.itdt-target.itd0)*1e6);
%caxis([-10 10]);

%% create struct
template.fs = fs;
template.fc = fc;
template.itd = template_itd;
template.ditd = template_ditd;
template.coef.itd = expcoef_itd;
template.coef.ditd = expcoef_ditd;
template.coef.H = expcoef_H;
template.H = template_hrtf;
template.coords = UNIcoords;

% if target required
if nargout > 1
    % target computation
    if(~isempty(kv.targ_az) || ~isempty(kv.targ_el))
        assert(numel(kv.targ_az)==numel(kv.targ_el))
        % interpolate to target coordinates
        Y_N = local_SH(kv.SHorder, [kv.targ_az*pi/180 kv.targ_el*pi/180]); 
        target_itd = (Y_N*expcoef_itd)';
        target_ditd = (Y_N*expcoef_ditd)';
        target_hrtf(:,:,1) = transpose(Y_N*squeeze(expcoef_H(:,:,1))');
        target_hrtf(:,:,2) = transpose(Y_N*squeeze(expcoef_H(:,:,2))');  
        [target_coords(:,1), target_coords(:,2), target_coords(:,3)] = ...
            sph2cart(kv.targ_az*pi/180,kv.targ_el*pi/180,ones(length(kv.targ_az),1));
        
    else % if no targets given, equal all template coordinates
        target_hrtf = template_hrtf;
        target_itd = template_itd;
        target_ditd = template_ditd;
        target_coords = UNIcoords;
    end

    target.fs = fs;
    target.fc = fc;
    target.itd = target_itd;
    target.ditd = target_ditd;
    target.S = S;
    target.H = target_hrtf;
    target.coords = target_coords;
end
end

function fc = fc_ERB(n_channels, freq_start, freq_end)
    % ERB computation according to Moore and Glasberg 1983
    c = 1;
    fc = zeros(n_channels, 1);
    fc(1) = freq_start;
    while (c < n_channels)
        c =c + 1; 
        fc(c) = fc(c-1) + 6.23 * (fc(c-1)/1000)^2 + 93.39 * (fc(c-1)/1000) + 28.52;
    end
    
    fc(fc > freq_end) = [];
end

function Y_N = local_SH(N, dirs)
% calculate spherical harmonics up to order N for directions dirs [azi ele;...] (in radiant)
% 
    N_dirs = size(dirs, 1);
    N_SH = (N+1)^2;
	dirs(:,2) = pi/2 - dirs(:,2); % convert to inclinations

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


