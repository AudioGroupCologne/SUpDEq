function [template, target] = reijniers2014_featureextraction(SOFAtemplate, varargin)
%REIJNIERS2014_FEATUREEXTRACTION - extract HRTF using gammatone frequency bands and ITDs from SOFA object
%
%   Usage: [template, target] = reijniers2014_featureextraction(SOFAobj)
%
%   Input parameters:
%       SOFAtemplate: Struct in SOFA format with HRTFs
%
%   Output parameters:
%     template : template struct with spectral components
%     target   : template struct with spectral components
%
%   REIJNIERS2014_FEATUREEXTRACTION(...) computes temporally integrated
%   spectral magnitude profiles and itd.
%
%   REIJNIERS2014_FEATUREEXTRACTION accepts the following optional parameters:
% 
%     'source_ir',source_ir   Set the sound source's impulse reponse. 
%                             Default value a broadband sound source 
%                             with 0dB amplitude.
%
%     'fs',fs           Set the sampling rate to fs. 
%                       Default value takes the fs of SOFA template object.
%
%     'fb_ch',fb_ch     Set the number of channels for the gammatone
%                       filterbank to fb_ch.
%                       Default value is 30.
%
%     'fb_low',fb_low   Set the lowest frequency in the filterbank to
%                       fb_low. Default value is 300 Hz.
%
%     'fb_high',fb_high    Set the highest frequency in the filterbank to
%                          fhigh. Default value is 15000 Hz.
%
%     'ir_pad',len      Define the padding length for the impulse responses
%                       before being convolved with gammatone filters. 
%                       Default value is 0.05 s.
%
%     'targ_az',targ_az    Set the azimuth of a set of sound sources
%                          to targ_el. It can be a scalar or a column vector
%                          Default value is []: all target azimuths are
%                          used. Must have the same size of targ_el.
%
%     'targ_el',targ_el    Set the elevation of a set of sound sources
%                          to targ_el. It can be a scalar or a column vector
%                          Default value is []: all target elevations are
%                          used. Must have the same size of targ_az.
%
%   See also: exp_reijniers2014 reijniers2014
%
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
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/reijniers2014_featureextraction.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB SOFA M-Stats M-Image
%   #Author: Michael Sattler 
%   #Author: Roberto Barumerli (2020)
%   #Author: Clara Hollomey (2021)
% (adapted from code provided by Jonas Reijniers)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


definput.import={'amt_cache'};
definput.flags.disp = {'no_debug','debug'};

% sampling rate
definput.keyvals.fs = SOFAtemplate.Data.SamplingRate;

% gammatone filterbank parameters
definput.keyvals.fb_ch = 30;
definput.keyvals.fb_low = 300;
definput.keyvals.fb_high = 15e3;

% 50ms' pad
definput.keyvals.ir_pad = 0.05*definput.keyvals.fs;

definput.keyvals.targ_az = [];
definput.keyvals.targ_el = [];

% sound source
% broad band sound source at 0dB
definput.keyvals.source_ir = 0;

[flags, kv]  = ltfatarghelper({'source_ir', 'fs', 'fb_ch','fb_low','fb_high', ...
                            'ir_pad', 'targ_az','targ_el'}, definput, varargin);

if(kv.fs ~= SOFAtemplate.Data.SamplingRate)
    assert(kv.fs ~= SOFAtemplate.Data.SamplingRate)
end
% Get directions from SOFA file and convert them in cartesian
% coordinates
coords = SOFAcalculateAPV(SOFAtemplate);

% assume position on a sphere with radius of 1 meter
coords(:, 3) = 1;

% convert polar to cartesian
template_coords = zeros(size(SOFAtemplate.SourcePosition));
[template_coords(:,1), template_coords(:,2), template_coords(:,3)] = ...
    sph2cart(coords(:,1)*pi/180,coords(:,2)*pi/180,coords(:,3));

%% ITD computation 
% do alignment as was performed by Katz 2014
template_itd = itdestimator(SOFAtemplate,'Threshold', 'lp', ...
    'upper_cutfreq', 3000, 'butterpoly', 10, 'threshlvl', -10, flags.disp);

%% Pad HRIR vector
time_idx = find(SOFAtemplate.API.Dimensions.Data.IR == 'N');
dir_idx = find(SOFAtemplate.API.Dimensions.Data.IR == 'M');
ear_idx = find(SOFAtemplate.API.Dimensions.Data.IR == 'R');

% permute in order to use ufilterbankz
hrir = permute(double(SOFAtemplate.Data.IR),[time_idx, dir_idx, ear_idx]);
% pad to account for longer filters in the filterbank
pad_mat = zeros(kv.ir_pad - SOFAtemplate.API.('N'), SOFAtemplate.API.('M'), SOFAtemplate.API.('R'));
hrir = cat(1, hrir, pad_mat);

%% Gammatone filterbank
fc = fc_ERB(kv.fb_ch, kv.fb_low, kv.fb_high);
% if the number of channels exceed the fb_high 
% the vector will be shorter than kv.fb_ch
kv.fb_ch = length(fc);  

[bgt,agt] = gammatone(fc,kv.fs,'complex');

%% H_L and H_R generation
template_hrtf = 2*real(ufilterbankz(bgt,agt,hrir(:,:))); 
hrtf_size = size(hrir);
template_hrtf = reshape(template_hrtf,[hrtf_size(1),kv.fb_ch,hrtf_size(2),hrtf_size(3)]);
clear hrir
% Averaging over time (RMS)
template_hrtf = 20*log10(squeeze(rms(template_hrtf, 'dim', 1))+eps);      % in dB

% normalize
template_hrtf = template_hrtf - max(template_hrtf(:));

%% S computation
if isequal(kv.source_ir, definput.keyvals.source_ir)
    S = zeros(kv.fb_ch, 1);
else
    % if S is not the default value compute its spectrum
    % pad to account for longer filters in the filterbank
    temp = zeros(kv.ir_pad + length(kv.source_ir), 2);
    temp(1:length(kv.source_ir),:) = kv.source_ir;
    kv.source_ir = temp;
%     kv.source_ir = padarray(kv.source_ir(:), ...
%                                 [abs(kv.ir_pad - length(kv.source_ir)) 0],'post');
    S = 2*real(ufilterbankz(bgt,agt,kv.source_ir)); 
    % Averaging over time (RMS)
    S = 20*log10(squeeze(rms(S, 'dim', 1))+eps); 
end


% create struct
template.fs = kv.fs;
template.fc = fc;
template.itd = template_itd;
template.H = template_hrtf;
template.coords = template_coords;

% if target required
if nargout > 1
    %% target computation
    if(~isempty(kv.targ_az) || ~isempty(kv.targ_el))
        assert(numel(kv.targ_az)==numel(kv.targ_el))
        target_idx = SOFAfind(SOFAtemplate, kv.targ_az, kv.targ_el);
        if(numel(target_idx) ~= numel(kv.targ_az))
            amt_disp(['Requested HRTF''s points: ',numel(kv.targ_az)*numel(kv.targ_el),' Found:',numel(target_idx)],flags.disp);
        end
        target_hrtf = template_hrtf(:,target_idx,:);
        target_itd = template_itd(target_idx);
        target_coords = template_coords(target_idx, :);
    else
        target_hrtf = template_hrtf;
        target_itd = template_itd;
        target_coords = template_coords;
    end

    target.fs = kv.fs;
    target.fc = fc;
    target.itd = target_itd;
    target.S = S;
    target.H = target_hrtf;
    target.coords = target_coords;
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


