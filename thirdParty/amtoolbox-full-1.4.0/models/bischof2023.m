function [DynBU_OUT,BMLD_pred,BE_pred] = bischof2023(INtarget,INinterf,params)
%BISCHOF2023 binaural masking level differences moving sound sources
%   Usage: [DynBU_OUT,BMLD_pred,BE_pred] = bischof2023(INtarget,INinterf,params)
%
%   Input parameters:
%     INtarget       : 2-channel waveform of target signal in columns
%     INinterf       : 2-channel waveform of interferer signal in columns
%     params         : parameter structure with the following fields
%
%                      '.fs'          sampling frequency [Hz] (default: 44.1e3)
%
%                      '.f_range'     1x2 vector defining the frequency range in [Hz] for auditory filters (default: [10 15.5e3])
%
%                      '.Bark_ord'    filter order for Bark filters (default: 4)
%
%                      '.Bark_len'    filter length for Bark filters (default: 512)
%
%                      '.t_st'        short time evaluation window [sec] (default: 0.012)
%
%                      '.t_SLUGGint'  time constant for binaural sluggishness integration [sec] (default: 0.225)
%
%                      '.t_INTint'    time constant for intensity integration of better-ear SNR [sec] (default: 0.090)
%
%   Output parameters:
%     DynBU_OUT          : overall binaural unmasking to detect the given
%                          target signal in the given noise signal [in dB]
%     BMLD_pred          : predicted binaural masking level difference [in dB]
%     BE_pred            : predicted better-ear SNR [in dB]
%
%
%
%   BISCHOF2023
%   DynBU_fast: DYNamic Binaural Unmasking model with "fast" cue extraction.
%   Calculates the better ear and binaural benefit for detecting a dynamic
%   sound source in noise. 
%
%   Target and interferer must be binaural signals (two channels only)
%
%   References:
%     J. F. Culling, M. L. Hawley, and R. Y. Litovsky. The role of
%     head-induced interaural time and level differences in the speech
%     reception threshold for multiple interfering sound sources. J. Acoust.
%     Soc. Am., 116(2):1057--1065, august 2004.
%     
%     N. Bischof, P. Aublin, and B. Seeber. Fast processing models effects of
%     reflections on binaural unmasking. Acta Acustica, 2023.
%     
%     N. Durlach. Binaural signal detection: Equalization and cancellation
%     theory. Academic Press, New York, 1972. Fundamentals of Modern Auditory
%     Theory, Volume II.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/bischof2023.php


%   #Author: Norbert F. Bischof (2023)
%   #Author: Pierre G. Aublin
%   #Author: Bernhard Seeber (2023)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.


%% check input signals
% working on column vectors for target and interferer
if size(INtarget,2)>size(INtarget,1)
    warning('binaural target signal was transposed to work on columns!');
    INtarget = INtarget';
end
if size(INinterf,2)>size(INinterf,1)
    warning('binaural interferer signal was transposed to work on columns!');
    INinterf = INinterf';
end
% check if target and interferer are binaural signals (two channels)
if size(INtarget,2)~=2
    error('Target must be a binaural signal (two channels)');
end
if size(INinterf,2)~=2
    error('Interferer must be a binaural signal (two channels)');
end
% check if target and interferer have the same signal length
if length(INtarget)~=length(INinterf)
    error('Target and interferer signal must have the same length!');
end

% if params are not defined use default parameters
if nargin < 3
    params = [];
end
if ~isfield(params,'fs'); params.fs = 44100; end
if ~isfield(params,'f_range'); params.f_range = [10 15.5e3]; end
if ~isfield(params,'Bark_ord'); params.Bark_ord = 4; end
if ~isfield(params,'Bark_len'); params.Bark_len = 512; end
if ~isfield(params,'t_st'); params.t_st = 0.024; end
if ~isfield(params,'t_SLUGGint'); params.t_SLUGGint = 0.225; end
if ~isfield(params,'t_INTint'); params.t_INTint = 0.090; end
%% define some additional fixed parameters needed in the model
% sample time of signal
params.ts = 1/params.fs;

% sigma_e & sigma_d as defined by Durlach (1972) and mentioned in Culling 
% et al. (2004) for EC process
params.sigma_e = 0.25;
params.sigma_d = 0.0000105;

% short time windows for given t_st
params.st_block_size = ceil(params.t_st*params.fs);

% short time Hanning window
params.st_win = hann(params.st_block_size);
% number of short time blocks with 50% overlap
params.st_nframes = 2*floor(length(INtarget)/params.st_block_size)-1;
% sampling frequency according to short time blocks with 50% overlap
params.st_fs = 2/params.t_st;

% get center frequencies of auditory filters according to Bark scale
[~,~,fc] = bischof2023_filterbank(0,params.Bark_ord,params.Bark_len,params.fs,params.f_range);
%% initialize vectors for loop
% output variables
Dyn_BMLD = zeros(params.st_nframes,length(fc));
Dyn_BE_SNR = zeros(params.st_nframes,length(fc));
%% filter target and interferer signal with gammatone filters according to Bark scale before further processing
target_sig_BARK(:,:,1) = bischof2023_filterbank(INtarget(:,1),params.Bark_ord,params.Bark_len,params.fs,params.f_range);
target_sig_BARK(:,:,2) = bischof2023_filterbank(INtarget(:,2),params.Bark_ord,params.Bark_len,params.fs,params.f_range);
interf_sig_BARK(:,:,1) = bischof2023_filterbank(INinterf(:,1),params.Bark_ord,params.Bark_len,params.fs,params.f_range);
interf_sig_BARK(:,:,2) = bischof2023_filterbank(INinterf(:,2),params.Bark_ord,params.Bark_len,params.fs,params.f_range);

%% Get better-ear SNR and BMLD
% iterate across all short-time analysis windows
for ii = 1:params.st_nframes
    % derive current target and interferer short time block
    target_BARK_shorttime = target_sig_BARK(floor((ii-1)*(params.st_block_size/2))+1 : floor((ii-1)*(params.st_block_size/2))+params.st_block_size,:,:).*params.st_win;
    interf_BARK_shorttime = interf_sig_BARK(floor((ii-1)*(params.st_block_size/2))+1 : floor((ii-1)*(params.st_block_size/2))+params.st_block_size,:,:).*params.st_win;
    
    % Short-time signal-to-noise ratio on both ear signals across filters
    % The maximum across ears is further used as better-ear SNR.
    for jj = 1:length(fc)
        % derive SNR on left and right ear
        SNR_left = 20*log10(sqrt(mean(target_BARK_shorttime(:,jj,1).^2))) - 20*log10(sqrt(mean(interf_BARK_shorttime(:,jj,1).^2)));
        SNR_right = 20*log10(sqrt(mean(target_BARK_shorttime(:,jj,2).^2))) - 20*log10(sqrt(mean(interf_BARK_shorttime(:,jj,2).^2)));

        Dyn_BE_SNR(ii,jj) = max(max([SNR_left SNR_right]),0);
    end
    
    % Short-time BMLD derived with the formula according to Culling et al.
    % (2004). Therefore, IPDs of target and interferer as well as
    % interaural coherence of the interferer signal are derived using a
    % normalized cross-correlation. Only positive BMLDs are used.
    for jj = 1:length(fc)
        k = (1 + params.sigma_e^2) * exp((2*pi*fc(jj)).^2 * params.sigma_d^2);
        
        [IPD_target,~] = calc_crosscorr(target_BARK_shorttime(:,jj,1),target_BARK_shorttime(:,jj,2),params.fs,fc(jj));
        [IPD_interf,IACC_interf] = calc_crosscorr(interf_BARK_shorttime(:,jj,1),interf_BARK_shorttime(:,jj,2),params.fs,fc(jj));
        
        % using the formula provided in Culling et al. (2004)
        Dyn_BMLD(ii,jj) = max((k - cos(IPD_target - IPD_interf)) ./ (k - IACC_interf),1);
    end
end

% define 1st order IIR exponential decay integration filters for
% sluggishness and intensity integration
num_sluggishness = 1 - exp(- (1/(params.t_SLUGGint*params.st_fs)));
denom_sluggishness = [1 num_sluggishness-1];

num_intensity = 1 - exp(- (1/(params.t_INTint*params.st_fs)));
denom_intensity = [1 num_intensity-1];

% apply 1st order IIR integration window to short time BMLD and
% better-ear SNR
pred_BMLD = 10*log10(filter(num_sluggishness,denom_sluggishness,Dyn_BMLD,[],1));
% prevent BMLD to be come negative
pred_BMLD = max(pred_BMLD,0);

% for intensity integration transform first back to intensities, then
% transform to decibels
Dyn_BE_SNR = (10.^(Dyn_BE_SNR./20)).^2;
pred_BE_SNR = filter(num_intensity,denom_intensity,Dyn_BE_SNR,[],1);
pred_BE_SNR = 10*log10(pred_BE_SNR);

% sum both contributions
PRED = pred_BMLD + pred_BE_SNR;

% apply MAX picking accross frequency and time for detection tasks
% select maximum contribution accross frequencies
[PRED_fmax,idx_fmax] = max(PRED,[],2,'omitnan');

% select maximum contribution over time
[PRED_tmax,idx_tmax] = max(PRED_fmax,[],1,'omitnan');

BMLD_pred = pred_BMLD(idx_tmax,idx_fmax(idx_tmax));
BE_pred = pred_BE_SNR(idx_tmax,idx_fmax(idx_tmax));
DynBU_OUT = PRED_tmax;

function [phase,coherence] = calc_crosscorr(left,right,fs,fc)
%   Usage: [phase,coherence] = calc_crosscorr(left,right,fs,fc)
%
% 
%   Input parameters:
%        left   : left ear signal
%        right  : right ear signal
%        fs     : sampling frequency [Hz]
%        fc     : center frequency to calculate the phase difference
%
%   Output parameters:
%        phase      : interaural phase difference at a given center frequency
%        coherence  : interaural coherence
%
%
%   CALC_CROSSCORR
%   calculates the interaural corss-correlation of a binaural signal at a
%   given center-frequency and returns the interaural phase difference at a
%   given ceneter-frequency and the interaural coherence.

% Version:
%           v1.0.2023-01-10
% History:

% (c) Bischof, N.F. & Seeber, B.U., AIP TUM Jan 2020 - 2023

% using a normalized cross-correlation of left and right ear signals to
% derive interaural coherence, as the maximum of the cross-correlation, and 
% the phase difference at a given center frequency fc.

iacc = xcorr(left,right,round(fs./(fc.*2)),'coeff');

[coherence, delay_samp] = max(iacc);

delay_samp=floor(delay_samp-length(iacc)/2);
phase = 2*pi*fc.*delay_samp/fs;
end
end

