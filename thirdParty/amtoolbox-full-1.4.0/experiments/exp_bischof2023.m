function [OUT_struct] = exp_bischof2023(varargin)
%EXP_BISCHOF2023 experiments from Bischof et al. 2023
%   Usage: [OUT_struct] = exp_bischof2023(flags)
%
%
%   Input parameters:
%     flags : string to reproduce specific figure from Bischof et al. (2023)
%             ('fig3','fig4','fig7' or 'fig8')
%
%   Output parameters:
%     OUT_struct    : structure with all predicted SNRs, BMLDs and
%                     better ear SNRs for the experiment reported in
%                     Bischof et al. (2023).
%
%
%   EXP_BISCHOF2023
%   runs the model bischof2023 for DYNamic Binaural Unmasking (DynBU)with
%   binaural recordings of the original stimuli used in the detection 
%   experiment and returns an output sturcture containing the model
%   predictions as well as the experimental data for further analysis or
%   plotting of the results.
%
%
%   The following flags can be specified
%
%     'fig3'    Reproduce Fig.3:
%               Medians and quartiles of the measured binaural detection
%               thresholds of a reverberant harmonic complex tone for
%               different truncations of the room impule response in the
%               presence of an anechoic bandpass noise with 60 dB SPL from
%               the front. Solid lines indicate thresholds for a collocated
%               target sound source at 0°, dashed lines for a target sound
%               source at 60°. Data in blue correspond to measured
%               thresholds with an absorption coefficient of 0.5, data in
%               red for an absorption coefficient of 0.1.
%
%     'fig4'    Reproduce Fig.4:
%               Medians and quartiles of the measured binaural detection
%               thresholds of a reverberant harmonic complex tone located 
%               at 0° for different time conditions of cut early 
%               reflections from the room impulse response in the presence 
%               of an anechoic bandpass noise with 60 dB SPL from the 
%               front. Data in blue correspond to measured thresholds with
%               an absorption coefficient of 0.8, data in red for an
%               absorption coefficient of 0.1.
%
%     'fig7'    Reproduce Fig.7:
%               Predictions of bischof2023 are shown with green squares 
%               connected with dashed lines along with measured thresholds
%               reported in Figure 3 and Figure 4. The left column shows
%               data for an absorpiton coefficient of 0.1, the right column
%               for an absorption coefficient of 0.5 respectively. The
%               first row refers to data with a collocated target and noise
%               at 0°, the second row for a target at 60° and a noise
%               masker at 0°, both for different truncations of the room 
%               impule response. The thirs row refers to data with a 
%               collocated target and noise at 0° for different time 
%               conditions of cut early reflections from the room impulse 
%               response.
%               
%     'fig8'    Reproduce Fig.8:
%               Contributions of better-ear SNR (dark green shaded area)
%               and BMLD (light green shaded area) to the overall
%               predicted binaural benefit using bischof2023. The overall
%               prediction is ploted as detection benefit for all different
%               experimental conditions shown in Figure 7. The division of
%               the individual panels corresponds to that in Figure 7.
%
%
%   Examples:
%   ---------
%
%   To display Fig.3 use :
%
%     exp_bischof2023('fig3');
%
%   To display Fig.4 use :
%
%     exp_bischof2023('fig4');
%
%   To display Fig.7 use :
%
%     exp_bischof2023('fig7');
%
%   To display Fig.8 use :
%
%     exp_bischof2023('fig8');
%
%   See also: bischof2023_filterbank data_bischof2023 plot_bischof2023
%             bischof2023
%
%
%   References:
%     N. Bischof, P. Aublin, and B. Seeber. Fast processing models effects of
%     reflections on binaural unmasking. Acta Acustica, 2023.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_bischof2023.php


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

%% check input variables
definput.flags.type = {'missingflag','fig3','fig4','fig7', 'fig8'};
definput.flags.plot = {'plot','no_plot'};


[flags,~]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%if nargin < 1; do_plot = 0; end
%% load experimental data for Bischof et al. 2023 if not yet in workspace
if ~exist('expdata','var') == 1
    [expdata,fs] = data_bischof2023;
end
%% define additional variables for the experiment in Bischof et al 2023
interf_length = length(expdata.bischof2023_interf); % overall length of interferer signal

target_length = length(expdata.bischof2023_target_cell{1}); % overall length of recorded target signal

% define start sample index for time centered target signal
start_sample = ceil((interf_length - target_length)/2);
interf_sig = expdata.bischof2023_interf(start_sample+1:start_sample+target_length,:);
%% Define Model parameters for bischof2023
Model_params.fs = fs;               % sampling frequency
Model_params.f_range = [300,770];   % frequency range to be evaluated
Model_params.Bark_ord = 4;          % filter order for gammatone filters
Model_params.Bark_len = 512;        % filter length for gammatone filters
Model_params.t_st = 0.012;          % short time analysis window in sec
Model_params.t_SLUGGint = 0.225;    % time constant for sluggishness integration in sec
Model_params.t_INTint = 0.09;       % time constant for intensity integration in sec
%% run bischof2023
% initialize output variables
OUT_struct.pred_SNR_fast_bischof2023 = zeros(size(expdata.bischof2023_target_cell));
OUT_struct.pred_BMLD_fast_bischof2023 = zeros(size(OUT_struct.pred_SNR_fast_bischof2023));
OUT_struct.pred_BE_fast_bischof2023 = zeros(size(OUT_struct.pred_SNR_fast_bischof2023));

% run bischof2023 for all target signals
for ii = 1:numel(expdata.bischof2023_target_cell)
    target_sig = expdata.bischof2023_target_cell{ii};
    [OUT_struct.pred_SNR_fast_bischof2023(ii),...
     OUT_struct.pred_BMLD_fast_bischof2023(ii),...
     OUT_struct.pred_BE_fast_bischof2023(ii)] = bischof2023(target_sig,...
                                                            interf_sig,...
                                                            Model_params);
end
%% add expdata to OUT_struct
OUT_struct.expdata = expdata;
%% do plot
if flags.do_plot
    plot_bischof2023(OUT_struct, flags);
end

