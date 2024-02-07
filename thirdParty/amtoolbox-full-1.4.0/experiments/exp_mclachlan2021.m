function varargout = exp_mclachlan2021(varargin)
%EXP_MCLACHLAN2021 Results from McLachlan et al. (2021)
%   Usage: data = exp_mclachlan2021(flag) 
%
%   EXP_MCLACHLAN2021(flag) reproduces figures of the study from 
%   McLachlan et al. (2021).
%
%   The following flags can be specified
%
%     'fig6'    Reproduce Fig.6
%
%     'fig7a'    Reproduce Fig.7a
%
%     'fig7b'    Reproduce Fig.7b
%
%     'fig7c'    Reproduce Fig.7c
%
%   Examples:
%   ---------
%
%   To display Fig.6 use :
%
%     exp_mclachlan2021('fig6');
%
%   To display Fig.7a use :
%
%     exp_mclachlan2021('fig7a');
%
%   To display Fig.7b use :
%
%     exp_mclachlan2021('fig7b');
%
%   To display Fig.7c use :
%
%     exp_mclachlan2021('fig7c');
%
%
%   See also: mclachlan2021
%
%   References:
%     G. McLachlan, P. Majdak, J. Reijniers, and H. Peremans. Towards
%     modelling active dynamic sound localisation based on Bayesian
%     inference. Acta Acustica, 2021.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_mclachlan2021.php


%   #Author: Glen McLachlan

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.

%% ------ Check input options ---------------------------------------------
definput.import = {'amt_cache'};
definput.keyvals.MarkerSize = 6;
definput.keyvals.FontSize = 12;

definput.flags.type = {'missingflag','fig6','fig7a',...
    'fig7b','fig7c'};
definput.flags.plot = {'plot', 'no_plot'};
definput.flags.plot_type = {'interp','scatter'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},...
             definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.', ...
      upper(mfilename),flagnames);
end

%% Get listener's data
SOFA_obj = amt_load('mclachlan2021', 'HRIR_L2354.sofa');

%% Run model
    fig = [];
    
    % Preprocessing
    [template, target] = ...
        mclachlan2021_preproc(SOFA_obj,flags.type);
    
    [fig.doa, fig.params] = ...
        mclachlan2021(template, target, flags.type);
    %amt_cache('set','fig',fig); % save model output in cache

    if flags.do_plot
        mclachlan2021_metrics(fig.doa,'polar');
        met = mclachlan2021_metrics(fig.doa, 'fbc');
    else
        met = mclachlan2021_metrics(fig.doa);
    end

    % Calcualte performance measures 
    amt_disp('------------------------')
    amt_disp('Performance Predictions:')
    amt_disp('------------------------')
    
    met.entropy = mean(fig.params.entropy);
    met.information = mean(fig.params.information);
    
    %amt_disp(sprintf('Lateral accuracy: %0.2fdeg', met.accL))
    amt_disp(sprintf('Lateral RMS: %0.2f', met.rmsL))
    %amt_disp(sprintf('Elevation accuracy: %0.2fdeg', met.accE))
    amt_disp(sprintf('Polar RMS: %0.2f', met.rmsP))
    %amt_disp(sprintf('Quadrant error: %0.6f%%',met.querr))
    amt_disp(sprintf('Mean entropy: %0.2fbits', met.entropy))
    amt_disp('------------------------')
