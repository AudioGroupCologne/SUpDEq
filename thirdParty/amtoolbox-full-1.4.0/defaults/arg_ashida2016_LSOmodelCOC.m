function  definput = arg_ashida2016_LSOmodelCOC(definput)
%ARG_ASHIDA2016_LSOMODELCOC
%
% [lso] = lso_default_param() loads mpar for lso
% 
% Original file name: lso_default_mpar.m
%
%   #Author: Go Ashida (2016)
%   #Author: Alejandro Osses (2023) : integrated to AMT 1.4.0
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_ashida2016_LSOmodelCOC.php


definput.flags.lso_type = {'coc_mid','not_specified'};
definput.keyvals.lso.fs = 100e3; % Hz
definput.keyvals.lso.max_rate = 600;
definput.keyvals.lso.data_name = 'klug_2020';

% structure parameters
definput.keyvals.lso.inputfiberLists_ipsi_ex =[];
definput.keyvals.lso.inputfiberLists_contra_in =[];
definput.keyvals.lso.inputfiberLists_ipsi_in =[];
definput.keyvals.lso.inputfiberLists_contra_ex =[];
definput.keyvals.lso.mapping_method = @uniform_within_cf_mapping;
definput.keyvals.lso.cf = [];
definput.keyvals.lso.neurons_per_cf = 1;

% parameters (ashida 2016)
definput.keyvals.lso.Tref = 1.6e-3; % refractory time
definput.keyvals.lso.ThEx = 3;   % threshold
definput.keyvals.lso.WiEx = 1.1e-3; % length coincidence window
definput.keyvals.lso.gIn = 2;   % threshold increase by inhibition
definput.keyvals.lso.WiIn = 3.1e-3; % length inhibition window
definput.keyvals.lso.fibersPerNeuron_ipsi_ex = 20;
definput.keyvals.lso.fibersPerNeuron_contra_in = 8;
definput.keyvals.lso.fibersPerNeuron_ipsi_in = 0;
definput.keyvals.lso.fibersPerNeuron_contra_ex = 0;
definput.keyvals.lso.best_ipd = 0;

% default_set = 'default';

% if ~isdir('param_store')%CH
%     mkdir('param_store')
% end
% save(['param_store/',default_set,'_lso.mat'],'lso')

