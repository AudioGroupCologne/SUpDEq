function definput = arg_bruce2018(definput)
% ARG_BRUCE2018
%
%   #License: GPL
%   #Author: Clara Hollomey (2021)
%   #Author: Piotr Majdak (2021): adapted to AMT 1.1
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_bruce2018.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% General
definput.flags.species = {'human','cat'};
definput.keyvals.fsmod   = 100e3; % processing rate of the model (Hz) 

%% Parameters of the IHC
definput.keyvals.nrep = 1;
definput.keyvals.reptime = 1.2;

  % fitaudiogram: derive the inner and outer hair cell impairment factors 
  %   (ohcs and ihcs) via audiogram fitting given by ag_fs and ag_dbloss
  % no_fitaudiogram: use user-defined ohcs and ohcs parameters
definput.flags.fitHL = {'no_fitaudiogram', 'fitaudiogram'}; 
  % Used only if fitaudiogram
definput.keyvals.ag_fs = [125 250 500 1e3 2e3 4e3 8e3]; 
definput.keyvals.ag_dbloss = [0 0 0 0 0 0 0];
  % Used only if no_fitaudiogram
definput.keyvals.cihcs = 1; % scalar or vector of fc
definput.keyvals.cohcs = 1; % scalar or vector of fc

%% parameters to generate an auditory nerve population
  % 'outputPerSynapse' model will output results per synapse (slower)
  % 'outputPerCF' model will output average results over all synapses (faster)
definput.flags.output = {'outputPerCF', 'outputPerSynapse'};

  %PSTH parameters
definput.keyvals.psthbinwidth_mr = 100e-6;
definput.keyvals.windur_ft=32;
definput.keyvals.windur_mr=128;

  %synapse parameters
definput.flags.noiseType = {'fixedFGn','varFGn'};
definput.flags.powerLawImp = {'approxPL','actualPL'};%converted to '0' and '1' later on

  %nerve fiber paramters
definput.flags.SR={'autoSR','specificSR','specificSRautoTiming'};
  % Used if autoSR or specificSRautoTiming
definput.keyvals.numL = 4;   % Number of low SR fibers at each CF
definput.keyvals.lossL = 1;
definput.keyvals.numM = 4;   % Number of medium SR fibers at each CF
definput.keyvals.lossM = 1;
definput.keyvals.numH = 12;   % Number of high SR fibers at each CF
definput.keyvals.lossH = 1;
  % Used only if specificSRautoTiming 
definput.keyvals.SRL  = 0.1; % spikes/s, spontaneous rate
definput.keyvals.SRM  = 4; % spikes/s, spontaneous rate
definput.keyvals.SRH  = 70; % spikes/s, spontaneous rate
  % Used only if specificSR (SR specified as below)
definput.keyvals.numsponts = 10; % number of fibers
definput.keyvals.spont = 50; % SR (spikes/s); can be vector of numsponts
definput.keyvals.tabs = 0.6e-3; % absolute timing; can be vector of numsponts
definput.keyvals.trel = 0.6e-3; % relative timing; can be vector of numsponts
  
  



