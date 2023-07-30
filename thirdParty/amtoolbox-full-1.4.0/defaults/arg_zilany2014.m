function definput = arg_zilany2014(definput)
% ARG_ZILANY2014
%
%   #License: GPL
%   #Author: Piotr Majdak (2021)
%% General
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_zilany2014.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
definput.flags.disp = {'no_debug','debug'};
definput.flags.species = {'human','cat'};
definput.keyvals.dboffset = dbspl(1);

%% Outer and middle ear

%% Auditory filterbank
definput.flags.BMtuning = {'shera2002','glasberg1990'}; % selects the BM tuning, either Shera et al. (2002) or Glasberg & Moore (1990)

%% IHC & OHC
definput.keyvals.cohc    = 1;
definput.keyvals.cihc    = 1;
definput.keyvals.nrep    = 1; % repeated calculations, 1: fast; 500: nice PSTH
definput.keyvals.reptime = 2; % repeat the stimulus in intervals of reptime*signallength

%% AN
definput.keyvals.fiberType = 2; % Simulate a neuron with the following SR: 1=Low; 2=Medium; 3=High; 4:Custom (see numH, numM and numL)
definput.flags.noiseType = {'fixedFGn','varFGn'};
definput.flags.powerLawImp = {'approxPL','actualPL'};
definput.keyvals.fsmod   = 100e3;

definput.keyvals.numH = 12; % # of high SR neurones, if fiberType=4
definput.keyvals.numM = 4; % # of high SR neurones, if fiberType=4
definput.keyvals.numL = 4; % # of high SR neurones, if fiberType=4

definput.keyvals.psth_binwidth = []; % width for PSTH calculations (in s). Set to [] for a raw PSTH (i.e., time resolution 1/fsmod)


