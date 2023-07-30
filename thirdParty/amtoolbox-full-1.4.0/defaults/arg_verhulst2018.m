function definput=arg_verhulst2018(definput)
% ARG_VERHULST2018
%
%   #License: GPL
%   #Author: Piotr Majdak (2021)
%% General
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_verhulst2018.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
definput.flags.disp = {'debug','no_debug'}; % flag to provide debugging information when model called with 'debug', see amt_disp
definput.flags.ihc  = {'ihc','no_ihc'}; % Calculate IHC and provide detailed IHC response in the output
definput.flags.an   = {'an','no_an'};   % Calculate IHC, AN and provide detailed AN response in the output
definput.flags.cn   = {'cn','no_cn'};   % Calculate IHC, AN, CN and provide detailed CN response in the output
definput.flags.ic   = {'ic','no_ic'};   % Calculate IHC, AN, CN, IC and provide detailed IC response in the output
definput.flags.mfb  = {'no_mfb','mfb'}; % Determines the type of W3 and W5 in the output

%% Outer and middle ear
definput.flags.outerear = {'no_outerear' ,'outerear'};
definput.flags.middleear= {'middleear','no_middleear','jepsen2008'};

%% cochlear filterbank parameters

definput.keyvals.hearing_profile = 'Flat00'; % hearing profile
definput.keyvals.subject= 1; % to control randomness in cochlear irregularities
definput.keyvals.irr_on = 1; % irregularities on (1) or off (0)
definput.keyvals.non_linear_type = 'vel'; % 'vel' or 'disp'
definput.keyvals.IrrPct = 0.05;
definput.keyvals.fs_up  = 100000;
definput.flags.y    = {'no_y','y'}; % Include the detailed BM displacement in the output
definput.flags.v    = {'no_v','v'}; % Include the detailed BM velocity in the output. 
definput.flags.oae  = {'no_oae','oae'}; % Include OAE in the output

%% IHC
definput.keyvals.ihc_scal_constant=0.118;

%% AN
definput.keyvals.numH = 13; % Number of high-SR neurones
definput.keyvals.kSR_H  = 68.5; % (spikes/s) spontaneous rate
definput.keyvals.kmax_H = 3000; % Peak exocytosis rate (spikes/s)

definput.keyvals.numM =  3; % Number of medium-SR neurones
definput.keyvals.kSR_M  = 10; % (spikes/s) spontaneous rate
definput.keyvals.kmax_M = 1000; % Peak exocytosis rate (spikes/s)

definput.keyvals.numL =  3; % Number of low-SR neurones
definput.keyvals.kSR_L  = 1; % (spikes/s), spontaneous rate
definput.keyvals.kmax_L = 800; % Peak exocytosis rate (spikes/s)

definput.flags.anfH = {'no_anfH','anfH'}; % Include the detailed high-SR data in the output
definput.flags.anfM = {'no_anfM','anfM'}; % Include the detailed medium-SR data in the output
definput.flags.anfL = {'no_anfL','anfL'}; % Include the detailed low-SR data in the output
%% Brainstem (CN and IC)
definput.keyvals.model_version = 1.20;
definput.keyvals.version_year = 2018;
definput.keyvals.M1 = 4.2767e-14;
definput.keyvals.M3 = 5.1435e-14;
definput.keyvals.M5 = 13.3093e-14;

definput.keyvals.Tex_cn = 0.5e-3; % Tau excitation
definput.keyvals.Tin_cn = 2e-3; % Tau inhibition
definput.keyvals.dly_cn = 1e-3; % dly ms
definput.keyvals.Acn    = 1.5;
definput.keyvals.Scn    = 0.6;

definput.keyvals.Tex_ic = 0.5e-3; % Tau excitation
definput.keyvals.Tin_ic = 2e-3; % Tau inhibition
definput.keyvals.dly_ic = 2e-3; % dly ms
definput.keyvals.Aic    = 1.0;
definput.keyvals.Sic    = 1.5;
definput.keyvals.subfs  = 20000;


