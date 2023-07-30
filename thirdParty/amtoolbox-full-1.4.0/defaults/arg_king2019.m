function definput=arg_king2019(definput)
% ARG_KING2019
%
%   #License: GPL
%   #Author: Alejandro Osses (2020): first adaptation for AMT
%   #Author: Piotr Majdak (2021): further integration in the AMT 1.0
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_king2019.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% General
definput.flags.disp = {'no_debug','debug'};
definput.flags.afb = {'afb', 'no_afb'};
definput.flags.ihc = {'ihc','no_ihc'}; 
definput.flags.adt = {'adt','no_adt'};
definput.flags.mfb = {'mfb','no_mfb'}; 

%% AFB Parameters
definput.keyvals.flow=[];
definput.keyvals.fhigh=[];
definput.keyvals.gamma_order = 4; % order of the Gammatone filters in the auditory filterbank

definput.flags.compression = {'compression_brokenstick','compression_power', 'no_compression'};
definput.keyvals.compression_n    = 0.3; % 0.6 for normal hearing: see Wallaert2017
definput.keyvals.compression_knee_dB = 30; % dB
definput.keyvals.compression_a = [];
definput.keyvals.compression_b = [];

%% IHC
definput.flags.ihc_type={'ihc_king2019','ihc_undefined','ihc_bernstein1999','ihc_breebaart2001','ihc_dau1996','hilbert', ...
                    'ihc_lindemann1986','ihc_meddis1990','ihc_relanoiborra2019'};

%% Adaptation loops
definput.keyvals.adt_HP_fc = 3; % Hz
definput.keyvals.adt_HP_order = 1;

%% Modulation filterbank
definput.flags.modfilter_150Hz_LP = {'no_LP_150_Hz','LP_150_Hz'}; % modbank_LPfilter
definput.keyvals.mflow  =   2; % Hz, modbank_fmin
definput.keyvals.mfhigh = 150; % Hz, modbank_fmax
definput.keyvals.modbank_Nmod    = []; % number of filters, for overalpped 
                               % filters choose 'modbank_Nmod'
definput.keyvals.modbank_Qfactor =  1; % Q factor for the filters
definput.flags.modfilter_phase = {'phase_insens_hilbert', 'no_phase_insens'};
definput.keyvals.phase_insens_cut = 10; % Hz 
definput.keyvals.subfs = 16000; % Hz 


