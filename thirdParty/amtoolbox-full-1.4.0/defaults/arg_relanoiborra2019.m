function definput=arg_relanoiborra2019(definput)
% ARG_RELANOIBORRA2019
%
%   #License: GPL
%   #Author: Piotr Majdak (2021): created for the AMT 1.0
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_relanoiborra2019.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% General
definput.flags.disp = {'no_debug','debug'};
definput.keyvals.dboffset = dbspl(1); % dB Full scale convention of this model

%% AFB Parameters
definput.keyvals.flow=100;
definput.keyvals.fhigh=8000;
definput.flags.afb = {'erbspace', 'erbspacebw'};
  % used when erbspace is used to calculat fc's
definput.keyvals.fcnt = 60; % number of fc
  % used when erbspacebw is used to calculate fc's
definput.keyvals.bwmul = 0.5; % bandwidth of each fc 
definput.keyvals.basef = 8000; % one of the fc's will be exactly basef

definput.keyvals.subject='NH';
definput.keyvals.N_org = [];
definput.flags.internalnoise= {'internalnoise','no_internalnoise'};


%% IHC
definput.flags.ihc = {'ihc','no_ihc'}; 

%% Adaptation loops
definput.flags.an = {'an','no_an'};
definput.keyvals.limit=10; % arbitrary units
definput.keyvals.minspl=dbspl(2e-7,[],100); % approx. -34 dB, lowest audible SPL of the signal (in dB)
definput.keyvals.tau=[0.005 0.050 0.129 0.253 0.500];

%% Modulation filterbank



