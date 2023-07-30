function definput=arg_auditoryfilterbank(definput)
% ARG_AUDITORYFILTERBANK
%
%   #License: GPL
%   #Author: Peter Soendergaard (2011): Initial version
%   #Author: Alejandro Osses (2020): Extensions
%   #Author: Piotr Majdak (2021): Adapted to AMT 1.0
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_auditoryfilterbank.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% General
definput.keyvals.dboffset = dbspl(1); % dB Full scale convention
definput.flags.outerear = {'no_outerear','outerear'};
definput.flags.middleear= {'no_middleear','middleear','jepsen2008'};
definput.flags.singlefc = {'default', 'lavandier2022'};

%% Auditory filterbank
definput.keyvals.flow=80;
definput.keyvals.fhigh=8000;
definput.keyvals.basef=[];
definput.keyvals.bwmul=1;

definput.keyvals.fs_up  = [];
definput.flags.internalnoise= {'no_internalnoise', 'internalnoise'};

%% Groups
definput.groups.afb_dau1997 = {'dboffset',100,'basef',1000};

% relanoiborra2019: basef=8000 Hz, to match this frequency
%                   bwmul=0.5 because their filter design uses 60 bands
%                      spaced at 0.5 ERBN between 100 and 8000 Hz
definput.groups.drnl_relanoiborra2019 = {'outerear','middleear','basef',8000, ...
    'hearing_profile','NH','internalnoise','bwmul',0.5};

definput.groups.afb_osses2021 = {'outerear','middleear','basef',[]}; 


