function definput=arg_adaptloop(definput)
%ARG_ADAPTLOOP
%   #License: GPL
%   #Author: Peter Soendergaard (2011): Initial version
%   #Author: Alejandro Osses (2020): Extensions
%   #Author: Piotr Majdak (2021): Adapted to AMT 1.0
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_adaptloop.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
 
  definput.flags.adt = {'adt','no_adt'};
  definput.keyvals.limit=10;
  definput.keyvals.minspl=0; % lowest audible SPL of the signal (in dB)
  definput.keyvals.tau=[0.005 0.050 0.129 0.253 0.500];
    
  definput.groups.adt_dau1996 = {'tau',[0.005 0.050 0.129 0.253 0.500],'limit',0};
  definput.groups.adt_dau1997 = {'tau',[0.005 0.050 0.129 0.253 0.500],'limit',10};
  definput.groups.adt_breebaart2001 = {'tau',linspace(0.005,0.5,5),'limit',0};
  definput.groups.adt_puschel1988 = {'tau',linspace(0.005,0.5,5),'limit',0};
  definput.groups.adt_osses2021 = {'tau',[0.005 0.050 0.129 0.253 0.500],'limit',5};
  definput.groups.adt_relanoiborra2019 = {'tau',[0.005 0.050 0.129 0.253 0.500],'limit',10,'minspl',dbspl(2e-7,[],100)};

