function [sponts,tabss,trels]=bruce2018_generateanpopulation(numcfs,numsponts)
%BRUCE2018_GENERATEANPOPULATION generates an AN population for a given number of nerve fibres and numsponts
%
%   Usage:
%     [sponts,tabss,trels]=bruce2018_generateanpopulation(numcfs,numsponts)
%
%   Input parameters:
%     numfcs    : number of frequencies
%     numsponts : number of low- mid- and highly-spontaneously firing nerve fibres
%                 [numlow nummid numhigh]
%
%   Output parameters:
%     sponts    : nerve fibres
%     tabs      : absolute timing info
%     trels     : relative timing info
%
%   BRUCE2018_GENERATEANPOPULATION generates an AN population for a 
%   given absolute number of nerve fibres and relative proportion of
%   low- mid- and highly-spontaneously firing nerve fibres
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/bruce2018_generateanpopulation.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB MEX M-Signal
%   #Author: Ian Bruce: basic code of the model
%   #Author: Alejandro Osses (2020): original implementation
%   #Author: Clara Hollomey (2021): adapted to the AMT 1.0, removed saving of ANpopulation.mat
%   #Author: Piotr Majdak (2021): adaptations to exp_osses2022; specificSRautoTiming added

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


tabsmax = 1.5*461e-6;
tabsmin = 1.5*139e-6;
trelmax = 894e-6;
trelmin = 131e-6;

% generate sponts, tabss & trels for LS fibers (fiberType = 1)
sponts.LS = min(max(0.1+0.1*randn(numcfs,numsponts(1)),1e-3),0.2);
refrand = rand(numcfs,numsponts(1));
tabss.LS = (tabsmax - tabsmin)*refrand + tabsmin;
trels.LS = (trelmax - trelmin)*refrand + trelmin;

% generate sponts, tabss & trels for MS fibers (fiberType = 2)
sponts.MS = min(max(4+4*randn(numcfs,numsponts(2)),0.2),18);
refrand = rand(numcfs,numsponts(2));
tabss.MS = (tabsmax - tabsmin)*refrand + tabsmin;
trels.MS = (trelmax - trelmin)*refrand + trelmin;

% generate sponts, tabss & trels for HS fibers (fiberType = 3)
sponts.HS = min(max(70+30*randn(numcfs,numsponts(3)),18),180);
refrand = rand(numcfs,numsponts(3));
tabss.HS = (tabsmax - tabsmin)*refrand + tabsmin;
trels.HS = (trelmax - trelmin)*refrand + trelmin;

%if exist ('OCTAVE_VERSION', 'builtin') ~= 0
%    save('-mat','ANpopulation.mat','sponts','tabss','trels')
%else
%    save('ANpopulation.mat','sponts','tabss','trels')
%end



