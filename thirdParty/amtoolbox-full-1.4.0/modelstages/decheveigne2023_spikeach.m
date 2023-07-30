function [h,bins]=decheveigne2023_spikeach(spikes,binwidth,maxinterval)
%DECHEVEIGNE2023_SPIKEACH
%
%   Usage:
%     [h,bins]=decheveigne2023_spikeach(spikes,binwidth,maxinterval) - auto-coincidence histogram
%
%   Input parameters:
%     spikes   : spike times
%     binwidth : width of bins to count intervals [default 0.0001 s]
%                maxinterval [default 0.02 s]
%
%   Output parameters:
%     h    : histogram
%     bins : (s) bin centers
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/decheveigne2023_spikeach.php


%   #StatusDoc: Unknown
%   #StatusCode: Unknown
%   #Verification: Unknown
%   #Requirements: Unknown
%   #Author: Alain de Cheveigne (2023)
%   #Authors: Alejandro Osses (2023): integration in AMT 1.4

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if nargin<2||isempty(binwidth); binwidth=0.0001; end
if nargin<3||isempty(maxinterval); maxinterval=0.02; end
    
[h,bins] = decheveigne2023_spikecch(spikes,spikes,binwidth,maxinterval);
