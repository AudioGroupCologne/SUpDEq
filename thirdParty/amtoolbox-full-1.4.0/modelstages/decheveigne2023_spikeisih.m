function [h,bins]=decheveigne2023_spikeisih(spikes,binwidth)
%DECHEVEIGNE2023_SPIKEISIH
%
%   Usage:
%     [h,bins]=decheveigne2023_spikeisih(spikes,binwidth) - first-order interspike interval histogram
%
%   Input parameters:
%     spikes   : array of spike times
%     binwidth : s, width of histogram bins
%
%   Output parameters:
%     h        : histogram
%     bins     : (s) bin centers
%
%  The size of the histogram may vary depending on the largest interval in 
%  the spike train. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/decheveigne2023_spikeisih.php


%   #StatusDoc: Unknown
%   #StatusCode: Unknown
%   #Verification: Unknown
%   #Requirements: Unknown
%   #Author: Alain de Cheveigne (2023)
%   #Authors: Alejandro Osses (2021): integration in AMT 1.4

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if nargin<2||isempty(binwidth); binwidth=0.001; end % 1ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isis=diff(spikes);
bins=linspace(binwidth/2, max(isis), ceil(max(isis)/binwidth));
h=hist(isis,bins);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout==0
    disp('spike_isi: no output requested, plot');
    stairs(bins,h); 
    xlim([bins(1) bins(end)])
    xlabel('inter-spike interval (s)'); ylabel('count per bin'); title('first order ISI histogram');
    clear h;
end
end % spike_isih
