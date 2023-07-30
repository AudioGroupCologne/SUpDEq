function [h,bins]=decheveigne2023_spikepsth(spikes,binwidth, stim_times)
%DECHEVEIGNE2023_SPIKEPSTH
%
%   Usage:
%     [h,bins]=decheveine2023_spikepsth(spikes,binwidth,stim_times) - peri-stimulus histogram
%
%   Input parameters:
%     spikes   : spike times
%     binwidth : width of bins to count intervals [default 0.0001 s]
%                stim_times (array or single number) [default 0.01 s]
%
%   Output parameters:
%     h    : histogram
%     bins : bin centers
%
% If stim_times is a single number, it is treated as the stimulus period.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/decheveigne2023_spikepsth.php


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
if nargin<3||isempty(stim_times); stim_times=0.01; end
    
if numel(stim_times)==1;
    stim_times=0:stim_times:max(spikes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stim_times=[stim_times,max(spikes)]; 
D=max(diff(stim_times));
bins=linspace(0,D,ceil(D/binwidth));
h=zeros(size(bins));
for iStim=1:numel(stim_times)-1
    segment=spikes(spikes >= stim_times(iStim) & spikes < stim_times(iStim+1));
    segment=segment-stim_times(iStim);
    h=h+hist(segment,bins);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout==0
    disp('spike_psth: no output requested, plot');
    stairs(bins,h); 
    xlim([bins(1) bins(end)])
    xlabel('time re stim (s)'); ylabel('count per bin'); title('PST histogram');
    clear h;
end

end % spike_psth
