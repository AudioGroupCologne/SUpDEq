function [h,bins]=decheveigne2023_spikecch(spikes1,spikes2,binwidth,maxinterval)
%DECHEVEIGNE2023_SPIKECCH
%
%   Usage:
%     [h,bins]=decheveigne2023_spikecch(spikes1,spikes2,binwidth,maxinterval) - cross-coincidence histogram
%
%   Input parameters:
%     spikes1,spikes2: spike times
%     binwidth: width of bins to count intervals [default 0.0001 s]
%     maxinterval [default 0.02s]
%
%   Output parameters:
%     h   : histogram
%     bins: (s) bin centers
%
% spike_cch counts only positive intervals (elements of spike2 later than
% those of spike1). To get negative intervals, swap spikes1 & spikes2
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/decheveigne2023_spikecch.php


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

if nargin<3||isempty(binwidth); binwidth=0.0001; end % s
if nargin<4||isempty(maxinterval); maxinterval=0.02; end  % s

spikes1 = spikes1(:); % to column array
spikes2 = spikes2(:); % to column array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shift=0;
isis=[];
while 1
    a=spikes2(shift+1:end)-spikes1(1:end-shift);
    a=a(a<maxinterval);
    if isempty(a); break; end
    a=a(a~=0); % avoids counting identical spikes when used by spike_ach
    isis=[isis; a];
    shift=shift+1;
end
bins=linspace(binwidth/2, maxinterval, ceil(maxinterval/binwidth));
h=hist(isis,bins);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout==0
    disp('spike_cch: no output requested, plot');
    stairs(bins,h); 
    xlim([bins(1) bins(end)])
    xlabel('inter-spike interval (s)'); ylabel('count per bin'); title('cross-coincidence histogram');
    clear h;
end
end % spike_cch
