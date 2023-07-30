function jittered=decheveigne2023_spikejitter(spikes,jitter)
%DECHEVEIGNE2023_SPIKEJITTER
%
%   Usage:
%     jittered=decheveigne2023_spikejitter(spikes,jitter)
%
%   Input parameters:
%     spikes  : spike times
%     jitter  : array of jitters or std of Gaussian jitter [default: 0.1 ms]
%
%   Output parameters:
%     jittered: jittered spike times
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/decheveigne2023_spikejitter.php


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

if nargin<2||isempty(jitter); jitter=0.0001; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(jitter)==1;
    jitter = jitter*randn(size(spikes));
end
if numel(jitter)~=numel(spikes); error('!'); end
jittered=spikes+jitter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout==0
    disp('spike_jitter: no output requested, plot');
    subplot 211;
    binwidth=0.00005;
    spike_ach(spikes,binwidth);
    subplot 212
    spike_ach(jittered,binwidth);
    title('jittered')
    clear jittered;
end

end % spike_jitter

