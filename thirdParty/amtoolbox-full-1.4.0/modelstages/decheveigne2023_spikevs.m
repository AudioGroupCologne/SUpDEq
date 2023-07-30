function v=decheveigne2023_spikevs(spikes,T)
% DECHEVEIGNE2023_SPIKEVS vector strength
%
%   Usage:
%     v=decheveigne2023_spikevs(spikes,T);
%
%   Input parameters:
%     spikes : (s) spike times
%     T      : (s) period
%
%   Output parameters:
%     v      : vector strength
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/decheveigne2023_spikevs.php


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

if nargin==0; test_code; return; end

if nargin<2; error('!'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phases=(spikes-T*floor(spikes/T))/T*2*pi;
v = sqrt(sum(cos(phases)).^2 + sum(sin(phases)).^2)/numel(phases);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~nargout
    disp('spike_vs: no output requested, print'); 
    disp(['vector strength: ', num2str(v)]);
end

end % spike_vs
