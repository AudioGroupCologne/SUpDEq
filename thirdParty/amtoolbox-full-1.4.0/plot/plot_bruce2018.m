function h = plot_bruce2018(t,CFs,neurogram,varargin)
%PLOT_BRUCE2018 plots the output of the bruce 2018 model
%
%   Usage:
%     h = plot_bruce2018(t,CFs,neurogram)
%
%   Input parameters:
%     t            : time index vector
%     CFs          : characteristic frequencies vector
%     neurogram    : neurogram coefficients [time CF]
%
%   This function provides a framework for plotting neurograms.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_bruce2018.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB MEX M-Signal
%   #Author: Ian Bruce: basic code of the model
%   #Author: Alejandro Osses (2020): original implementation
%   #Author: Clara Hollomey (2021): adapted to the AMT 1.0
%   #Author: Piotr Majdak (2021): adaptations to exp_osses2022; specificSRautoTiming added

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if (nargin > 4)
    error('Too many input arguments')
elseif (nargin > 3)
    hin = varargin{1};
    
    if strncmp(get(hin,'Type'),'figu',4)
        
        figure(hin)
        h = axes;
        
    elseif strncmp(get(hin,'Type'),'axes',4)
        
        h = hin;
        
    else
        error(['Handle type of ' get(hin,'Type') ' not supported'])
    end
    
else
    
    h = axes;
    
end

axes(h)
imagesc(t,log10(CFs/1e3),neurogram')
axis xy
yticks = [0.125 0.5 2 8 16];
set(h,'ytick',log10(yticks))
set(h,'yticklabel',yticks)
ylabel('CF (kHz)')
xlabel('Time (s)')
hcb = colorbar;
set(get(hcb,'ylabel'),'string','spikes')

end



