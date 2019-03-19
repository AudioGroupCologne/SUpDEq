%DEMO_ZILANY2014  Demo of the Zilany et al., (2014) model
%
%   This demos generates a simple figure that shows the behaviour of the Zilany et al. (2014) model
% 
%   Figure 1: Figure from Zilany et al. (2014) model
%
%
%   References:
%     M. S. A. Zilany, I. C. Bruce, and L. H. Carney. Updated parameters and
%     expanded simulation options for a model of the auditory periphery. The
%     Journal of the Acoustical Society of America, 135(1):283-286, Jan.
%     2014.
%     
%     M. Zilany, I. Bruce, P. Nelson, and L. Carney. A phenomenological model
%     of the synapse between the inner hair cell and auditory nerve:
%     Long-term adaptation with power-law dynamics. J. Acoust. Soc. Am.,
%     126(5):2390 - 2412, 2009.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/demos/demo_zilany2014.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%Parameter Settings

% model fiber parameters
CF    = 1.5e3;   % CF in Hz;   
fiberType = 1;  % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High

% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
fsstim = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 50e-3;  % stimulus duration in seconds
rt = 2.5e-3; % rise/fall time in seconds
stimdb = 65; % stimulus intensity in dB SPL

% peri-stimulus time histogram (PSTH) parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3;  % binwidth in seconds;


%% Computations

% Stimulus generation
t = 0:1/fsstim:T-1/fsstim; % time vector
mxpts = length(t);
irpts = rt*fsstim;
stim = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
stim(1:irpts)= stim(1:irpts).*(0:(irpts-1))/irpts; 
stim((mxpts-irpts):mxpts)=stim((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

% AN modeling
[ANresp,fc,vihc,psth] = zilany2014(...
  stimdb,stim,fsstim,...
  'flow',CF','fhigh',CF,'nfibers',1,'fiberType',fiberType);

% PSTH conversion
timeout = (1:length(psth))*1/fsstim;
psthbins = round(psthbinwidth*fsstim);  % number of psth bins per psth bin
psthtime = timeout(1:psthbins:end); % time vector for psth
pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
Psth = pr/psthbinwidth; % psth in units of spikes/s


%% Plots

figure
subplot(4,1,1)
plot(timeout,[stim zeros(1,length(timeout)-length(stim))])
title('Input Stimulus')
ylabel('Pascal')

subplot(4,1,2)
plot(timeout,vihc(1:length(timeout)))
title('IHC Output')
ylabel('Volts')

subplot(4,1,3)
plot(timeout,ANresp);
xl = xlim;
title('Mean Rate Output')
ylabel('spikes/s')

subplot(4,1,4)
bar(psthtime,Psth)
xlim(xl)
title('Peri-stimulus Time Histogram')
xlabel('Time (s)')
ylabel('spikes/s')

