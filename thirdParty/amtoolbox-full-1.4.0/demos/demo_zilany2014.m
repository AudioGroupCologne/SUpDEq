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
%     Journal of the Acoustical Society of America, 135(1):283--286, Jan.
%     2014.
%     
%     M. Zilany, I. Bruce, P. Nelson, and L. Carney. A phenomenological model
%     of the synapse between the inner hair cell and auditory nerve:
%     Long-term adaptation with power-law dynamics. J. Acoust. Soc. Am.,
%     126(5):2390 -- 2412, 2009.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_zilany2014.php


%   #Author: Peter L. Soendergaard (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%Parameter Settings
display_level = 'no_debug'; % set to 'debug' to see more information, set to 'no_debug' to have less mess on your display

% model fiber parameters
CF    = 1.5e3;   % CF in Hz;   
fiberType=4; % Simulate a neuron with the following SR: 1=Low; 2=Medium; 3=High; 4: custom, see numH, numM and numL
numH=12; % # of high SR neurones, if fiberType=4
numM=4; % # of medium SR neurones, if fiberType=4
numL=4; % # of low SR neurones, if fiberType=4

% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
fsstim = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 50e-3;  % stimulus duration in seconds
rt = 2.5e-3; % rise/fall time in seconds
stimdb = 65; % stimulus intensity in dB SPL

% peri-stimulus time histogram (PSTH) parameters
nrep = 50;               % number of stimulus repetitions (e.g., 50);
psth_binwidth = 0.0005;  % binwidth in seconds, set to [] for raw PSTH


%% Computations

% Stimulus generation
t = 0:1/fsstim:T-1/fsstim; % time vector
mxpts = length(t);
irpts = rt*fsstim;
stim = scaletodbspl(sin(2*pi*F0*t),stimdb); % unramped stimulus
stim(1:irpts)= stim(1:irpts).*(0:(irpts-1))/irpts; 
stim((mxpts-irpts):mxpts)=stim((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

% AN modeling
par={'fiberType',fiberType,'numH',numH,'numM',numM,'numL',numL,'nrep',nrep,'psth_binwidth',psth_binwidth};
[r_mean,psth,ihc] = zilany2014(stim,fsstim,CF,par{:},display_level);


%% Plots
timestim = (1:length(r_mean))*1/fsstim;
if isempty(psth_binwidth), psth_binwidth=1/fsstim; end
psthbins = round(psth_binwidth*fsstim);  % number of PSTH bins per PSTH bin
timebins = timestim(1:psthbins:end); % time vector for psth

figure
subplot(4,1,1)
plot(timestim,[stim zeros(1,length(timestim)-length(stim))])
title('Input Stimulus')
ylabel('Pascal')

subplot(4,1,2)
plot(timestim,ihc(1:length(timestim)))
title('IHC Output')
ylabel('Volts')

subplot(4,1,3)
plot(timestim,r_mean);
xl = xlim;
title('Mean Rate Output')
ylabel('spikes/s')

subplot(4,1,4)
bar(timebins,psth)
xlim(xl)
title('Peri-stimulus Time Histogram')
xlabel('Time (s)')
ylabel('spikes/s')



