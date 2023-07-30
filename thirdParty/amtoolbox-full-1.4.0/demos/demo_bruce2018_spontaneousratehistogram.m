%DEMO_BRUCE2018_SPONTANEOUSRATEHISTOGRAM  plot the spontaneous firing rate histogram
%
%   DEMO_BRUCE2018_SPONTANEOUSRATEHISTOGRAM yields the histogram of the
%   synapse firing rates
%
%   Figure 1: Histogram of the spontaneous firing rates
% 
%
%   See also: bruce2018 exp_bruce2018 demo_bruce2018_bestthresholdcurve
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_bruce2018_spontaneousratehistogram.php


%   AUTHOR : Ian Bruce

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% model fiber parameters
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
species = 1;   % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
noiseType = 1; % 1 for variable fGn; 0 for fixed (frozen) fGn
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

trials = 100; 738;

numcfs = 30;

CFs   = logspace(log10(300),log10(20e3),numcfs);  % CF in Hz;

numsponts = [10 10 30];


[sponts,tabss,trels] = bruce2018_generateanpopulation(numcfs,numsponts);



% stimulus parameters
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 30;  % stimulus duration in seconds

% PSTH parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

t = 0:1/Fs:T-1/Fs; % time vector

vihc = zeros(1,length(t));

trial=1;
m = zeros(1,trials);
p = zeros(1,trials);
spontval = zeros(1,trials);

meanrate = [];

while trial<=trials
    
    fprintf(1,'trial = %i/%i',trial,trials);
    
    % flush the output for the display of the coutput in Octave
    if exist ('OCTAVE_VERSION', 'builtin') ~= 0
        fflush(stdout);
    end
    
    cfind = ceil(numcfs*rand(1));
    
    CF = CFs(cfind);
    
    sponts_concat = [sponts.LS(cfind,1:numsponts(1)) sponts.MS(cfind,1:numsponts(2)) sponts.HS(cfind,1:numsponts(3))];
    tabss_concat = [tabss.LS(cfind,1:numsponts(1)) tabss.MS(cfind,1:numsponts(2)) tabss.HS(cfind,1:numsponts(3))];
    trels_concat = [trels.LS(cfind,1:numsponts(1)) trels.MS(cfind,1:numsponts(2)) trels.HS(cfind,1:numsponts(3))];
    
    spontind = ceil(sum(numsponts)*rand(1));
    
    spont = sponts_concat(spontind);
    tabs = tabss_concat(spontind);
    trel = trels_concat(spontind);
    
    psth = bruce2018_synapse(vihc,CF,nrep,1/Fs,noiseType,implnt,spont,tabs,trel);
    
    sptimes= find (psth==1)/Fs;
    nspikes=length(sptimes);
    
    ISI = diff(sptimes); % Compute ISIs from spike times
    N = length(ISI);
    
    display([':  ' num2str(N) ' spikes in this trial']);
    if exist ('OCTAVE_VERSION', 'builtin') ~= 0
        fflush(stdout);
    end
    
    m(trial) = mean(ISI);
    trial= trial+1;
    spontval(trial) = spont;
    
end

edges = 0:120;
n = histc(1./m,edges);

figure
bar(edges,n)
xlabel('Spontaneous Rate (/s)')
ylabel('Number of units')
xlim([-0.5 120])


