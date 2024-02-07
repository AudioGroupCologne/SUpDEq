%DEMO_BRUCE2018_BESTTHRESHOLDCURVE estimates the firing threshold [dB SPL]
%
%   DEMO_BRUCE2018_BESTTHRESHOLDCURVE finds the firing threshold for a given
%   center frequency by incrementing the stimulus intensity by 1 dB until the 
%   firing rate increases the fibres' spontaneous rate plus a margin of 10.
%
%   Figure 1: Composite best threshold curve (CBTC) based on new curve (NBTC) from Miller et al. 1997 above 1 kHz and the Liberman 1978 curve (LBTC)
% 
%
%   See also: bruce2018 exp_bruce2018 demo_bruce2018
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_bruce2018_bestthresholdcurve.php


%   #Author : Ian Bruce

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


numcfs = 20;
% numcfs = 40; % increase to get a better sampling of CFs

CFs   = logspace(log10(125),log10(15e3),numcfs);  % CF in Hz (range to check BTC vs CF)
numsponts = [0 0 15]; % Just use high spont fibers to check BTC vs CF

[sponts,tabss,trels] = bruce2018_generateanpopulation(numcfs,numsponts);


cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
species = 1;   % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
noiseType = 1; % 1 for variable fGn; 0 for fixed (frozen) fGn
Fs = 100e3;    % sampling rate in Hz (must be 100, 200 or 500 kHz)

thrsh = zeros(numcfs,sum(numsponts));
thrsh_meanHS = zeros(numcfs,1);
thrsh_stdHS = zeros(numcfs,1);

for cflp = 1:numcfs
    
    CF = CFs(cflp);
    
    F0 = CF;  % stimulus frequency in Hz
    
    for spontlp = 1:sum(numsponts)
        
        disp(['CFlp = ' int2str(cflp) '/' int2str(numcfs) '; spontlp = ' int2str(spontlp) '/' int2str(sum(numsponts))])
        
        % flush the output for the display of the coutput in Octave
        if exist ('OCTAVE_VERSION', 'builtin') ~= 0
            fflush(stdout);
        end
        
        sponts_concat = [sponts.LS(cflp,1:numsponts(1)) sponts.MS(cflp,1:numsponts(2)) sponts.HS(cflp,1:numsponts(3))];
        tabss_concat = [tabss.LS(cflp,1:numsponts(1)) tabss.MS(cflp,1:numsponts(2)) tabss.HS(cflp,1:numsponts(3))];
        trels_concat = [trels.LS(cflp,1:numsponts(1)) trels.MS(cflp,1:numsponts(2)) trels.HS(cflp,1:numsponts(3))];
        
        spont = sponts_concat(spontlp);
        tabs = tabss_concat(spontlp);
        trel = trels_concat(spontlp);
        
 %-------------------------------------------------------------------------
 %findcfthreshold
 %-------------------------------------------------------------------------
        stimdb_in = -10;
        F0_in=CF;
        psthbinwidth_in = 0.5e-3; % binwidth in seconds;
        nrep_in = 200;  % number of stimulus repetitions - Liberman (1978) used 10;
        T_in  = 50e-3;  % stimulus duration in seconds
        rt_in = 2.5e-3; % rise/fall time in seconds

        t_in = 0:1/Fs:T_in-1/Fs; % time vector
        mxpts_in = length(t_in);
        irpts_in = rt_in*Fs;

        SpontRate_in = spont;

        firingRate_Icreased_to_in  = SpontRate_in;

        while ((firingRate_Icreased_to_in  <(SpontRate_in + 10)) && (stimdb_in < 50))
            stimdb_in =  stimdb_in+1;
            if exist ('OCTAVE_VERSION', 'builtin') ~= 0
                fflush(stdout);
            end
            pin_in = sqrt(2)*20e-6*10^(stimdb_in/20)*sin(2*pi*F0_in*t_in); % unramped stimulus
            pin_in(1:irpts_in) = pin_in(1:irpts_in).*(0:(irpts_in-1))/irpts_in;
            pin_in((mxpts_in-irpts_in):mxpts_in) = pin_in((mxpts_in-irpts_in):mxpts_in).*(irpts_in:-1:0)/irpts_in;
    

            vihc_in = bruce2018_innerhaircells(pin_in,CF,nrep_in,1/Fs,T_in*2,cohc,cihc,species);
            psth_in = bruce2018_synapse(vihc_in,CF,nrep_in,1/Fs,noiseType,implnt,spont,tabs,trel);
    
            psthbins_in = round(psthbinwidth_in*Fs);  % number of psth bins per psth bin
            pr_in = sum(reshape(psth_in,psthbins_in,length(psth_in)/psthbins_in))/nrep_in; % pr of spike in each bin
            psTH_in = pr_in/psthbinwidth_in; % psth in units of spikes/s
    
            ronset_in =  round(1.5e-3/psthbinwidth_in)+1;
            roffset_in = round(T_in/psthbinwidth_in);
    
            SpontRate_in = mean(psTH_in(roffset_in+1:end));
    
            firingRate_Icreased_to_in = mean(psTH_in(ronset_in:ronset_in+roffset_in));
    
        end
        thrsh(cflp,spontlp) = stimdb_in;
 %-------------------------------------------------------------------------
 %-------------------------------------------------------------------------
    end
    
    thrsh_meanHS(cflp) = mean(thrsh(cflp,:));
    thrsh_stdHS(cflp) = std(thrsh(cflp,:));
    
end

thrsh_CBTC = [20.0000 -4.5000 -4.5000 -4.5000 -4.5000 6.3000 -2.9000 -4.6000 -3.1000];
CFs_CBTC = 1e3*[0.1800 0.8400 1.0000 1.2200 2.4100 4.0500 4.8000 7.2000 10.0000];

figure
semilogx(CFs/1e3,thrsh,'kx')
hold on
semilogx(CFs_CBTC/1e3,thrsh_CBTC,'k--')
ylabel('Threshold (dB SPL)')
xlabel('CF (kHz)')
xlim([0.1 20])
set(gca,'xtick',[0.1 1 10])
set(gca,'xticklabel',[0.1 1 10])



