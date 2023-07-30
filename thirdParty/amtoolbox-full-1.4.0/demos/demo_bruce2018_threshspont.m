%DEMO_BRUCE2018_THRESHSPONT plot the relative threshold as a function of the spontaneous firing rate
%
%   DEMO_BRUCE2018_THRESHSPONT depicts the firing behaviour close to the best
%   threshold curve
%
%   Figure 1: The relative firing threshold as a function of characteristic frequency
%
%   Figure 2: The relative firing threshold as a function of the adjusted spontaneous firing rate
% 
%
%   See also: bruce2018 exp_bruce2018 demo_bruce2018
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_bruce2018_threshspont.php


%   #AUTHOR : Ian Bruce

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
amt_disp('Warning: running this demo may take more than one hour.');
% numcfs = 40; % Number of CF bins from Liberman 1978
numcfs = 10; % Use a smaller number of CF bins for a quicker simulation
CFs   = logspace(log10(400),log10(15e3),numcfs);  % CF in Hz (range for Fig. 10 of Liberman 1978)

numsponts = [4 4 12];


[sponts,tabss,trels] = bruce2018_generateanpopulation(numcfs,numsponts);


cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
species = 1;   % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
noiseType = 1; % 1 for variable fGn; 0 for fixed (frozen) fGn
Fs = 100e3;    % sampling rate in Hz (must be 100, 200 or 500 kHz)

% T  = 10;  % duration in seconds over which the SR is calculated
T  = 1;  % use a shorter duration for a quicker simulation
nrep = 1;               % number of stimulus repetitions (e.g., 50);

t = 0:1/Fs:T-1/Fs; % time vector

vihc = zeros(1,length(t));

thrsh = zeros(numcfs,sum(numsponts));
thrsh_meanHS = zeros(numcfs,1);
SR = zeros(numcfs,sum(numsponts));
spontval = zeros(numcfs,sum(numsponts));

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
        
        psth = bruce2018_synapse(vihc,CF,nrep,1/Fs,noiseType,implnt,spont,tabs,trel);
        
        SR(cflp,spontlp) = sum(psth)/T;
        
        spontval(cflp,spontlp) = spont;
        
        %thrsh(cflp,spontlp) = bruce2018_findcfthreshold(CF,Fs,cohc,cihc,species,noiseType,implnt,spont,tabs,trel);
 %-------------------------------------------------------------------------
 %findcfthreshold
 % find CF Threshold using the STB tone and incremental intensity by 1dB until when the firing rate for
% a specific fiber passed 10 plus the spontrate of the fiber
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
    
    thrsh_meanHS(cflp) = mean(thrsh(cflp,SR(cflp,:)>18));
    
end

figure
semilogx(CFs/1e3,thrsh,'kx')
hold on
semilogx(CFs/1e3,thrsh_meanHS,'r-')
ylabel('Threshold (dB SPL)')
xlabel('CF (kHz)')
xlim([0.1 20])
set(gca,'xtick',[0.1 1 10])
set(gca,'xticklabel',[0.1 1 10])

relthrsh = thrsh - repmat(thrsh_meanHS,1,sum(numsponts));

p = polyfit(log10(max(0.1,SR(SR<=18))),relthrsh(SR<=18),1);

p_allANFs = polyfit(log10(max(0.1,SR)),relthrsh,1);

figure
semilogx(max(0.1,SR(SR<=18)),relthrsh(SR<=18),'b^')
hold on
semilogx(SR(SR>18),relthrsh(SR>18),'r^')
hold on
semilogx(logspace(-1,2,100), p(1)*log10(logspace(-1,2,100))+ p(2),'b-','linewidth',2.0)
text(0.15,-5,['thrsh = ' num2str(p(1),3) '*log10(spont)+' num2str(p(2),3)])
legend('Low & medium spont fibers','High spont fibers','Fit to low & medium spont fibers')
xlabel('Adjusted Spont Rate (/s)')
ylabel('Relative Threshold (dB)')
xlim([0.1 150])
set(gca,'xtick',[0.1 1 10 100])
set(gca,'xticklabel',[0.1 1 10 100])


