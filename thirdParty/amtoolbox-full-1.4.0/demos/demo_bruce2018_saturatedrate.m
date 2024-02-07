%DEMO_BRUCE2018_SATURATEDRATE plots the fibre discharge rate at saturation
%
%   demo_bruce2018_SATURATEDRATE depicts the firing behaviour close to the best
%   threshold curve
%
%   Figure 1: The discharge rate at saturation as a function of characteristic frequency
% 
%
%   See also: bruce2018 exp_bruce2018 demo_bruce2018_bestthresholdcurve
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_bruce2018_saturatedrate.php


%   AUTHOR : Ian Bruce

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

numcfs = 30;
CFs   = logspace(log10(300),log10(20e3),numcfs);  % CF in Hz;

numsponts_healthy = [0 0 4];
numsponts = numsponts_healthy(3);

[sponts,tabss,trels] = bruce2018_generateanpopulation(numcfs,numsponts_healthy);


numstims = 55;

rates = zeros(numcfs,numsponts,numstims);

cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
species = 1;   % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
noiseType = 1; % 1 for variable fGn; 0 for fixed (frozen) fGn
Fs = 100e3;    % sampling rate in Hz (must be 100, 200 or 500 kHz)

% stimulus parameters
T  = 50e-3;  % stimulus duration in seconds
rt = 2.5e-3; % rise/fall time in seconds

% PSTH parameters
psthbinwidth = 0.5e-3; % binwidth in seconds;
nrep = 10;  % number of stimulus repetitions - Liberman (1978) used 10;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;


for cflp = 1:numcfs
    
    CF = CFs(cflp);
    
    F0 = CF;  % stimulus frequency in Hz
    
    for spontlp = 1:numsponts
        
        spont = sponts.HS(cflp,spontlp);
        tabs  = tabss.HS(cflp,spontlp);
        trel  = trels.HS(cflp,spontlp);
        
        %thrsh = bruce2018_findcfthreshold(CF,Fs,cohc,cihc,species,noiseType,implnt,spont,tabs,trel);
 %-------------------------------------------------------------------------
 % find CF Threshold using the STB tone and incremental intensity by 1dB until when the firing rate for
% a specific fiber passed 10 plus the spontrate of the fiber
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
        thrsh = stimdb_in;
 %-------------------------------------------------------------------------     
        stimdbs = thrsh:thrsh+54;
        
        for stimlp = 1:numstims
            
            disp(['CFlp = ' int2str(cflp) '/' int2str(numcfs) '; spontlp = ' int2str(spontlp) '/' int2str(sum(numsponts)) '; stimlp = ' int2str(stimlp) '/' int2str(sum(numstims))])
            
            % flush the output for the display of the coutput in Octave
            if exist ('OCTAVE_VERSION', 'builtin') ~= 0
                fflush(stdout);
            end
            
            stimdb_in = stimdbs(stimlp);
            pin = sqrt(2)*20e-6*10^(stimdb_in/20)*sin(2*pi*F0*t); % unramped stimulus
            pin(1:irpts) = pin(1:irpts).*(0:(irpts-1))/irpts;
            pin((mxpts-irpts):mxpts) = pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
            
            vihc = bruce2018_innerhaircells(pin,CF,nrep,1/Fs,T*2,cohc,cihc,species);
            psth = bruce2018_synapse(vihc,CF,nrep,1/Fs,noiseType,implnt,spont,tabs,trel);
            
            timeout = (0:length(psth)-1)*1/Fs;
            psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
            psthtime = timeout(1:psthbins:end); % time vector for psth
            pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
            psTH = pr/psthbinwidth; % psth in units of /s
            
            ronset = round(1.5e-3/psthbinwidth)+1;
            roffset = round(T_in/psthbinwidth);
            
            rates(cflp,spontlp,stimlp)= mean(psTH(ronset:ronset+roffset));
            
        end
        
    end
    
end

figure
semilogx(CFs/1e3,max(rates,[],3),'ko')
xlim([0.1 40])
ylim([100 350])
xlabel('Characteristic Frequency (kHz)')
ylabel('Discharge Rate at Saturation (/s)')




