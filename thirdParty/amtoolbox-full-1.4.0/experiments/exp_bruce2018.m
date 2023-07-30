function varargout = exp_bruce2018(varargin)
%EXP_BRUCE2018 Experiments from Bruce et al. 2018
%
%   exp_bruce2011(fig) reproduce Fig. no. fig from the Bruce et
%   al. 2018 paper.
%
%   The following flags can be specified;
%
%     'fig7b'    Reproduce Fig. 7 panel b.
%
%     'fig8b'    Reproduce Fig. 8, panel b.
%
%     'fig10'    Reproduce Fig. 10.
%
%
%   Examples:
%   ---------
%
%   To display Fig. 7b use :
%
%     exp_bruce2018('fig7b');
%
%   To display Fig. 8b use :
%
%     exp_bruce2018('fig8b');
%
%   To display Fig. 10 use :
%
%     exp_bruce2018('fig10');
%
%   See also: bruce2018
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_bruce2018.php


%   #Author: Ian Bruce
%   #Author: Clara Hollomey (2021): adaptations for AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

  
definput.import={'amt_cache'};
definput.flags.type={'missingflag', ...
    'fig7b', ... % 'popratelevel'
    'fig8b', ... % 'fanofactor'
    'fig10'}; 

definput.flags.plot={'plot','no_plot'};

flags  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%%% Common parameters:
Fs = 100e3;    % sampling rate in Hz (must be 100, 200 or 500 kHz)
dt=1/Fs; %  each time in length

%% demo_bruce2018_popratelevel.m
if flags.do_fig7b
    
    numcfs = 1;
    CFs   = 8e3;  % CF in Hz;

    numsponts = [1 1 3]; % reduce the number of ANFs to have a shorter simulation time
    %numsponts = [20 20 60];
    stimdbs = -10:5:100;
    
    numstims = length(stimdbs);
    rates = zeros(numcfs,sum(numsponts),numstims);

    cohc  = 1.0;   % normal ohc function
    cihc  = 1.0;   % normal ihc function
    species = 'cat';   % cat or human
    implnt = 'approxPL';    % 'approxPL', or'actualPL'
    noiseType = 'varFGn'; % 'fixedFGn', or'varFGn'
    
    % stimulus parameters
    T  = 50e-3;  % stimulus duration in seconds
    rt = 2.5e-3; % rise/fall time in seconds

    % PSTH parameters
    psthbinwidth = 0.5e-3; % binwidth in seconds;
    nrep = 100;  % number of stimulus repetitions - Liberman (1978) used 10;

    t = 0:dt:T-dt; % time vector
    mxpts = length(t);
    irpts = rt*Fs;
    
    pin = sqrt(2)*sin(2*pi*CFs*t); % unramped stimulus
    pin(1:irpts) = pin(1:irpts).*(0:(irpts-1))/irpts;
    pin((mxpts-irpts):mxpts) = pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
    
    figure
    for ii = 1:numstims
       pin=scaletodbspl(pin,stimdbs(ii));
       output = bruce2018(pin, round(1/dt), CFs, 'nrep', nrep, ...
           'numL', numsponts(1), 'numM',numsponts(2) , 'numH',numsponts(3) ,...
           'psthbinwidth_mr', psthbinwidth, 'reptime', 2, ...
           noiseType, implnt, 'cohcs', cohc, 'cihcs', cihc, species, 'outputPerSynapse');
    
       psth = output.psth_ft;
                
       psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
       pr = zeros(size(psth, 3), length(psth)/psthbins);
       psTH = zeros(size(pr));
       
       ronset = round(1.5e-3/psthbinwidth)+1;
       roffset = round(T/psthbinwidth);
    
       for jj = 1: size(psth, 3)
         % pr of spike in each bin
         pr(jj,:) = sum(reshape(psth(:,:,jj),psthbins,length(psth)/psthbins))/nrep; 
         psTH(jj,:) = pr(jj,:)/psthbinwidth; % psth in units of /s
         rates(1,jj, ii)= mean(psTH(jj,ronset:ronset+roffset));
       end
    end

 for spontlp = 1:sum(numsponts)
    if (spontlp<=numsponts(1))
      plot(stimdbs,squeeze(rates(1,spontlp,:)),'r')
      hold on
    elseif (spontlp<=sum(numsponts([1 2])))
      plot(stimdbs,squeeze(rates(1,spontlp,:)),'b')
      hold on
    else
      plot(stimdbs,squeeze(rates(1,spontlp,:)),'m')
      hold on
    end
 end
    xlabel('Stimulus Level (dB SPL)')
    ylabel('Firing Rate (/s)')
end    

%% demo_bruce2018_fanofactor.m
if flags.do_fig8b
    
    % model fiber parameters

    %default nerve fiber parameters of the models, therefore not explicitly
    %passed to bruce2018
    spont = 50;      % spontaneous firing rate
    tabs   = 0.6e-3; % absolute refractory period
    trel   = 0.6e-3; % baseline mean relative refractory period
    CF    = 1.5e3;   % CF in Hz;
    cohc  = 1.0;     % normal ohc function
    cihc  = 1.0;     % normal ihc function
    species = 'cat';     % cat or human (with Shera et al. tuning)
    noiseType = 'varFGn';   % variable fGn; for fixed (frozen) fGn: 'fixedFGn'
    implnt = 'approxPL';      % approximate 
                              %for actual implementation of the power-law functions in the Synapse: 'actualPL'

    % stimulus parameters
    F0 = CF;     % stimulus frequency in Hz
    T  = 25;  % stimulus duration in seconds
    rt = 2.5e-3; % rise/fall time in seconds

    stimdb = -inf; % stimulus intensity in dB SPL; set to -inf to get spont activity

    numsponts = 10;

    numTs = 14;
    Ts = logspace(log10(1e-3),log10(10),numTs);
    Ts = round(Ts/dt)*dt;

    Ft = zeros(numsponts,numTs);
    Ft_shuf = zeros(numsponts,numTs);
    meanrate = zeros(numsponts,numTs);

    % PSTH parameters

    nrep = 1;               % number of stimulus repetitions (e.g., 50);
    t = 0:1/Fs:T-1/Fs; % time vector
    mxpts = length(t);
    irpts = rt*Fs;

    pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
    pin(1:irpts)= pin(1:irpts).*(0:(irpts-1))/irpts;
    pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
    pin=scaletodbspl(pin,stimdb);
    output = bruce2018(pin, round(1/dt),  CF, 'nrep', nrep, 'fsmod', Fs, ...
            noiseType, implnt, 'specificSR','numsponts', numsponts, 'spont',spont,'tabs',tabs,'trel',trel,...
          'cohcs', cohc, 'cihcs', cihc, 'reptime', 2, species, 'outputPerSynapse');

for trial = 1:numsponts
        psth = squeeze(output.psth_ft(:,1, trial));        
        simtime = length(psth)/Fs;
        tvect = 0:1/Fs:simtime-1/Fs;

        ISIs = diff(tvect(logical(psth)));
        ISIs_shuf = ISIs(randperm(length(ISIs)));
        spiketimes_shuf = cumsum(ISIs_shuf);
        psth_shuf = histc(spiketimes_shuf,tvect);

        for Tlp = 1:numTs

            psthbinwidth = Ts(Tlp);
            psthbins = round(psthbinwidth*Fs);  % number of time bins per Psth bin
            numPsthbins = floor(length(psth)/psthbins);

            Psth = sum(reshape(psth(1:psthbins*numPsthbins),psthbins,numPsthbins)); %

            Psth_shuf = sum(reshape(psth_shuf(1:psthbins*numPsthbins),psthbins,numPsthbins));

            Ft(trial,Tlp) = std(Psth)^2/(mean(Psth)+eps);
            Ft_shuf(trial,Tlp) = std(Psth_shuf)^2/(mean(Psth_shuf)+eps);
            meanrate(trial,Tlp) = mean(Psth)/Ts(Tlp);

        end

end

    figure
    % Plot the calculated Fano factor curve for each trial
    loglog(Ts*1e3,Ft) 
    hold on
    % Plot the mean Fano factor curve
    loglog(Ts*1e3,mean(Ft),'k-','linewidth',2) 
    % Plot the calculated Fano factor curves for the shuffled ISIs for each trial
    loglog(Ts*1e3,Ft_shuf,'--') 
    % Plot the calculated Fano factor curves for the shuffled ISIs for each trial
    loglog(Ts*1e3,mean(Ft_shuf),'k--','linewidth',2) 
    ylabel('F(T)')
    xlabel('T (ms)')
    xlim([1e0 1e4])
    ylim([0.2 10])    
end

%% FIG 10
if flags.do_fig10 % flags.do_analyticalmeanvarcounts

    % Description:
    %%% 1. Stimuli:
    % Pure tone centred at 8 kHz (CF)
    % Duration of 0.25 s (T)
    % Presentation level: 20 dB SPL (stimdb)
    % Sampling frequency of 100 kHz (Fs)
    % Stimuli with 2.5-ms linear ramps (rt)
    % The stimuli start after 25 ms (ondelay)
    %
    %%% 2. Model configuration:
    % Normal-hearing configuration: cohc = cihc = 1
    % Cat model (species = 1)
    % Spontaneous firing rate of 100 spikes/s (spont)
    % tabs = trel = 0.6 ms
    % Variable fractional Gaussian noise (noiseType = 1)
    % Approximate implementation (implnt = 0)
    % PSTHs using bin sizes of 0.5 ms, 5 ms, and 50 ms
    % Stimuli repeated once each time the model is called (nrep = 1), but the 
    %   model is run 1000 trials (ntrials). PSTHs after the ntrials are assessed
    
    % model parameters
    CF    = 8e3;   % CF in Hz;
    cohc  = 1.0;    % normal ohc function
    cihc  = 1.0;    % normal ihc function
    species = 'cat';    % cat or human
    noiseType = 1;   % variable fGn; or fixed (frozen) fGn
    spont = 100; % spontaneous firing rate 
    tabs   = 0.6e-3; % Absolute refractory period
    trel   = 0.6e-3; % Baseline mean relative refractory period
    implnt = 0;     % "0" for approximate or "1" for actual power-law implementation

    % stimulus parameters
    F0 = CF;     % stimulus frequency in Hz
    T  = 0.25;  % stimulus duration in seconds
    rt = 2.5e-3; % rise/fall time in seconds
    stimdb = 20; % stimulus intensity in dB SPL

    % PSTH parameters
    nrep = 1;               % number of stimulus repetitions (e.g., 50);
    psthbinwidths = [5e-4 5e-3 5e-2]; % binwidth in seconds;
    numpsthbinwidths = length(psthbinwidths);
    trials = 1e3;
    % trials = 10e3; % higher number of trials gives more accurate estimates but takes longer to run

    t = 0:dt:T-dt; % time vector
    mxpts = length(t);
    irpts = rt*Fs;

    ondelay = 25e-3;

    onbin = round(ondelay*Fs);

    pin = zeros(1,onbin+mxpts);

    pin(onbin+1:onbin+mxpts) = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
    pin(onbin+1:onbin+irpts)= pin(onbin+1:onbin+irpts).*(0:(irpts-1))/irpts;
    pin(onbin+(mxpts-irpts):onbin+mxpts)=pin(onbin+(mxpts-irpts):onbin+mxpts).*(irpts:-1:0)/irpts;
    pin=scaletodbspl(pin,stimdb); 

    output = bruce2018(pin, round(1/dt), CF, 'nrep', nrep, 'fsmod', Fs, ... % signal parameters
           'varFGn', 'approxPL', species, 'outputPerSynapse', ... % general parameters
           'specificSR','numsponts', 1, 'spont',spont, 'tabs',tabs, 'trel',trel, ... % spontaneous rate parameters
           'cohcs', cohc, 'cihcs', cihc, 'reptime', 2); % hearing loss parameters

       
        
    psth = squeeze(output.psth_ft(:, 1,:));
    vihc = output.vihc;
    meanrate = output.meanrate;
    varrate = output.varrate;
    
    timeout = (0:length(psth)-1)*1/Fs;

    psthbins = zeros(numpsthbinwidths,1);
    psthtime = cell(numpsthbinwidths);
    psths = cell(numpsthbinwidths);
    mrs = cell(numpsthbinwidths);
    vrs = cell(numpsthbinwidths);

    for binlp = 1:numpsthbinwidths
    
        psthbins(binlp) = round(psthbinwidths(binlp)*Fs);  % number of psth bins per psth bin
        psthtime{binlp} = timeout(1:psthbins(binlp):end); % time vector for psth
        cnt = sum(reshape(psth,psthbins(binlp),length(psth)/psthbins(binlp))); % spike cnt in each psth bin
        mr = mean(reshape(meanrate,psthbins(binlp),length(psth)/psthbins(binlp))); % mean average of theor spike cnt in each psth bin
        vr = mean(reshape(varrate,psthbins(binlp),length(psth)/psthbins(binlp))); % mean var of theor spike cnt in each psth bin
    
        psths{binlp}(1,:) = cnt;
        mrs{binlp}(1,:) = mr;
        vrs{binlp}(1,:) = vr;
    
    end

    for lp = 2:trials
    
        amt_disp(['lp = ' int2str(lp) '/' int2str(trials)], 'volatile');
    
        [psth,meanrate,varrate,~,~,~] = bruce2018_synapse(vihc,CF,nrep,dt,noiseType,implnt,spont,tabs,trel);
    
        timeout = (0:length(psth)-1)*1/Fs;
    
        for binlp = 1:numpsthbinwidths
        
            psthbins(binlp) = round(psthbinwidths(binlp)*Fs);  % number of psth bins per psth bin
            psthtime{binlp} = timeout(1:psthbins(binlp):end); % time vector for psth
            cnt = sum(reshape(psth,psthbins(binlp),length(psth)/psthbins(binlp))); % spike cnt in each psth bin
            mr = mean(reshape(meanrate,psthbins(binlp),length(psth)/psthbins(binlp))); % mean average of theor spike cnt in each psth bin
            vr = mean(reshape(varrate,psthbins(binlp),length(psth)/psthbins(binlp))); % mean var of theor spike cnt in each psth bin
        
            psths{binlp}(lp,:) = cnt;
            mrs{binlp}(lp,:) = mr;
            vrs{binlp}(lp,:) = vr;
        
        end
    
    end
    amt_disp();
    
    for binlp = 1:numpsthbinwidths
    
        if psthbinwidths(binlp)>10e-3
            mrksize = 6;
        else
            mrksize = 2;
        end
        figure
        subplot(2,1,1)
        h1 = bar(psthtime{binlp},mean(psths{binlp}),'histc');
        set(h1,'edgecolor','k','facecolor',0.8*ones(1,3))
        ylabel('E[count]')
        xlabel('Time (s)')
        hold on
        plot(psthtime{binlp}+psthbinwidths(binlp)/2,mean(mrs{binlp})*psthbinwidths(binlp),'ro','markerfacecolor','r','markersize',mrksize,'linewidth',1)
        xlim(psthtime{binlp}([1 end]))
        title(['PSTH bin width = ' num2str(psthbinwidths(binlp)*1e3,2) 'ms'])
        subplot(2,1,2)
        h2 = bar(psthtime{binlp},var(psths{binlp}),'histc');
        set(h2,'edgecolor','k','facecolor',0.8*ones(1,3))
        xlabel('Time (s)')
        ylabel('var[count]')
        hold on
        plot(psthtime{binlp}+psthbinwidths(binlp)/2,mean(vrs{binlp})*psthbinwidths(binlp),'go','markerfacecolor','g','markersize',mrksize,'linewidth',1)
        plot(psthtime{binlp}+psthbinwidths(binlp)/2,mean(mrs{binlp})*psthbinwidths(binlp),'bo','markerfacecolor','b','markersize',mrksize,'linewidth',1)
        xlim(psthtime{binlp}([1 end]))
    
    end  
end




