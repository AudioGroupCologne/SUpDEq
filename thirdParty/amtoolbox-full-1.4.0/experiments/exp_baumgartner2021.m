function data = exp_baumgartner2021(varargin)
%EXP_BAUMGARTNER2021 Simulations of Baumgartner and Majdak (2021)
%   Usage: data = exp_baumgartner2021(flag) 
%
%   EXP_BAUMGARTNER2021(flag) reproduces figures of the study from 
%   Baumgartner and Majdak (2021).
%
%   The following flags can be specified
%
%     'fig2'   Externalization ratings: actual data from psychoacoustic 
%              experiments (closed circles) and simulations of the single-cue 
%              models (open symbols). 
%              Effects of low-frequency modifications tested by Hartmann and 
%              Wittenberg (1996). Exp. I: IID set to zero (bilateral average 
%              magnitude); actual data from their Fig.7, average from N=2. 
%              Exp. II: ipsilateral magnitudes flattened (IID compensated by 
%              contralateral magnitudes); actual data from their Fig.8, average 
%              from N=4. Simulated results for various cues, average from 
%              N=21. 
%              Exp.III: effect of spectral smoothing of low-frequency sounds 
%              presented from various azimuths (left: 0; right: 50); 
%              actual data represents direct-sound condition from 
%              Hassager et al. (2016), average from N=7. Simulated N=21. 
%              Exp.IV: effect of spectral smoothing in high frequencies as a 
%              function of spectral contrast (C=1: natural listener-specific 
%              spectral profile; C=0: flat spectra); actual data calculated from 
%              the paired-comparison data from Baumgartner et al. (2017), N=10 
%              (actual and simulated). 
%              Exp.V: effects of stimulus bandwidth and microphone casing for 
%              various mixes between simulations based on listener-specific 
%              BRIRs (100%) and time-delay stereophony (0%); actual data from 
%              Boyd et al. (2012), N=3 (actual and simulated). 
%              ITE: in-the-ear casing; BTE: behind-the-ear casing; 
%              BB: broadband stimulus; LP: low-pass filtered at 6.5kHz; 
%              Error bars denote standard errors of the mean.
%              Add flag supplements to also show the predictions for the
%              different combination strategies.
%
%     'fig3'   Optimization and performance of single-cue models. 
%              Cue-specific sensitivities used as optimization parameters. 
%              Higher values denote steeper mapping functions. 
%              Simulation errors as the RMS difference between the actual and 
%              simulated externalization ratings. Per experiment, the smallest 
%              error indicates the most informative cue.
%
%     'fig4'   Simulation errors for different decision strategies and cue 
%              combinations show that static combination (WSM) based on monaural 
%              and interaural spectral shape cues (MSS, ISS) performs best. 
%              RMS simulation errors for different strategies and pooled 
%              experimental data, N=54. Error bars denote 95% confidence 
%              intervals estimated via bootstrapping (1000 resamples). 
%              WSM: weighted-sum model; L/M/WTA: looser/median/winner takes all. 
%              Individual cue contributions. Cue abbreviations defined in Tab.1. 
%              Top: Simulation errors for pairwise combinations of considered cues. 
%              Dashed line shows the border of significant difference to the best 
%              pair (MSS and ISS). Bottom: Considered cue pairs with their 
%              respective weights (encoded by brightness). 
%              Add flag supplements to also show the results for evaluations
%              based on rank correlations.
%
%     'hartmann1996'   models experiments from Hartmann & Wittenberg (1996; Fig.7-8)
%                      1st panel: Synthesis of zero-ILD signals. Only the harmonics 
%                      from 1 to nprime had zero interaural level difference; 
%                      harmonics above nprime retained the amplitudes of the baseline  
%                      synthesis. Externalization scores as a function of the boundary  
%                      harmonic number nprime. Fundamental frequency of 125 Hz.
%                      2nd panel: Synthesis of signals to test the ISLD hypothesis. 
%                      Harmonics at and below the boundary retained only the interaural 
%                      spectral level differences of the baseline synthesis. Higher 
%                      harmonics retained left and right baseline harmonic levels. 
%                      Externalization scores as a function of the boundary
%                      frequency.
%
%     'hassager2016'   models experiments from Hassager et al. (2016; Fig.6). 
%                      The mean of the seven listeners perceived sound source 
%                      location (black) as a function of the bandwidth factor 
%                      and the corresponding model predictions (colored). 
%                      The model predictions have been shifted slightly to the right 
%                      for a better visual interpretation. The error bars are one 
%                      standard error of the mean.
%
%     'baumgartner2017'  models experiments from Baumgartner et al. (2017). 
%                        Effect of HRTF spectral contrast manipulations on sound externalization. 
%                        Externalization scores were derived from paired comparisons via 
%                        Bradley-Terry-Luce modeling.
%
%     'boyd2012'   models experiments from Boyd et al. (2012; Fig.1, top).
%                  Average externalization ratings of 1 talker for NH participants 
%                  against mix point as a function of microphone position (ITE/BTE) 
%                  and frequency response (BB/LP). The reference condition (ref) is 
%                  the same as ITE/BB. Error bars show SEM. 
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v0.4.3 or higher from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2017
%
%   3) Statistics Toolbox for Matlab (for some of the figures)
%
%   Examples:
%   ---------
%
%   To display indivdual predictions for all experiments use :
%
%     exp_baumgartner2021('fig2');
%
%   To display summary results for all experiments use :
%
%     exp_baumgartner2021('fig3');
%
%   To display results for different decision strategies use :
%
%     exp_baumgartner2021('fig4');
%
%   References:
%     R. Baumgartner, D. K. Reed, B. Tóth, V. Best, P. Majdak, H. S. Colburn,
%     and B. Shinn-Cunningham. Asymmetries in behavioral and neural responses
%     to spectral cues demonstrate the generality of auditory looming bias.
%     Proceedings of the National Academy of Sciences, 2017. [1]http ]
%     
%     A. W. Boyd, W. M. Whitmer, J. J. Soraghan, and M. A. Akeroyd. Auditory
%     externalization in hearing-impaired listeners: The effect of pinna cues
%     and number of talkers. J. Acoust. Soc. Am., 131(3):EL268--EL274, 2012.
%     [2]www: ]
%     
%     W. M. Hartmann and A. Wittenberg. On the externalization of sound
%     images. J. Acoust. Soc. Am., 99(6):3678--88, June 1996.
%     
%     H. G. Hassager, F. Gran, and T. Dau. The role of spectral detail in the
%     binaural transfer function on perceived externalization in a
%     reverberant environment. J. Acoust. Soc. Am., 139(5):2992--3000, 2016.
%     
%     References
%     
%     1. http://www.pnas.org/content/early/2017/08/16/1703247114.abstract
%     2. http://dx.doi.org/10.1121/1.3687015
%     
%
%   See also: baumgartner2021
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_baumgartner2021.php


%   #Author: Robert Baumgartner (2021), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.import={'amt_cache','baumgartner2021'};
definput.flags.type = {'all','fig2','fig3','fig4','hartmann1996','hassager2016','baumgartner2017','boyd2012','tabSingleCue','tabStrats','diotic'};
definput.flags.gradients={'positive','negative','both'};
definput.flags.quickCheck = {'','quickCheck'};
definput.flags.supplements = {'no_supps','supplements'}; % for plotting supplementary figures
definput.flags.disp = {'no_debug','debug'};
definput.keyvals.Sintra = 2;
definput.keyvals.Sinter = 2;
definput.keyvals.numCue = 2;
definput.keyvals.ignore = {'ITSD'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

Ncues = 8;

% Normalized externalization rating scale
Erange = 1;
Eoffset = 0;

% For cue combination analysis
cueSelLbl = {{'MSS','ISS','MSSD','ISSD','ITIT','IC','MI'};...
             {'MSS','ISS','ISSD','ITIT','IC','MI'};...
             {'MSS','ISS','ISSD','ITIT','IC'};...
             {'MSS','ISS','ISSD','ITIT'};...
             {'MSS','ISS','ISSD'};...
             {'MSS','ISS'};...
             {'MSS'};...
            };
cueSelLblTxt = {'MSS, ISS, MSSD, ISSD, ITIT, IC, MI';...
                'MSS, ISS, ISSD, ITIT, IC, MI';...
                'MSS, ISS, ISSD, ITIT, IC';...
                'MSS, ISS, ISSD, ITIT';...
                'MSS, ISS, ISSD';...
                'MSS, ISS';...
                'MSS'};
              
strategies = {'Oracle','WSM','LTA','MTA','WTA'};

% Symbol order for plotting
symb = {'-o','-s','-h','-d','-<','->','-^','-x','-p',':*','--+',':v'};
colors = flipud([228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;247,129,191;166,86,40;0,0,0]/256);

%% Hartmann & Wittenberg (1996)
if flags.do_hartmann1996
  
  cond = {'ILD','ISLD'};
  nprime = [0,1,8,14,19,22,25,38];
  f0 = 125; %Hz
  
  fncache = flags.type; % hartmann1996
  if not(flags.do_positive)
    fncache = [fncache,'_',flags.gradients,'Grad'];
  end
  [Pext,cues,cueLbl] = amt_cache('get',fncache,flags.cachemode);
  if isempty(Pext)
    azi = -37;
    data = data_baumgartner2017;
    if flags.do_quickCheck
      data = data(1:5);
    end
    Pext = repmat({nan(length(data),length(nprime))},2,2);
    cues = repmat({nan(length(data),length(nprime),Ncues)},2,1);
    for isub = 1:length(data)
      Obj = data(isub).Obj;
      template = sig_hartmann1996(0,'Obj',Obj,'dur',0.1);
      for ee = 1:length(cond)
        for nn = 1:length(nprime)
          target = sig_hartmann1996(nprime(nn),'Obj',Obj,'dur',0.1,cond{ee});
          [~,cues{ee}(isub,nn,:),cueLbl] = baumgartner2021(target,template,'flow',100,'fhigh',5000,'lat',azi,flags.gradients);
        end
      end
      amt_disp([num2str(isub),' of ',num2str(length(data)),' subjects completed.']);
    end
    if not(flags.do_quickCheck)
      amt_cache('set',fncache,Pext,cues,cueLbl);
    end
  end
  
  %% Cue weighting optimization procedure
  
  Pext = cell(Ncues-length(kv.ignore),length(cond));
  PextNumCue = cell(length(strategies),length(cond));
  %figure
  hax{1} = local_tight_subplot(1,2,0,[.15,.1],[.1,.05]);
  %figure;
  %hax{2} = local_tight_subplot(1,2,0,[.15,.1],[.1,.05]);
  for ee = 1:length(cond)
    act = data_hartmann1996(cond{ee});
    actual = act.avg.Escore(:)/3; % normalized response scale
    
    % align actual and simulated data format
    tmp = permute(cues{ee},[2,1,3]);
    tmp = interp1(nprime,tmp,act.avg.nprime);
    predCue = squeeze(num2cell(tmp,1:2));
    
    fncache2 = [fncache,'_exp',num2str(ee)];
    amt_disp([flags.type,', Exp. ',num2str(ee)],'verbose');
    amt_cache('set',[fncache2,'_fitting'],predCue,actual);
    
    % cue analysis: first optimal then fixed sensitivities
    [Pext(:,ee),PextNumCue(:,ee),cueTab{ee},cueSelTab{ee}] = local_cueAnalysis(...
      cueLbl,cueSelLbl,cueSelLblTxt,predCue,actual,Erange,Eoffset,fncache2,kv);
  
  %% Plot
    dataLbl = [{'Actual'};cueTab{ee}.Properties.RowNames];
    Nsubj = size(Pext{1},2);
    act = data_hartmann1996(cond{ee});
    x = 1e-3*f0*act.avg.nprime;
    if flags.do_supplements
      yall = {Pext,PextNumCue};
    else
      yall = {Pext};
    end
    for iy = 1:length(yall)
      axes(hax{iy}(ee)); clear h
      y = yall{iy};
      for cc = 1:size(y,1)
          dx = 0.03*(cc - size(y,1)/2);
          h(cc+1) = errorbar(x+dx,mean(y{cc,ee},2),std(y{cc,ee},0,2)/sqrt(Nsubj),...
            symb{cc+1},'Color',colors(cc+1,:));
          hold on
      end
      h(1) = plot(x,act.avg.Escore/3,symb{1},'Color',colors(1,:));
      set(h(2:end),'MarkerFaceColor','w')
      set(h(1),'MarkerFaceColor','k')
      if ee == 1
        xlabel('Highest frequency altered in kHz') % n^{\prime}
        ylabel('Externalization')
      else
        xlabel('Highest frequency altered in kHz')
        set(gca,'YTickLabel',{})
      end
      axis([x(1)-.5,x(end)+.5,Eoffset-0.1*Erange,Eoffset+1.1*Erange])
      box off
      if ee == 1
        if iy == 1
          leg= [{'Actual'};cueTab{ee}.Properties.RowNames];
          legend(h,leg(1:end-1),'Location','southwest','box','off');
        else
          leg= [{'Actual'};strategies(:)];
          legend(h,leg,'Location','southwest','box','off');
        end
      end
    end
  end

end

%% Hassager et al. (2016)
if flags.do_hassager2016
  azi = [0,50];
  flp = 6000; % lowpass filtered noise stimuli
  
  Pext_A = data_hassager2016;
  B = Pext_A.B;
  
  fncache = flags.type;
  if not(flags.do_positive)
    fncache = [fncache,'_',flags.gradients,'Grad'];
  end
  [Pext,cues,cueLbl] = amt_cache('get',fncache,flags.cachemode);
  if isempty(cues)
    
    data = data_baumgartner2017;
    if flags.do_quickCheck
      data = data(1:5);
    end

    cues = nan(length(B),length(azi),length(data),Ncues);
    for isubj = 1:length(data)
      Obj = data(isubj).Obj;
      for iazi = 1:length(azi)
        idazi = Obj.SourcePosition(:,1) == azi(iazi) & Obj.SourcePosition(:,2) == 0;
        template = squeeze(shiftdim(Obj.Data.IR(idazi,:,:),2));
        for iB = 1:length(B)
          amt_disp(num2str(iB),'volatile');
          if isnan(B(iB))
            target = template;
          else
            Obj_tar = sig_hassager2016(Obj,B(iB));
            target = squeeze(shiftdim(Obj_tar.Data.IR(idazi,:,:),2));
          end
          [~,cues(iB,iazi,isubj,:),cueLbl] = baumgartner2021(target,template,'flow',100,'fhigh',flp,'lat',azi(iazi),flags.gradients);
        end
      end
      amt_disp([num2str(isubj),' of ',num2str(length(data)),' subjects completed.']);
    end
     
    if not(flags.do_quickCheck)
      amt_cache('set',fncache,Pext,cues,cueLbl);
    end
  end
  
  %% Cue weighting optimization procedure
  dimSubj = 3;
  Nsubj = size(cues,dimSubj);

  predCue = squeeze(num2cell(reshape(cues,[],Nsubj,Ncues),1:2));
  % optimize sensitivities
  actual = Pext_A.rating(:);%repmat(Pext_A.rating(:),Nsubj,1); 
  actual = (actual-1)/4; % normalized response scale
  amt_disp(flags.type,'verbose');
  
  amt_cache('set',[fncache,'_fitting'],predCue,actual);
  
  % cue analysis
  [Pext,PextNumCue,cueTab,cueSelTab] = local_cueAnalysis(...
    cueLbl,cueSelLbl,cueSelLblTxt,predCue,actual,Erange,Eoffset,fncache,kv);

  
  %% Plot
  BplotTicks = logspace(log10(1),log10(64),3);
  BplotTicks = round(BplotTicks*100)/100;
  BplotStr = ['Ref.';[repmat(' ',3,1),num2str(BplotTicks(:)),repmat(' ',3,1)]];
  BplotTicks = [0.25/2,BplotTicks];
  B(1) = B(2)/2;
  if flags.do_supplements
    yall = {Pext,PextNumCue};
  else
    yall = {Pext};
  end
  for iy = 1:length(yall)
    y = yall{iy};
    figure 
    hax = local_tight_subplot(1,2,0,[.15,.1],[.1,.05]);
    for iazi = 1:length(azi)
      axes(hax(iazi)); clear h
      h(1) = plot(B,(Pext_A.rating(:,iazi)-1)/4,symb{1},'Color',colors(1,:));
      hold on
      for m = 1:length(y)
        dx = 1+0.01*(m - (Ncues+1)/2);
        pred = reshape(y{m},length(B),[],Nsubj);
        h(m+1) = errorbar(dx*B,mean(pred(:,iazi,:),dimSubj),std(pred(:,iazi,:),0,dimSubj)/sqrt(Nsubj),...
          symb{m+1},'Color',colors(m+1,:));
      end
      set(h,'MarkerFaceColor','w')
      set(h(1),'MarkerFaceColor','k')
      set(gca,'XTick',BplotTicks,'XTickLabel',BplotStr,'XScale','log')
      box off
      axis([BplotTicks(1)/1.5,BplotTicks(end)*1.5,Eoffset-0.1*Erange,Eoffset+1.1*Erange])
      xlabel('Bandwidth of smoothing filters in ERB')
      if iazi==1
        ylabel('Externalization')
      else
        set(gca,'YTickLabel',{})
      end
      text(0.125,Eoffset,[num2str(azi(iazi)),'\circ'])
    end
    if iy == 1
      leg=[{'Actual'};cueTab.Properties.RowNames];
      leg = legend(h,leg(1:end-1),'Location','southwest');
    else
      leg = legend(h,[{'Actual'};strategies(:)],'Location','southwest');
    end
    set(leg,'Box','off')
  end

  
end

%% Baumgartner et al. (2017) ----------------------------------------------
if flags.do_baumgartner2017
  
  % Behavioral results
  fncache = [flags.type,'_actual']; %'baumgartner2017_E_BTL'
  data = amt_cache('get',fncache,flags.cachemode);
  if isempty(data)
    exp2 = data_baumgartner2017looming('exp2');
    Nsubj = length(exp2.meta.subject);
    iISI = exp2.meta.ISI == 0.1;
    iresp = strcmp(exp2.meta.response,'approaching');
    data.C = 0:0.5:1;
    data.subj = [];%exp2.meta.subject;
    data.Escore = [];%nan(Nsubj,length(A));
    data.azi = [exp2.rawData.azimuth];
    M = zeros(3,3);
    iM = [7,4,8,3,2,6]; % indices to map judgements to preference matrix
    iss = 0;
    for ss = 1:Nsubj
      M(iM) = exp2.data(ss,1:6,iresp,iISI)+eps;
      Ns = 200; % max # hits  
      E = nansum(M,2)/Ns;
      E = E/E(3);
      if E(1)<=1 % if C = 0 is more externalized than the reference (C = 1) either the subject misunderstood the task or something was wrong with the reference 
        iss = iss+1;
        data.subj{iss} = exp2.meta.subject{ss};
        data.Escore(iss,:) = E;
      end
    end 
    amt_cache('set',fncache,data)
  end
  Nsubj = length(data.subj);
  
  hrtf = data_baumgartner2017looming('exp2','hrtf');
  
  % Modeling
  fncache = flags.type;
  if not(flags.do_positive)
    fncache = [fncache,'_',flags.gradients,'Grad'];
  end
  [Pext,cues,cueLbl] = amt_cache('get',fncache,flags.cachemode);
  if isempty(cues)
    
    Pext = nan(length(data.C),Nsubj);
    cues = nan(length(data.C),Nsubj,Ncues);
    for isubj = 1:Nsubj
      template = hrtf(ismember({hrtf.id},data.subj(isubj))).Obj;
      for iC = 1:length(data.C)
        target = sig_baumgartner2017looming(template,data.C(iC));
        [Pext(iC,isubj),cues(iC,isubj,:),cueLbl] = baumgartner2021(target,template,...
                          'lat',data.azi(isubj),'flow',1000,'fhigh',18e3,flags.gradients);
        
      end
      amt_disp([num2str(isubj),' of ',num2str(Nsubj),' subjects completed.']);
    end
    amt_cache('set',fncache,Pext,cues,cueLbl);
  end
 
  %% Cue weighting optimization procedure
  dimSubj = 2;
  dimCues = 3;
  
  Nsubj = size(cues,dimSubj);
  actual = permute(data.Escore,[2,1]);
  predCue = squeeze(num2cell(cues,1:dimCues-1));
  amt_disp(flags.type,'verbose');
  amt_cache('set',[fncache,'_fitting'],predCue,actual);
  [Pext,PextNumCue,cueTab,cueSelTab] = local_cueAnalysis(...
    cueLbl,cueSelLbl,cueSelLblTxt,predCue,actual,Erange,Eoffset,fncache,kv);
  
  %% Figure
  flags.do_plot_individual = false;
  if flags.do_supplements
    yall = {[{data.Escore'};Pext],[{data.Escore'};PextNumCue]};
  else
    yall = {[{data.Escore'};Pext]};
  end
  for iy = 1:length(yall)
    y = yall{iy};
    figure; clear h;
    if not(flags.do_plot_individual)
      for iE = 1:length(y)
        dx = .01*(iE - length(y)/2);
        h(iE) = errorbar(data.C+dx,mean(y{iE},2),std(y{iE},0,2)/sqrt(Nsubj),...
          symb{iE},'Color',colors(iE,:));
        set(h,'MarkerFaceColor','w')
        hold on
      end
      set(h(1),'MarkerFaceColor','k')
      set(gca,'XTick',data.C)
      axis([-0.2,1.2,Eoffset-0.1*Erange,Eoffset+1.1*Erange])
      box off
      ylabel('Externalization')
      xlabel('Spectral contrast, C')
    else % flags.do_plot_individual
      for isubj = 1:Nsubj
        subplot(3,ceil(Nsubj/3),isubj); clear h;
        for iE = 1:length(y)
          dx = .005*(iE - length(y)/2);
          h(iE) = plot(data.C+dx,y{iE}(:,isubj),...
            symb{iE},'Color',colors(iE,:));
          set(h,'MarkerFaceColor','w')
          hold on
        end
        set(h(1),'MarkerFaceColor','k')
        set(gca,'XTick',data.C)
        axis([-0.2,1.2,Eoffset-0.1*Erange,Eoffset+1.1*Erange])
        box off
        ylabel('Externalization')
        xlabel('Spectral contrast, C')
      end
    end
    if iy == 1
      leg=[{'Actual'};cueTab.Properties.RowNames];
      legend(h,leg(1:end-1),'Location','EastOutside','box','off')
    else
      leg=[{'Actual'};strategies(:)];
      legend(h,leg,'Location','EastOutside','box','off')
    end
  end
  
end

%% Boyd et al. (2012) -----------------------------------------------------
if flags.do_boyd2012
  
  flags.do_BRIR = true; % directly based on BRIRs
  flags.do_noise = false; % or based on BRIR-convoluted noise
  % otherwise: original speech samples
  
  flp = [nan,6500]; % Low-pass cut-off frequency
  azi = -30;
  
  subjects = data_boyd2012;
  if flags.do_noise || flags.do_BRIR
    subjects = subjects([3,6,7]); % only the ones with BRIRs
  else
    subjects = subjects([1,3,4,6,7]); % only the ones with stimuli
  end
  mix = subjects(1).Resp.mix/100;
  idn = not(isnan(mix));
  mix = mix(idn); % omit reference
  
  % determine cache name
  fncache = flags.type;%['boyd2012'];
  if not(flags.do_positive)
    fncache = [fncache,'_',flags.gradients,'Grad'];
  end
  if flags.do_noise
    fncache = [fncache,'_noise'];
  elseif flags.do_BRIR
    % standard
  else
    fncache = [fncache,'_speech'];
  end
    
  [E,cues,cueLbl] = amt_cache('get',fncache,flags.cachemode);
  if isempty(cues)
    
    if flags.do_noise || flags.do_BRIR
      sig = noise(50e3,1,'pink');
      fs = subjects(1).BRIR.fs;
      [b,a]=butter(10,2*flp(2)/fs,'low');
    end

    E.all = repmat({nan([length(mix),2,2,length(subjects)])},1,1);
    cues = nan(length(mix),2,2,length(subjects),Ncues);
    for isub = 1:length(subjects)
      E.all{1}(:,:,:,isub) = cat(3,...
        [subjects(isub).Resp.ITE_BB_1T(idn),subjects(isub).Resp.BTE_BB_1T(idn)],...
        [subjects(isub).Resp.ITE_LP_1T(idn),subjects(isub).Resp.BTE_LP_1T(idn)]);
      if flags.do_noise
        stim{1} = lconv(sig,subjects(isub).BRIR.ITE(:,:,1));
        stim{2} = lconv(sig,subjects(isub).BRIR.BTE(:,:,1));
        stim{3} = lconv(sig,subjects(isub).BRIR.noHead(:,:,1));
        template = stim{1};
        targetSet = cell(length(mix),2,2);
      elseif flags.do_BRIR
        template = subjects(isub).BRIR.ITE(:,:,1);
        stim{1} = subjects(isub).BRIR.ITE(:,:,1);
        stim{2} = subjects(isub).BRIR.BTE(:,:,1);
        stim{3} = subjects(isub).BRIR.noHead(:,:,1);
        targetSet = cell(length(mix),2,2);
      else % speech samples
        template = subjects(isub).Reference_1T;
        targetSet = cat(3,...
          [subjects(isub).Target.ITE_BB_1T(:),subjects(isub).Target.BTE_BB_1T(:)],...
          [subjects(isub).Target.ITE_LP_1T(:),subjects(isub).Target.BTE_LP_1T(:)]);
      end

      %% Model simulations
      flow = 100;
      fhigh = [16e3,16e3];  % set fhigh(2) to flp(2) for manual bandwidth adjustment
            
      for c = 1:2
        for lp = 1:2
          for m = 1:length(mix)
            % mixing
            if flags.do_noise || flags.do_BRIR
              targetSet{m,c,lp} = mix(m)*stim{c} + (1-mix(m))*stim{3};
              % low-pass filtering
              if lp == 2
                targetSet{m,c,lp} = filter(b,a,targetSet{m,c,lp});
              end
            end
            target = targetSet{m,c,lp};
            [~,cues(m,c,lp,isub,:),cueLbl] = baumgartner2021(target,template...
              ,'argimport',flags,kv,'flow',flow,'fhigh',fhigh(lp),'lat',azi,flags.gradients);%,'reflectionOnsetTime',5e-3);
             
          end
        end
      end
      amt_disp([num2str(isub),' of ',num2str(length(subjects)),' subjects completed.']);

    end
    
    if not(flags.do_quickCheck)
      amt_cache('set',fncache,E,cues,cueLbl);
    end
  end

  %% Cue weighting optimization procedure
  dimSubj = 4;
  
  Nsubj = size(cues,dimSubj);
  if flags.do_BRIR % comparison on individual basis
    actual = E.all{1};
    actual = reshape(actual,[],Nsubj);
  else % average comparison
    actual = mean(E.all{1},4);
    actual = actual(:);
  end
  
  actual = actual/100; % normalized response range
  E.all{1} = E.all{1}/100;
  
  predCue = squeeze(num2cell(reshape(cues,[],Nsubj,Ncues),1:2));
  amt_disp(flags.type,'verbose');
  amt_cache('set',[fncache,'_fitting'],predCue,actual);
  [Pext,PextNumCue,cueTab,cueSelTab] = local_cueAnalysis(...
    cueLbl,cueSelLbl,cueSelLblTxt,predCue,actual,Erange,Eoffset,fncache,kv);
  
  dims = [length(mix),size(cues,2),size(cues,3),size(cues,4)];
  for ii = 1:length(Pext)
    E.all{ii+1} = reshape(Pext{ii},dims);
  end
  EnumCue = E.all(1);
  for ii = 1:length(PextNumCue)
    EnumCue{ii+1} = reshape(PextNumCue{ii},dims);
  end
  
  %% Plot
  flags.do_plot_individual = false;
  
  condLbl = {'ITE | BB','BTE | BB','ITE | LP','BTE | LP'};
  mix = mix*100;
  
  if not(flags.do_plot_individual)
    if flags.do_supplements
      yall = {E.all,EnumCue};
    else
      yall = {E.all};
    end
    for iy = 1:length(yall)
      y = yall{iy};
      figure
      hax = local_tight_subplot(2,2,0,[.15,.1],[.1,.05]);
      for cc = 1:length(condLbl)
        axes(hax(cc))
        for ee = 1:length(y)
          dx = 1.2*(ee - length(y)/2);
          m = mean(y{ee},4);
          s = std(y{ee},0,4);
          h = errorbar(mix+dx,m(:,cc),s(:,cc),...
            symb{ee},'Color',colors(ee,:));
          if ee == 1
            set(h,'MarkerFaceColor','k')
          else
            set(h,'MarkerFaceColor','w')
          end
          hold on
        end
        set(gca,'XDir','reverse')

        if cc == 3 || cc == 4
          xlabel('IID presence in %')
        else
          set(gca,'XTickLabel',[])
        end

        if cc == 1 || cc == 3
          ylabel('Externalization')
        else
          set(gca,'YTickLabel',[])
        end
        text(110,0,condLbl{cc},'HorizontalAlignment','left')
        axis([-20,120,Eoffset-0.1*Erange,Eoffset+1.1*Erange])
        box off
        if cc == 4
        end
      end
    end
  else % flags.do_plot_individual
    
    for isub = 1:size(E.all{1},4)
      figure
      hax = local_tight_subplot(1,4,0,[.15,.1],[.1,.05]);
      for cc = 1:length(condLbl)
        axes(hax(cc))
        for ee = 1:Nee
          Eisub = E.all{ee}(:,:,:,isub);
          h = plot(mix,Eisub(2:end,cc),symb{ee},'Color',colors(ee,:));
          set(h,'MarkerFaceColor','w')
          hold on
        end
        set(gca,'XDir','reverse')
        xlabel('IID presence in %')
        if cc == 1
          ylabel('Externalization')
        else
          set(gca,'YTickLabel',[])
        end
        title(condLbl{cc})
        axis([-20,120,Eoffset-0.1*Erange,Eoffset+1.1*Erange])
        if cc == 4
          leg = legend(dataLbl);
          set(leg,'Box','off','Location','north')
        end
      end
    end
    
  end


Pext = E;

end

%% Output for individual experiments
if exist('cueTab','var')
  data = {cueTab,cueSelTab};
end

%% Figure 2
if flags.do_all || flags.do_fig2
  exp = {'hartmann1996','hassager2016','baumgartner2017','boyd2012'};
  for ii = 1:length(exp)
    exp_baumgartner2021(exp{ii},flags.cachemode);
  end
  %close figure 2
end%-do_all

%% Figure 4
if flags.do_fig4
  exp = {'hartmann1996_exp1','hartmann1996_exp2','hassager2016','baumgartner2017','boyd2012'};
  Nexp = length(exp);
  [~,~,cueLbl] = amt_cache('get','hassager2016');
  predCue = [];
  actual = [];
  expnum = []; % experiment number
  for ii = 1:Nexp
    if flags.do_positive
      fn = [exp{ii},'_fitting'];
    else
      if ii <=2
        fn = [exp{ii}(1:12),'_',flags.gradients,'Grad',exp{ii}(13:end),'_fitting'];
      else
        fn = [exp{ii},'_',flags.gradients,'Grad','_fitting'];
      end
    end
    [p,a] = amt_cache('get',fn);
    if isempty(p)
      amt_disp('Execution stopped because cached data not found. Please run calculations for fig2 first.');
      return
    end
    p = mean(cat(3,p{:}),2);
    predCue = cat(1,predCue,p);
    actual = cat(1,actual,mean(a,2));
    expnum = cat(1,expnum,ii*ones(size(a,1),1));
  end
  predCue = squeeze(num2cell(predCue,1:2));
  
  amt_disp('Evaluating performance of strategies including all cues...')
  Pext = cell(length(cueSelLbl),1);
  Nstrat = length(strategies)-1; % # strategies w/o Oracle
  Wsel = nan(length(cueSelLbl),length(cueSelLbl{1}),Nstrat);
  Ssel = Wsel;
  RMSE = nan(length(cueSelLbl),Nstrat);
  RMSE_CI = nan(length(cueSelLbl),Nstrat,2);
  BIC = RMSE;
  RC = RMSE;
  RC_CI = RMSE_CI;
  for ics = 1:length(cueSelLbl)
    cueSelId = ismember(cueLbl(:,1),cueSelLbl{ics});
    [Pext{ics},Stmp,Wtmp,MSEtmp,BICtmp,RCtmp] = local_mapAndDecide(...
      predCue(cueSelId),actual,Erange,Eoffset,[],[],expnum);
    Wsel(ics , ismember(cueSelLbl{1},cueSelLbl{ics}),:) = Wtmp(:,2:end);
    Ssel(ics , ismember(cueSelLbl{1},cueSelLbl{ics}),:) = Stmp(:,2:end);
    RMSE(ics,:) = sqrt(MSEtmp.m(end-Nstrat+1:end));
    RMSE_CI(ics,:,:) = sqrt(MSEtmp.ci(end-Nstrat+1:end,:));
    BIC(ics,:) = BICtmp(end-Nstrat+1:end);
    RC(ics,:) = RCtmp.m(end-Nstrat+1:end);
    RC_CI(ics,:,:) = RCtmp.ci(end-Nstrat+1:end,:);
  end
  % Performance stats for all combinations
  cueSelTab{1} = table(RMSE(:,1),RMSE(:,2),RMSE(:,3),RMSE(:,4),...
            RC(:,1),RC(:,2),RC(:,3),RC(:,4),...%BIC(:,1),BIC(:,2),BIC(:,3),BIC(:,4),...
      'RowNames',cueSelLblTxt,...
      'VariableNames',...
      {'RMSE_WSM','RMSE_LTA','RMSE_MTA','RMSE_WTA',...
      'RC_WSM','RC_LTA','RC_MTA','RC_WTA'...
%        'BIC_WSM','BIC_LTA','BIC_MTA','BIC_WTA',...
       });
  amt_disp('Prediciton errors for tested decisions:','verbose');
  amt_disp(cueSelTab{1},'verbose');
  % Weights or percentages for all combinations and strategies
  for ss = 1:Nstrat
    amt_disp([strategies{1+ss} ' weights/percentages'],'verbose');
    cueSelTab{1+ss} = array2table([squeeze(Wsel(:,:,ss)),RMSE(:,ss),RC(:,ss)],...
      'VariableNames',[cueSelLbl{1},{'RMSE','RC'}],...
      'RowNames',cueSelLblTxt);
     amt_disp(cueSelTab{1+ss},'verbose');
     amt_disp([strategies{1+ss} ' mapping sensitivity'],'verbose');
     SselTab{ss} = array2table(squeeze(Ssel(:,:,ss)),...
      'VariableNames',cueSelLbl{1},...
      'RowNames',cueSelLblTxt);
     amt_disp(SselTab{ss},'verbose');
  end

  %% Plot performance (RMSE) as a function of # considered cues
  y = {RMSE,RC};
  ci = {RMSE_CI,RC_CI};
  ylab = {'RMS prediction error','Rank correlation'};
  if ~flags.do_supplements
    y = y(1);
  end
  
  colors = [0.8500    0.3250    0.0980;... % red
            0.9290    0.6940    0.1250;... % yellow
            0.4940    0.1840    0.5560;... % purple
            0.4660    0.6740    0.1880]; % green
  for iy = 1:length(y)
    
    errlow = y{iy}-ci{iy}(:,:,1);
    errhigh = ci{iy}(:,:,2)-y{iy};
    
    if flags.do_supplements
      figure;
      x = length(cueSelLbl{1}):-1:1;
      b = bar(x,y{iy});
      for i = 1:length(b)
        b(i).FaceColor = colors(i,:);
      end
      hold on
      ngroups = size(y{iy}, 1);
      nbars = size(y{iy}, 2);
      % Calculating the width for each bar group
      groupwidth = min(0.8, nbars/(nbars + 1.5));
      for i = 1:nbars
          x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
          errorbar(fliplr(x), y{iy}(:,i), errlow(:,i), errhigh(:,i), 'k.');
      end
      hold off
      leg = legend(strategies(2:end),'Location','northwest');
      leg.Box = 'off';
      xlabel('Number of cues')
      ylabel(ylab{iy})
    end

    % Performance (RMSE and RC) only for all 7 cues
    figure
    for ii = 1:4
      b(ii) = bar(ii,y{iy}(1,ii));
      b(ii).FaceColor = colors(ii,:);
      hold on
      errorbar(ii, y{iy}(1,ii), errlow(1,ii), errhigh(1,ii), 'k.');
    end
    hold off
    set(gca,'XLim',[.4,4.6],'XTick',1:4,'XTickLabel',strategies(2:end))
    set(gca,'XTickLabelRotation',45)
    xlabel('Strategy')
    ylabel(ylab{iy})
    box off
  end

  %% Cue weight or selection stats
  figure
  x = 1:length(cueSelLbl{1});
  Wsel(:,:,1) = 100*Wsel(:,:,1); % also weights in percent
  b = barh(fliplr(x),squeeze(Wsel(1,:,:)));
  for i = 1:length(b)
    b(i).FaceColor = colors(i,:);
  end
  set(gca,'YTick',x,'YTickLabel',fliplr(cueSelLbl{1}))
  xlabel('Weight or selection rate in %')
  box off

  %% Check all pairwise combinations
  amt_disp('Checking pairwise combinations...')
  Ncues = length(cueSelLbl{1});
  combs = nchoosek(1:Ncues,2);
  Ncombs = size(combs,1);
  Pext = cell(Ncombs,1);
  Wsel = nan(Ncombs,Ncues);
  S = nan(Ncombs,Ncues);
  Nstrat = length(strategies)-1; % no Oracle
  RMSE = nan(Ncombs,Nstrat);
  RMSE_CI = nan(Ncombs,Nstrat,2);
  BIC = RMSE;
  RC = RMSE;
  RC_CI = RMSE_CI;
  for ics = 1:Ncombs
    cueSelId = false(length(predCue),1);
    cueSelId(combs(ics,:)) = true;
    [Pext{ics},Stmp,Wtmp,MSEtmp,BICtmp,RCtmp] = local_mapAndDecide(...
      predCue(cueSelId),actual,Erange,Eoffset,[],[],expnum);
    Wsel(ics , cueSelId(1:Ncues)) = Wtmp(:,2);
    S(ics , cueSelId(1:Ncues)) = Stmp(:,2);
    RMSE(ics,:) = sqrt(MSEtmp.m(end-Nstrat+1:end));
    RMSE_CI(ics,:,:) = sqrt(MSEtmp.ci(end-Nstrat+1:end,:));
    BIC(ics,:) = BICtmp(end-Nstrat+1:end);
    RC(ics,:) = RCtmp.m(end-Nstrat+1:end);
    RC_CI(ics,:,:) = RCtmp.ci(end-Nstrat+1:end,:);
  end
  Stab = array2table(S,'VariableNames',cueSelLbl{1});
  data = array2table(cat(2,Wsel,RMSE(:,1),squeeze(RMSE_CI(:,1,:)),RC(:,1),squeeze(RC_CI(:,1,:))),...
          'VariableNames',cat(2,cueSelLbl{1},{'WSM_RMSE','RMSE_CIlow','RMSE_CIhigh','RC','RC_CIlow','RC_CIhigh'}));
  [data,isort]=sortrows(data,8);
  Stab = Stab(isort,:);
  amt_disp('Contributions and Prediction errors:','verbose');
  amt_disp(data,'verbose');
  amt_disp('Sensitivities:','verbose');
  amt_disp(Stab,'verbose');
  
  %% Plot errors and weights for all paired WSM combinations
  y = {RMSE(isort,1),RC(isort,1)};
  ci = {RMSE_CI(isort,1,:),RC_CI(isort,1,:)};
  ylab = {'RMS prediction error','Rank correlation'};
  ylim = {[.01,.5],[.01,1]};
  if ~flags.do_supplements
    y = y(1);
  end
  for iy = 1:length(y)
    figure
    ha = local_tight_subplot(2,1,0,[.12 .03],[.15 .22]);
    axes(ha(1))
    x = 1:length(isort);
    xlim = [x(1)-.5,x(end)+.5];
    try
      h = errorbar(x,y{iy},y{iy}-ci{iy}(:,1,1),ci{iy}(:,1,2)-y{iy},'k.','CapSize',0);
    catch
      h = errorbar(x,y{iy},y{iy}-ci{iy}(:,1,1),ci{iy}(:,1,2)-y{iy},'k.');
    end
    hold on
    plot(xlim,ones(2,1)*ci{iy}(iy,1,2),'k--')
    ylabel(ylab{iy})
    axis([xlim,ylim{iy}])
    box off
    set(gca,'XTickLabel',{},'XTick',x)
    axes(ha(2))
    imagesc(Wsel(isort,:)')
    colormap(flipud(hot(100)));
    c = colorbar('eastoutside');
    c.Position(1) = c.Position(1)+.1;
    c.Label.String = 'Weight';
    c.Box = 'off';
    XTickLabel = num2cell(x);
    XTickLabel(2:2:end) = repmat({[]},floor(length(x)/2),1);
    set(gca,'YTick',1:length(cueSelLbl{1}),'YTickLabel',cueSelLbl{1},'XTick',x,'XTickLabel',XTickLabel)
    xlabel('Error rank')
    box off
  end
end%-do_fig4  

%% Figure 3
if flags.do_tabSingleCue || flags.do_fig3

  if flags.do_redo
    warning('Redo is ignored for fig3 because results are taken from calculations for fig2. To re-run all calculations, redo fig2 instead.')
  end
  
  exp = {'hartmann1996_exp1','hartmann1996_exp2','hassager2016','baumgartner2017','boyd2012'};
  Nexp = length(exp);
  for ii = 1:Nexp
    
    if flags.do_positive
      fn = [exp{ii},'_tab2'];
    else
      if ii <=2 % two experiments from hartmann1996
        fn = [exp{ii}(1:12),'_',flags.gradients,'Grad',exp{ii}(13:end),'_tab2'];
      else
        fn = [exp{ii},'_',flags.gradients,'Grad','_tab2'];
      end
    end
    
    tab = amt_cache('get',fn,'localonly');
    
    tabT = array2table(transpose(table2array(tab(:,1:3))));
    for rr = 1:3
      tabT.Properties.RowNames{rr} = ['Exp',num2str(ii),'_',tab.Properties.VariableNames{rr}];
    end
    tab.Properties.RowNames = strrep(tab.Properties.RowNames,'Oracle','Oracle');
    tabT.Properties.VariableNames = tab.Properties.RowNames;
    
    if ii == 1
      tab2 = tabT;
    else
      tab2 = [tab2;tabT];
    end
    
  end
  Error = tab2{2:3:end,:};
  Sensitivity = tab2{1:3:end,:};
  Weight = tab2{3:3:end,:};
  
  MeanError = mean(Error);
  WavgS_cue = sum(Sensitivity.*Error.^-2./repmat(sum(Error.^-2),5,1));
  MeanWeight = mean(Weight);
  
  sumTab = array2table([WavgS_cue;MeanError;MeanWeight]);
  sumTab.Properties.VariableNames = tab2.Properties.VariableNames;
  sumTab.Properties.RowNames = {'W. Avg. S_cue','Mean Error','Mean Weight'};
  tab2 = [tab2;sumTab];

  if ~flags.do_fig2, amt_disp(tab2,'verbose'); end
  % Output
  data = tab2;
  
  if flags.do_fig3
    
    figure
    x = 1:5;
    imagesc(Error(:,1:7)')
    colormap(flipud(copper(100)));
    c = colorbar('eastoutside');
    c.Label.String = 'Error';
    c.Box = 'off';
    c.Limits = round(c.Limits*10)/10;
    c.Ticks = c.Limits(1):.1:c.Limits(2);
    XTickLabel = {'I','II','III','IV','V'};
    set(gca,'YTick',1:length(cueSelLbl{1}),'YTickLabel',cueSelLbl{1},'XTick',x,'XTickLabel',XTickLabel)
    xlabel('Experiment')
    box off
    
    figure
    imagesc(log10(1./Sensitivity(:,1:7)'))
    colormap(flipud(summer(100)));
    c = colorbar('eastoutside');
    c.Label.String = 'Log Sensitivity';
    c.Box = 'off';
    XTickLabel = {'I','II','III','IV','V'};
    set(gca,'YTick',1:length(cueSelLbl{1}),'YTickLabel',cueSelLbl{1},'XTick',x,'XTickLabel',XTickLabel)
    xlabel('Experiment')
    box off
  end
end%-do_fig3

%% Results for iteratively restricted set of cues 
if flags.do_tabStrats 
    
  experiments = {'Exp.I','Exp.II','Exp.III','Exp.IV','Exp.V'};
  exp = {'hartmann1996_exp1','hartmann1996_exp2','hassager2016','baumgartner2017','boyd2012'};
  Nexp = length(exp);
  MSE = [];
  BIC = [];
  for ii = 1:Nexp
    
    if flags.do_positive
      fn = [exp{ii},'_tab3'];
    else
      if ii <=2
        fn = [exp{ii}(1:12),'_',flags.gradients,'Grad',exp{ii}(13:end),'_tab3'];
      else
        fn = [exp{ii},'_',flags.gradients,'Grad','_tab3'];
      end
    end
    tab = amt_cache('get',fn,flags.cachemode);
    
    if isempty(tab)
      if ii == 1
        exp{ii} = 'hartmann1996';
      end
      if ii == 2
        tab = amt_cache('get',[exp{ii},'_tab3']); % has been (re)computed together with exp. 1
      else
        data = exp_baumgartner2021(exp{ii},flags.cachemode);
        if ii < 3
          tab = data{2}{2};
        else
          tab = data{2};
        end
      end
    end
    
    if ii == 1
      meanTab = tab;
    else
      meanTab{:,:} = meanTab{:,:} + tab{:,:};
    end
    
    MSE(:,:,ii) = tab{:,7+(1:length(strategies))};
    BIC(:,:,ii) = tab{:,8+length(strategies):end};
    
  end
  meanTab{:,:} = 1/Nexp * meanTab{:,:};
  
  if flags.do_supplements
    amt_disp(meanTab,'verbose');
    
    % plot MSE and BIC for every strategy and experiment
    dBIC = BIC-repmat(BIC(:,1,:),1,length(strategies),1);
    y = {MSE,dBIC};
    ylbl = {'Mean squared prediction error','BIC'};
    for yy = 1:length(y)
        figure
        Ncues = size(y{yy},1):-1:1;
        b = bar(Ncues,mean(y{yy},3));
        hold on
        symbols = {'o','s','<','p','h'};
        hleg = nan(Nexp,1);
        for ii = 1:Nexp
          for jj = 1:length(strategies)
              h = plot(Ncues+ .9*(jj/(1+length(strategies))-0.5),y{yy}(:,jj,ii),['k',symbols{ii}]);
              set(h,'MarkerFaceColor',b(jj).FaceColor)
          end
          hleg(ii) = h(1);
        end
        xlabel('Number of cues')
        ylabel(ylbl{yy})
        l = legend([b';hleg],[strategies,experiments]);
        set(l,'Box','off','Location','eastoutside')
        box off
    end
  end
  
  % Output
  data = meanTab;
 
end%-do_tabStrats

%% Diotic condition
if flags.do_diotic
    subj = data_baumgartner2017;
    target = repmat(noise(1000,1),1,1,2);
    strat = {'WSM','LTA','MTA','WTA'};
    % sensitivities obtained for strategy-specific optimizations (do_fig4)
    W = [0.6     0.4]; % weights
    S{1} = [0.20029    0.51084]; % WSM
    S{2} = [1.4042    0.55273]; % LTA
    S{3} = [0.38715    0.26584]; % MTA
    S{4} = [0.212    0.081351]; % WTA
    azi = 0;
    for ss = 1:length(subj)
        idazi = subj(ss).Obj.SourcePosition(:,1) == azi & subj(ss).Obj.SourcePosition(:,2) == 0;
        template = squeeze(shiftdim(subj(ss).Obj.Data.IR(idazi,:,:),2));
        for ii = 1:length(strat)
            E(ii,ss) = baumgartner2021(target,template,'S',S{ii},'cueWeights',W,strat{ii});
        end
    end
    amt_disp('Predicted externalization scores for diotic stimuli:','verbose')
    tab = table(mean(E,2),'RowNames',strat,'VariableNames',{'Externalization'});
    amt_disp(tab,'verbose')
    % Output
    data = tab;
end%-do_diotic

end%-exp_baumgartner2021


%%%%%%%%%%%%%%%%%%%% INTERNAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%
function ha = local_tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% local_tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% ha = local_tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   going row-wise as in
%
%  Example: ha = local_tight_subplot(3,2,[.01 .03],[.1 .1],[.1 .1])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 20.6.2010   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1; 
    gap = [gap gap];
end
if numel(marg_w)==1; 
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1; 
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh; 

ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
end

function [Pext,PextNumCue,cueTab,cueSelTab] = local_cueAnalysis(cueLbl,cueSelLbl,cueSelLblTxt,predCue,actual,Erange,Eoffset,fncache,kv)
% [Pext,PextNumCue,cueTab,cueSelTab] = local_cueAnalysis(cueLbl,cueSelLbl,cueSelLblTxt,predCue,actual,Erange,Eoffset,fncache,kv)
%
% Description:
% local_cueAnalysis first runs local_mapAndDecide to obtain predicted 
% externalization (Pext) scores for all individual cues. 
% Then, it runs local_mapAndDecide again to obtain WSM 
% predictions (PextNumCue) for an iteratively restricted set of cues and
% for the different combination strategies. 

% remove to be ignored cues
idkeep = not(ismember(cueLbl(:,1),kv.ignore));
cueLbl = cueLbl(idkeep,:);
predCue = predCue(idkeep);
  
% Predictions for individual cues -> Pext
[Pext,S,W,MSE] = local_mapAndDecide(predCue,actual,Erange,Eoffset);

PextLbl = [cueLbl(:,1);{'Oracle'}];
cueTab = table([S(:,1);nan],MSE.m(1:length(S)+1),[W(:,1);nan],...
  [cueLbl(:,2);...
  {'Weighted combination'}],...
  'RowNames',PextLbl,...
  'VariableNames',{'S_cue','Error','Weight','Description'});

amt_disp(['Showing modelling details for ' fncache],'verbose');
amt_disp(cueTab,'verbose'); 

amt_cache('set',[fncache,'_tab2'],cueTab);  

% Predictions for combination strategies -> PextNumCue
try
  tab2 = exp_baumgartner2021('tabSingleCue','cached');
catch
  amt_disp('Table 2 cannot be created. Re-run individual experiments.');
end 
if exist('tab2','var') && not(isempty(tab2))
  Sfix = tab2{16,:};
  MSE = nan(length(cueSelLbl),5);
  BIC = nan(length(cueSelLbl),5);
  Wsel = nan(length(cueSelLbl),length(cueSelLbl{1}));
  try 
    tab3 = exp_baumgartner2021('tabStrats','cached'); 
  catch
    tab3 = array2table(nan(length(cueSelLbl),length(cueSelLbl{1})));  
    amt_disp('Re-run to calculate WSM');
  end
  PextNumCue = cell(length(cueSelLbl),1);
  for ics = 1:length(cueSelLbl)
    Wfix = tab3{ics,ismember(cueSelLbl{1},cueSelLbl{ics})};
    cueSelId = find(ismember(cueLbl(:,1),cueSelLbl{ics}));
    [PextNumCue{ics},~,Wtmp,MSEtmp,BICtmp] = local_mapAndDecide(...
        predCue(cueSelId),actual,Erange,Eoffset,Sfix(cueSelId),Wfix);
    Wsel(ics , ismember(cueSelLbl{1},cueSelLbl{ics})) = Wtmp(:,1);
    MSE(ics,:) = MSEtmp.m(end-4:end);
    BIC(ics,:) = BICtmp(end-4:end);
  end
  cueSelTab = [array2table(Wsel,'VariableNames',cueSelLbl{1}),...
      table(MSE(:,1),MSE(:,2),MSE(:,3),MSE(:,4),MSE(:,5),...
            BIC(:,1),BIC(:,2),BIC(:,3),BIC(:,4),BIC(:,5),...
      'RowNames',cueSelLblTxt,...
      'VariableNames',...
      {'Error_Oracle','Error_WSM','Error_LTA','Error_MTA','Error_WTA',...
       'BIC_Oracle','BIC_WSM','BIC_LTA','BIC_MTA','BIC_WTA'})];
   amt_disp(['Showing cue contributions for ' fncache],'verbose');
   amt_disp(cueSelTab,'verbose');
  amt_cache('set',[fncache,'_tab3'],cueSelTab);
else
  cueSelTab = nan;
end

% only output to be plotted data
Pext = Pext(1:length(S)); % only single-cue predictions
if exist('PextNumCue','var')
  PextNumCue = PextNumCue{end-kv.numCue+1}; % only specific # cues
  PextNumCue = PextNumCue(kv.numCue+1:end); % only combination strategies
else
  PextNumCue = Pext(kv.numCue+1:end);
  amt_disp('Warning: Used *Pext* instead of *PextNumCue*.');
end
end

function [Pext,S,W,MSE,BIC,RC] = local_mapAndDecide(predCue,actual,Erange,Eoffset,S,Wfix,expnum)
% [Pext,S,W,MSE,BIC,RC] = local_mapAndDecide(predCue,actual,Erange,Eoffset,S,Wfix)
%
% Input parameters:
% predCue ... cue-specific deviation metrics; 
%             matrix dimensions: conditions x listeners
% actual ... listeners' externalization scores; matrix dimensions must
%            match *predCue*
% Erange ... range of externalization rating scale, typically 1
% Eoffset ... offset of externalization rating scale, typically 0
% S ... sensitivity parameters, one for each cue [optional, otherwise optimized]
% Wfix ... fixed weights for cues [optional, otherwise optimized]
% expnum ... experiment number for separate rank evaluations
% 
% Output parameters:
% Pext ... predicted externalization scores (first) for individual cues and 
%          (then) for cue combinations (see decision strategies)
% S ... sensitivity parameters, one for each cue
% W ... weights to combine selected cues and percentages of cue selection 
%       for optimally predicting *actual*
% MSE ... mean squared error
% BIC ... Bayesian information criterion
% RC ... rank correlation
%
% Decision strategies:
% 1) optimal weighting ("oracle")
% 2) fixed weighting ("weighted sum model")
% 3) minimum externalization selection ("loser takes all")
% 4) median externalization selection ("median takes all")
% 5) maximum externalization selection ("winner takes all")

Ncues = length(predCue);

% Definitions of error functions
f_cue2E = @(cue,s) baumgartner2021_mapping(cue,s,Erange,Eoffset);
if size(actual,2) > 1 % Individual data
  f_predErr = @(p) nanmean((p(:)-actual(:)).^2);%nanrms(p(:)-actual(:)) / (Erange-Eoffset);
  f_predErrA = @(p,a) nanmean((p(:)-a(:)).^2); % used for bootstrapping
else % Average data
  f_predErr = @(p) nanmean((mean(p,2) - actual).^2);%nanrms(p(:)-actual(:)) / (Erange-Eoffset);
  f_predErrA = @(p,a) nanmean((mean(p,2) - a).^2);
end

% Optimize sensitivities for single-cue externalization predictions
if not(exist('S','var')) || isempty(S)
  S.WSM = ones(Ncues,1); % sensitivity to cue
  for cc = 1:Ncues
    S.WSM(cc) = fminsearch(@(s) f_predErr(f_cue2E(predCue{cc},s)),S.WSM(cc));
  end
end

if ~isfield(S,'WSM')
  Stmp = S;
  clear S
  S.WSM = Stmp;
end

% Compute single-cue predictions
predExtSingle = nan([size(predCue{1}),Ncues]);
for cc = 1:Ncues
  if all(isnan(predCue{cc}(:)))
    predExtSingle(:,:,cc) = zeros(size(predCue{cc}));
  else
    predExtSingle(:,:,cc) = f_cue2E(predCue{cc},S.WSM(cc));
  end
end

% Optimize weights (S adopted from single-cue predictions)
W = 1/Ncues*ones(size(S.WSM)); % weighting of cue
f_weight = @(p,w) reshape( reshape(p,[],length(w))*w(:) ,size(p,1),[]);
f_opt = @(p,w) f_predErr( f_weight(p,abs(w)./sum(abs(w))) );
o = optimset('Display','off');
W = fminsearch(@(w) f_opt(predExtSingle,w) , W,o); 
W = abs( W );
W = W/sum(W);
W(W<.005) = 0; % threshold the weights
W = W/sum(W);

if not(exist('Wfix','var')) || isempty(Wfix) || sum(isnan(Wfix))>0
  Wfix = W;
end

% Selection strategies
function [E,cueStat] = local_f_select(cue,s,type) 
    Eall = reshape(...
      f_cue2E(cue(:),repmat(s(:),length(cue(:))/length(s(:)),1)),...
      size(cue,1),size(cue,2),size(cue,3));
    switch type
      case 'LTA'
        [E,iCue] = nanmin(Eall,[],3);
      case 'MTA'
        E = median(Eall,3,'omitnan');
        [~,iCue] = nanmin((Eall-repmat(E,1,1,size(Eall,3))).^2,[],3);
      case 'WTA'
        [E,iCue] = nanmax(Eall,[],3);
      otherwise
        error('Undefined type of decision strategy')
    end
    cueStat = 100*histcounts(iCue,1:Ncues+1)/length(iCue);
end

% Optimize S for selective strategies
predCueMtx = cat(3,predCue{:});
cueStat = nan(Ncues,3);
o = optimset('MaxFunEvals',300*length(S),'Display','off');
S.LTA=fminsearch(@(s) f_predErr(local_f_select(predCueMtx,s,'LTA')),S.WSM,o);
S.MTA=fminsearch(@(s) f_predErr(local_f_select(predCueMtx,s,'MTA')),S.WSM,o);
S.WTA=fminsearch(@(s) f_predErr(local_f_select(predCueMtx,s,'WTA')),S.WSM,o);

% Externalization predictions
predExtComb(:,:,1) = f_weight(predExtSingle,W); % Oracle
predExtComb(:,:,2) = f_weight(predExtSingle,Wfix); % WSM
[predExtComb(:,:,3),cueStat(:,1)] = local_f_select(predCueMtx,S.LTA,'LTA');
[predExtComb(:,:,4),cueStat(:,2)] = local_f_select(predCueMtx,S.MTA,'MTA');
[predExtComb(:,:,5),cueStat(:,3)] = local_f_select(predCueMtx,S.WTA,'WTA');
predExt = cat(3,predExtSingle,predExtComb); % concatenate single-cue and combinations
Pext = squeeze(num2cell(predExt,1:2));

% Sensitivities
S = [S.WSM(:),S.WSM(:),S.LTA(:),S.MTA(:),S.WTA(:)];

% Cue weights and selection frequencies
W = [W(:),Wfix(:),cueStat];

%% Performance metrics
% Mean squared errors
Nstrat = 5;
MSE.m = nan(Ncues+Nstrat,1);
MSE.ci = nan(Ncues+Nstrat,2);
for ii = 1:Ncues+Nstrat
  if ii <= Ncues
    p = predExtSingle(:,:,ii);
  else
    p = predExtComb(:,:,ii-Ncues);
  end
  MSE.m(ii) = f_predErr(p);
  MSE.ci(ii,:) = bootci(1000,f_predErrA,p,actual);
end

% Bayesian Information Criteria
n = size(actual,1); % # conditions (no individual parameters for # listeners)
mc = ones(size(MSE.m)); % model complexity = # parameters; 1 sensitivity parameter per cue
mc(Ncues+(1:2)) = Ncues; % weighted combinations -> adopted sensitivity parameters but optimized weights
mc(Ncues+3:end) = Ncues; % selection startegies -> optimized sensitivity parameters
BIC = n*log(MSE.m) + mc*log(n);

% Rank correlations evaluated within experiments and combined afterwards
function RC = local_f_cumrankcorr(p,actual,expnum)
  n = length(unique(expnum));
  rc = zeros(n,1);
  nee = zeros(n,1);
  for ee = 1:n
    idee = expnum==ee;
    nee(ee) = sum(idee);
    if nee(ee) > 1
      rc(ee) = corr(p(idee),actual(idee),'type','spearman');
    else
      rc(ee) = 0;
      nee(ee) = 0;
    end
  end
  RC = nee'./sum(nee)*rc;
end
if exist('expnum','var')
  RC.m = nan(Ncues+Nstrat,1);
  RC.ci = nan(Ncues+Nstrat,2);
  for ii = 1:Ncues+Nstrat
    if ii <= Ncues
      p = predExtSingle(:,:,ii);
    else
      p = predExtComb(:,:,ii-Ncues);
    end
    RC.m(ii) = local_f_cumrankcorr(p,actual,expnum);
    RC.ci(ii,:) = bootci(300,@local_f_cumrankcorr,p,actual,expnum);
  end
end

end



