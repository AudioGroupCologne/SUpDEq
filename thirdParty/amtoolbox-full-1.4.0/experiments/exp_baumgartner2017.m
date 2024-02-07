function data = exp_baumgartner2017(varargin)
%EXP_BAUMGARTNER2017 Experiments of Baumgartner et al. (2017)
%   Usage: data = exp_baumgartner2017(flag) 
%
%   EXP_BAUMGARTNER2017(flag) reproduces figures of the study from 
%   Baumgartner et al. (2017).
%
%   The following flags can be specified
%
%     'boyd2012'    models experiments from Boyd et al. (2012; Fig.1, top).
%                   Average externalization ratings of 1 talker for NH participants 
%                   against mix point as a function of microphone position (ITE/BTE) 
%                   and frequency response (BB/LP). The reference condition (ref) is 
%                   the same as ITE/BB. Error bars show SEM. 
%
%     'hartmann1996'  models experiments from Hartmann & Wittenberg (1996; Fig.7-8)
%                     1st panel: Synthesis of zero-ILD signals. Only the harmonics 
%                     from 1 to nprime had zero interaural level difference; 
%                     harmonics above nprime retained the amplitudes of the baseline  
%                     synthesis. Externalization scores as a function of the boundary  
%                     harmonic number nprime. Fundamental frequency of 125 Hz.
%                     2nd panel: Synthesis of signals to test the ISLD hypothesis. 
%                     Harmonics at and below the boundary retained only the interaural 
%                     spectral level differences of the baseline synthesis. Higher 
%                     harmonics retained left and right baseline harmonic levels. 
%                     Externalization scores as a function of the boundary
%                     frequency.
%
%     'hassager2016'  models experiments from Hassager et al. (2016; Fig.6). 
%                     The mean of the seven listeners perceived sound source 
%                     location (black) as a function of the bandwidth factor 
%                     and the corresponding model predictions (colored). 
%                     The model predictions have been shifted slightly to the right 
%                     for a better visual interpretation. The error bars are one 
%                     standard error of the mean.
%
%     'baumgartner2017'   models experiments from Baumgartner et al. (2017). 
%                         Effect of HRTF spectral contrast manipulations on sound externalization. 
%                         Externalization scores were derived from paired comparisons via 
%                         Bradley-Terry-Luce modeling.
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
%   To display results for Fig.1 from Boyd et al. (2012) use :
%
%     exp_baumgartner2017('boyd2012');
%
%   To display results for Fig.7-8 from Hartmann & Wittenberg (1996) use :
%
%     exp_baumgartner2017('hartmann1996');
%
%   To display results for Fig.6 from Hassager et al. (2016) use :
%
%     exp_baumgartner2017('hassager2016');
%
%   References:
%     R. Baumgartner, P. Majdak, H. Colburn, and B. Shinn-Cunningham.
%     Modeling sound externalization based on listener-specific spectral
%     cues. In Acoustics ‘17 Boston: The 3rd Joint Meeting of the Acoustical
%     Society of America and the European Acoustics Association, Boston, MA,
%     Jun 2017.
%     
%     A. W. Boyd, W. M. Whitmer, J. J. Soraghan, and M. A. Akeroyd. Auditory
%     externalization in hearing-impaired listeners: The effect of pinna cues
%     and number of talkers. J. Acoust. Soc. Am., 131(3):EL268--EL274, 2012.
%     [1]www: ]
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
%     1. http://dx.doi.org/10.1121/1.3687015
%     
%
%   See also: hassager2016 baumgartner2017
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_baumgartner2017.php


%   #Author: Robert Baumgartner (2017), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


definput.import={'amt_cache','baumgartner2017'};
definput.flags.type = {'missingflag','boyd2012','hartmann1996','hassager2016','baumgartner2017'};
definput.flags.quickCheck = {'','quickCheck'};
definput.keyvals.Sintra = 2;
definput.keyvals.Sinter = 2;
definput.keyvals.Interaction = {'MSG','ITIT'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

Ncues = 7;
% Symbol order for plotting
symb = {'-o','-s','-d','-<','->','-^','-v','-p',':*','--x'};
colors = 0.8*[zeros(1,3);hsv(Ncues+2)];
InteractionLbl = [kv.Interaction{1},' x ',kv.Interaction{2}];

%% Hassager et al. (2016)
if flags.do_hassager2016
  azi = [0,50];
  flp = 6000; % lowpass filtered noise stimuli
  
  Pext_A = data_hassager2016;
  B = Pext_A.B;
  
  fncache = ['hassager2016'];
  [Pext,cues,cueLbl] = amt_cache('get',fncache,flags.cachemode);
  if isempty(cues)
    
    data = data_baumgartner2017;
    if flags.do_quickCheck
      data = data(1:5);
    end

%     Pext = nan(length(B),length(data),length(azi));
%     Pext = {Pext,Pext};
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
          [~,cues(iB,iazi,isubj,:),cueLbl] = baumgartner2017(target,template,'flow',100,'fhigh',flp,'lat',azi(iazi));
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
  dimCues = 4;
  Nsubj = size(cues,dimSubj);

  cues = num2cell(cues,1:dimCues-1);
  % optimize sensitivities
%   predCue = reshape(cues,length(B)*length(azi)*Nsubj,Ncues);
  actual = repmat(Pext_A.rating(:),Nsubj,1);
  Erange = 4;
  Eoffset = 1;
  InteractionIds = find(ismember(cueLbl(:,1),kv.Interaction));
  [S,W,Pext,dE] = local_cueWeightOpti(cues,actual,Erange,Eoffset,InteractionIds);
  
  PextLbl = [cueLbl(:,1);{'Comb.'};{InteractionLbl}];
  
  cueTab = table([S;nan(2,1)],dE,[W;nan(2,1)],[cueLbl(:,2);{'Weighted combination'};{[InteractionLbl, ' interaction']}],...
    'RowNames',PextLbl,...
    'VariableNames',{'Sensitivity','Error','Weight','Description'});
  amt_disp(flags.type,'documentation');
  amt_disp(cueTab,'documentation');
  
  %% Plot
  BplotTicks = logspace(log10(0.25),log10(64),9);
  BplotTicks = round(BplotTicks*100)/100;
  BplotStr = ['Ref.';num2str(BplotTicks(:))];
  BplotTicks = [BplotTicks(1)/2,BplotTicks];
  B(1) = B(2)/2;
  dataLbl = [{'Actual'};PextLbl];
  idleg = not(ismember(dataLbl,'ITSD'));
  figure 
  hax = local_tight_subplot(1,2,0,[.15,.1],[.1,.05]);
  for iazi = 1:length(azi)
%     subplot(1,2,iazi)
    axes(hax(iazi))
    h(1) = plot(B,Pext_A.rating(:,iazi),symb{1},'Color',colors(1,:));
    hold on
    for m = 1:length(Pext)
      dx = 1+0.02*(m - (Ncues+1)/2);
      h(m+1) = errorbar(dx*B,mean(Pext{m}(:,iazi,:),dimSubj),std(Pext{m}(:,iazi,:),0,dimSubj)/sqrt(Nsubj),...
        symb{m+1},'Color',colors(m+1,:));
    end
    set(h,'MarkerFaceColor','w')
    set(gca,'XTick',BplotTicks,'XTickLabel',BplotStr,'XScale','log')
    axis([BplotTicks(1)/1.5,BplotTicks(end)*1.5,0.8,5.2])
    xlabel('Bandwidth Factor [ERB]')
    if iazi==1
      ylabel('Mean Externalization Rating')
    else
      set(gca,'YTickLabel',{})
    end
    title([num2str(azi(iazi)),'\circ'])
  end
  leg = legend(h(idleg),dataLbl(idleg),'Location','southwest');
  set(leg,'Box','off')
  
  %% Output
  data = {cueTab,Pext};
  
end

%% Hartmann & Wittenberg (1996)
if flags.do_hartmann1996
  
  cond = {'ILD','ISLD'};
  condPlotLbl = {'ILD to 0','ILD const.'};
  nprime = [0,1,8,14,19,22,25,38];
  
  fncache = ['hartmann1996'];
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
          [~,cues{ee}(isub,nn,:),cueLbl] = baumgartner2017(target,template,'flow',100,'fhigh',5000,'lat',azi);
        end
      end
      amt_disp([num2str(isub),' of ',num2str(length(data)),' subjects completed.']);
    end
    if not(flags.do_quickCheck)
      amt_cache('set',fncache,Pext,cues,cueLbl);
    end
  end
  
  %% Cue weighting optimization procedure
  
  % align actual and modeled data format
  actual = [];
  predCue = [];
  for ee = 1:length(cond)
    act = data_hartmann1996(cond{ee});
    actual = cat(1,actual,act.avg.Escore(:));
    tmp = squeeze(mean(cues{ee},1));
    tmp = interp1(nprime,tmp,act.avg.nprime);
    predCue = cat(1,predCue,tmp);
  end
  predCue = num2cell(predCue,1);
  
  
  Erange = 3;
  Eoffset = 0;
  InteractionIds = find(ismember(cueLbl(:,1),kv.Interaction));
  if length(InteractionIds) ~=2
    error('RB: Check interaction labels.')
  end
  [S,W,~,dE] = local_cueWeightOpti(predCue,actual,Erange,Eoffset,InteractionIds);
  % calculate separate errors for Exp 1 and 2
  dE1 = nan(size(dE));
  dE2 = nan(size(dE));
  for ss = 1:length(S)
    E = local_cue2E(S(ss),predCue{ss},Erange,Eoffset);
    dE1(ss) = local_nanrms(E(1:5)-actual(1:5)) / (Erange-Eoffset);
    dE2(ss) = local_nanrms(E(6:9)-actual(6:9)) / (Erange-Eoffset);
  end
  
  PextLbl = [cueLbl(:,1);{'Comb.'};{InteractionLbl}];
  
  cueTab = table([S;nan(2,1)],dE1,dE2,dE,[W;nan(2,1)],[cueLbl(:,2);{'Weighted combination'};{[InteractionLbl, ' interaction']}],...
    'RowNames',PextLbl,...
    'VariableNames',{'Sensitivity','ErrorExp1','ErrorExp2','Error','Weight','Description'});
  amt_disp(flags.type,'documentation');
  amt_disp(cueTab,'documentation');

  %% Individual data
  dimCues = 3;
%   Ncues = size(cues{1},dimCues);
  Pext = cell(Ncues+1,length(cond));
  for ee = 1:length(cond)
    for cc = 1:Ncues
      Pext{cc,ee} = local_cue2E(S(cc),cues{ee}(:,:,cc),Erange,Eoffset);
    end
    Pext{Ncues+1,ee} = W(1)*Pext{1,ee};
    for cc = 2:Ncues
      Pext{Ncues+1,ee} = nansum(cat(dimCues,Pext{Ncues+1,ee},W(cc)*Pext{cc,ee}),dimCues);
    end
    Pext{Ncues+2,ee} = sqrt(Pext{InteractionIds(1),ee}.*Pext{InteractionIds(2),ee}); % interaction term
  end
  
  %% Plot
  dataLbl = [{'Actual'};PextLbl];
  idleg = not(ismember(dataLbl,'ITSD'));
  Ns = size(Pext{1},1);
  figure
  hax = local_tight_subplot(1,2,0,[.15,.1],[.1,.05]);
  for ee = 1:length(cond)
    act = data_hartmann1996(cond{ee});
%     subplot(1,2,ee)
    axes(hax(ee))
    h(1) = plot(act.avg.nprime,act.avg.Escore,symb{1},'Color',colors(1,:));
    hold on
    for cc = 1:size(Pext,1)
      if idleg(cc+1)
        dx = 0.1*(cc - size(Pext,1)/2);
        h(cc+1) = errorbar(nprime+dx,mean(Pext{cc,ee}),std(Pext{cc,ee})/sqrt(Ns),...
          symb{cc+1},'Color',colors(cc+1,:));
      end
    end
    set(h(idleg),'MarkerFaceColor','w')
    if ee == 1
      xlabel('n^{\prime} (Highest harmonic with ILD = 0)')
      ylabel('Externalization score')
    else
      xlabel('n^{\prime} (Highest harmonic with altered amplitudes)')
      set(gca,'YTickLabel',{})
    end
    title(condPlotLbl{ee})
    axis([-2.5,42.5,-0.4,3.9])
    idx = false(size(idleg));
    if ee == 1
      idx(1:Ncues) = true;
    else
      idx(Ncues+1:end) = true;
    end
    legend(h(and(idleg,idx)),dataLbl(and(idleg,idx)),'Location','southwest','box','off');
  end
  
  % Output
  data = {cueTab,Pext};
end

%% Boyd et al. (2012) -----------------------------------------------------
if flags.do_boyd2012
  
  flags.do_noise = false;
  flags.do_georganti2013 = false;
  flags.do_BRIR = true;
  
  flp = [nan,6500]; % Low-pass cut-off frequency
  azi = -30;
  
  subjects = data_boyd2012;
  if flags.do_noise || flags.do_BRIR
    subjects = subjects([3,6,7]); % only the ones with BRIRs
  else
    subjects = subjects([1,3,4,6,7]); % only the ones with stimuli
  end
%   fprintf(['Note that only 3 of originally 7 listeners are simulated and compared because\n',...
%     'the listener-specific BRIRs for the other 4 listeners are not available.\n'])
  mix = subjects(1).Resp.mix/100;
    
  fncache = ['boyd2012_lp'];
  [E,cues,cueLbl] = amt_cache('get',fncache,flags.cachemode);
  if isempty(cues)
    
    mix(1) = 1;

    if flags.do_noise || flags.do_BRIR
      sig = noise(50e3,1,'pink');
      fs = subjects(1).BRIR.fs;
      [b,a]=butter(10,2*flp(2)/fs,'low');
    end

    E.meta = {'actual','interaural','monaural'};
    if flags.do_georganti2013
      E.meta(4) = {'BSMD'};
    end
    E.all = repmat({nan([length(mix),2,2,length(subjects)])},1,length(E.meta));
    cues = nan(length(mix),2,2,length(subjects),Ncues);
    for isub = 1:length(subjects)
      E.all{1}(:,:,:,isub) = cat(3,...
        [subjects(isub).Resp.ITE_BB_1T(:),subjects(isub).Resp.BTE_BB_1T(:)],...
        [subjects(isub).Resp.ITE_LP_1T(:),subjects(isub).Resp.BTE_LP_1T(:)]);
%       ITE = subjects(isub).Obj;
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
        targetSet = [repmat({template},[1,2,2]);targetSet];
      end

      %% Model simulations
      flow = 100;
      fhigh = [16e3,flp(2)];
            
      for c = 1:2
        for lp = 1:2
          for m = 1:length(mix)
            if m==1 && (c~=1 || lp~=1) % same reference condition
              for ee = 2:length(E.all)
                E.all{ee}(m,c,lp,isub) =  E.all{ee}(m,1,1,isub);
                targetSet{m,c,lp} = targetSet{m,1,1};
              end
            else
  %             target{m,c,lp} = sig_boyd2012(stim{c},sig,mix(m),flp(lp),fs);
              % mixing
              if flags.do_noise || flags.do_BRIR
                targetSet{m,c,lp} = mix(m)*stim{c} + (1-mix(m))*stim{3};
                % low-pass filtering
                if lp == 2
                  targetSet{m,c,lp} = filter(b,a,targetSet{m,c,lp});
                end
              end
              target = targetSet{m,c,lp};
  % % Description in paper is ambiguous about whether broadband ILDs were
  % % adjusted to the template:
  %           temSPL = dbspl(template); 
  %             for ch = 1:2
  %               target{m,c,lp}(:,ch) = scaletodbspl(target{m,c,lp}(:,ch),temSPL(ch), 100);
  %             end
              amt_disp(num2str(mix(m)),'volatile');
              if flags.do_BRIR
                flags.do_spectralCueEchoSuppression = true;
              end
              [~,cues(m,c,lp,isub,:),cueLbl] = baumgartner2017(target,template...
                ,'argimport',flags,kv,'flow',flow,'fhigh',fhigh(1),'lat',azi,'reflectionOnsetTime',5e-3); % replace fhigh(1) with fhigh(lp) for manual bandwidth adjustment

              if flags.do_georganti2013
                geor.fs = subjects(1).fs;
                geor.timeFr = 1;
                geor.fmin = flow;
                geor.fmax = fhigh(lp);
                E.all{4}(m,c,lp,isub) = median(georganti2013(target,geor)); % median over time
              end
              
            end
          end
        end
      end
      amt_disp([num2str(isub),' of ',num2str(length(subjects)),' subjects completed.']);

    end

    if flags.do_georganti2013
      E.all{4} = 100*E.all{4}./repmat(E.all{4}(1,1,1,:),[length(mix),2,2,1]);
    end
    
    if not(flags.do_quickCheck)
      amt_cache('set',fncache,E,cues,cueLbl);
    end
  end

  
  %% Cue weighting optimization procedure
  dimSubj = 4;
  dimCues = 5;
  Erange = 100;
  Eoffset = 0;
  
  Nsubj = size(cues,dimSubj);
  if flags.do_BRIR % comparison on individual basis
    actual = E.all{1};
  else % average comparison
    actual = repmat(mean(E.all{1},4),[ones(1,Nsubj-1),Nsubj]);
  end
  
  cuesC = squeeze(num2cell(cues,1:dimCues-1));
  InteractionIds = find(ismember(cueLbl(:,1),kv.Interaction));
  [S,W,Pext,dE] = local_cueWeightOpti(cuesC,actual,Erange,Eoffset,InteractionIds);
  
  E.all(2:Ncues+3) = Pext;
  
  PextLbl = [cueLbl(:,1);{'Comb.'};{InteractionLbl}];
  
  cueTab = table([S;nan(2,1)],dE,[W;nan(2,1)],[cueLbl(:,2);{'Weighted sum'};{[InteractionLbl, ' interaction']}],...
    'RowNames',PextLbl,...
    'VariableNames',{'Sensitivity','Error','Weight','Description'});
  amt_disp(flags.type,'documentation');
  amt_disp(cueTab,'documentation');
  
  %% Average data
  for ee = 1:length(E.all)
    E.m{ee} = mean(E.all{ee},4);
    E.se{ee} = std(E.all{ee},0,4);
  end
  
  %% Plot
  flags.do_plot_individual = false;
  
  condLbl = {'ITE / BB','BTE / BB','ITE / LP','BTE / LP'};
  mix = mix(2:end)*100;
  
  leglbl = [{'Actual'};PextLbl];
  idleg = not(ismember(leglbl,'ITSD'));
  if not(flags.do_plot_individual)
    figure
    hax = local_tight_subplot(2,2,0,[.15,.1],[.1,.05]);
    for cc = 1:length(condLbl)
      axes(hax(cc))
      for ee = 1:length(E.m)
        if idleg(ee)
          dx = 2*(ee - length(E.m)/2);
          h = errorbar(mix+dx,E.m{ee}(2:end,cc),E.se{ee}(2:end,cc),...
            symb{ee},'Color',colors(ee,:));
          set(h,'MarkerFaceColor','w')
          hold on
        end
      end
      set(gca,'XDir','reverse')

      if cc == 3 || cc == 4
        xlabel('Mix (%)')
      else
        set(gca,'XTickLabel',[])
      end

      if cc == 1 || cc == 3
        ylabel('Externalization score')
      else
        set(gca,'YTickLabel',[])
      end
%       title(condLbl{cc})
      text(110,0,condLbl{cc})
      axis([-20,120,-15,115])
      if cc == 4
%         leg = legend(leglbl(idleg));
%         set(leg,'Box','off','Location','eastoutside')
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
        xlabel('Mix (%)')
        if cc == 1
          ylabel('Externalization score')
        else
          set(gca,'YTickLabel',[])
        end
        title(condLbl{cc})
        axis([-20,120,-5,105])
        if cc == 4
          leg = legend(dataLbl);
          set(leg,'Box','off','Location','north')
        end
      end
    end
    
  end
%   set(leg,'Location','eastoutside','Position',get(leg,'Position')+[.1,.2,0,0])

%% Output
data = {cueTab,E};

end

if flags.do_baumgartner2017
  
  % Behavioral results
  fncache = 'baumgartner2017_E_BTL';
  data = amt_cache('get',fncache,flags.cachemode);
  if isempty(data)
    exp2 = data_baumgartner2017looming('exp2');
    Nsubj = length(exp2.meta.subject);
    iISI = exp2.meta.ISI == 0.1;
    iresp = strcmp(exp2.meta.response,'approaching');
    A = {1;2;3}; % arbitrary indices to define BTL model
  %   Alabel = {'C = 0','C = 0.5','C = 1'};
    data.C = 0:0.5:1;
    data.subj = [];%exp2.meta.subject;
    data.Escore = [];%nan(Nsubj,length(A));
    data.BTL_chistat = [];%nan(Nsubj,1);
    data.azi = [exp2.rawData.azimuth];
    M = zeros(3,3);
    iM = [7,4,8,3,2,6];
    iss = 0;
    for ss = 1:Nsubj
      M(iM) = exp2.data(ss,1:6,iresp,iISI)+eps;
      [E,chistat,~,lL] = local_OptiPt(M,A);
      if chistat(1) < 10 && chistat(1) > 0.001
        iss = iss+1;
        data.subj{iss} = exp2.meta.subject{ss};
        data.BTL_chistat(iss) = chistat(1);
        data.BTL_logLikeli(iss) = lL;
        data.Escore(iss,:) = E/E(3);
      end
    end 
    amt_cache('set',fncache,data)
  end
  Nsubj = length(data.subj);
  
  hrtf = data_baumgartner2017looming('exp2','hrtf');
  
  % Modeling
  fncache = ['baumgartner2017'];
  [Pext,cues,cueLbl] = amt_cache('get',fncache,flags.cachemode);
  if isempty(cues)
    
    Pext = nan(length(data.C),Nsubj);
    cues = nan(length(data.C),Nsubj,Ncues);
    for isubj = 1:Nsubj
      template = hrtf(ismember({hrtf.id},data.subj(isubj))).Obj;
      for iC = 1:length(data.C)
        target = sig_baumgartner2017looming(template,data.C(iC));
        [~,cues(iC,isubj,:),cueLbl] = baumgartner2017(target,template,...
                          'lat',data.azi(isubj),'flow',1000,'fhigh',18e3);
        
      end
      amt_disp([num2str(isubj),' of ',num2str(Nsubj),' subjects completed.']);
    end
    amt_cache('set',fncache,Pext,cues,cueLbl);
  end
 
  %% Cue weighting optimization procedure
  dimSubj = 2;
  dimCues = 3;
  Erange = 1;
  Eoffset = 0;
  
  Nsubj = size(cues,dimSubj);
  actual = permute(data.Escore,[2,1]);
  
  cuesC = squeeze(num2cell(cues,1:dimCues-1));
  InteractionIds = find(ismember(cueLbl(:,1),kv.Interaction));
  [S,W,Pext,dE] = local_cueWeightOpti(cuesC,actual,Erange,Eoffset,InteractionIds);
  
  PextLbl = [cueLbl(:,1);{'Comb.'};{InteractionLbl}];
  
  cueTab = table([S;nan(2,1)],dE,[W;nan(2,1)],[cueLbl(:,2);{'Weighted sum'};{[InteractionLbl, ' interaction']}],...
    'RowNames',PextLbl,...
    'VariableNames',{'Sensitivity','Error','Weight','Description'});
  amt_disp(flags.type,'documentation');
  amt_disp(cueTab,'documentation');
  
  %% Figure
  flags.do_plot_individual = false;
  E_all = [{data.Escore'};Pext];
  dataLbl = [{'Actual'};PextLbl];
  figure
  if not(flags.do_plot_individual)
    for iE = 1:length(E_all)
      dx = 0.01*(iE - length(E_all)/2);
      h(iE) = errorbar(data.C+dx,mean(E_all{iE},2),std(E_all{iE},0,2)/sqrt(Nsubj),...
        symb{iE},'Color',colors(iE,:));
      set(h,'MarkerFaceColor','w')
      hold on
    end
    set(gca,'XTick',data.C)
    axis([-0.2,1.2,-0.1,1.1])
    ylabel('Externalization')
    xlabel('Spectral contrast, C')
  else % flags.do_plot_individual
    for isubj = 1:Nsubj
      subplot(3,ceil(Nsubj/3),isubj)
      for iE = 1:length(E_all)
        dx = 0.01*(iE - length(E_all)/2);
        h(iE) = plot(data.C+dx,E_all{iE}(:,isubj),...
          symb{iE},'Color',colors(iE,:));
        set(h,'MarkerFaceColor','w')
        hold on
      end
      set(gca,'XTick',data.C)
      axis([-0.2,1.2,-0.1,1.1])
      ylabel('Externalization')
      xlabel('Spectral contrast, C')
    end
  end
  idleg = not(ismember(dataLbl,'ITSD'));
  legend(h(idleg),dataLbl(idleg),'Location','EastOutside','box','off')

  %% Output
  data = {cueTab,Pext};
  
end




end

%%%%%%%%%% INTERNAL FUNCTIONS FOR VISUALIZATION %%%%%%%%%%
function ha = local_tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
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
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .1],[.1 .1])
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

function [p,chistat,u,lL_eba,lL_sat,fit,cova] = local_OptiPt(M,A,s)
% OptiPt parameter estimation for BTL/Pretree/EBA models
%   p = OptiPt(M,A) estimates the parameters of a model specified
%   in A for the paired-comparison matrix M. M is a matrix with
%   absolute frequencies. A is a cell array.
%
%   [p,chistat,u] = OptiPt(M,A) estimates parameters and reports
%   the chi2 statistic as a measure of goodness of fit. The vector
%   of scale values is stored in u.
%
%   [p,chistat,u,lL_eba,lL_sat,fit,cova] = OptiPt(M,A,s) estimates
%   parameters, checks the goodness of fit, computes the scale values,
%   reports the log-likelihoods of the model specified in A and of the
%   saturated model, returns the fitted values and the covariance
%   matrix of the parameter estimates. If defined, s is the starting
%   vector for the estimation procedure. Otherwise each starting value
%   is set to 1/length(p).
%   The minimization algorithm used is FMINSEARCH.
%
%   Examples
%     Given the matrix M = 
%                            0    36    35    44    25
%                           19     0    31    37    20
%                           20    24     0    46    24
%                           11    18     9     0    13
%                           30    35    31    42     0
%
%     A BTL model is specified by A = {[1];[2];[3];[4];[5]}
%     Parameter estimates and the chi2 statistic are obtained by
%       [p,chistat] = OptiPt(M,A)
%
%     A Pretree model is specified by A = {[1 6];[2 6];[3 7];[4 7];[5]} 
%     A starting vector is defined by s = [2 2 3 4 4 .5 .5]
%     Parameter estimates, the chi2 statistic, the scale values, the
%     log-likelihoods of the Pretree model and of the saturated model,
%     the fitted values, and the covariance matrix are obtained by
%       [p,chistat,u,lL_eba,lL_sat,fit,cova] = OptiPt(M,A,s)
%
% Authors: Florian Wickelmaier (wickelmaier@web.de) and Sylvain Choisel
% Last mod: 03/JUL/2003
% For detailed information see Wickelmaier, F. & Schmid, C. (2004). A Matlab
% function to estimate choice model parameters from paired-comparison data.
% Behavior Research Methods, Instruments, and Computers, 36(1), 29-40.

I = length(M);  % number of stimuli
mmm = 0;
for i = 1:I
  mmm = [mmm max(A{i})];
end
J = max(mmm);  % number of pt parameters
if(nargin == 2)
  p = ones(1,J)*(1/J);  % starting values
elseif(nargin == 3)
  p = s;
end

for i = 1:I
  for j = 1:I
    diff{i,j} = setdiff(A{i},A{j});  % set difference
  end
end

p = fminsearch(@ebalik,p,optimset('Display','iter','MaxFunEvals',10000,...
    'MaxIter',10000),M,diff,I);  % optimized parameters
lL_eba = -ebalik(p,M,diff,I);  % likelihood of the specified model

lL_sat = 0;  % likelihood of the saturated model
for i = 1:I-1
  for j = i+1:I
    lL_sat = lL_sat + M(i,j)*log(M(i,j)/(M(i,j)+M(j,i)))...
                    + M(j,i)*log(M(j,i)/(M(i,j)+M(j,i)));
  end
end

fit = zeros(I);  % fitted PCM
for i = 1:I-1
  for j = i+1:I
    fit(i,j) = (M(i,j)+M(j,i))/(1+sum(p(diff{j,i}))/sum(p(diff{i,j})));
    fit(j,i) = (M(i,j)+M(j,i))/(1+sum(p(diff{i,j}))/sum(p(diff{j,i})));
  end
end

chi = 2*(lL_sat-lL_eba);
df =  I*(I-1)/2 - (J-1);
chistat = [chi df];  % 1-chi2cdf(chi,df)];  % goodness-of-fit statistic

u = sum(p(A{1}));  % scale values
for i = 2:I
  u = [u sum(p(A{i}))];
end

H = local_hessian('ebalik',p',M,diff,I);
C = inv([H ones(J,1); ones(1,J) 0]);
cova = C(1:J,1:J);
end

function lL_eba = ebalik(p,M,diff,I)  % computes the likelihood

if min(p)<=0  % bound search space
  lL_eba = inf;
  return
end

thesum = 0;
for i = 1:I-1
  for j = i+1:I
    thesum = thesum + M(i,j)*log(1+sum(p(diff{j,i}))/sum(p(diff{i,j})))...
                    + M(j,i)*log(1+sum(p(diff{i,j}))/sum(p(diff{j,i})));
  end
end
lL_eba = thesum;
end

function H = local_hessian(f,x,varargin)  % computes numerical Hessian
% based on a solution posted on Matlab Central by Paul L. Fackler

k = size(x,1);
fx = feval(f,x,varargin{:});
h = eps.^(1/3)*max(abs(x),1e-2);
xh = x+h;
h = xh-x;
ee = sparse(1:k,1:k,h,k,k);

g = zeros(k,1);
for i = 1:k
  g(i) = feval(f,x+ee(:,i),varargin{:});
end

H = h*h';
for i = 1:k
  for j = i:k
    H(i,j) = (feval(f,x+ee(:,i)+ee(:,j),varargin{:})-g(i)-g(j)+fx)...
                 / H(i,j);
    H(j,i) = H(i,j);
  end
end
end

function [S,W,Pext,dE] = local_cueWeightOpti(predCue,actual,Erange,Eoffset,InteractionIds)
% optimization routine for cue-specific sensitivities S and optimal
% weighting W of cues. Make sure that actual matches the dimension of each
% cue-specific cell of predCue.

% Definitions of functions
f_cue2E = @(s,cue) local_cue2E(s,cue,Erange,Eoffset);
f_predErr = @(p) local_nanrms(p(:)-actual(:)) / (Erange-Eoffset);

% optimize sensitivities
Ncues = length(predCue);
S = ones(Ncues,1); % sensitivity to cue
Pext = cell(Ncues+1,1);
dE = nan(Ncues+1,1);
predE = nan(length(predCue{1}(:)),Ncues);
for cc = 1:Ncues
  if not(all(isnan(predCue{cc}(:))))
    S(cc) = fminsearch(@(s) f_predErr(f_cue2E(s,predCue{cc})) ,S(cc));
% splot = logspace(log10(S(cc))-1,log10(S(cc))+1,100);
% for ii = 1:length(splot)
%   Eplot(ii) = f_predErr(f_cue2E(splot(ii),predCue{cc}));
% end
% figure; semilogx(splot,Eplot)
  end
  
  Pext{cc} = f_cue2E(S(cc),predCue{cc});
  predE(:,cc) = Pext{cc}(:);
  dE(cc) = f_predErr(predE(:,cc)); % mean(abs( predE(:,cc) - actual )) / (Erange-Eoffset);
end

% optimize weighting
W = zeros(size(S)); % weighting of cue
isan = not(all(isnan(predE)));
Wtmp = abs( fminsearch(@(w) f_predErr(predE(:,isan) * abs(w)./sum(abs(w))),W(isan)) ); %rms( predE * abs(w)./sum(abs(w)) - actual ),W) );
Wtmp = Wtmp/sum(Wtmp);
dE(cc+1) = f_predErr(predE(:,isan)*Wtmp);
W(isan) = Wtmp;

% Combined externalization score
dimCues = ndims(Pext{1})+1;
Pext{Ncues+1} = W(1)*Pext{1};
for cc = 2:Ncues
  Pext{Ncues+1} = nansum(cat(dimCues,Pext{Ncues+1},W(cc)*Pext{cc}),dimCues);
end

% Interaction term
dE(cc+2) = f_predErr(sqrt(predE(:,InteractionIds(1)).*predE(:,InteractionIds(2))));
Pext{Ncues+2} = sqrt(Pext{InteractionIds(1)}.*Pext{InteractionIds(2)});

end

function E = local_cue2E(s,cue,Erange,Eoffset)
E = Erange * exp(-cue/s) + Eoffset;
end

function y = local_nanrms(x)
id = not(isnan(x));
y = rms(x(id));
end


