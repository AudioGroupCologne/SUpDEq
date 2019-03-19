function varargout = data_baumgartner2017looming( varargin )
%data_baumgartner2017looming - Results from Baumgartner et al. (2017)
%
%   Usage: data = data_baumgartner2017looming(dataFlag,measure)
%          data = data_baumgartner2017looming(fig)
%
%   DATA_BAUMGARTNER2017LOOMING provides individually measured HRTFs and
%   experimental results from Baumgartner et al. (2017). Use the fig
%   flag to obtain data shown in figures from Baumgartner et al. (2017).
%
%   The dataFlag flag may be used to choose between HRTFs and various
%   experimental results:
%
%     'hrtf'   HRTFs used in all experiments.
%     'pretest' Behavioral results of pre-test.
%     'exp1'    Behavioral or ERP results of Exp. I. Default.
%     'exp2'    Behavioral results of Exp. II.
%
%   The measure flag may be one of:
%
%     'judgment'  Judgment of relative distance change (motion direction).
%                 Default.
%     'rt'        Response time.
%     'erp'       ERP magnitude measures.
%
%   The fig flag may be one of:
%
%     'fig1b'     Effect of spectral contrast manipulation according to
%                 factor C on magnitude responses of listener-specific
%                 stimuli of Exp. I (B, Top) as well as their frequency-specific
%                 and overall loudness changes relative to C = 1. Shaded
%                 areas denote ±1 standard error of the mean (SEM; N = 15).
%                 Note that changes in overall loudness oppose the intended
%                 effect of contrast switch.
%     'fig2'      Behavioral responses were more consistent for sounds
%                 perceived as approaching compared to sounds perceived as
%                 receding if instantaneous spectral changes were presented
%                 in continuous stimuli. Mean behavioral responses in the
%                 3?AFC motion discrimination task of Exp. I (fist figure; N = 15)
%                 and the 2?AFC motion discrimination task of Exp. II
%                 (second figure; N = 13). Results of Exp. II are separated
%                 between trials presenting instantaneous but continuous
%                 (Left; as in Exp. I) and discontinuous (Right) spectral
%                 contrast switches. Decreasing spectral contrast switches
%                 (orange lines) were predominantly perceived as approaching
%                 (orange triangles), increasing spectral contrast switches
%                 (green lines) as receding (green triangles), and constant
%                 spectral contrasts (no lines) as static (gray squares).
%                 Statistical analyses focused on these predominant response
%                 associations. Values reflect mean ±1 SEM.
%     'fig3a'     Extracted N1 and P2 amplitudes evoked by stimulus onset.
%                 Error bars reflect SEM.
%     'fig3b'     Extracted N1 and P2 amplitudes evoked by stimulus switch.
%                 Error bars reflect SEM.
%
%   Additional flags may be:
%
%     'plot'    Plot results as published.
%     'noplot'  No plots. Default.
%     'stat'   Analyze and display inferential statistics.
%     'nostat' No statistics. Default.
%     'onset'   To use onset ERPs.
%     'switch'  To use switch-ERPs. Default.
%
%   Output parameters:
%     data    : structure that contains either HRTFs (id, and Obj) or
%               experimental results including raw and averaged data
%               (rawData,data, opional meta information).
%               Statisitcs (stat) and figure handles (fig) are
%               provided if requested.
%
%   Requirements:
%   -------------
%
%   1) SOFA API v0.4.3 or higher from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
%
%   2) Data in hrtf/baumgartner2017looming and auxdata/baumgartner2017looming (downloaded on the fly)
%
%   3) Statistics Toolbox for Matlab (for some fig)
%
%   Examples:
%   ---------
%
%   To display results of Fig.1B :
%
%     data_baumgartner2017looming('fig1b','plot');
%
%   To display results of Fig.2 :
%
%     data_baumgartner2017looming('fig2','plot');
%
%   To display results of Fig.3A :
%
%     data_baumgartner2017looming('fig3a','plot');
%
%   To display results of Fig.3B :
%
%     data_baumgartner2017looming('fig3b','plot');
%
%   References:
%     R. Baumgartner, D. K. Reed, B. TÃ³th, V. Best, P. Majdak, H. S.
%     Colburn, and B. Shinn-Cunningham. Asymmetries in behavioral and neural
%     responses to spectral cues demonstrate the generality of auditory
%     looming bias. Proceedings of the National Academy of Sciences, 2017.
%     [1]arXiv | [2]http ]
%     
%     References
%     
%     1. http://arxiv.org/abs/http://www.pnas.org/content/early/2017/08/16/1703247114.full.pdf
%     2. http://www.pnas.org/content/early/2017/08/16/1703247114.abstract
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_baumgartner2017looming.php

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

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

% definput.import={'amt_cache'};
definput.flags.expirement = {'exp1','exp2','pretest','hrtf'};
definput.flags.measure = {'judgment','rt','erp'};
definput.flags.erp = {'switch','onset'};
definput.flags.plot = {'noplot','plot'};
definput.flags.stat = {'nostat','stat'};
definput.flags.fig = {'nofig','fig1b','fig2','fig3a','fig3b','fig3c'};

[flags]=ltfatarghelper({},definput,lower(varargin));

if flags.do_nofig
  %% HRTFs
  if flags.do_hrtf
    id = amt_load('baumgartner2017looming','subjects.mat');
    out = struct('id',{},'Obj',{});
    for ii = 1:length(id.subjects)
      out(ii).id = id.subjects{ii};
      out(ii).Obj = SOFAload(fullfile(SOFAdbPath,'baumgartner2017looming',...
        [id.subjects{ii},'_eq.sofa']));
    end
    varargout{1} = out;
    return
  end

  %% Behavioral Results
  if flags.do_judgment || flags.do_rt

    raw=amt_load('baumgartner2017looming',[flags.expirement,'.mat']);
    Nsubj = length(raw.data);

    ISI = unique(raw.data(1).ISI);
    nISI = length(ISI);
    ISIlbl = {'continuous','discont.'};

    RTprctile = 25;
    conditions = {[0,1];[0,0.5];[0.5,1];[1,0];[0.5,0];[1,0.5];[0,0];[0.5,0.5];[1,1]};
    condLabelPlot = { '0\leftrightarrow1';'0\leftrightarrow.5';'.5\leftrightarrow1';...
                  'C_1=C_2'};
    response = {'receding','approaching','constant'};
    
    %% Evaluate Data
    pResp = nan(Nsubj,length(conditions),length(response),nISI);
    RTall = pResp;
    for ss = 1:Nsubj
      E = raw.data(ss).judgment;
      C12 = raw.data(ss).C12;
      sISI = raw.data(ss).ISI;
      RT = raw.data(ss).RT;
      for iISI = 1:nISI
        for cc = 1:length(conditions)
          idc = C12(:,1) == conditions{cc}(1) & C12(:,2) == conditions{cc}(2) & sISI(:) == ISI(iISI);
          N = sum(idc);

          Irecede = E(idc) > 0;
          pResp(ss,cc,1,iISI) = sum(Irecede) / N;
          RTall(ss,cc,1,iISI) = percentile(RT(Irecede),RTprctile);

          Iappr = E(idc) < 0;
          pResp(ss,cc,2,iISI) = sum(Iappr) / N;
          RTall(ss,cc,2,iISI) = percentile(RT(Iappr),RTprctile);

          Iconst = E(idc) == 0;
          pResp(ss,cc,3,iISI) = sum(Iconst) / N;
          RTall(ss,cc,3,iISI) = percentile(RT(Iconst),RTprctile);
        end
      end
    end


    if flags.do_judgment
      meas = 100*pResp; % responses in percent
      YLabel = 'Response (%)';
      if flags.do_exp2
        YLim = [43,105];
      else
        YLim = [-5,105];
      end
    elseif flags.do_rt
      meas = 1e3*RTall; % response times in ms
      YLabel = 'Response time (ms)';
      YLim = 1e3*[0.451,1.149];
    end

    % outlier removal
    rawData = struct2cell(raw.data);
    ID = rawData(1,:);
    if flags.do_exp2
      pc = reshape(cat(2,pResp(:,1:3,1,:),pResp(:,4:6,2,:)),[Nsubj,6,nISI]); % percent correct
      bias = squeeze(mean(pc(:,4:6,:) - pc(:,1:3,:),2));
      [~,outlier] = max(bias(:,2));
      amt_disp(['Subject ' raw.data(outlier).ID ' identified as outlier and removed from further analyses.'])
      iNew = (1:Nsubj)~=outlier;
      ID = ID(iNew);
      meas = meas(iNew,:,:,:);
      Nsubj = Nsubj-1;
    end

    % Average constant conditions
    meas = cat(2,meas(:,1:6,:,:),nan_mean(meas(:,7:9,:,:),2));
    conditions = [conditions(1:6);'constant'];

    % Standard errors of the means
    sem = nan_std(meas)/sqrt(Nsubj);

    % Output
    out.data = meas;
    out.rawData = raw.data;
    out.meta.dim = 'subject_C_response_ISI';
    out.meta.subject = ID;
    out.meta.C = conditions;
    out.meta.response = response;
    out.meta.ISI = ISI;

    if flags.do_stat && ~verLessThan('matlab','8.2')
      % ANOVA
      pc = reshape(cat(2,meas(:,1:3,1,:),meas(:,4:6,2,:)),Nsubj,[]); % percent correct
      DV = array2table(pc);
      contrast = strrep(condLabelPlot(1:3),'\leftrightarrow','-');
      contrast = repmat(contrast,[2*nISI,1]);
      direction = cell(6,1);
      direction(1:3) = {'receding'};
      direction(4:6) = {'approaching'};
      direction = repmat(direction,[nISI,1]);
      idISI = ceil((1:6*nISI)/6);
      ISInom = nominal(ISI(idISI))';
      IVs = table(contrast,direction,ISInom);
      rm = fitrm(DV,['pc1-pc',num2str(length(contrast)),' ~ 1'],'WithinDesign',IVs);
      if length(ISI) == 1
        [ranovaResult,~,C,~] = ranova(rm,'WithinModel','direction*contrast');
      else
        [ranovaResult,~,C,~] = ranova(rm,'WithinModel','direction*contrast*ISInom');
      end

      ranovaResult.Properties.RowNames = strrep(ranovaResult.Properties.RowNames,'(Intercept):','');

      % Sphericity corrections
      spherCorr = epsilon(rm,C);

      % Add corrected DFs to ranova table
      idrep = round(0.5:0.5:length(spherCorr.GreenhouseGeisser)); % repeat iteratively
      ranovaResult.DFGG = ranovaResult.DF .* ...
        reshape(spherCorr.GreenhouseGeisser(idrep),size(ranovaResult.DF));

      % Add effect sizes to ranova table
      SSeffect = ranovaResult.SumSq(1:2:end);
      SSerror = ranovaResult.SumSq(2:2:end);
      eta_pSq = nan(2*length(SSerror),1);
      eta_pSq(1:2:end) = SSeffect./(SSeffect+SSerror); % effect size per (eta_partial)^2
      ranovaResult.eta_pSq = eta_pSq;

      amt_disp(ranovaResult(:,[4,6,9,10]))

      mc = multcompare(rm,'contrast');
      amt_disp(mc)

      out.stat = ranovaResult;
    end

    %% Plot
    if flags.do_plot
      out.fig = figure;
      x = [1:3,1:3,4];
      XTick = [x(1:3),x(end)];
      if flags.do_exp2
        x = x(1:6);
        XTick = x(1:3);
        response = response(1:2);
      end
      XLim = [XTick(1)-.6,XTick(end)+.6];
      dx2 = .05*[-1,0,1]; % between responses
      symb = '^vs';
      lineStyle = {':';'-';'-.'}; % increase, decrease
      colorR = [5,120,100;250,30,0;149,123,109]/255;
      colorC = [5,120,100;250,30,0;255,255,255]/255;
      for iisi = 1:nISI
        subplot(1,nISI,iisi)
        for rr = 1:length(response)
          y = nan_mean(meas(:,:,rr,iisi));
          l = sem(1,:,rr,iisi);
          u = sem(1,:,rr,iisi);
          idx = {1:3;4:6;7};
          if flags.do_exp2
            idx = idx(1:2);
          end
          for ii = 1:length(idx)
            hC(rr,ii) = plot(x(idx{ii})+dx2(rr),y(idx{ii}),lineStyle{ii});
            hold on
            hR(rr,ii) = errorbar(x(idx{ii})+dx2(rr),y(idx{ii}),l(idx{ii}),u(idx{ii}),symb(rr));
            set(hC(rr,ii),'Color',colorC(ii,:))
          end
          set(hR(rr,:),'MarkerFaceColor',colorR(rr,:),'Color',colorR(rr,:))
        end
        set(hR(1:2:size(hR,1),:),'MarkerFaceColor','w')
        if flags.do_exp2
          set([hC(2:3),hR(2:3)],'Visible','off')
        end

        set(gca,'XTick',XTick,'XTickLabel',condLabelPlot(1:length(XTick)))
        axis([XLim,YLim])

        ylabel(YLabel)
        xlabel('Spectral contrast pair')
        if flags.do_exp2
          title(['Exp. II: ',ISIlbl{iisi}])
        end
      end
      stimulus = {'C increase','C decrease','C constat'};
      legend([hC(1,:)';hR(:,1)],[stimulus(1:size(hC,2)),response],'Location','eastoutside')
    end
  end

  %% Physiological Results
  if flags.do_erp

    if not(flags.do_exp1)
      error('ERPs only available for exp1.')
    end

    if flags.do_onset % onset
      erp = amt_load('baumgartner2017looming','onsetERP.mat');
    else % flags.do_switch
      erp = amt_load('baumgartner2017looming','switchERP.mat');
    end
    erp.compLbl = {'N1','P2'};
    out.rawData = erp;

    if flags.do_stat
      for rr = 1:length(erp.compLbl)
        amt_disp(erp.compLbl{rr})
        amt_disp(erp.compStats{rr}.ranova(:,[4,6,9,10]))
        if isfield(erp.compStats{rr},'posthoc')
          if isfield(erp.compStats{rr}.posthoc,'combination')
            amt_disp(erp.compStats{rr}.posthoc.combination)
          else
            amt_disp(erp.compStats{rr}.posthoc)
          end
        end
      end
    end

    if flags.do_plot

      if flags.do_onset % onset

        condLabelData = {'0','0.5','1'};
        condLabel = condLabelData;
        legLbl = erp.compLbl;
        XLabel = 'Spectral contrast';
        YLim = [-3.8,4.9];

        condOrder = [1,3,2]; % to reorder erp.condLbl acc. to condLabel
        resp = permute(erp.compAmp(condOrder,:,:),[2,1,3]); % subjects in first dimension
        idx = {1:3};
        dx = 0;

      else % flags.do_switch

        condLabel = { '0\leftrightarrow1';'0\leftrightarrow.5';'.5\leftrightarrow1'};
        legLbl = {[erp.compLbl{1},', C increase'],[erp.compLbl{1},', C decrease'],...
                  [erp.compLbl{2},', C increase'],[erp.compLbl{2},', C decrease']};
        XLabel = 'Spectral contrast pair';
        YLim = [-3.5,3.5];

        condOrder = [5,6,4,2,3,1]; % to reorder erp.condLbl acc. to condLabel
        resp = permute(erp.compAmp(condOrder,:,:),[2,1,3]);

        idx = {1:3;4:6};
        dx = .1*[-1,1];

      end

      % Standard errors
      seResp = std(resp)/sqrt(size(erp.compAmp,2));

      out.fig = figure;
      x = 1:3;
      symb = 'oo';
      lineStyle = {':','-'};
      color = 0.8*[.65,.35,1;1,.5,0];
      for rr = 1:length(erp.compLbl)
        y = mean(resp(:,:,rr));
        l = seResp(1,:,rr);
        u = seResp(1,:,rr);
        for ii = 1:length(idx)
          h(rr,ii) = errorbar(x+dx(ii),y(idx{ii}),l(idx{ii}),u(idx{ii}),...
            [symb(rr),lineStyle{ii}]);
          hold on
        end
        set(h(rr,1),'MarkerFaceColor','w','Color',color(rr,:))
        set(h(rr,length(idx)),'MarkerFaceColor',color(rr,:),'Color',color(rr,:))

        set(gca,'XTick',x,'XLim',x([1,3])+[-1,1])
        set(gca,'XTickLabel',condLabel)
        ylabel('Cz potential (uV)')
        xlabel(XLabel)
      end
      set(gca,'YLim',YLim,'YMinorTick','on')
      legend(legLbl,'Location','eastoutside')
    end
  end
end

%% Fig. 1B
if flags.do_fig1b

  stim = sig_baumgartner2017looming( 'exp1');

  %% Top panel: Transfer characterisitcs
  fs = stim(1).fs;
  Nfft = 2^10;
  freq = 0:fs/Nfft:fs/2; % frequency vector
  mag = nan(length(freq),length(stim(1).C_IR),2,length(stim));
  for ss = 1:length(stim)
    if stim(ss).azi == -90
      ipsiContra = [2,1];
    else
      ipsiContra = [1,2];
    end
    for cc = 1:length(stim(1).C_IR)
      Sig = stim(ss).IR{cc};
      for ich = 1:2
        ch = ipsiContra(ich);
        mag(:,cc,ich,ss) = db(abs(fftreal(Sig(:,ch),Nfft)));
      end
    end
  end
  mag = mag-3; % arbitrary adjustment to set ipsi. C=0 at 0 dB

  if flags.do_plot
    XLim = [800,17e3];
    YLim = [-34,12];
    blue = [0,0,0.7];
    green = [0,0.7,0];
    red = [0.7,0,0];
    color = {blue,1.4*blue;green,1.4*green;red,1.4*red};
    lineStyle = {'-','--'};
    out.fig(1) = figure;
    ii = 1;
    for cc = 1:3
      for ch = 1:2
        lMEAN = mean(mag(:,cc,ch,:),4);
        lSEM = std(mag(:,cc,ch,:),0,4)/sqrt(15);
        h(ii) = shadedErrorBar(freq,lMEAN,lSEM,{'LineStyle',lineStyle{ch},'Color',color{cc,ch}},1);
        hold on
        ii = ii+1;
      end
    end
    set(gca,'XScale','log','XLim',XLim,'YLim',YLim)
    leg = legend([h.mainLine],'C=0 (ipsi)','C=0 (contra)','C=0.5 (ipsi)','C=0.5 (contra)','C=1 (ipsi)','C=1 (contra)');
    set(leg,'Location','eastoutside','box','off')
    ylabel('Magnitude (dB)')
    xlabel('Frequency (Hz)')
  end

  out.magnitudeResponse.data = permute(mag,[1,3,4,2]);
  out.magnitudeResponse.meta.dim = 'freq_channel_subject_C';
  out.magnitudeResponse.meta.freq = freq;
  out.magnitudeResponse.meta.channel = {'ipsi','contra'};
  out.magnitudeResponse.meta.subject = cat(1,stim.ID);
  out.magnitudeResponse.meta.C = 0:0.5:1;

  %% Bottom panels: Loudness predictions
  % The following loudness predictions were performed with the
  % LoudnessToolbox 1.2 provided by Genesis (http://genesis-acoustics.com);
  % models used:
  %   M1: Loudness_ISO532B_from_sound (calculated but not shown)
  %   M2: Loudness_ANSI_S34_2007 (used for publication);
  % Data dimensions:
  %   subject (1:15) x C (0:.5:1) x model (M1,M2) [x channel (ipsi,contra)];
  %   frequency (M1:BarkAxis; M2:fc) x channel (ipsi,contra)
  L = amt_load('baumgartner2017looming','specificLoudness.mat');

  % Difference to reference
  dLL_specif = cell(3,1);
  for ss = 1:length(stim)
    LLC1 = L.loudnessLevel_specif{ss,3,2};
    for m = 1:3
      dLL_specif{m}(:,:,ss) =  L.loudnessLevel_specif{ss,m,2} - LLC1;
    end
  end
  dLoudnessLevel = L.loudnessLevel(:,1:2,:,:) - repmat(L.loudnessLevel(:,3,:,:),[1,2,1,1]);

  if flags.do_plot
    out.fig(2) = figure;
    YLim = [-9,13];
    ii = 1;
    for m = 1:3
      for ch = 1:2
        lMEAN = mean(dLL_specif{m}(:,ch,:),3);
        lSEM = std(dLL_specif{m}(:,ch,:),0,3)/sqrt(15);
        h(ii) = shadedErrorBar(L.fc,lMEAN,lSEM,{'LineStyle',lineStyle{ch},'Color',color{m,ch}},1);
        hold on
        ii = ii+1;
      end
    end
    set(gca,'XScale','log','XLim',XLim,'YLim',YLim)
    leg = legend([h.mainLine],'C=0 (ipsi)','C=0 (contra)','C=.5 (ipsi)','C=.5 (contra)');
    set(leg,'Location','northwest','box','off')
    ylabel('Loudness level difference to C=1 (phon)')
    xlabel('Frequency (Hz)')
  end

  out.specificLoudnessLevelDiff.data = cat(4,dLL_specif{:});
  out.specificLoudnessLevelDiff.meta.dim = 'freq_channel_subject_C';
  out.specificLoudnessLevelDiff.meta.freq = L.fc;
  out.specificLoudnessLevelDiff.meta.channel = {'ipsi','contra'};
  out.specificLoudnessLevelDiff.meta.subject = cat(1,stim.ID);
  out.specificLoudnessLevelDiff.meta.C = 0:0.5:1;

  % overall loudness level
  dLoudnessLevelP = dLoudnessLevel(:,:,2,:);
  if flags.do_plot
    try
        out.fig(3) = figure;
        boxplot(dLoudnessLevelP(:,:),... %,{{'M1';'M1';'M2';'M2'},[0;0.5;0;.5]}
          'Factorgap',10,'FactorSeparator',1,'Whisker',Inf,...
          'Colors',cat(1,color{[1,2,4,5]}))
        set(gca,'YLim',YLim)
        ylabel('Loudness level difference to C=1 (phon)')
        xlabel('Spectral contrast (C)')
    catch
    	warning('Statistics Toolbox not available, omitting figure 3.')
    end
    
  end

  out.loudnessLevelDiff.data = permute(dLoudnessLevelP,[4,1,2,3]);
  out.loudnessLevelDiff.meta.dim = 'channel_subject_C';
  out.loudnessLevelDiff.meta.channel = {'ipsi','contra'};
  out.loudnessLevelDiff.meta.subject = cat(1,stim.ID);
  out.loudnessLevelDiff.meta.C = [0,0.5];

end

%% Fig. 2
if flags.do_fig2
  amt_disp('Exp. I:')
  out.exp1 = data_baumgartner2017looming('exp1',flags.plot,'stat');
  title('Exp. I')
  amt_disp('Exp. II:')
  out.exp2 = data_baumgartner2017looming('exp2',flags.plot,'stat');
  legend off
end

%% Fig. 3A
if flags.do_fig3a
  out = data_baumgartner2017looming('erp','onset',flags.plot,flags.stat);
end

%% Fig. 3B
if flags.do_fig3b
  out = data_baumgartner2017looming('erp','switch',flags.plot,flags.stat);
end

%% Fig. 3C
if flags.do_fig3c
  out = amt_load('baumgartner2017looming','ERPclusterAnalysis.mat');
end

%% Output
if nargout == 1
  varargout{1} = out;
end

end

%% Internal plotting functions
function varargout=shadedErrorBar(x,y,errBar,lineProps,transparent)
% function H=shadedErrorBar(x,y,errBar,lineProps,transparent)
%
% Purpose
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar color is chosen automatically.
%
% Inputs
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of n observations by m cases
%     where m has length(x);
% errBar - if a vector we draw symmetric errorbars. If it has a size
%          of [2,length(x)] then we draw asymmetric error bars with
%          row 1 being the upper bar and row 2 being the lower bar
%          (with respect to y). ** alternatively ** errBar can be a
%          cellArray of two function handles. The first defines which
%          statistic the line should be and the second defines the
%          error bar.
% lineProps - [optional,'-k' by default] defines the properties of
%             the data line. e.g.:
%             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
% transparent - [optional, 0 by default] if ==1 the shaded error
%               bar is made transparent, which forces the renderer
%               to be openGl. However, if this is saved as .eps the
%               resulting file will contain a raster not a vector
%               image.
%
% Outputs
% H - a structure of handles to the generated plot objects.
%
%
% Examples
% y=randn(30,80); x=1:size(y,2);
% shadedErrorBar(x,mean(y,1),std(y),'g');
% shadedErrorBar(x,y,{@median,@std},{'r-o','markerfacecolor','r'});
% shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'});
%
% Overlay two transparent lines
% y=randn(30,80)*10; x=(1:size(y,2))-40;
% shadedErrorBar(x,y,{@mean,@std},'-r',1);
% hold on
% y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
% shadedErrorBar(x,y,{@mean,@std},'-b',1);
% hold off
%
%
% Rob Campbell - November 2009



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error checking
narginchk(3,5)


%Process y using function handles if needed to make the error bar
%dynamically
if iscell(errBar)
    fun1=errBar{1};
    fun2=errBar{2};
    errBar=fun2(y);
    y=fun1(y);
else
    y=y(:)';
end

if isempty(x)
    x=1:length(y);
else
    x=x(:)';
end


%Make upper and lower error bars if only one was specified
if length(errBar)==length(errBar(:))
    errBar=repmat(errBar(:)',2,1);
else
    s=size(errBar);
    f=find(s==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

if length(x) ~= length(errBar)
    error('length(x) must equal length(errBar)')
end

%Set default options
defaultProps={'-k'};
if nargin<4, lineProps=defaultProps; end
if isempty(lineProps), lineProps=defaultProps; end
if ~iscell(lineProps), lineProps={lineProps}; end

if nargin<5, transparent=0; end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot to get the parameters of the line
H.mainLine=plot(x,y,lineProps{:});


% Work out the color of the shaded region and associated lines
% Using alpha requires the render to be openGL and so you can't
% save a vector image. On the other hand, you need alpha if you're
% overlaying lines. There we have the option of choosing alpha or a
% de-saturated solid colour for the patch surface .

col=get(H.mainLine,'color');
edgeColor=col+(1-col)*0.55;
patchSaturation=0.15; %How de-saturated or transparent to make patch
if transparent
    faceAlpha=patchSaturation;
    patchColor=col;
    set(gcf,'renderer','openGL')
else
    faceAlpha=1;
    patchColor=col+(1-col)*(1-patchSaturation);
    set(gcf,'renderer','painters')
end


%Calculate the error bars
uE=y+errBar(1,:);
lE=y-errBar(2,:);


%Add the patch error bar
holdStatus=ishold;
if ~holdStatus, hold on,  end


%Make the patch
yP=[lE,fliplr(uE)];
xP=[x,fliplr(x)];

%remove nans otherwise patch won't work
xP(isnan(yP))=[];
yP(isnan(yP))=[];


H.patch=patch(xP,yP,1,'facecolor',patchColor,...
              'edgecolor','none',...
              'facealpha',faceAlpha);


%Make pretty edges around the patch.
H.edge(1)=plot(x,lE,'-','color',edgeColor);
H.edge(2)=plot(x,uE,'-','color',edgeColor);

%Now replace the line (this avoids having to bugger about with z coordinates)
uistack(H.mainLine,'top')


if ~holdStatus, hold off, end


if nargout==1
    varargout{1}=H;
end
end

function prc = percentile(x,k)
    % percentile function to replace prctile in statistics toolbox
    % x .. data vector
    % k .. percentage in % (k >= 1)
    % if k is outside the range the min or max value of x gets assigned
    
    len = length(x);
    if isempty(x)
        prc = NaN;
        return
    end
    if len == 1
        prc = x; return
    end
    
    y = sort(x);
    z = 100*(0.5:1:(len-0.5))/len;
    
    if k<z(1)
        prc=k(1); return
    end
    if k>z(end)
        prc=z(end); return
    end
    
    if isempty(find(z==k, 1))
        prc = interp1(z,y,k);
    else
        prc = y(find(z==k, 1));
    end
end

function y = nan_mean(x,dim)
% FORMAT: Y = NANMEAN(X,DIM)
% 
%    Average or mean value ignoring NaNs
%
%    This function enhances the functionality of NANMEAN as distributed in
%    the MATLAB Statistics Toolbox and is meant as a replacement (hence the
%    identical name).  
%
%    NANMEAN(X,DIM) calculates the mean along any dimension of the N-D
%    array X ignoring NaNs.  If DIM is omitted NANMEAN averages along the
%    first non-singleton dimension of X.
%
%    Similar replacements exist for NANSTD, NANMEDIAN, NANMIN, NANMAX, and
%    NANSUM which are all part of the NaN-suite.
%
%    See also MEAN

    if isempty(x)
    	y = NaN;
    	return
    end

    if nargin < 2
        dim = min(find(size(x)~=1));
        if isempty(dim)
            dim = 1;
        end
    end

    % Replace NaNs with zeros.
    nans = isnan(x);
    x(isnan(x)) = 0; 

    % denominator
    count = size(x,dim) - sum(nans,dim);

    % Protect against a all NaNs in one dimension
    i = find(count==0);
    count(i) = ones(size(i));
    y = sum(x,dim)./count;
    y(i) = i + NaN;
end

function y = nan_std(x,dim,flag)
% FORMAT: Y = NANSTD(X,DIM,FLAG)
% 
%    Standard deviation ignoring NaNs
%
%    This function enhances the functionality of NANSTD as distributed in
%    the MATLAB Statistics Toolbox and is meant as a replacement (hence the
%    identical name).  
%
%    NANSTD(X,DIM) calculates the standard deviation along any dimension of
%    the N-D array X ignoring NaNs.  
%
%    NANSTD(X,DIM,0) normalizes by (N-1) where N is SIZE(X,DIM).  This make
%    NANSTD(X,DIM).^2 the best unbiased estimate of the variance if X is
%    a sample of a normal distribution. If omitted FLAG is set to zero.
%    
%    NANSTD(X,DIM,1) normalizes by N and produces the square root of the
%    second moment of the sample about the mean.
%
%    If DIM is omitted NANSTD calculates the standard deviation along first
%    non-singleton dimension of X.
%
%    Similar replacements exist for NANMEAN, NANMEDIAN, NANMIN, NANMAX, and
%    NANSUM which are all part of the NaN-suite.
%
%    See also STD


    if isempty(x)
        y = NaN;
        return
    end

    if nargin < 3
        flag = 0;
    end

    if nargin < 2
        dim = min(find(size(x)~=1));
    	if isempty(dim)
            dim = 1; 
        end	  
    end


    % Find NaNs in x and nanmean(x)
    nans = isnan(x);
    avg = nan_mean(x,dim);

    % create array indicating number of element 
    % of x in dimension DIM (needed for subtraction of mean)
    tile = ones(1,max(ndims(x),dim));
    tile(dim) = size(x,dim);

    % remove mean
    x = x - repmat(avg,tile);

    count = size(x,dim) - sum(nans,dim);

    % Replace NaNs with zeros.
    x(isnan(x)) = 0; 


    % Protect against a  all NaNs in one dimension
    i = find(count==0);

    if flag == 0
    	y = sqrt(sum(x.*x,dim)./max(count-1,1));
    else
    	y = sqrt(sum(x.*x,dim)./max(count,1));
    end
    y(i) = i + NaN;
end
