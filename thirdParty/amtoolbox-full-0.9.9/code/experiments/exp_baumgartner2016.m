function varargout = exp_baumgartner2016(varargin)
%EXP_BAUMGARTNER2016 evaluation of baumgartner2016 model
%   Usage: data = exp_baumgartner2016(flag)
%
%   EXP_BAUMGARTNER2016(flag) reproduces figures of the study from 
%   Baumgartner et al. (2014).
%
%   The following flags can be specified
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%     'fig2' or 'ratelevelcurves'
%               Rate-level curves of the three different fiber types
%               represented in the auditory-periphery model. Firing rates were
%               evaluated at a CF of 4 kHz in response to Gaussian white noise at
%               various SPLs. Note that high-SR fibers saturate already at low SPLs,
%               medium-SR fibers at moderate SPLs, and low-SR fibers not at all.
%
%     'fig3' or 'baseline'
%               Correspondence between actual and predicted
%               baseline performance for the 23 normal-hearing listeners after
%               listener-specific calibration of the models sensitivity parameter (S).
%
%     'fig4' or 'tab1' or 'hearingthreshold'
%               Hearing thresholds estimated for simulated OHC gains
%               (COHC) within the range of 1 (normal active cochlea) to 0 (passive
%               cochlea). The selected set of OHC gains results in approximately
%               equal increments of high-frequency thresholds.
%               Table I: Simulated Conditions of OHC Dysfunction,  
%               Estimated PTAs, and Corresponding Hearing Loss Categories.
%
%     'fig5' or 'numchan' and 'spatstrat'
%               Model evaluation for normal-hearing listeners tested
%               on the effects of spectral resolution (by number of vocoder
%               channels in Goupell et al., 2010) and spectral warping (Majdak,
%               Walder, et al., 2013). Model data (filled circles) are compared
%               with actual data (open circles) from the two studies. Error bars
%               represent SDs. Symbols are slightly shifted along the abscissa for
%               better visibility. BB = broadband noise burst; CL = broadband click
%               train (infinite number of channels); LP = low-pass filtered at
%               8.5 kHz; W = HRTFs spectrally warped from 2.8 to 16 kHz to
%               2.8 to 8.5 kHz.
%
%     'fig6' or 'evalSPLtem'
%               Effect of template SPL on predictive power of the
%               model for the two studies (Goupell et al., 2010; Majdak, Walder,
%               et al., 2013) shown in Figure 5. Predictions based on a single
%               template SPL equivalent to the actual SPL of the target sounds of
%               60 dB result in similar prediction residues as based on templates
%               mixed across a broad range of SPLs. Higher plausibility of the
%               mixed-SPL templates was the reason to choose this representation
%               for all further simulations (including predictions shown in Figure 5).
%
%     'fig7' or 'impairment'
%               Effects of OHC dysfunctions and selective activity of
%               AN fibers on predicted quadrant error rates (top) and local RMS
%               errors (bottom). Thick bar: interquartile range (IQR). Thin
%               bar: data range within 1.5 IQR. Horizontal line within thick
%               bar: average. Dashed horizontal line: chance performance.
%
%     'fig8' or 'sensitivity'
%               Sensitivity (dprime) of AN fibers in level discrimination
%               as function of SPL predicted for different fiber types and OHC
%               dysfunctions. Sensitivities were evaluated for SPL increments of
%               10 dB and averaged across 28 CFs from 0.7 to 18 kHz. Gray area:
%               stimulus range of target sounds at 60 dB SPL.
%
%     'fig9' or 'effectOnCues'
%               Effect of OHC dysfunction on positive spectral gradients.  
%               Exemplary median-plane HRTFs from one listener (NH46). Note the
%               distinct direction-specific patterns for the normal and moderate 
%               OHC dysfunctions (COHC>=0.4), which are almost absent in the 
%               cases of the severe and complete OHC dysfunctions (COHC<=0.1).
%
%     'baseline_ex' 
%               Prediction examples. Actual responses and response predictions 
%               for three exemplary listeners when listening to median-plane 
%               targets in the baseline condition. 
%               Actual response angles are shown as open circles. 
%               Probabilistic response predictions are encoded by brightness 
%               according to the color bar to the right. Actual (A:) and 
%               predicted (P:) quadrant error rates (QE) and local polar 
%               RMS errors (PE) are listed above each panel.
%
%     'spatstrat_ex' 
%               Effect of band limitation and spectral warping. Actual 
%               responses and response predictions for listener NH12 when 
%               listening to broadband (BB), low-pass filtered (LP), or 
%               spectrally warped (W) DTFs of the median plane. Data were 
%               pooled within pm15^circ of lateral angle.
%
%     'numchan_ex'
%               Effect of spectral resolution in terms of varying the number 
%               of spectral channels of a channel vocoder. Actual responses 
%               and response predictions for exemplary listener NH12. 
%               Results for 24, 9, and 3 channels are shown. All other 
%               conventions are as in Fig.3.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v0.4.3 or higher from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2014
%
%   3) Statistics Toolbox for Matlab (for some of the figures)
%
%   Examples:
%   ---------
%
%   To display the rate-level curves use :
%
%     exp_baumgartner2016('fig2');
%
%   To display the baseline prediction use :
%
%     exp_baumgartner2016('fig3');
%
%   To display estimations of hearing thresholds use :
%
%     exp_baumgartner2016('fig4');
%
%   To display model evaluation for normal-hearing listeners use :
%
%     exp_baumgartner2016('fig5');
%
%   To display evaluation results for template SPL use :
%
%     exp_baumgartner2016('fig6');
%
%   To display predicted effects of sensorineural hearing loss use (requires Matlab 2013b or higher) :
%
%     exp_baumgartner2016('fig7');
%
%   To display sensitivity evaluation for different fiber types use :
%
%     exp_baumgartner2016('fig8');
%
%   To display effect of OHC damage on exemplary spectral cue representation use :
%
%     exp_baumgartner2016('fig9');
%
%   See also: baumgartner2016 data_baumgartner2016
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling the effects of
%     sensorineural hearing loss on auditory localization in the median
%     plane. Trends in Hearing, 20:1-11, 2016.
%     
%     P. Majdak, T. Walder, and B. Laback. Effect of long-term training on
%     sound localization performance with spectrally warped and band-limited
%     head-related transfer functions. J. Acoust. Soc. Am., 134:2148-2159,
%     2013.
%     
%     M. J. Goupell, P. Majdak, and B. Laback. Median-plane sound
%     localization as a function of the number of spectral channels using a
%     channel vocoder. J. Acoust. Soc. Am., 127:990-1001, 2010.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/experiments/exp_baumgartner2016.php

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


% AUTHOR: Robert Baumgartner


%% ------ Check input options --------------------------------------------

definput.import={'amt_cache','localizationerror','baumgartner2014_pmv2ppp'};

definput.flags.experiment = {'missingflag',...
  'fig2','fig3','fig4','fig5','fig6','fig7','fig8','fig9','tab1',...
   'sensitivity','baseline','baseline_ex',...
   'baseline_lat','parametrization',...
   'spatstrat_ex','spatstrat','numchan_ex','numchan',...
   'evalSPLtem','cOHCvsSens',...
   'sabin2005','impairment','effectOnCues','dynrangecheck',...
   'localevel','hearingthreshold','ratelevelcurves','evalSpectralContrast'};
     
definput.keyvals.ModelSettings = {};
 
% Figure Settings
definput.flags.plot = {'plot','noplot'};
definput.keyvals.FontSize = 10;
definput.keyvals.MarkerSize = 6;
definput.keyvals.gap = 0;
definput.keyvals.marg_h = [.08,.05];
definput.keyvals.marg_w = .05;
definput.keyvals.TickLength = [0.02,0.04];

% For: sabin2005
definput.keyvals.SL2SPL = 10;
definput.flags.gainInDeg = {'','gainInDeg'};
definput.keyvals.sabin2005_Sdivisor = 1;
definput.flags.sabin2005nomrs = {'','nomrs'};

% For: impairment
definput.keyvals.subjects = [];
definput.keyvals.SPLset = 60;
definput.keyvals.FTset = {1:3,1,2,3};
definput.keyvals.cOHCset = [1,0.4,0.1,0];
definput.flags.impairment={'FTlabel','noFTlabel'};
definput.flags.SPLseparation={'splitSPL','compriseSPL'};

% For: localevel
definput.flags.localevelplot = {'performance','pmv'};

% For: parametrization
definput.flags.parametrize = {'gamma','mrsandgamma'};

% For: dynrangecheck
definput.flags.dynrangecheck = {'ratelevel','dynrangeDiff','dprime'};
definput.flags.dynrangecheck_comb = {'separate','combined'};

definput.flags.effectOnCues={'OHC','FT'};

% General: availability of Statistics Toolbox
definput.flags.statistics = {'stat','nostat'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

model.definput.import={'baumgartner2016'};
[model.flags,model.kv] = ltfatarghelper({},model.definput,kv.ModelSettings);

errorflag = [flags.errorflag,flags.ppp];

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.experiment{2:end-2}),...
             sprintf('%s or %s',definput.flags.experiment{end-1},definput.flags.experiment{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% Define cache name according to settings for auditory periphery model

cachename = ['g' num2str(model.kv.gamma,'%u') ...
      '_mrs' num2str(model.kv.mrsmsp,'%u') ...
      '_do' num2str(model.kv.do,'%u') ...
      '_tem' num2str(model.kv.SPLtem,'%u') 'dB_' model.flags.fbank];
if model.flags.do_gammatone
  cachename = [cachename '_'  num2str(1/model.kv.space,'%u') 'bpERB'];
  if model.flags.do_middleear; cachename = [cachename '_middleear']; end
  if model.flags.do_ihc; cachename = [cachename '_ihc']; end
else % zilany
  cachename = [cachename '_' model.flags.fibertypeseparation];
end
if model.kv.prior > 0 
  cachename = [cachename '_prior' num2str(model.kv.prior,'%u')];
end
if model.kv.tiwin < 0.5
  cachename = [cachename '_tiwin' num2str(model.kv.tiwin*1e3) 'ms']; 
end
cachename = [cachename '_mgs' num2str(model.kv.mgs)]; 


%% Hearing thresholds following OHC dysfunction
if flags.do_hearingthreshold || flags.do_fig4

  cOHC = kv.cOHCset;
  flow = 125;%700; % Hz
  fhigh = 18000; % Hz
  spl = -20:100; % dB
  Ncf = 40; % # CF

  fs = 100e3; % Hz
  t = 0:1/fs:0.1; % s
  cf = audspace(flow,fhigh,Ncf); % CFs under test

  cachename = [model.flags.fbank '_' model.flags.fibertypeseparation];
  cachename = [cachename '_cohc' num2str(cOHC)];
  cachename = strrep(cachename,' ','');cachename = strrep(cachename,'.','p');
  cachename = ['hearingthreshold_' cachename];
  afr = amt_cache('get',cachename,flags.cachemode);
  if isempty(afr)
    afr = nan(length(cf),length(cOHC),length(spl),3);
    for ff = 1:length(cf)
      sig = sin(2*pi*cf(ff)*t);
      for iiOHC = 1:length(cOHC)
        for iispl = 1:length(spl)
          for ft = 1:3
            ANoutTem = zilany2014(spl(iispl),sig,fs,...
                      'flow',cf(ff),'fhigh',cf(ff),...
                      'nfibers',1,'fiberType',ft,...
                      'cohc',cOHC(iiOHC),'cihc',1);
            afr(ff,iiOHC,iispl,ft) = mean(ANoutTem,2);
          end
        end
      end
      amt_disp([num2str(ff) ' of ' num2str(length(cf)) '  done.'],'progress')
    end
  amt_cache('set',cachename,afr)
  end


  ftd = [0.16,0.23,0.61]; % fiber type distribution from Liberman (1978)

  %% Evaluate avg. firing rate at normal hearing threshold (afrLMHc1)
  sizeafr = size(afr);
  afrLMH = sum(afr.*repmat(shiftdim(ftd,-2),[sizeafr(1:3),1]),4);
  afrLMHc1 = nan(length(cf),1);
  HT0 = nan(length(cf),1);
  for ff = 1:length(cf)
    [tmp,id] = min(abs(spl-absolutethreshold(cf(ff),'hda200')));
    HT0(ff) = spl(id);
    afrLMHc1(ff) = afrLMH(ff,1,id);
  end

  %% 
  % ftd(1:2) = 0;
  afr = sum(afr.*repmat(shiftdim(ftd,-2),[sizeafr(1:3),1]),4);
  HT = nan(length(cf),length(cOHC));
  afr0 = nan(length(cf),1);
  for ff = 1:length(cf)
    for iiOHC = 1:length(cOHC)
      id = find(afr(ff,iiOHC,:)>=afrLMHc1(ff),1,'first');
      if isempty(id)
        id = length(spl);
      end
      HT(ff,iiOHC) = spl(id);
    end
  end

  % Pure-tone average HL
  PTAf{1} = [500,1000,2000]; % Rakerd et al. (1998) 
  PTAf{2} = [3150,5000,8000]; % Rakerd et al. (1998) 
  PTAf{3} = [4000,8000,11000]; % Otte et al. (2013) 
  for cc = 1:length(cOHC)
    for ii = 1:length(PTAf)
      for ff = 1:length(PTAf{ii})
        HL(ff,cc,ii) = interp1(cf,HT(:,cc),PTAf{ii}(ff));
      end
    end
  end
  HL = HL(:,2:end,:) - repmat(HL(:,1,:),[1,length(cOHC)-1,1]);
  PTA = shiftdim(mean(HL));
  legendentries = cat(2,repmat('C_{OHC} = ',length(cOHC),1),num2str(cOHC(:),'%2.1f'));
  for cc = 1:length(cOHC)-1
    RowNames{cc} = legendentries(cc+1,:);
  end
  PTA = round(PTA);
  if verLessThan('matlab','8.2'),
    disp('PTAlow_Rakerd98 PTAhigh_Rakerd98 PTAhigh_Otte13');
    PTA
  else
    table(PTA(:,1),PTA(:,2),PTA(:,3),'RowNames',RowNames,'VariableNames',{'PTAlow_Rakerd98','PTAhigh_Rakerd98','PTAhigh_Otte13'})
  end

  varargout{1} = HT;
  varargout{2} = cf;
  
  if flags.do_plot

    figure
    symb = {'-k*','-kd','-k>','-kp'};
    for ii = 1:size(HT,2)
      h(ii) = semilogx(cf,HT(:,ii),symb{ii});
      hold on
    end
    set(h,'MarkerFaceColor','w','MarkerSize',kv.MarkerSize)
    set(gca,'XLim',[cf(1),cf(end)],'YDir','reverse','FontSize',kv.FontSize,'TickLength',kv.TickLength)
    % set(gca,'XTickLabel',{'0.2','0.3','','0.5','','0.7','','','1','2','3','','5','','7','','','10'})
    set(gca,'Layer', 'top')
    xlabel('Frequency (Hz)','FontSize',kv.FontSize)
    ylabel('Hearing threshold (dB SPL)','FontSize',kv.FontSize)

    leg = legend(legendentries,'Location','southwest');
    set(leg,'FontSize',kv.FontSize)

  end

end

%% ------ BASELINE EXAMPLES --------------------------------------------------------
if flags.do_baseline_ex
  
  SL = 50; % presentation level of stimuli
  model.kv.SPL = SL + kv.SL2SPL;
  
  latseg = 0;ii=1;%[-20,0,20]; ii = 2; % centers of lateral segments
%   dlat =  10;         % lateral range (+-) of each segment

  s = data_baumgartner2016('argimport',model.flags,model.kv);
  
%   idselect = ismember({s.id},{'NH15','NH22','NH62','NH12','NH39','NH18'});
  idselect = 19:23;
  s = s(idselect);

  %% LocaMo
  qe = zeros(length(s),length(latseg));
  pe = zeros(length(s),length(latseg));
  for ll = 1:length(s)

      s(ll).sphrtfs{ii} = 0;     % init
      s(ll).p{ii} = 0;        % init

      [s(ll).sphrtfs{ii},polang] = extractsp( latseg(ii),s(ll).Obj );
      [s(ll).p{ii},respangs] = baumgartner2016(...
          s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},'argimport',model.flags,model.kv,...
          'ID',s(ll).id,'fs',s(ll).fs,...
          'mrsmsp',s(ll).mrs,'S',s(ll).S,'lat',latseg(ii),'polsamp',polang,...
          'priordist',s(ll).priordist); 

      [ qe(ll,ii),pe(ll,ii) ] = baumgartner2014_pmv2ppp( ...
          s(ll).p{ii} , polang , respangs , s(ll).target{ii});

      if flags.do_plot
        if ll ==1; figure; end
        subplot(2,3,ll)
        Nmax = min(150,s(ll).Ntar{ii});
        idplot = round(1:s(ll).Ntar{ii}/Nmax:s(ll).Ntar{ii});
        plot_baumgartner2014(s(ll).p{ii},polang,respangs,...
                  s(ll).target{ii}(idplot),s(ll).response{ii}(idplot),...
                  'MarkerSize',kv.MarkerSize,'cmax',0.05,'nocolorbar');
        title({['A: PE = ' num2str(s(ll).pe_exp_lat(ii),2) '\circ, QE = ' num2str(s(ll).qe_exp_lat(ii),2) '%'];['P: PE = ' num2str(pe(ll,ii),2) '\circ, QE = ' num2str(qe(ll,ii),2) '%']},...
          'FontSize',kv.FontSize-1)
        text(90,240,s(ll).id,'FontSize',kv.FontSize,...
          'Color','w','HorizontalAlignment','center')
        xlabel('Target Angle (deg)','FontSize',kv.FontSize)
        ylabel('Response Angle (deg)','FontSize',kv.FontSize)
        set(gca,'FontSize',kv.FontSize-1)
        set(gca,'XTickLabel',{[];[];0;[];60;[];120;[];180;[];[]})
        set(gca,'YTickLabel',{-60;[];0;[];60;[];120;[];180;[];240})
      end

  end
  
  s = rmfield(s,{'Obj','itemlist','sphrtfs'}); % reduce file size 
  
  varargout{1} = s;
  
end

%% ------ PARAMETRIZATION -----------------------------------------------------------
if flags.do_parametrization
  
  if flags.do_mrsandgamma
  
    [gamma,mrs] = amt_cache('get','parametrization', flags.cachemode);
    if isempty(gamma)
      amt_disp('Note that this procedure lasts at least 2 hours!','progress')

      tempfn = fullfile(amt_basepath,'experiments','exp_baumgartner2014_parametrization'); % temporary folder
      mkdir(tempfn)

      s = data_baumgartner2016('argimport',model.flags,model.kv);
      [gamma,mrs] = baumgartner2016_parametrization(s,'SPLtem',kv.SPLtem,...
        flags.fbank,flags.fibertypeseparation,'mgs',kv.mgs);

      amt_cache('set','parametrization',gamma,mrs);
    end
    varargout{1} = {gamma,mrs};
  
  else
    
    gamma = 7:2:30;%[3:2:7,8:13,15:2:30];
    
    dtot = nan(size(gamma));
    for g = 1:length(gamma)
      d = exp_baumgartner2016('argimport',model.flags,model.kv,'baseline','gamma',gamma(g));
      dtot(g) = d.total;
    end
    [d_opt,id_opt] = min(dtot);
    gamma_opt = gamma(id_opt);
    varargout{1} = gamma_opt;
    
    if flags.do_plot
      figure; 
      plot(gamma,dtot)
      xlabel('Gamma')
      ylabel('Prediction deviation')
    end
    
  end
  
end

%% ------ BASELINE ----------------------------------------------------------
if flags.do_baseline || flags.do_fig3
  
  SL = 50; % presentation level of stimuli
  model.kv.SPL = SL + kv.SL2SPL;
  
  cachename = ['baseline_tar' num2str(model.kv.SPL,'%u') 'dB_' cachename];
  if not(isempty(model.flags.errorflag))
    cachename = [cachename '_' model.flags.errorflag];
  end
  [Pcorr,d,s] = amt_cache('get',cachename,flags.cachemode);
  
  if isempty(Pcorr)

    s = data_baumgartner2016('argimport',model.flags,model.kv);
    
    % # of targets for evaluation of prediction residues
    Ntargets = cat(1,s.Ntar); % # of targets
    Ntargets = cat(1,Ntargets{:});
    relfreq = Ntargets/sum(Ntargets(:));

    if isempty(model.flags.errorflag)
      errorflag = 'querrMiddlebrooks';
    else
      errorflag = model.flags.errorflag;
    end
    for ii = 1:length(s)
      [s(ii).err,pred,m] = baumgartner2016(s(ii).Obj,s(ii).Obj,...
        'argimport',model.flags,model.kv,'ID',s(ii).id,'Condition','baseline',...
        'mrsmsp',s(ii).mrs,'S',s(ii).S,'priordist',s(ii).priordist,errorflag);
      [s(ii).qe_pred,s(ii).pe_pred] = baumgartner2014_pmv2ppp(pred,'exptang',s(ii).itemlist(:,6));
    end
    
    if isempty(model.flags.errorflag) % QE and PE
    
      qe_exp = cat(1,s.qe_exp);
      pe_exp = cat(1,s.pe_exp);
      qe_pred = cat(1,s.qe_pred);
      pe_pred = cat(1,s.pe_pred);
      
      % correlation
      [Pcorr.qe.r,Pcorr.qe.p] = corrcoeff(qe_exp,qe_pred);
      [Pcorr.pe.r,Pcorr.pe.p] = corrcoeff(pe_exp,pe_pred);

      % prediction residues
      sd_pe = (pe_pred-pe_exp).^2; % squared differences
      d.pe = sqrt(relfreq(:)' * sd_pe(:));    % weighted RMS diff.
      sd_qe = (qe_pred-qe_exp).^2;
      d.qe = sqrt(relfreq(:)' * sd_qe(:));
      
    else % arbitrary localization measure
      
      for ii = 1:length(s)
        s(ii).err_exp = localizationerror(s(ii).itemlist,model.flags.errorflag);
      end
      err_exp = cat(1,s.err_exp);
      err_pred = cat(1,s.err);
      [Pcorr.r,Pcorr.p] = corrcoeff(err_exp,err_pred,2);
      sd = (err_pred-err_exp).^2; % squared differences
      d = sqrt(relfreq(:)' * sd(:));    % weighted RMS diff.
      
    end

    s = rmfield(s,{'Obj'});
    
    amt_cache('set',cachename,Pcorr,d,s)
  end
  
  if isempty(model.flags.errorflag)
    d.total = (d.pe/90 + d.qe/100) /2;
    amt_disp(['Corr. QE: ' num2str(Pcorr.qe.r,'%2.2f') ' (p = ',num2str(Pcorr.qe.p,'%0.3f'),'), PE: ' num2str(Pcorr.pe.r,'%2.2f') ' (p = ',num2str(Pcorr.pe.p,'%0.3f'),'), QE+PE: ' num2str((Pcorr.qe.r+Pcorr.pe.r)/2,'%2.2f')])
  end
  
  if nargout >0; varargout{1} = d;	end
  if nargout >1; varargout{2} = r;	end
  if nargout >2; varargout{3} = s;	end
  
  if flags.do_plot
    
    if isempty(model.flags.errorflag)
    
      qe_exp = cat(1,s.qe_exp);
      pe_exp = cat(1,s.pe_exp);
      qe_pred = cat(1,s.qe_pred);
      pe_pred = cat(1,s.pe_pred);

      figure
      subplot(122)
      minqe = min([qe_exp(:);qe_pred(:)])-2;
      maxqe = max([qe_exp(:);qe_pred(:)])+2;
      limqe = [minqe,maxqe];
      plot(limqe,limqe,'k--')
      hold on
      h(1) = plot(qe_exp,qe_pred,'kd');
      axis equal
      axis([minqe maxqe minqe maxqe])

      xlabel('Actual QE','FontSize',kv.FontSize)
      ylabel('Predicted QE','FontSize',kv.FontSize)
      title(['e_{QE} = ' num2str(d.qe,'%0.1f') '% , r_{QE} = ' num2str(Pcorr.qe.r,'%0.2f')],...
          'FontSize',kv.FontSize)

      subplot(121)
      minpe = min([pe_exp(:);pe_pred(:)])-2;
      maxpe = max([pe_exp(:);pe_pred(:)])+2;
      limpe = [minpe,maxpe];
      plot(limpe,limpe,'k--')
      hold on
      h(2) = plot(pe_exp,pe_pred,'ks');
      set(h,'MarkerFaceColor','w')
      axis equal
      axis([minpe maxpe minpe maxpe])
      xlabel('Actual PE','FontSize',kv.FontSize)
      ylabel('Predicted PE','FontSize',kv.FontSize)
      title(['e_{PE} = ' num2str(d.pe,'%0.1f') '\circ , r_{PE} = ' num2str(Pcorr.pe.r,'%0.2f')],...
          'FontSize',kv.FontSize)

    else
      
      figure
      plot(cat(1,s.err_exp),cat(1,s.err),'ko')
      axis equal square
      xlim = get(gca,'XLim');
      hold on
      plot(xlim,xlim,'k:')
      xlabel('Actual','FontSize',kv.FontSize)
      ylabel('Predicted','FontSize',kv.FontSize)
      title({model.flags.errorflag;...
        ['(e = ' num2str(d,'%0.1f') ' , r = ' num2str(r,'%0.2f') ')']},...
        'FontSize',kv.FontSize)
      
    end
  end
  
end


%% ------ SPATSTRAT EXAMPLES -----------------------------------------------------------
if flags.do_spatstrat_ex
  
  SL = 50; % presentation level of stimuli
  model.kv.SPL = SL + kv.SL2SPL;
  
  latdivision = 0;  % lateral angle
  dlat = 15;

  % Experimental Settings
  Conditions = {'BB','LP','W'};


  %% Computations
  s = data_baumgartner2016('argimport',model.flags,model.kv);  
  s = s(ismember({s.id},'NH58'));
  amt_disp(['Listener: ' s.id])
  chance = [];
  for C = 1:length(Conditions)

    Cond = Conditions{C};

    %% Data

    % Experimental data
    data = data_majdak2013(Cond);
    for ll = 1:length(s)
      if sum(ismember({data.id},s(ll).id)) % participant ?
        s(ll).itemlist=data(ismember({data.id},s(ll).id)).mtx;
        for ii = 1:length(latdivision)
          latresp = s(ll).itemlist(:,7);
          idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
          mm2 = s(ll).itemlist(idlat,:);
          chance = [chance;mm2];
          s(ll).target{ii} = mm2(:,6); % polar angle of target
          s(ll).response{ii} = mm2(:,8); % polar angle of response
        end
      end
    end


    for ll = 1:length(s)
        for ii = 1:length(latdivision)
            s(ll).spdtfs{ii} = 0;     % init
            s(ll).polang{ii} = 0;   % init
            [s(ll).spdtfs{ii},s(ll).polang{ii}] = extractsp(...
              latdivision(ii),s(ll).Obj);

            if C == 1       % Learn 
                s(ll).spdtfs_c{ii} = s(ll).spdtfs{ii};
            elseif C == 2   % Dummy
                temp=amt_load('baumgartner2014','spatstrat_lpfilter.mat');
                s(ll).spdtfs_c{ii} = filter(temp.blp,temp.alp,s(ll).spdtfs{ii});
            elseif C == 3   % Warped
                s(ll).spdtfs_c{ii} = warp_hrtf(s(ll).spdtfs{ii},s(ll).fs);
            end

        end
    end


    %% Run Model

    for ll = 1:length(s)
      qe = zeros(1,length(latdivision));
      pe = zeros(1,length(latdivision));
      qe_t = zeros(1,length(latdivision));
      pe_t = zeros(1,length(latdivision));
      for ii = 1:length(latdivision)

        switch Cond
          case 'BB'
            clabel = 'baseline';
          case 'LP'
            clabel = 'lowpassed';
          case 'W'
            clabel = 'warped';
        end
        [s(ll).p{ii},rang] = baumgartner2016(...
              s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},...
              'argimport',model.flags,model.kv,...
              'ID',s(ll).id,'Condition',clabel,'fs',s(ll).fs,...
              'mrsmsp',s(ll).mrs,'S',s(ll).S,'lat',latdivision(ii),...
              'polsamp',s(ll).polang{ii},'priordist',s(ll).priordist);
        respangs{ii} = rang;

        [ qe(ii),pe(ii) ] = baumgartner2014_pmv2ppp(s(ll).p{ii} , s(ll).polang{ii} , rang);

        if sum(ismember({data.id},s(ll).id)) % if participant then actual targets
          [ qe_t(ii),pe_t(ii) ] = baumgartner2014_pmv2ppp( ...
              s(ll).p{ii} , s(ll).polang{ii} , rang , s(ll).target{ii} );

        end

      end

      % Model results of pool
      s(ll).qe_pool(C,1) = mean(qe); 
      s(ll).pe_pool(C,1) = mean(pe);

      if sum(ismember({data.id},s(ll).id)) % participant ?
        % Actual experimental results
        s(ll).qe_exp(C,1) = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
        s(ll).pe_exp(C,1) = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
        s(ll).Nt(C,1) = size(s(ll).itemlist,1);
        % Model results of participants (actual target angles)
        if length(latdivision) == 3
          s(ll).qe_part(C,1) = (qe_t(1)*length(s(ll).target{1}) + ...
              qe_t(2)*length(s(ll).target{2}) + ...
              qe_t(3)*length(s(ll).target{3}))/...
              (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
          s(ll).pe_part(C,1) = (pe_t(1)*length(s(ll).target{1}) + ...
              pe_t(2)*length(s(ll).target{2}) + ...
              pe_t(3)*length(s(ll).target{3}))/...
              (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
        else 
          s(ll).qe_part(C,1) = mean(qe_t);
          s(ll).pe_part(C,1) = mean(pe_t);
        end

        if flags.do_plot

          if C == 1; figure; end
          subplot(1,3,C)
          ii = find(latdivision==0);
          responses = [];
          targets = [];
          for jj = ii
            responses = [responses;s(ll).response{jj}];
            targets = [targets;s(ll).target{jj}];
          end
          plot_baumgartner2014(s(ll).p{ii},s(ll).polang{ii},rang,...
                targets,responses,'MarkerSize',kv.MarkerSize,'cmax',0.05,'nocolorbar')
          text(90,240,Cond,...
            'FontSize',kv.FontSize,'Color','w','HorizontalAlignment','center')
          Nt = length(targets);
          tmp.m = [zeros(Nt,5) targets(:) zeros(Nt,1) responses(:)];
          tmp.qe = localizationerror(tmp.m,'querrMiddlebrooks');
          tmp.pe = localizationerror(tmp.m,'rmsPmedianlocal');
          title({['A: PE = ' num2str(tmp.pe,2) '\circ, QE = ' num2str(tmp.qe,2) '%'];...
            ['P: PE = ' num2str(pe_t(ii),2) '\circ, QE = ' num2str(qe_t(ii),2) '%']},'FontSize',kv.FontSize-1)
          xlabel('Target Angle (deg)','FontSize',kv.FontSize)
          ylabel('Response Angle (deg)','FontSize',kv.FontSize)
          set(gca,'FontSize',kv.FontSize-1)
          set(gca,'XTickLabel',{[];[];0;[];60;[];120;[];180;[];[]})
          set(gca,'YTickLabel',{-60;[];0;[];60;[];120;[];180;[];240})

        end
        
      end

    end

  end
  
  varargout{1} = s;

end


%% ------ SPATSTRAT -----------------------------------------------------------
if flags.do_spatstrat
  
  SL = 50; % presentation level of stimuli
  model.kv.SPL = SL + kv.SL2SPL;
  
  cachename = ['spatstrat_tar' num2str(model.kv.SPL,'%u') 'dB_' cachename];
  [r,d,s,act,pred] = amt_cache('get',cachename,flags.cachemode);
  
  if isempty(r)
    
    latdivision = 0;%[-20,0,20];            % lateral angle
    dlat = 30;%10;

    % Experimental Settings
    Conditions = {'BB','LP','W'};

    %% Computations
    s = data_baumgartner2016('argimport',model.flags,model.kv);
    for C = 1:length(Conditions)

      Cond = Conditions{C};

      %% Data

      % Experimental data
      data = data_majdak2013(Cond);
      for ll = 1:length(s)
        if sum(ismember({data.id},s(ll).id)) % if actual participant
          s(ll).itemlist=data(ismember({data.id},s(ll).id)).mtx; 
          for ii = 1:length(latdivision)
            latresp = s(ll).itemlist(:,7);
            idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
            mm2 = s(ll).itemlist(idlat,:);
            s(ll).target{ii} = mm2(:,6); % polar angle of target
            s(ll).response{ii} = mm2(:,8); % polar angle of response
          end
        end
      end


      for ll = 1:length(s)
          for ii = 1:length(latdivision)
              s(ll).spdtfs{ii} = 0;     % init
              s(ll).polang{ii} = 0;   % init
              [s(ll).spdtfs{ii},s(ll).polang{ii}] = extractsp(...
                latdivision(ii),s(ll).Obj);

              if C == 1       % Learn 
                  s(ll).spdtfs_c{ii} = s(ll).spdtfs{ii};
              elseif C == 2   % Dummy
                temp=amt_load('baumgartner2014','spatstrat_lpfilter.mat');
                s(ll).spdtfs_c{ii} = filter(temp.blp,temp.alp,s(ll).spdtfs{ii});
              elseif C == 3   % Warped
                  s(ll).spdtfs_c{ii} = warp_hrtf(s(ll).spdtfs{ii},s(ll).fs);
              end

          end
      end


      %% Run Model
      idpart = [];
      for ll = 1:length(s)
        qe = zeros(1,length(latdivision));
        pe = zeros(1,length(latdivision));
        qe_t = zeros(1,length(latdivision));
        pe_t = zeros(1,length(latdivision));
        for ii = 1:length(latdivision)

          switch Cond
            case 'BB'
              clabel = 'baseline';
            case 'LP'
              clabel = 'lowpassed';
            case 'W'
              clabel = 'warped';
          end

          if sum(ismember({data.id},s(ll).id)) % if actual participant actual targets
            idpart = [idpart,ll];
            [s(ll).p{ii},rang] = baumgartner2016(...
              s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},'argimport',model.flags,model.kv,...
              'ID',s(ll).id,'Condition',clabel,'fs',s(ll).fs,...
              'mrsmsp',s(ll).mrs,'S',s(ll).S,'lat',latdivision(ii),...
              'polsamp',s(ll).polang{ii},'priordist',s(ll).priordist);
            
            [ qe_t(ii),pe_t(ii) ] = baumgartner2014_pmv2ppp( ...
                s(ll).p{ii} , s(ll).polang{ii} , rang , s(ll).target{ii} );

          end

        end

        if sum(ismember({data.id},s(ll).id)) % if actual participant 
          % Actual experimental results
          s(ll).qe_exp(C,1) = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
          s(ll).pe_exp(C,1) = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
          s(ll).Nt(C,1) = size(s(ll).itemlist,1);
          % Model results of participants (actual target angles)
          if length(latdivision) == 3
            s(ll).qe_part(C,1) = (qe_t(1)*length(s(ll).target{1}) + ...
                qe_t(2)*length(s(ll).target{2}) + ...
                qe_t(3)*length(s(ll).target{3}))/...
                (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
            s(ll).pe_part(C,1) = (pe_t(1)*length(s(ll).target{1}) + ...
                pe_t(2)*length(s(ll).target{2}) + ...
                pe_t(3)*length(s(ll).target{3}))/...
                (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
          else 
            s(ll).qe_part(C,1) = mean(qe_t);
            s(ll).pe_part(C,1) = mean(pe_t);
          end

        end

      end
    end
    s = rmfield(s,{'Obj','spdtfs_c','spdtfs'});% reduce file size

    act.qe = [s(idpart).qe_exp]';
    act.pe = [s(idpart).pe_exp]';
    pred.qe = [s(idpart).qe_part]';
    pred.pe = [s(idpart).pe_part]';
    
    % Correlation coefficients
    r.qe =  corrcoeff(act.qe,pred.qe);
    disp(['QE: r = ' num2str(r.qe,'%0.2f')]);

    r.pe =  corrcoeff(act.pe,pred.pe);
    disp(['PE: r = ' num2str(r.pe,'%0.2f')]);


    % RMS Differences
    % individual:
    Ntargets = [s.Nt]'; % # of targets
    relfreq = Ntargets/sum(Ntargets(:));
    sd_pe = (pred.pe-act.pe).^2; % squared differences
    d.pe = sqrt(relfreq(:)' * sd_pe(:));    % weighted RMS diff.
    sd_qe = (pred.qe-act.qe).^2;
    d.qe = sqrt(relfreq(:)' * sd_qe(:));
  
    amt_cache('set',cachename,r,d,s,act,pred);
  else
    
    data = data_majdak2013;
    idpart = ismember({s.id},{data.id});
  
  end
  
  s = s(idpart);
  
  d.total = (d.pe/90 + d.qe/100) /2;
  
  if nargout >0; varargout{1} = d;	end
  if nargout >1; varargout{2} = r;	end
  if nargout >2; varargout{3} = s;	end
  
  if flags.do_plot
    
    %% Measures

    % Quartiles
    quart_pe_part = quantile(pred.pe,[.25 .50 .75]);
    quart_qe_part = quantile(pred.qe,[.25 .50 .75]);

    quart_pe_exp = quantile(act.pe,[.25 .50 .75]);
    quart_qe_exp = quantile(act.qe,[.25 .50 .75]);

    % Chance performance
    [qe0,pe0] = baumgartner2014_pmv2ppp('chance',...
      'rang',-30-model.kv.mrsmsp:1:210+model.kv.mrsmsp);
    
    
    %% Plots
    dx = 0.15;
    
    figure 

    subplot(121)
    errorbar((1:3)+dx,quart_pe_part(2,:),...
        quart_pe_part(2,:) - quart_pe_part(1,:),...
        quart_pe_part(3,:) - quart_pe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar((1:3),quart_pe_exp(2,:),...
        quart_pe_exp(2,:) - quart_pe_exp(1,:),...
        quart_pe_exp(3,:) - quart_pe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');

    plot([0,4],[pe0,pe0],'k--')

    title(['e_{PE} = ' num2str(d.pe,'%0.1f') '\circ , r_{PE} = ' num2str(r.pe,'%0.2f')],...
      'FontSize',kv.FontSize)
    ylabel('Local Polar RMS Error (deg)','FontSize',kv.FontSize)
    set(gca,...
        'XLim',[0.5 3.5],...
        'XTick',1:3,...
        'YLim',[27 57.5],...
        'XTickLabel',{'BB';'LP';'W'},...
        'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))
    l = legend('Model','Actual');
    set(l,'Location','north','FontSize',kv.FontSize-1)

    subplot(122)
    errorbar((1:3)+dx,quart_qe_part(2,:),...
        quart_qe_part(2,:) - quart_qe_part(1,:),...
        quart_qe_part(3,:) - quart_qe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar((1:3),quart_qe_exp(2,:),...
        quart_qe_exp(2,:) - quart_qe_exp(1,:),...
        quart_qe_exp(3,:) - quart_qe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');
    
    plot([0,4],[qe0 qe0],'k--')

    title(['e_{QE} = ' num2str(d.qe,'%0.1f') '% , r_{QE} = ' num2str(r.qe,'%0.2f')],...
      'FontSize',kv.FontSize)
    ylabel('Quadrant Error (%)','FontSize',kv.FontSize)
    set(gca,...
        'XLim',[0.5 3.5],...
        'XTick',1:3,...
        'YLim',[0.1 49],...
        'XTickLabel',{'BB';'LP';'W'},...
        'YAxisLocation','left',...
        'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))

    set(gcf,'PaperPosition',[1,1,10,3.5])
    
  end
end

%% ------ NUMCHAN EXAMPLES -----------------------------------------------------------
if flags.do_numchan_ex
  
  % Stimulus Level
  SL = 50; % presentation level of stimuli
  
  % Model Settings
  model.kv.SPL = SL + kv.SL2SPL;
  latdivision = 0;            % lateral angle
  dlat = 10;

  % Experimental Settings
  Conditions = {'BB','CL','N24','N9','N3'};

  % Vocoder Settings 
  N = [inf,inf,24,9,3];
  flow = 300;     % lowest corner frequency
  fhigh = 16000;  % highest corner frequency

  %% Computations
  s = data_baumgartner2016('argimport',model.flags,model.kv);
  s = s(ismember({s.id},'NH33')); 
  disp(['Listener: ' s.id])
  chance = [];
  for C = 1:length(Conditions)

    Cond = Conditions{C};

    %% Data

    % Experimental data
    data = data_goupell2010(Cond);
    for ll = 1:length(s)
      if sum(ismember({data.id},s(ll).id)) % if actual participant
        s(ll).itemlist=data(ismember({data.id},s(ll).id)).mtx; 
        for ii = 1:length(latdivision)
          latresp = s(ll).itemlist(:,7);
          idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
          mm2 = s(ll).itemlist(idlat,:);
          chance = [chance;mm2];
          s(ll).target{ii} = mm2(:,6); % polar angle of target
          s(ll).response{ii} = mm2(:,8); % polar angle of response
        end
      end
    end

    % SP-DTFs
    for ll = 1:length(s)
        for ii = 1:length(latdivision)
            s(ll).spdtfs{ii} = 0;   % init
            s(ll).polang{ii} = 0;   % init
            [s(ll).spdtfs{ii},s(ll).polang{ii}] = extractsp(latdivision(ii),s(ll).Obj);
        end
    end


    %% Genereate conditional HRIRs

    stimPar.SamplingRate = s(ll).fs;
    imp = [1;zeros(2^12-1,1)]; % smooth results for 2^12
    for ll = 1:length(s)
      for ii = 1:length(latdivision)

        if strcmp(Cond,'BB') || strcmp(Cond,'CL')
            s(ll).spdtfs_c{ii} = s(ll).spdtfs{ii};

        else
          n = N(C);

          [syncrnfreq, GETtrain] = GETVocoder('',imp,n,flow,fhigh,0,100,stimPar);
          corners = [syncrnfreq(1);syncrnfreq(:,2)];

          ref = s(ll).spdtfs{ii};

          cond = zeros(length(imp),size(ref,2),2);

          for ch = 1:size(ref,3)
              for ang = 1:size(ref,2)
                  cond(:,ang,ch) = channelize('', 0.5*ref(:,ang,ch), ref(:,1), imp, n, corners, [], ...
                                  GETtrain, stimPar, 1, 0.01*s(ll).fs, 0.01*s(ll).fs);
              end
          end

          s(ll).spdtfs_c{ii} = cond;

        end
      end
    end

    if flags.do_plot 
     if C==1; fp = figure; end
     subplot(2,length(Conditions),C)
      if strcmp(Cond,'CL')
        stim = repmat([1;zeros(s(ll).fs/100-1,1)],17,1); % pulse train with 100pps
      else
        stim = noise(8e3,1,'white');
      end
      sig = lconv(stim,s(ll).spdtfs_c{ii});
     [mp,fc] = baumgartner2016_spectralanalysis(sig,kv.SPL,...
       'argimport',model.flags,model.kv,'target','ID',s(ll).id,'Condition',Cond);
     pcolor(fc,s(ll).polang{ii},mp(:,:,1)')
     shading flat
     xlabel('Frequency (Hz)')
     ylabel('Discharge rate')
     title(Cond)
    end
    

    %% Run Model

    for ll = 1:length(s)
      clear qe pe qe_t pe_t
      for ii = 1:length(latdivision)
        if strcmp(Cond,'CL')
          stim = repmat([1;zeros(s(ll).fs/100-1,1)],10,1); % pulse train with 100pps
        else
          stim = [];
        end
        [p,rang] = baumgartner2016(...
              s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},'argimport',model.flags,model.kv,...
              'ID',s(ll).id,'Condition',Cond,...
              'stim',stim,'fsstim',s(ll).fs,...
              'mrsmsp',s(ll).mrs,'S',s(ll).S,'lat',latdivision(ii),...
              'polsamp',s(ll).polang{ii},'priordist',s(ll).priordist);

        [ qe(ii),pe(ii) ] = baumgartner2014_pmv2ppp(p , s(ll).polang{ii} , rang);

        if sum(ismember({data.id},s(ll).id)) % if participant then actual targets
          [ qe_t(ii),pe_t(ii) ] = baumgartner2014_pmv2ppp( ...
              p , s(ll).polang{ii} , rang , s(ll).target{ii} );
        end

      end

      % Model results of pool
      s(ll).qe_pool(C,1) = mean(qe); 
      s(ll).pe_pool(C,1) = mean(pe);

      if sum(ismember({data.id},s(ll).id)) % if actual participant
        % Actual experimental results
        s(ll).qe_exp(C,1) = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
        s(ll).pe_exp(C,1) = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
        s(ll).Nt(C,1) = size(s(ll).itemlist,1);
        % Model results of participants (actual target angles)
        if length(latdivision) == 3
          s(ll).qe_part(C,1) = (qe_t(1)*length(s(ll).target{1}) + ...
              qe_t(2)*length(s(ll).target{2}) + ...
              qe_t(3)*length(s(ll).target{3}))/...
              (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
          s(ll).pe_part(C,1) = (pe_t(1)*length(s(ll).target{1}) + ...
              pe_t(2)*length(s(ll).target{2}) + ...
              pe_t(3)*length(s(ll).target{3}))/...
              (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
        else 
          s(ll).qe_part(C,1) = mean(qe_t);
          s(ll).pe_part(C,1) = mean(pe_t);
        end

        if flags.do_plot && latdivision(ii) == 0

          subplot(2,length(Conditions),C+length(Conditions))

          plot_baumgartner2014(p,s(ll).polang{ii},rang,...
                s(ll).target{ii},s(ll).response{ii},...
                    'MarkerSize',kv.MarkerSize,'cmax',0.05,'nocolorbar');
          title({['A: PE = ' num2str(s(ll).pe_exp(C,1),2) '\circ, QE = ' num2str(s(ll).qe_exp(C,1),2) '%'];...
            ['P: PE = ' num2str(s(ll).pe_part(C,1),2) '\circ, QE = ' num2str(s(ll).qe_part(C,1),2) '%']},'FontSize',kv.FontSize-1)
          text(90,240,Cond,...
            'FontSize',kv.FontSize,'Color','w','HorizontalAlignment','center')
          xlabel('Target Angle (deg)','FontSize',kv.FontSize)
          ylabel('Response Angle (deg)','FontSize',kv.FontSize)
          set(gca,'FontSize',kv.FontSize-1)
          set(gca,'XTickLabel',{[];[];0;[];60;[];120;[];180;[];[]})
          set(gca,'YTickLabel',{-60;[];0;[];60;[];120;[];180;[];240})

        end

      end

    end
  end
  
  varargout{1} = s;
  
end

%% ------ FIG NUMCHAN ----------------------------------------------------------
if flags.do_numchan
  
  % Stimulus Level
  SL = 50; % presentation level of stimuli
  model.kv.SPL = SL + kv.SL2SPL;
  
  cachename = ['numchan_tar' num2str(model.kv.SPL,'%u') 'dB_' cachename];
  [N,r,d,s] = amt_cache('get',cachename,flags.cachemode);
  if isempty(N)
    
    % Model Settings
    latdivision = 0; % lateral angle
    dlat = 30;

    % Experimental Settings
    Conditions = {'BB','CL','N24','N18','N12','N9','N6','N3'};

    % Vocoder Settings 
    N = fliplr([3,6,9,12,18,24,30,36]);	% # of vocoder channels
    flow = 300;     % lowest corner frequency
    fhigh = 16000;  % highest corner frequency


    %% Computations
    s = data_baumgartner2016('argimport',model.flags,model.kv);
    for C = 1:length(Conditions)

      Cond = Conditions{C};

      %% Data

      % Experimental data
      data = data_goupell2010(Cond);
      for ll = 1:length(s)
        if sum(ismember({data.id},s(ll).id)) % if actual participant
          s(ll).itemlist=data(ismember({data.id},s(ll).id)).mtx; 
          for ii = 1:length(latdivision)
            latresp = s(ll).itemlist(:,7);
            idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
            mm2 = s(ll).itemlist(idlat,:); 
%             chance = [chance;mm2];
            s(ll).target{ii} = mm2(:,6); % polar angle of target
            s(ll).response{ii} = mm2(:,8); % polar angle of response
          end
        end
      end

      % SP-DTFs
      for ll = 1:length(s)
        for ii = 1:length(latdivision)
          s(ll).spdtfs{ii} = 0;   % init
          s(ll).polang{ii} = 0;   % init
          [s(ll).spdtfs{ii},s(ll).polang{ii}] = extractsp(latdivision(ii),s(ll).Obj);
        end
      end


      %% Genereate conditional HRIRs

      stimPar.SamplingRate = s(ll).fs;
      imp = [1;zeros(2^12-1,1)]; % smooth results for 2^12
      for ll = 1:length(s)
        for ii = 1:length(latdivision)

          if strcmp(Cond,'BB') || strcmp(Cond,'CL')
            s(ll).spdtfs_c{ii} = s(ll).spdtfs{ii};

          else
            n = N(C);
            cachenameGET = ['numchan_GET_N' num2str(n) '_' s(ll).id];
            cond = amt_cache('get',cachenameGET);
            if isempty(cond)
              
              [syncrnfreq, GETtrain] = GETVocoder('',imp,n,flow,fhigh,0,100,stimPar);
              corners = [syncrnfreq(1);syncrnfreq(:,2)];

              ref = s(ll).spdtfs{ii};

              cond = zeros(length(imp),size(ref,2),2);

              for ch = 1:size(ref,3)
                for ang = 1:size(ref,2)
                    cond(:,ang,ch) = channelize('', 0.5*ref(:,ang,ch), ref(:,1), imp, n, corners, [], ...
                                    GETtrain, stimPar, 1, 0.01*s(ll).fs, 0.01*s(ll).fs);
                end
              end
              
              amt_cache('set',cachenameGET,cond);
            end
            s(ll).spdtfs_c{ii} = cond;

          end
        end
      end

      %% Run Model

      idpart = [];
      for ll = 1:length(s)
        clear qe pe qe_t pe_t
        if sum(ismember({data.id},s(ll).id)) % if actual participant actual targets
          
          for ii = 1:length(latdivision)
            if strcmp(Cond,'CL')
              stim = repmat([1;zeros(s(ll).fs/100-1,1)],17,1); % pulse train with 100pps
            else
              stim = [];
            end
            
            [p,rang] = baumgartner2016(...
                s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},...
                'argimport',model.flags,model.kv,...
                'ID',s(ll).id,'Condition',Cond,...
                'stim',stim,'fsstim',s(ll).fs,...
                'mrsmsp',s(ll).mrs,'S',s(ll).S,'lat',latdivision(ii),...
                'polsamp',s(ll).polang{ii},'priordist',s(ll).priordist);
              
            [ qe_t(ii),pe_t(ii) ] = baumgartner2014_pmv2ppp( ...
                p , s(ll).polang{ii} , rang , s(ll).target{ii} );
          end
          
          idpart = [idpart,ll];
          % Actual experimental results
          s(ll).qe_exp(C,1) = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
          s(ll).pe_exp(C,1) = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
          s(ll).Nt(C,1) = size(s(ll).itemlist,1);
          % Model results of participants (actual target angles)
          if length(latdivision) == 3
            s(ll).qe_part(C,1) = (qe_t(1)*length(s(ll).target{1}) + ...
                qe_t(2)*length(s(ll).target{2}) + ...
                qe_t(3)*length(s(ll).target{3}))/...
                (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
            s(ll).pe_part(C,1) = (pe_t(1)*length(s(ll).target{1}) + ...
                pe_t(2)*length(s(ll).target{2}) + ...
                pe_t(3)*length(s(ll).target{3}))/...
                (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
          else 
            s(ll).qe_part(C,1) = mean(qe_t);
            s(ll).pe_part(C,1) = mean(pe_t);
          end
          
        end

      end
      disp(['Condition ' Cond ' completed.'])
    end
    
    act.qe = [s(idpart).qe_exp];
    act.pe = [s(idpart).pe_exp];
    pred.qe = [s(idpart).qe_part];
    pred.pe = [s(idpart).pe_part];
    
    r.qe =  corrcoeff(act.qe,pred.qe);
    disp(['QE: r = ' num2str(r.qe,'%0.2f')]);

    r.pe =  corrcoeff(act.pe,pred.pe);
    disp(['PE: r = ' num2str(r.pe,'%0.2f')]);
    
    % RMS Differences
    % individual:
    Ntargets = [s.Nt]; % # of targets
    relfreq = Ntargets/sum(Ntargets(:));
    sd_pe = (pred.pe-act.pe).^2; % squared differences
    d.pe = sqrt(relfreq(:)' * sd_pe(:));    % weighted RMS diff.
    sd_qe = (pred.qe-act.qe).^2;
    d.qe = sqrt(relfreq(:)' * sd_qe(:));
    
    s = rmfield(s,{'spdtfs','spdtfs_c','Obj','itemlist'});
    
    amt_cache('set',cachename,N,r,d,s)
  end
  
  data = data_goupell2010;
  idpart = ismember({s.id},{data.id});
  s = s(idpart);
  
  d.total = (d.pe/90 + d.qe/100) /2;
  
  if nargout >0; varargout{1} = d;	end
  if nargout >1; varargout{2} = r;	end
  if nargout >2; varargout{3} = s;	end
  if nargout >3; varargout{4} = N;	end
  

  if flags.do_plot
    
    %% Quartiles
    quart_pe_part = fliplr(quantile([s(idpart).pe_part]',[.25 .50 .75]));
    quart_qe_part = fliplr(quantile([s(idpart).qe_part]',[.25 .50 .75]));

    quart_pe_exp = fliplr(quantile([s(idpart).pe_exp]',[.25 .50 .75]));
    quart_qe_exp = fliplr(quantile([s(idpart).qe_exp]',[.25 .50 .75]));

    [qe0,pe0] = baumgartner2014_pmv2ppp('chance',...
      'rang',-30-model.kv.mrsmsp:1:210+model.kv.mrsmsp);
    
    
    %% Plot
    dx = 0.7;
    figure

    %% PE
    subplot(121)
    errorbar(fliplr(N)+dx,quart_pe_part(2,:),...
        quart_pe_part(2,:) - quart_pe_part(1,:),...
        quart_pe_part(3,:) - quart_pe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar(fliplr(N),quart_pe_exp(2,:),...
        quart_pe_exp(2,:) - quart_pe_exp(1,:),...
        quart_pe_exp(3,:) - quart_pe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');
    plot([0,2*max(N)],[pe0,pe0],'k--')
    xlabel('Num. of channels','FontSize',kv.FontSize)
    ylabel('RMS of local errors (deg)','FontSize',kv.FontSize)

    title(['e = ' num2str(d.pe,'%0.1f') '\circ , r = ' num2str(r.pe,'%0.2f')],...
      'FontSize',kv.FontSize,'FontWeight','normal')
    set(gca,'XLim',[1 max(N)+2],'XTick',fliplr(N),...%[3 6 9 12 18 24 30],...
        'XTickLabel',{3;6;'';12;18;24;'CL';'BB'},...
        'YLim',[27 57.5],...
        'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))
      
    l = legend('Model','Actual');
    set(l,'Location','northeast','FontSize',kv.FontSize)

    %% QE
    subplot(122)
    errorbar(fliplr(N)+dx,quart_qe_part(2,:),...
        quart_qe_part(2,:) - quart_qe_part(1,:),...
        quart_qe_part(3,:) - quart_qe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar(fliplr(N),quart_qe_exp(2,:),...
        quart_qe_exp(2,:) - quart_qe_exp(1,:),...
        quart_qe_exp(3,:) - quart_qe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');
      
    plot([0,2*max(N)],[qe0,qe0],'k--')
      
    title(['e = ' num2str(d.qe,'%0.1f') '% , r = ' num2str(r.qe,'%0.2f')],...
      'FontSize',kv.FontSize,'FontWeight','normal')
    xlabel('Num. of channels','FontSize',kv.FontSize)
    ylabel('% Quadrant errors','FontSize',kv.FontSize)
    set(gca,'XLim',[1 max(N)+2],'XTick',fliplr(N),...%[3 6 9 12 18 24 30],...
        'XTickLabel',{3;6;'';12;18;24;'CL';'BB'},...
        'YLim',[0.1 49],...
        'YMinorTick','on',...
        'YAxisLocation','left','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))
      

    set(gcf,'PaperPosition',[1,1,10,3.5])

  end
end

%% SPLtem evaluation based on SpatStrat and NumChan
if flags.do_evalSPLtem  || flags.do_fig6
  
  SL2SPL = kv.SL2SPL;

  SPLtem = 40:10:80;
  SPLtemR = [40,80];

  for ll = 1:length(SPLtem)
%     [bl(ll).d,bl(ll).r,bl(ll).s] = exp_baumgartner2016('baseline','noplot',...
%       'ModelSettings',{'SPLtem',SPLtem(ll)},'SL2SPL',SL2SPL);
    [nc(ll).d,nc(ll).r,nc(ll).s] = exp_baumgartner2016('numchan','noplot',...
      'ModelSettings',{'argimport',model.flags,model.kv,'SPLtem',SPLtem(ll)},...
      'SL2SPL',SL2SPL);
    [ss(ll).d,ss(ll).r,ss(ll).s] = exp_baumgartner2016('spatstrat','noplot',...
      'ModelSettings',{'argimport',model.flags,model.kv,'SPLtem',SPLtem(ll)},...
      'SL2SPL',SL2SPL);
  end

  ll = length(SPLtem)+1;
  [nc(ll).d,nc(ll).r,nc(ll).s] = exp_baumgartner2016('numchan','noplot',...
    'ModelSettings',{'argimport',model.flags,model.kv,'SPLtem',SPLtemR},'SL2SPL',SL2SPL);
  [ss(ll).d,ss(ll).r,ss(ll).s] = exp_baumgartner2016('spatstrat','noplot',...
    'ModelSettings',{'argimport',model.flags,model.kv,'SPLtem',SPLtemR},'SL2SPL',SL2SPL);

  %% Number of listener-specific data points (#subjects * #conditions)
  N.nc = length(nc(1).s)*length(nc(1).s(1).pe_exp);
  N.ss = length(ss(1).s)*length(ss(1).s(1).pe_exp);
  N.all = N.nc + N.ss;

  %% Pool residues across experiments (averaged acc. to #data)
  for ll = 1:length(SPLtem)+1
    dtotal(ll) = (N.nc*nc(ll).d.total + N.ss*ss(ll).d.total) / N.all;
    dpe(ll) = (N.nc*nc(ll).d.pe + N.ss*ss(ll).d.pe) / N.all;
    dqe(ll) = (N.nc*nc(ll).d.qe + N.ss*ss(ll).d.qe) / N.all;
  end

  %% Prediction residuum obtained by chance prediction
  chance.dpe = 1;
  chance.dqe = 1;
  chance.dtotal = 1;

  %% Plots
  XLim = [SPLtem(1),SPLtem(end)]+[-1,1]*3;
  
  figure
  hax = tight_subplot(1,1,kv.gap,kv.marg_h,kv.marg_w);
  axes(hax(1))
  % single-SPL
  hi(1,1) = plot(SPLtem,dqe(1:length(SPLtem))/chance.dqe,'d-');
  hold on
  hi(2,1) = plot(SPLtem,dpe(1:length(SPLtem))/chance.dpe,'s-');
  % multiple-SPL
  hi(3,1) = plot(SPLtemR,dqe(end)/chance.dqe*[1,1],'--');
  hi(5,1) = plot(mean(SPLtemR),dqe(end)/chance.dqe,'d');
  hi(4,1) = plot(SPLtemR,dpe(end)/chance.dpe*[1,1],'--');
  hi(6,1) = plot(mean(SPLtemR),dpe(end)/chance.dpe,'s');
  % general
  set(gca,'YLim',[3.1,10.9],'XLim',XLim)
%   xlabel('Template SPL (dB)')
  ylabel('Prediction residuum')

  amt_disp('Predictive power for SpatStrat and NumChan.')

  xlabel('Template SPL (dB)')
  
  htmp = plot([0,0],[0,0],'k-');

  leg = legend([hi(5:6,1);htmp;hi(3,1)],...
    '% Quadrant errors','Local RMS error (deg)','Single SPL','Multiple SPLs');
  set(leg,'Location','northoutside','Orientation','vertical')
  
  set(hi,'LineWidth',1,'MarkerSize',kv.MarkerSize,'Color',zeros(1,3))
  set(hi,'MarkerFaceColor',zeros(1,3))
  set(hi(3:4,:),'LineStyle','--')
  
end

%% Modeling Sabin et al. (2005)
if flags.do_sabin2005
  
  % Presentation level of baseline stimuli (for calibration)
  SL = 50; 
  model.kv.SPL = SL + kv.SL2SPL;
  
  cachename = ['sabin2005_calibtar' num2str(model.kv.SPL,'%u') 'dB_' cachename];
  cachename = [cachename '_SL2SPL' num2str(kv.SL2SPL) 'dB'];
  if model.kv.lat ~= 0
    cachename = [cachename '_lat' num2str(model.kv.lat)];
  end
  if model.flags.do_SPLtemAdapt
    cachename = [cachename '_SPLtemAdapt'];
  end
  if model.flags.do_gammatone
    cachename = [cachename '_minSPL' num2str(model.kv.GT_minSPL)];
  end
  if kv.sabin2005_Sdivisor ~= 1
    cachename = [cachename '_Sdiv' num2str(kv.sabin2005_Sdivisor*10)];
  end
  if flags.do_nomrs
    cachename = [cachename '_nomrs'];
  end
  pred = amt_cache('get',cachename,flags.cachemode);
  if isempty(pred)
  
    s = data_baumgartner2016('argimport',model.flags,model.kv);
    
    SPL = [0:5:20,30:10:70]+kv.SL2SPL;

    if flags.do_nomrs
      model.kv.mrsmsp = 0;
    end
    
    pred.pvfront = nan(length(SPL),length(s));
    pred.pvrear = nan(length(SPL),length(s));
    pred.gfront = nan(length(SPL),length(s));
    pred.grear = nan(length(SPL),length(s));
    pred.precfront = nan(length(SPL),length(s));
    pred.precrear = nan(length(SPL),length(s));
    pred.qe = nan(length(SPL),length(s));
    pred.pe = nan(length(SPL),length(s));
    pred.prob = cell(length(SPL),length(s));
    for ii = 1:length(s)
      for ll = 1:length(SPL)
        if flags.do_nostat
          [pred.qe(ll,ii),pred.prob{ll,ii},m] = baumgartner2016(...
            s(ii).Obj,s(ii).Obj,...
            'argimport',model.flags,model.kv,...
            'ID',s(ii).id,'S',s(ii).S/kv.sabin2005_Sdivisor,'SPL',SPL(ll),...
            'priordist',s(ii).priordist,'QE');
          pred.pvfront(ll,ii) = nan;
          pred.pvrear(ll,ii) = nan;
          pred.gfront(ll,ii) = nan;
          pred.grear(ll,ii) = nan;
          pred.precfront(ll,ii) = nan;
          pred.precrear(ll,ii) = nan;
        else
          [pred.pvfront(ll,ii),pred.prob{ll,ii},m] = baumgartner2016(...
            s(ii).Obj,s(ii).Obj,...
            'argimport',model.flags,model.kv,...
            'ID',s(ii).id,'S',s(ii).S/kv.sabin2005_Sdivisor,'SPL',SPL(ll),...
            'priordist',s(ii).priordist,'pVeridicalPfront');
          pred.pvrear(ll,ii) = localizationerror(m,'pVeridicalPrear');
          pred.gfront(ll,ii) = localizationerror(m,'gainPfront');
          pred.grear(ll,ii) = localizationerror(m,'gainPrear');
          pred.precfront(ll,ii) = localizationerror(m,'precPregressFront');
          pred.precrear(ll,ii) = localizationerror(m,'precPregressRear');
          pred.qe(ll,ii) = localizationerror(m,'querrMiddlebrooks');
          pred.pe(ll,ii) = localizationerror(m,'rmsPmedianlocal');
        end
      end
      amt_disp([num2str(ii,'%u') ' of ' num2str(length(s),'%u') ' done'],'progress')
    end

    pred.SPL = SPL;
    amt_cache('set',cachename,pred)
  end
  
  data = data_sabin2005;
  data.SPL = data.SL + kv.SL2SPL; % assumption on SL
  
  %% Gain conversion: slope in deg -> limited range and equidistant
  if flags.do_gainInDeg
    pred.gfront = gain2slope(pred.gfront);
    pred.grear = gain2slope(pred.grear);
    data.gain.f.m = gain2slope(data.gain.f.m);
    data.gain.r.m = gain2slope(data.gain.r.m);
    data.gain.f.sd = gain2slope(data.gain.f.sd);
    data.gain.r.sd = gain2slope(data.gain.r.sd);
  end
  
  %% Restriction to reliable data (at least 75% audible trials in Sabin2005)
  dvar = {'gain','pqv','var'}; % data variable names
  mvar = {'g','pv','prec'}; % modeled variable names
  SL = pred.SPL-kv.SL2SPL; % assumption on SL
  SPL = pred.SPL;
  
  minSL = 15;
  iddata = data.SL >= minSL;
  idpred = SL >= minSL;

  % Remove
  for ii = 1:length(mvar)
    eval(['data.' dvar{ii} '.f.m = data.' dvar{ii} '.f.m(iddata);'])
    eval(['data.' dvar{ii} '.f.sd = data.' dvar{ii} '.f.sd(iddata);'])
    eval(['pred.' mvar{ii} 'front = pred.' mvar{ii} 'front(idpred,:);'])
    
    eval(['data.' dvar{ii} '.r.m = data.' dvar{ii} '.r.m(iddata);'])
    eval(['data.' dvar{ii} '.r.sd = data.' dvar{ii} '.r.sd(iddata);'])
    eval(['pred.' mvar{ii} 'rear = pred.' mvar{ii} 'rear(idpred,:);'])
  end
  data.SL = data.SL(iddata);
  data.SPL = data.SPL(iddata);
  SL = SL(idpred);
  SPL = SPL(idpred);
  pred.SPL = SPL;
  pred.SL = SL;
  
  %% Pool data for 20 dB
  id20 = find(data.SL == 20,2,'first');
  data20p = data; % init
  for ii = 1:length(mvar)
    eval(['data20p.' dvar{ii} '.f.m = [data.' dvar{ii} '.f.m(1:' num2str(id20(1)-1) '),mean(data.' dvar{ii} '.f.m(' num2str(id20(1)) ':'  num2str(id20(2)) ')),data.' dvar{ii} '.f.m(' num2str(id20(2)+1) ':end)];']) % 20dB SL pooled
    eval(['data20p.' dvar{ii} '.r.m = [data.' dvar{ii} '.r.m(1:' num2str(id20(1)-1) '),mean(data.' dvar{ii} '.r.m(' num2str(id20(1)) ':'  num2str(id20(2)) ')),data.' dvar{ii} '.r.m(' num2str(id20(2)+1) ':end)];'])
  end
  
  %% Prediction deviation score
  limits = {'90','100','45'}; % for normalization
  idc = ismember(SL,data.SL);
  
  d = zeros(length(mvar),2);
  for ii = 1:length(mvar)
    if flags.do_nostat
      eval(['d(ii,1) = sqrt(mean((mean(pred.' mvar{ii} 'front(idc,:)),2) - transpose(data20p.' dvar{ii} '.f.m)).^2))/' limits{ii} ';'])
      eval(['d(ii,2) = sqrt(mean((mean(pred.' mvar{ii} 'rear(idc,:),2) - transpose(data20p.' dvar{ii} '.r.m)).^2))/' limits{ii} ';'])
    else
      eval(['d(ii,1) = sqrt(nanmean((nanmean(pred.' mvar{ii} 'front(idc,:),2) - transpose(data20p.' dvar{ii} '.f.m)).^2))/' limits{ii} ';'])
      eval(['d(ii,2) = sqrt(nanmean((nanmean(pred.' mvar{ii} 'rear(idc,:),2) - transpose(data20p.' dvar{ii} '.r.m)).^2))/' limits{ii} ';'])
    end
  end
  d = mean(d(:));
  amt_disp(['Prediction deviation score: ' num2str(d)])
  
  %% Correlation coefficient
  
  r = zeros(length(mvar),2);
  for ii = 1:length(mvar)
    if flags.do_stat
      eval(['r(ii,1) = corrcoeff( nanmean(pred.' mvar{ii} 'front(idc,:),2) , transpose(data20p.' dvar{ii} '.f.m) );'])
      eval(['r(ii,2) = corrcoeff( nanmean(pred.' mvar{ii} 'rear(idc,:),2) , transpose(data20p.' dvar{ii} '.r.m) );'])
    end
  end
  r = mean(r(:));
  amt_disp(['Correlation coefficient: ' num2str(r)])
  
  %% Output
  varargout{1}=d;
  varargout{2}=r;
  varargout{3}=pred;
  varargout{4}=data;

  %% Plot
  if flags.do_plot
    
    minSPLf = minSL + kv.SL2SPL; % minSLf
    minSPLr = minSL + kv.SL2SPL; % minSLr
    maxSPL = 70 + kv.SL2SPL;
    marSPL = 4.9; % margin of SL
    
    if flags.do_gainInDeg
      ylim = {[-2,60],[-5,105],[8,42]}; % ylimits acc. to Sabin et al. (2005)
      elabel = {'Slope (deg)','% Quasi-veridical','Variability (deg)'}; % plot ylabels
    else
      ylim = {[-0.2,1.6],[-5,105],[8,42]}; % ylimits acc. to Sabin et al. (2005)
      elabel = {'Gain','% Quasi-veridical','Variability (deg)'}; % plot ylabels
    end
    
    figure;
    ha = tight_subplot(3,2,kv.gap,kv.marg_h,kv.marg_w);
    for ii = 1:length(mvar)
      axes(ha(1+(ii-1)*2))
      if ii == 1
        plot([-10,100],[45,45],'k--') % ideal slope
        hold on
      end
      eval(['p1 = errorbar(SPL-0.3,nanmean(pred.' mvar{ii} 'front,2),nanstd(pred.' mvar{ii} 'front,1,2));'])
      hold on
      eval(['p2 = errorbar(data.SPL+0.3,data.' dvar{ii} '.f.m,data.' dvar{ii} '.f.sd);'])
      set(p1,'MarkerFaceColor','k','Marker','^','Color','k','MarkerSize',kv.MarkerSize)
      set(p2,'MarkerFaceColor','w','Marker','^','Color','k','MarkerSize',kv.MarkerSize)
      axis([minSPLf-marSPL,maxSPL+marSPL,ylim{ii}])
      ylabel(elabel{ii},'FontSize',kv.FontSize)
      set(gca,'XTick',round(minSPLf/10)*10:10:maxSPL)
      set(gca,'TickLength',2*get(gca,'TickLength'),'FontSize',kv.FontSize)

      if ii == 1
        title('Front','FontSize',kv.FontSize)
      elseif ii == 3
        xlabel('SPL (dB)','FontSize',kv.FontSize)
      end
      if ii<3
        set(gca,'XTickLabel',[])
      end
      
      axes(ha(2+(ii-1)*2))  
      if ii == 1
        plot([-10,100],[45,45],'k--') % ideal slope
        hold on
      end    
      eval(['p1 = errorbar(SPL-0.3,nanmean(pred.' mvar{ii} 'rear,2),nanstd(pred.' mvar{ii} 'rear,1,2));'])
      hold on
      eval(['p2 = errorbar(data.SPL+0.3,data.' dvar{ii} '.r.m,data.' dvar{ii} '.r.sd);'])
      set(p1,'MarkerFaceColor','k','Marker','o','Color','k','MarkerSize',kv.MarkerSize)
      set(p2,'MarkerFaceColor','w','Marker','o','Color','k','MarkerSize',kv.MarkerSize)
      axis([minSPLr-marSPL,maxSPL+marSPL,ylim{ii}])
      set(gca,'XTick',round(minSPLr/10)*10:10:maxSPL,'YTickLabel',[])
      set(gca,'TickLength',2*get(gca,'TickLength'),'FontSize',kv.FontSize)
      if ii == 1
        title('Rear','FontSize',kv.FontSize)
      elseif ii == 3
        xlabel('SPL (dB)','FontSize',kv.FontSize)
      end
      
    end
  end
  
end

if flags.do_impairment
  
  if isempty(errorflag)
    errorflag = 'QE';
    amt_disp('Localization performance measure not chosen -> QE used.')
  end
  
  cohc = sort(kv.cOHCset,'descend'); % default: [1,0.4,0.1,0];
  ft = kv.FTset; % default: {1:3;1;2;3};
  SPL = sort(kv.SPLset,'ascend'); % default: [80, 50]
  
  cachename = ['impairment_' cachename '_' errorflag];
  if length(SPL) == 1
    cachename = [cachename '_SPL' num2str(SPL)];
  elseif length(SPL) == 2
    cachename = [cachename '_SPLs' num2str(SPL(1)) 'and' num2str(SPL(2))];
  end
  if model.flags.do_SPLtemAdapt
    cachename = [cachename '_SPLtemAdapt'];
  end
  if model.flags.do_NHtem
    cachename = [cachename '_NHtem'];
  end
  cachename = strrep([cachename '_cohc' num2str(cohc)],' ','');cachename = strrep(cachename,'.','p');
  s = amt_cache('get',cachename,flags.cachemode);
  if isempty(s)
    
    s = data_baumgartner2016('argimport',model.flags,model.kv);
      
    if not(isempty(kv.subjects))
      s = s(kv.subjects);
    end

    for ii=1:length(s)
      for cc=1:length(cohc)
        for ll=1:length(SPL)
          for ff=1:length(ft)
            err = baumgartner2016(s(ii).Obj,s(ii).Obj,...
              'argimport',model.flags,model.kv,'ID',s(ii).id,errorflag,'fiberTypes',ft{ff},...
              'S',s(ii).S,'cohc',cohc(cc),'SPL',SPL(ll),'priordist',s(ii).priordist);
            s(ii).err(cc,ff,ll) = err;
            s(ii).cohc(cc,ff,ll) = cohc(cc);
            s(ii).SPL(cc,ff,ll) = SPL(ll);
            s(ii).ft{cc,ff,ll} = ft{ff};
          end
        end
      end
      amt_disp([num2str(ii,'%u') ' of ' num2str(length(s),'%u') ' done'],'progress')
    end
    
    s = rmfield(s,{'Obj','itemlist'});
    amt_cache('set',cachename,s);
  end
  varargout{1} = s;
  
  err = cat(4,s.err);
  
  Ncohc = size(err,1);
  Nft = size(err,2);
  Nspl = size(err,3);
  Ncond = length(s(1).cohc(:)); % #conditions
  Nsub = length(s); % #subjects
  mtx.err = reshape(shiftdim(err,3),[Nsub,Ncond]);
  
  % labels
  cohcstrA = num2str(s(1).cohc(:),'%2.1f');
  cohcstr = cell(length(cohcstrA),1);
  for ii = 1:length(cohcstrA)
    cohcstr{ii} = strrep(cohcstrA(ii,:),'1.0','1');
    cohcstr{ii} = strrep(cohcstr{ii},'0.0','0');
  end
  SPLstr = num2str(s(1).SPL(:),'%2.0f');
  ftnum = s(1).ft(:);
  ftstr = cell(length(ftnum),1);
  XTickLabel = ftstr;
  for ii = 1:length(ftnum)
    lab = ['cohc: ' cohcstr(ii,:) ', ' SPLstr(ii,:) 'dB, '];
    if ftnum{ii} == 1
      ftstr{ii} = 'low-SR';
    elseif ftnum{ii} == 2
      ftstr{ii} = 'med-SR';
    elseif ftnum{ii} == 3
      ftstr{ii} = 'high-SR';
    else % ft{ii} == 1:3
      ftstr{ii} = 'all SRs';
    end
    XTickLabel{ii} = [lab ftstr{ii}];
  end
    
  % sort data acc. to ascending SPL
  [tmp,idSPLsort] = sort(s(1).SPL(:));
  ftstr = ftstr(idSPLsort);
  cohcstr = cellstr(cohcstr(idSPLsort,:));
  SPLstr = cellstr(SPLstr(idSPLsort,:));
  mtx.err = mtx.err(:,idSPLsort,:);
    
  % Output meta data
  meta(1).name = 'Listener';
  meta(2).name = 'Condition';
  meta(2).data = {ftstr(:),cohcstr(:),SPLstr(:)};
  
  if flags.do_plot
    
    % interaction plots
%     merr = mean(err,4); % average across subjects; dims: [cohc,ft,SPL]
%     fig(1) = figure('Name','Interaction Plots');
%     if length(kv.SPLset) > 1
%       % OHC-SPL
%       subplot(1,3,1)
%       err_OHC_SPL = squeeze(mean(merr,2));
%       plot(cohc,err_OHC_SPL)
%       legend(num2str(SPL(:)))
%       xlabel('C_{OHC}')
%       ylabel(errorflag)
%       % OHC-FT
%       subplot(1,3,2)
%       err_OHC_FT = squeeze(mean(merr,3));
%       plot(cohc,err_OHC_FT)
%       legend('All','LSR','MSR','HSR')
%       xlabel('C_{OHC}')
%       % SPL-FT
%       subplot(1,3,3)
%       err_FT_SPL = squeeze(mean(merr,1));
%       plot(SPL,err_FT_SPL)
%       legend('All','LSR','MSR','HSR')
%       xlabel('SPL')
%     else % only OHC-FT
%       plot(cohc,merr)
%       legend('All','LSR','MSR','HSR')
%       xlabel('C_{OHC}')
%     end
    
    emax = max(mtx.err(:))+0.5;
    emin = min(mtx.err(:))-0.5;
    De = emax-emin;
    
    [qe0,pe0] = baumgartner2014_pmv2ppp('chance',...
      'rang', -30-model.kv.mrsmsp : 1 : 210+model.kv.mrsmsp );
    
    if flags.do_splitSPL
      
      colors = {(1-1/length(kv.SPLset))*ones(1,3),zeros(1,3)};
      for ll = length(kv.SPLset):-1:1
        
        NcondPlot = Ncond/length(kv.SPLset);
         
        fig(ll+1) = figure;
        id = SPL(ll)==str2double(SPLstr);
        b = boxplot(mtx.err(:,id),{ftstr(id),cohcstr(id),SPLstr(id)},...
          'plotstyle','compact','medianstyle','line','colors',colors{ll},...
          'factorgap',[],'labelverbosity','all','symbol','');

        hold on
        % chance performance
        if strcmp(errorflag,'QE')
          echance = qe0;
        elseif strcmp(errorflag,'PE')
          echance = pe0;
        end
        h = plot([0.5,NcondPlot+0.5],[echance,echance],'k--');uistack(h,'bottom')
        
        if flags.do_FTlabel
          for ii = 1:Nft
            jj = 1+(ii-1)*Ncohc;
            yftlbl = 1.05*max(emax,echance);
            text(jj,yftlbl,ftstr{ii*Ncohc},'FontSize',kv.FontSize) % inside panel
          end
        end

        axis([0.5,NcondPlot+0.5,emin-De/20,emax+De/8])
        set(gca,'XTick',1:NcondPlot,'XTickLabel',cohcstr(1:NcondPlot),'FontSize',kv.FontSize);
        xlabel({'OHC gain, C_{OHC}'},'FontSize',kv.FontSize,'FontWeight','bold')
        ylabel(errorflag,'FontSize',kv.FontSize,'FontWeight','bold') 
        
      end
      
    else % compriseSPL
      
      fig(2) = figure;
      b = boxplot(mtx.err,{ftstr(:),cohcstr(:),SPLstr(:)},...
          'plotstyle','compact','colors',repmat([.5,.5,.5;0,0,0],Ncond/Nspl,1),'medianstyle','line',...
          'factorgap',[],'labelverbosity','all','symbol','');
      hold on
      if flags.do_FTlabel
        for ii = 1:Nft
          jj = 1+(ii-1)*Nspl*Ncohc;
          yftlbl = emax;
          text(jj,yftlbl,ftstr{ii*Ncohc},'FontSize',kv.FontSize) % inside panel
        end
      end

      % chance performance
      if strcmp(errorflag,'QE')
        h = plot([0.5,Ncond+0.5],[qe0,qe0],'k--');uistack(h, 'bottom')
      elseif strcmp(errorflag,'PE')
        h = plot([0.5,Ncond+0.5],[pe0,pe0],'k--');uistack(h, 'bottom')
      end

      axis([0.5,Ncond+0.5,emin-De/20,emax+De/8])
      set(gca,'XTick',1.5:2:Ncond,'XTickLabel',cohcstr(1:Ncohc),'FontSize',kv.FontSize);
      xlabel({'OHC gain, C_{OHC}'},'FontSize',kv.FontSize,'FontWeight','bold')
      ylabel(errorflag,'FontSize',kv.FontSize,'FontWeight','bold') 
    end
  else
    fig = [];
  end
  
  if flags.do_stat && ~verLessThan('matlab','8.2')
    s = rmfield(s,{'pe_exp','qe_exp','S','pe_exp_lat','qe_exp_lat','target','response','Ntar','mrs','fs'});
    
    t = array2table(mtx.err);
    if length(kv.SPLset) > 1 % 3-way repeated-measures ANOVA
      within = table(SPLstr,cohcstr,ftstr,'VariableNames',{'SPL','Cohc','FT'});
      rm = fitrm(t,['Var1-Var' num2str(Ncond) ' ~ 1'],'WithinDesign',within); % no between-subjects factors -> only intercept    
      [tbl.ranova,A,C,D] = ranova(rm,'WithinModel','Cohc*SPL*FT');
    else  % 2-way repeated-measures ANOVA
      within = table(cohcstr,ftstr,'VariableNames',{'Cohc','FT'});
      rm = fitrm(t,['Var1-Var' num2str(Ncond) ' ~ 1'],'WithinDesign',within); % no between-subjects factors -> only intercept    
      [tbl.ranova,A,C,D] = ranova(rm,'WithinModel','Cohc*FT');
    end
    tbl.ranova.Properties.RowNames = strrep(tbl.ranova.Properties.RowNames,'(Intercept):','');
    
    % Mauchly's test for sphericity
    tbl.mauchly = mauchly(rm,C);
    tbl.mauchly.Properties.RowNames = tbl.ranova.Properties.RowNames(1:2:end);
    
    % Sphericity corrections
    tbl.eps = epsilon(rm,C);
    tbl.eps.Properties.RowNames = tbl.ranova.Properties.RowNames(1:2:end);
    
    % Add corrected DFs to ranova table
    idrep = round(0.5:0.5:length(tbl.eps.GreenhouseGeisser)); % repeat iteratively
    tbl.ranova.DFGG = tbl.ranova.DF .* tbl.eps.GreenhouseGeisser(idrep);
    
    % Add effect sizes to ranova table
    SSeffect = tbl.ranova.SumSq(1:2:end);
    SSerror = tbl.ranova.SumSq(2:2:end);
    eta_pSq = nan(2*length(SSerror),1);
    eta_pSq(1:2:end) = SSeffect./(SSeffect+SSerror); % effect size per (eta_partial)^2
    tbl.ranova.eta_pSq = eta_pSq;
    
    % Post-hoc analysis
    tbl.posthoc.Cohc = multcompare(rm,'Cohc');
    tbl.posthoc.FT = multcompare(rm,'FT');
    
    % Display results
    amt_disp(['Repeated-measures ANOVA for ' errorflag])
    amt_disp(tbl.ranova)
    amt_disp('Mauchly test and sphericity corrections')
    amt_disp([tbl.mauchly,tbl.eps])
    amt_disp('Posthoc analysis')
    amt_disp(tbl.posthoc.Cohc)
    amt_disp(tbl.posthoc.FT)
    amt_disp('Reported in publication:')
    amt_disp(tbl.ranova(3:end,[9,4,6,10]))
    
  else
    tbl = [];
  end
  
  varargout = {tbl;fig;mtx.err;meta};
  
end

if flags.do_fig7
  
  errorflag = {...
    'QE','% Quadrant errors';...
    'PE','Local RMS error (deg)';...
    };
  
  SPLset = 60;
  
  NHtemflag = '';%'NHtem';

  for ii = 1:length(errorflag)
    tbl{ii} = exp_baumgartner2016('impairment',errorflag{ii,1},'SPLset',SPLset,...
      'noFTlabel','FontSize',kv.FontSize,flags.cachemode,'ModelSettings',{NHtemflag});
    ylabel(errorflag{ii,2},'FontSize',kv.FontSize)
    fig(ii) = gcf;
    ax(ii) = gca;
  end

  % Combine panels and add FT labels
  N.cohc = 4;
  N.ft = 4;
  labels = {'all SRs','low-SR','med-SR','high-SR'};
  
  figC = figure;
  marg = [.11,.06;.11,.03];
  ha = tight_subplot(length(errorflag),length(SPLset),0,marg(1,:),marg(2,:));
  for ii = 1:length(ha)
    
    axes(ha(ii));
    copyobj(allchild(ax(ii)),ha(ii))
    set(ha(ii),'XTick',get(ax(ii),'XTick'))
    set(ha(ii),'XTickLabel',get(ax(ii),'XTickLabel'))
    set(ha(ii),'YTick',get(ax(ii),'YTick'))
    set(ha(ii),'YTickLabel',get(ax(ii),'YTickLabel'))
    if ii <= length(SPLset) % add SPL and FT labels at top panels
      if length(SPLset) > 1
        title([num2str(SPLset(ii)) ' dB SPL'])
      end
      if length(SPLset) == 1
        yy = 49; % like title
      else
        yy = 55;
      end
      set(gca,'YLim',[-2,yy])
      for ff = 1:N.ft
        jj = (ff-0.5)*N.cohc+.5;
        text(jj,52,labels{ff},'FontWeight','bold','HorizontalAlignment','center')
      end
    else
      if length(SPLset) == 1
        yy = 51.5; % like title
      else
        yy = 54;
      end
      set(gca,'YLim',[17,yy])
    end
    if ii == 1 % show ylabel at left panels
      ylabel(errorflag{1,2},'FontWeight','bold')
    elseif ii == length(SPLset)+1
      ylabel(errorflag{2,2},'FontWeight','bold')
    else
      set(gca,'YTickLabel','')
      ylabel('')
    end
    if ii > length(SPLset) % show ylabel at bottom panels
      xlabel({' ';'OHC gain, C_{OHC}'},'FontWeight','bold')
    end
  end
  set(findall(figC,'-property','FontSize'),'FontSize',kv.FontSize)
  
  close(fig)
  
  varargout = {tbl,errorflag};
  
end

if flags.do_cOHCvsSens
  
  [sens,meta.sens] = exp_baumgartner2016('dynrangecheck','dprime','noplot',...
    'ModelSettings',{'argimport',model.flags,model.kv});
  idSPL = meta.sens(1).data == model.kv.SPL;
  Sensitivity = squeeze(sens(idSPL,:,:))'; % cOHC x  SR
  Sensitivity = circshift(Sensitivity,1,2); % all SRs first
  errorlabel = {'QE';'PE'};
  for ee = 1:length(errorlabel)
    
    [~,~,errorPrediction,meta.err] = exp_baumgartner2016('impairment',errorlabel{ee},...
      'noplot','nostat','ModelSettings',{'argimport',model.flags,model.kv});
    
    [r,p]=corrcoef(mean(errorPrediction,1)',Sensitivity(:));
    if not(isscalar(r))
      r = r(2);
    end
    
    amt_disp(['Correlation between average ',errorlabel{ee},...
      ' and intensity discriminability:'])
    amt_disp(['  r = ',num2str(r),' ( p = ',num2str(p),' )']);
    
  end
  
end

if flags.do_effectOnCues || flags.do_fig9
  
  sid = 10;    % listener No.
  
  s = data_baumgartner2016('argimport',model.flags,model.kv);
  [dtf,polang] = extractsp(0,s(sid).Obj);
  tmp = lconv(noise(8e3,1,'white'),dtf);
  sig = reshape(tmp,[size(tmp,1),size(dtf,2),size(dtf,3)]);
      
  amt_disp(['Exemplary listener: ' s(sid).id])
  
  ymin = 0;
  ymax = model.kv.mgs*pi;
  spl = kv.SPLset;
  
  if flags.do_FT
    % Effect of FT  
    FT = {1:3,1,2,3};
    mp_ft = cell(length(FT),length(spl));
    gp_ft = cell(length(FT),length(spl));
    for ll = 1:length(spl)
      for ft = 1:length(FT)
        [mp_ft{ft,ll},fc] = baumgartner2016_spectralanalysis(sig,spl(ll),...
          'target','ID',s(sid).id ,'Condition','baseline','fiberTypes',FT{ft},flags.cachemode);
        [gp_ft{ft,ll},gfc] = baumgartner2016_gradientextraction(mp_ft{ft,ll},fc,'mgs',kv.mgs);
      end
    end

    if flags.do_plot
      figure
      ha = tight_subplot(length(FT),length(cohc),kv.gap,kv.marg_h,kv.marg_w);
      colormap gray
      colormap(flipud(colormap))
      labels = {'LMH','L','M','H'};
      for ll = 1:length(spl)
        for ft = 1:length(FT)
          axes(ha(ft+(ll-1)*length(FT)))
          pcolor(gfc,polang,gp_ft{ft,ll}.m(:,:,1)')
          shading flat
          caxis([ymin,ymax])
          title(labels{ft})
          xlabel('Frequency (kHz)','FontWeight','bold')
          ylabel('Polar angle (deg)','FontWeight','bold')
          set(gca,'XScale','log','FontSize',kv.FontSize)
          set(gca,'layer','top',...
            'XTick',[1,2,4,8,16]*1e3,...
            'XTickLabel',[1,2,4,8,16],...
            'YTick',-30:30:180,...
            'FontSize',kv.FontSize)
        end
      end
    end

  else
    % Effect of C_OHC
    cohc = kv.cOHCset;
    gp_cohc = cell(length(cohc),length(spl));
    for ll = 1:length(spl)
      for cc = 1:length(cohc)
        [mp,fc] = baumgartner2016_spectralanalysis(sig,spl(ll),'target',...
          'argimport',flags,kv,'ID',s(sid).id,'Condition','baseline','cohc',cohc(cc));
        [gp_cohc{cc,ll},gfc] = baumgartner2016_gradientextraction(mp,fc,'mgs',model.kv.mgs);
      end
    end

    if flags.do_plot
      figure
      ha = tight_subplot(length(spl),length(cohc),[.02,.01],kv.marg_h,kv.marg_w);
      colormap gray
      colormap(flipud(colormap))
      labels = {'C_{OHC} = 1',['C_{OHC} = ' num2str(cohc(2),'%1.1f')],['C_{OHC} = ' num2str(cohc(3),'%1.1f')],'C_{OHC} = 0'};
      for ll = 1:length(spl)
        for cc = 1:length(cohc)
          axes(ha(cc+(ll-1)*length(cohc)))
          pcolor(gfc,polang,gp_cohc{cc,ll}.m(:,:,1)')
          shading flat
          caxis([ymin,ymax])
          % xlabel and COHC
          set(gca,'XTick',[1,2,4,8,16]*1e3,'XTickLabel',[1,2,4,8,16])
          if ll == 1
            title(labels{cc},'FontSize',kv.FontSize)
          end
          if ll == length(spl)
            if cc == 2
              xlabel(['                                            ',...
                'Frequency (kHz)'],'FontSize',kv.FontSize,'FontWeight','bold')
            end
          else
            set(gca,'XTickLabel',[])
          end
          % ylabel and SPL
          if cc == 1
            ylabel('Polar angle (deg)','FontSize',kv.FontSize,'FontWeight','bold')
            if length(spl) > 1
              text(180,200,[num2str(spl(ll)) ' dB'],'FontWeight','bold','FontSize',kv.FontSize)
            end
          else
            set(gca,'YTickLabel',[])
          end
          set(gca,'XScale','log','FontSize',kv.FontSize)
          set(gca,'layer','top',...
            'TickLength',2*get(gca,'TickLength'),... 
            'YTick',-30:30:180,...
            'FontSize',kv.FontSize)
        end
      end
      
      c = colorbar;
      set(c,'Position',[.93,.2,.02,.6])
      set(get(c,'Title'),'String','Spikes/s/ERB','FontSize',kv.FontSize)
      
    end
  end
  
end

if flags.do_evalSpectralContrast
  
  sid = 10;    % listener No.
  
  s = data_baumgartner2016('argimport',model.flags,model.kv);
  [dtf,polang] = extractsp(0,s(sid).Obj);
  tmp = lconv(noise(8e3,1,'white'),dtf);
  sig = reshape(tmp,[size(tmp,1),size(dtf,2),size(dtf,3)]);
      
  amt_disp(['Exemplary listener: ' s(sid).id])
  
  spl = kv.SPLset;
  
  FT = {1:3,1,2,3};
  mp_ft = cell(length(FT),length(spl));
  for ll = 1:length(spl)
    for ft = 1:length(FT)
      [mp_ft{ft,ll},fc] = baumgartner2016_spectralanalysis(sig,spl(ll),...
        'target','ID',s(sid).id ,'Condition','baseline','fiberTypes',FT{ft},flags.cachemode);
      spectRange(ft,ll) = range(mp_ft{ft,ll}(:));
    end
    spectralContrast(:,ll) = spectRange(2:end,ll)/spectRange(1,ll);
  end
  FTlabels = {'low-SR','med-SR','high-SR'};
  tab = table(spectralContrast,'RowNames',FTlabels);
  amt_disp(tab)
  varargout{1} = tab;
end

if flags.do_ratelevelcurves || flags.do_fig2
  
  splmax = 130; % dB
  splminplot = 5;
  fc = 4000; % Hz
  
  sig = noise(0.5*48e3,1,'white'); % 100-ms Gaussian white noise burst
  spl = 0:10:splmax;
  
  mp = zeros(1,length(spl),2,3);
  for ii = 1:length(spl)
    mp(:,ii,:,:) = baumgartner2016_spectralanalysis(cat(3,sig,sig),spl(ii),...
        'target','ID','','Condition','ratelevelcurves','fiberTypes',1:3,...
        'ftopt','cohc',1,'flow',fc,'fhigh',fc,'nf',1);
  end
  
  % Plot Rate-intensity curves
  if flags.do_plot
    figure
    sty = {'kv-','ko-','k^-'};
    for ft = 1:3
      h(ft) = plot(spl,mean(mp(:,:,1,ft),1),sty{ft});
      hold on
      xlabel('SPL (dB)','FontSize',kv.FontSize)
      ylabel('Firing rate (spikes/s)','FontSize',kv.FontSize)
    end
    set(h,'MarkerFaceColor','k')
    leg = legend(h,{'low-SR','med-SR','high-SR'});
    set(leg,'Location','northoutside','FontSize',kv.FontSize,'Orientation','horizontal')
    
    axis([splminplot,splmax-5,-30,369])
    XTick = round(splminplot/10)*10:10:splmax;
    XTickLabel = num2cell(XTick);
    XTickLabel(2:2:end) = {' '};
    set(gca,'XTick',XTick,'XTickLabel',XTickLabel,'FontSize',kv.FontSize)
  
  end
  
end

if flags.do_fig5
  exp_baumgartner2016('numchan','FontSize',kv.FontSize,'MarkerSize',4,flags.cachemode);
  singleFig(1) = gcf;
  axNumChan = get(gcf,'Children');
  exp_baumgartner2016('spatstrat','FontSize',kv.FontSize,'MarkerSize',4,flags.cachemode);
  singleFig(2) = gcf;
  axSpatStrat = get(gcf,'Children');

  % Adjust marker symbols
%   RB: set([allchild(axNumChan(1)),allchild(axSpatStrat(1))],'Marker','d') % QE
%   RB: set([allchild(axNumChan(3)),allchild(axSpatStrat(3))],'Marker','s') % PE

  % Combined plot
  fig = figure;
  ha = tight_subplot(2,2,0,[.1,.05],[.11,.02]);

  % Data
  copyobj(allchild(axNumChan(1)),ha(1))
  copyobj(allchild(axNumChan(3)),ha(3))
  copyobj(allchild(axSpatStrat(1)),ha(2))
  copyobj(allchild(axSpatStrat(3)),ha(4))

  % Labels
  QElabel = '% Quadrant errors';
  PElabel = 'Local RMS error (deg)';
  title(ha(1),'Goupell et al. (2010)')
  title(ha(2),'Majdak et al. (2013)')
  ylabel(ha(1),QElabel,'FontWeight','bold')
  ylabel(ha(3),PElabel,'FontWeight','bold')
  xlabel(ha(3),'Num. of channels','FontWeight','bold')
  xlabel(ha(4),'Spectral modification','FontWeight','bold')
  legend(ha(3),{'Model','Actual'},'Position',[0.5 0.80 0.1147 0.0440]);

  % Limits
  set(ha([1,3]),'XLim',get(axNumChan(1),'XLim'));
  set(ha([2,4]),'XLim',get(axSpatStrat(1),'XLim'));
  set(ha([1,2]),'YLim',[1,45.9])
  set(ha([3,4]),'YLim',[26,54])

  % Ticks
  set(ha,'TickLength',[0.02,.01],'Box','on')
  set(ha(1:2),'YTick',get(axNumChan(1),'YTick'))
  set(ha(3:4),'YTick',get(axNumChan(3),'YTick'))
  set(ha(1),'YTickLabel',get(axNumChan(1),'YTickLabel'))
  set(ha(3),'YTickLabel',get(axNumChan(3),'YTickLabel'))
  set(ha([1,3]),'XTick',get(axNumChan(1),'XTick'))
  set(ha([2,4]),'XTick',get(axSpatStrat(1),'XTick'))
  set(ha(3),'XTickLabel',get(axNumChan(1),'XTickLabel'))
  set(ha(4),'XTickLabel',get(axSpatStrat(1),'XTickLabel'))
  
  close(singleFig)
end

if flags.do_sensitivity || flags.do_fig8
  
  [data,meta] = exp_baumgartner2016('dynrangecheck','dprime','SPLset',60,...
    'FontSize',kv.FontSize,flags.cachemode,'marg_w',[.15,.01],'marg_h',[.05,.05]);
  dprime60dB = squeeze(data(meta(1).data == 60,[4,1:3],:))';
  
  [~,~,qe,errmeta] = exp_baumgartner2016('impairment','QE','SPLset',60,'noplot','nostat');
  [correlation_QE.r,correlation_QE.p] = corrcoeff(mean(qe,1)',dprime60dB(:));
  amt_disp('Correlation between dprime and quadrant errors:')
  amt_disp(correlation_QE)
  
  [~,~,pe] = exp_baumgartner2016('impairment','PE','SPLset',60,'noplot','nostat');
  [correlation_PE.r,correlation_PE.p] = corrcoeff(mean(pe,1)',dprime60dB(:));
  amt_disp('Correlation between dprime and local RMS errors:')
  amt_disp(correlation_PE)
  
end

if flags.do_dynrangecheck
  
  splmax = 130; % dB
  splHRTF = kv.SPLset; % dB (SPL range of targets)
  splminplot = 15;
  
  sig = noise(0.1*48e3,1,'white'); % 100-ms Gaussian white noise burst
  spl = 0:10:splmax;
  cohc = kv.cOHCset;
  
  panels = {'C_{OHC} = 1',['C_{OHC} = ' num2str(cohc(2),'%1.1f')],['C_{OHC} = ' num2str(cohc(3),'%1.1f')],'C_{OHC} = 0'};
  if flags.do_separate
    labels = {'low-SR','med-SR','high-SR','all SRs'};
  else
    labels = {'LMH','MH','H'};
  end
  
  figure
  ha = tight_subplot(4,1,kv.gap,kv.marg_h,kv.marg_w);
  if flags.do_ratelevel
    Y = nan(length(spl),4,length(cohc));
  else
    Y = nan(length(spl)-1,4,length(cohc));
  end
  for cc = 1:length(cohc)
    mp = zeros(model.kv.nf,length(spl),2,3);
    if flags.do_separate
      for ii = 1:length(spl)
        [mp(:,ii,:,:),fc] = baumgartner2016_spectralanalysis(cat(3,sig,sig),spl(ii),...
            'target','ID','','Condition','dynrangecheck','fiberTypes',1:3,'ftopt','cohc',cohc(cc));
      end
    else
      fiberTypes = {1:3,2:3,3};
      for ft = 1:length(fiberTypes)
        for ii = 1:length(spl)
          [mp(:,ii,:,ft),fc] = baumgartner2016_spectralanalysis(cat(3,sig,sig),spl(ii),...
              'target','dynrangecheck','fiberTypes',fiberTypes{ft},'cohc',cohc(cc));
        end
      end
    end

    splplot = spl;
    if flags.do_dynrangeDiff
      mp = diff(mp,1,2)/mean(diff(spl));
      splplot = spl(2:end);
    end
    
    if flags.do_dprime
      sd = 2.6*mp.^0.34;
      sdDiff = sqrt(sd(:,1:length(spl)-1,:,:).^2 + sd(:,2:length(spl),:,:).^2);
      mp = diff(mp,1,2)./sdDiff;
      splplot = spl(2:end);
    end

    axes(ha(cc))
    
    % target HRTF range
    color = {.8*ones(1,3)};
    if length(splHRTF) == 2
      color = {.9*ones(1,3),color{1}};
    end
    for ll = 1:length(splHRTF)
      a(1) = area(splHRTF(ll)+[-10,10],[1e3,1e3],'EdgeColor',ones(1,3));
      hold on
      a(2) = area(splHRTF(ll)+[-10,10],-[1e3,1e3],'EdgeColor',ones(1,3));
      set(a,'FaceColor',color{ll})
    end
    set(gca,'Layer','top')
    
    % Rate-intensity curves
    sty = {'kv-','ko-','k^-'};
    for ft = 1:3
      Y(:,ft,cc) = mean(mp(:,:,1,ft),1);
      h(ft) = plot(splplot,Y(:,ft,cc),sty{ft});
      hold on
      if cc == length(cohc)
        xlabel('SPL (dB)','FontSize',kv.FontSize,'FontWeight','bold')
      end
      if flags.do_dynrangeDiff
        text(splminplot+5,9,panels{cc},'FontSize',kv.FontSize)
        ylabel({'Rate difference (spikes/s/dB)'},'FontSize',kv.FontSize)
      elseif flags.do_dprime
        text(splminplot+5,5,panels{cc},'FontSize',kv.FontSize)
        if cc==3
          text(0,6,{'d'''},'FontSize',kv.FontSize,'FontWeight','bold')
        end
      else
        text(splminplot+5,330,panels{cc},'FontSize',kv.FontSize)
        ylabel('Firing rate (spikes/s)','FontSize',kv.FontSize)
      end
    end
    % All SRs combined
    ftd = [0.16,0.23,0.61]; % Liberman (1978)
    Y(:,4,cc) = Y(:,1:3,cc) * ftd';
    h(4) = plot(splplot,Y(:,4,cc),'k--');
    
    set(h,'MarkerFaceColor','k')
    if cc == 1%length(cohc)
      leg = legend(h,labels);
      set(leg,'Location','northoutside','FontSize',kv.FontSize,'Orientation','vertical')
      set(leg,'Position',[.4,.96,.33,.03])
    end
    
    if flags.do_dynrangeDiff
      plot([0,splmax],[0,0],'k--') % sensitivity threshold
      axis([splminplot,splmax-5,-1,10.5])
    elseif flags.do_dprime
      axis([splminplot,splmax-5,-1.5,5.9])  
    else
      axis([splminplot,splmax-5,-20,399])
    end
    XTick = 20:10:130;
    XTickLabel = num2cell(XTick);
    XTickLabel(2:2:end) = {' '};
    set(gca,'XTick',XTick,'XTickLabel',XTickLabel,'FontSize',kv.FontSize)
  
  end
  meta(1).name = 'SPL';
  meta(1).data = splplot;
  meta(2).name = labels;
  meta(2).data = {1,2,3,1:3};
  meta(3).name = 'OHC gain';
  meta(3).data = cohc;
  varargout = {Y;meta};
  
end
  
if flags.do_localevel
  
  dlat = 30; % deg
  condition = data_baumgartner2015('ConditionNames');
  
  cachename = ['localevel_' cachename];
  if not(model.kv.gammashortfact == 1)
    cachename = [cachename '_gsf' num2str(model.kv.gammashortfact,'%1.1f')];
  end
  if not(model.kv.Sshortfact == 1)
    cachename = [cachename '_Ssf' num2str(model.kv.Sshortfact,'%1.1f')];
  end
  if not(model.kv.psgeshort == 1)
    cachename = [cachename '_psgec' num2str(model.kv.psgeshort,'%1.1f')];
  end
  if model.flags.do_SPLtemAdapt
    cachename = [cachename '_SPLtemAdapt'];
  end
  if model.flags.do_gammatone
    cachename = [cachename '_minSPL' num2str(model.kv.GT_minSPL) '_maxSPL' num2str(model.kv.GT_maxSPL)];
  end
  [pred,ref] = amt_cache('get', cachename, flags.cachemode);
  if isempty(pred)
    pred.qe = nan(7,length(condition));
    pred.pe = pred.qe;
    ref.qe = pred.qe;
    ref.pe = pred.qe;
    for cc = 1:length(condition)

      data = data_baumgartner2016(condition{cc},'model','argimport',model.flags,model.kv);

      for ll = 1:length(data)

        latresp = data(ll).itemlist(:,7);
        idlat = latresp <= dlat & latresp > -dlat;
        itemlist = data(ll).itemlist(idlat,:);
        exptang = itemlist(:,6);
        exprang = itemlist(:,8);

        ref.pe(ll,cc) = real(localizationerror(itemlist,'rmsPmedianlocal'));
        ref.qe(ll,cc) = real(localizationerror(itemlist,'querrMiddlebrooks'));

        if strcmp(condition{cc},'Long')
          clbl = condition{cc};
        else % short
          clbl = '';
        end
        [err,pmv] = baumgartner2016(data(ll).Obj,data(ll).Obj,...
            'argimport',model.flags,model.kv,'ID',data(ll).id,'Condition',clbl,'S',data(ll).S,...
            'stim',data(ll).stim,'fsstim',data(ll).fsstim,'SPL',data(ll).SPL,...
            'QE_PE_EB','exptang',exptang,'priordist',data(ll).priordist);
        pred.pmv{ll,cc} = pmv;
        pred.qe(ll,cc) = err.qe;
        pred.pe(ll,cc) = err.pe;

      end
    end
    amt_cache('set',cachename,pred,ref);
  else
    data = data_baumgartner2016('all');
  end
  
  perrmtx = ref.pe';
  pred.perrmtx = pred.pe';
  qerrmtx = ref.qe';
  pred.qerrmtx = pred.qe';
  
  Nsub = length(data);
  
  %% Prediction residues 
  mm = 1;
  r_perr(mm) = corrcoeff([pred.pe],[ref.pe]);
  e_perr(mm) = mean(rms([pred.pe]-[ref.pe]));

  r_qerr(mm) = corrcoeff([pred.qe],[ref.qe]);
  e_qerr(mm) = mean(rms([pred.qe]-[ref.qe]));
    
  amt_disp(' e_PE  r_PE  e_QE   r_QE')
  amt_disp([num2str(e_perr(mm),'%2.1f') '\deg  ' num2str(r_perr(mm),'%2.2f') '  ' num2str(e_qerr(mm),'%2.1f') '%  ' num2str(r_qerr(mm),'%2.2f')])
  
  varargout{1} = 0.5* (e_perr(mm)/90 + e_qerr(mm)/100);


  %% Plots
  if flags.do_plot
    
    if flags.do_performance
    
      LineWidth = 1;

      symb = {...
            'rs-';...
            'bd-';...
            'gh-';...
            'mv-';...
            'c*-';...
            };
      symbExp = 'ko--';

      name{mm} = 'Pred.';
      % Legend
      legendentry = {'Actual, long';'Actual, short'};
      legendentry{2+2*mm-1} = [name{mm} ', long'];
      legendentry{2+2*mm} = [name{mm} ', short'];

      xval = 10:10:70;
      xtext = 10;

      % Local central RMS error
      hfig = figure;
      ha = tight_subplot(8,2,kv.gap,kv.marg_h,kv.marg_w);
      for ii=1:Nsub

        axes(ha(2*ii));
        hlong = plot(50,perrmtx(1,ii),symbExp(1:2));
        set(hlong,'MarkerFaceColor',symbExp(1),'LineWidth',LineWidth)

        set(gca,'YAxisLocation','right','YMinorTick','on')
        set(gca,'TickLength',kv.TickLength)

        hold on
        hshort = plot(xval,perrmtx(2:end,ii),symbExp);
        set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)

        for mm = 1:length(pred)
          hlong = plot(50,pred(mm).perrmtx(1,ii),symb{mm}(1:2));
          set(hlong,'MarkerFaceColor',symb{mm}(1),'LineWidth',LineWidth)
          hshort = plot(xval,pred(mm).perrmtx(2:end,ii),symb{mm});
          set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
        end

        if ii==1
          title('Local error (deg)','FontSize',kv.FontSize)
        end

        axis([5,75,26,54])
      end
      % Pooled
      axes(ha(2*ii+2));
      hlong = errorbar(50,mean(perrmtx(1,:)),std(perrmtx(1,:)),symbExp(1:2));
      set(hlong,'MarkerFaceColor',symbExp(1),'LineWidth',LineWidth)
      set(gca,'YAxisLocation','right','YMinorTick','on','TickLength',kv.TickLength)
      hold on
      hshort = errorbar(xval,mean(perrmtx(2:end,:),2),std(perrmtx(2:end,:),1,2),symbExp);
      set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
      hlong = errorbar(50,mean(pred(mm).perrmtx(1,:),2),std(pred(mm).perrmtx(1,:),1,2),symb{mm}(1:2));
      set(hlong,'MarkerFaceColor',symb{mm}(1),'LineWidth',LineWidth)
      hshort = errorbar(xval,mean(pred(mm).perrmtx(2:end,:),2),std(pred(mm).perrmtx(2:end,:),1,2),symb{mm});
      set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
      axis([5,75,26,54])

      xlabel('SL (dB)','FontSize',kv.FontSize)


      % Quadrant error
      for ii=1:Nsub
        axes(ha(2*ii-1));

        hlong = plot(50,qerrmtx(1,ii),symbExp(1:2));

        set(gca,'YMinorTick','on')
        set(gca,'TickLength',kv.TickLength)

        set(hlong,'MarkerFaceColor',symbExp(1),'LineWidth',LineWidth)
        hold on
        hshort = plot(xval,qerrmtx(2:end,ii),symbExp);
        set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)

        for mm = 1:length(pred)
          hlong = plot(50,pred(mm).qerrmtx(1,ii),symb{mm}(1:2));
          set(hlong,'MarkerFaceColor',symb{mm}(1),'LineWidth',LineWidth)
          hshort = plot(xval,pred(mm).qerrmtx(2:end,ii),symb{mm});
          set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
        end

        if ii==1
          title('Quadrant error (%)','FontSize',kv.FontSize)
        end

        axis([5,75,-4,54])

        % Listener ID
        ylabel([data(ii).id '            '],'FontSize',kv.FontSize,'Rotation',0,'FontWeight','bold');
      end
      % Pooled
      axes(ha(2*ii+1));
      hlong = errorbar(50,mean(qerrmtx(1,:)),std(qerrmtx(1,:)),symbExp(1:2));
      set(hlong,'MarkerFaceColor',symbExp(1),'LineWidth',LineWidth)
      set(gca,'YMinorTick','on','TickLength',kv.TickLength)
      hold on
      hshort = errorbar(xval,mean(qerrmtx(2:end,:),2),std(qerrmtx(2:end,:),1,2),symbExp);
      set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
      hlong = errorbar(50,mean(pred(mm).qerrmtx(1,:),2),std(pred(mm).qerrmtx(1,:),1,2),symb{mm}(1:2));
      set(hlong,'MarkerFaceColor',symb{mm}(1),'LineWidth',LineWidth)
      hshort = errorbar(xval,mean(pred(mm).qerrmtx(2:end,:),2),std(pred(mm).qerrmtx(2:end,:),1,2),symb{mm});
      set(hshort,'MarkerFaceColor','w','LineWidth',LineWidth)
      axis([5,75,-4,54])
      ylabel('Pooled            ','FontSize',kv.FontSize,'Rotation',0,'FontWeight','bold');
      xlabel('SL (dB)','FontSize',kv.FontSize)

    end
    
    if flags.do_pmv
      
      conditions = data_baumgartner2015('ConditionNames');
      idcond = [1,2,5,8]; %{'Long','10dB','30dB','50dB','70dB'};
      Nc = length(idcond);
      
      hfig = figure;
      ha = tight_subplot(Nsub,Nc,kv.gap,kv.marg_h,kv.marg_w);
      
      for cc=1:Nc
        idc = idcond(cc);
        data = data_baumgartner2015(conditions{idc});
        
        for ll=1:Nsub
          axes( ha(cc + Nc*(ll-1)) )
          plot_baumgartner2014(pred.pmv{ll,idc}.p,pred.pmv{ll,idc}.tang,pred.pmv{ll,idc}.rang,...
              data(ll).itemlist(:,6),data(ll).itemlist(:,8),...
              'MarkerSize',kv.MarkerSize/2,'nocolorbar','cmax',0.05)
            
          % Labels
          set(gca,'FontSize',kv.FontSize-1)
          if not(cc==1)
            set(gca,'YTickLabel',[])
            ylabel('')
          elseif not(ll==4)
              ylabel('')
          end
          if not(ll==Nsub)
            set(gca,'XTickLabel',[])
            xlabel('')
          else
            set(gca,'XTickLabelRotation',90)
            xlabel('')
          	if cc==3
              xlabel({' ';'Target Angle (deg)                                 '},'FontSize',kv.FontSize-1)
            end
          end
          if ll==1
            title(conditions{idc},'FontSize',kv.FontSize)
          end
          if cc==Nc
            set(gca,'YAxisLocation','right')
            ylabel(['            ' data(ll).id],'FontSize',kv.FontSize,'FontWeight','bold','Rotation',0);
          end
            
        end
        
      end
      
    end
  
  end
end

end


%% ------------------------------------------------------------------------
%  ---- INTERNAL FUNCTIONS ------------------------------------------------
%  ------------------------------------------------------------------------
function hM_warped = warp_hrtf(hM,fs)
% warps HRTFs acc. to Walder (2010)
% Usage: hM_warped = warp_hrtf(hM,fs)

N = fs;
fu = 2800;
fowarped = 8500;
fo = 16000;

fscala = [0:fs/N:fs-fs/N]';
hM_warped = zeros(512,size(hM,2),size(hM,3));
fuindex = max(find(fscala <= fu)); % 2800
fowindex = min(find(fscala >= fowarped));
foindex = min(find(fscala >= fo));

for canal = 1:size(hM,3)
   for el = 1:size(hM,2)
        yi = ones(fs/2+1,1)*(10^-(70/20));
        flin1 = [fscala(1:fuindex-1)];
        flin2 = [linspace(fscala(fuindex),fscala(foindex),fowindex-fuindex+1)]';
        fscalawarped = [flin1 ; flin2];

        % interpolate
        x = fscala(1:foindex);
        H = fft(hM(:,el,canal),N);
        Y = H(1:foindex);
        xi = fscalawarped;
        yi(1:length(xi),1) = interp1(x, Y, xi,'linear');

        yges=([yi; conj(flipud(yi(2:end-1)))]);
        hges = ifft([yges(1:end)],length(yges));
        hges=fftshift(ifft(yges));
        hwin=hges(fs/2-256:fs/2+768);
        hwinfade = FW_fade(hwin,512,24,96,192);
        hM_warped(1:end,el,canal)=hwinfade;
            
    end
end
end

function [syncrnfreq, GETtrain] = GETVocoder(filename,in,channum,lower,upper,alpha,GaussRate,stimpar)
warning('off')
% channel/timing parameters
srate=stimpar.SamplingRate;
nsamples=length(in); % length of sound record
duration=nsamples/srate; % duration of the signal (in s)
t=0:1/srate:duration;
t=t(1:nsamples);
extendedrange = 0;
if alpha == -1 % Log12ER case
    extendedrange = 1;
    alpha = 0.28;
elseif alpha == 0 % use predefined alphas
  switch channum
    case 3
      alpha=1.42;
    case 6
      alpha=0.67;
    case 9
      alpha=0.45;
    case 12
      alpha=0.33;
    case 18
      alpha=0.22;
    case 24
      alpha=0.17;
    otherwise 
      error(['Alpha not predefined for channels number of ' num2str(channu)]);
  end
end

% Synthesis: These are the crossover frequencies that the output signal is mapped to
crossoverfreqs = logspace( log10(lower), log10(upper), channum + 1);
if extendedrange == 1
    crossoverfreqs = [300,396,524,692,915,1209,1597,2110,2788,4200,6400,10000,16000];
end
syncrnfreq(:,1)=crossoverfreqs(1:end-1);
syncrnfreq(:,2)=crossoverfreqs(2:end);
for i=1:channum
    cf(i) = sqrt( syncrnfreq(i,1)*syncrnfreq(i,2) );
end

% Pulse train parameters
Gamma = alpha*cf;  
if extendedrange == 1
    Gamma(9) = 1412;
    Gamma(10) = 2200;
    Gamma(11) = 3600;
    Gamma(12) = 6000;
end
N = ceil(duration*GaussRate); % number of pulses
Genv=zeros(N,nsamples);
GETtrain=zeros(channum,nsamples);

% Make pulse trains
for i = 1:channum
    Teff = 1000/Gamma(i);
    if Teff > 3.75
    % if modulation depth is not 100%, make pulse train then modulate
        for n = 1:N
            % delay pulses by half a period so first Gaussian pulse doesn't
            % start at a max
            T = (n-0.5)/N*duration;
            Genv(n,:) = sqrt(Gamma(i)) * exp(-pi*(Gamma(i)*(t-T)).^2);
        end
        Genv_train(i,:) = sum(Genv);
        %modulate carrier
        GETtrain(i,:) = Genv_train(i,:) .* sin(2*pi*cf(i)*t);
        %normalize energy
        Energy(i) = norm(GETtrain(i,:))/sqrt(length(t));
%         Energy(i) = rms(GETtrain(i,:)); % !!!!!!!!!!!!
        GETtrain(i,:) = GETtrain(i,:)/Energy(i);
    else
    % if modulation depth is 100%, make modulated pulses and replicate
        T=(0.5)/N*duration;
        Genv=zeros(N,nsamples);
        Genv(1,:) = sqrt(Gamma(i)) * exp(-pi*(Gamma(i)*(t-T)).^2) .* sin(2*pi*cf(i)*t - T + pi/4); %!!! (t-T)
        Genv=repmat(Genv(1,:),[N 1]);
        for n=1:N
            T = round((n)/N*nsamples);
            Genv(n,:)=circshift(Genv(n,:),[1 T-1]);
        end
        GETtrain(i,:) = sum(Genv);
        %normalize energy
        Energy(i) = norm(GETtrain(i,:))/sqrt(length(t));
%         Energy(i) = rms(GETtrain(i,:)); % !!!!!!!!!!!!
        GETtrain(i,:) = GETtrain(i,:)/Energy(i);
    end
end

end

function out=channelize(fwavout, h, h0, in, channum, corners, syncrnfreq, ...
                        GETtrain, stimpar, amp, fadein, fadeout)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *** v1.2.0                                                           %
% GET Vocoder scaled by energy, not envelope                           %
%                                                                      %
% *** v1.1.0                                                           %
% Added Gaussian Envelope Tone (GET) Vocoder                           %
% GET pulse train is generated and passed to this function             %
% M. Goupell, May 2008                                                 %
%                                                                      %
% *** v1.0.0                                                           %
% Modified from ElecRang/matlab/makewav.m  v1.3.1                      %
% To be used with Loca by M. Goupell, Nov 2007                         % 
% Now program receives a sound rather than reading a speech file       %
% Rewrote according to the specifications (PM, Jan. 2008)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = hrtf
% h0 = hrtf for reference position
% in = reference noise
% noise = number of channels of noise vectors

srate=stimpar.SamplingRate;
N=length(in);               % length of sound record
d=0.5*srate;                % frequency scalar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read corner frequencies for channels                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis: These are the crossover frequencies that the output signal is mapped to
% calculated now in GETVocoder, passed to this function

% Analysis: These are the crossover frequencies that the input signal is subdivided by
if length(corners)<=channum
  error('You need at least one more corner frequency than number of channels');
end
anacrnfreq(:,1)=corners(1:end-1);
anacrnfreq(:,2)=corners(2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis: Filter Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inX=fftfilt(h,in);

order=4;     %order of butterworth filter
% out=zeros(1,N);
% filtX=zeros(channum,N);

for i=1:channum
   [b, a]=butter(order, anacrnfreq(i,:)/d);
   out=filter(b, a, inX);  % bandpass-filtering
%    E(i)=norm(out,1)/sqrt(N);
   E(i) = rms(out);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% each channel
for i=1:channum
    outX(i,:) = GETtrain(i,:) * E(i);    
end
% sum me up scotty
out=sum(outX,1);

% for i = 1:channum
%     subplot(channum/3,3,i)
%     plot(outX(i,(1:480)))
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save file                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out=out*(10^(amp/20))/sqrt(sum(out.^2))*sqrt(sum(in.^2))*sqrt(sum(h.^2))/sqrt(sum(h0.^2));
% disp(20*log10(sqrt(sum(out.^2))));
ii=max(max(abs(out)));
if ii>=1
  error(['Maximum amplitude value is ' num2str(20*log10(ii)) 'dB. Set the HRTF scaling factor lower to avoid clipping']);
end
out=FW_fade(out,0,fadein,fadeout);
% wavwrite(out,srate,stimpar.Resolution,fwavout);
end

function out = FW_fade(inp, len, fadein, fadeout, offset)
% FW_FADE crop/extend and fade in/out a vector.
%
% OUT = FW_FADE(INP, LEN, FADEIN, FADEOUT, OFFSET) crops or extends with zeros the signal INP
% up to length LEN. Additionally, the result is faded in/out using HANN window with 
% the length FADEIN/FADEOUT, respectively. If given, an offset can be added to show
% where the real signal begins.
%
% When used to crop signal, INP is cropped first, then faded out. 
% When used to extend signal, INP is faded out first, then extended too provide fading.
% 
% INP:     vector with signal (1xN or Nx1)
% LEN:     length of signal OUT (without OFFSET)
% FADEIN:  number of samples to fade in, beginning from OFFSET
% FADEOUT: number of samples to fade out, ending at the end of OUT (without OFFSET)
% OFFSET:  number of offset samples before FADEIN, optional
% OUT:     cropped/extended and faded vector
% 
% Setting a parameter to 0 disables corresponding functionality.

% ExpSuite - software framework for applications to perform experiments (related but not limited to psychoacoustics).
% Copyright (C) 2003-2010 Acoustics Research Institute - Austrian Academy of Sciences; Piotr Majdak and Michael Mihocic
% Licensed under the EUPL, Version 1.1 or ? as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence. 
% You may obtain a copy of the Licence at: http://ec.europa.eu/idabc/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% 7.11.2003
% 22.08.2005: improvement: INP may be 1xN or Nx1 now.
% Piotr Majdak (piotr@majdak.com)

ss=length(inp);    % get length of inp
	% offset
if ~exist('offset','var')
	offset = 0;
end
offset=round(offset);
len=round(len);
fadein=round(fadein);
fadeout=round(fadeout);
if offset>ss
	error('OFFSET is greater than signal length');
end
	% create new input signal discarding offset
inp2=inp(1+offset:end);
ss=length(inp2);

  % fade in
if fadein ~= 0 
  if fadein > ss
    error('FADEIN is greater than signal length');
  end
  han=hanning(fadein*2);
  if size(inp2,1)==1
    han=han';
  end
  inp2(1:fadein) = inp2(1:fadein).*han(1:fadein);
end
  % fade out window
  
if len == 0    
  len = ss;
end
if len <= ss
    % crop and fade out
  out=inp2(1:len);
    % fade out
  if fadeout ~= 0
    if fadeout > len
      error('FADEOUT is greater than cropped signal length');
    end
    han = hanning(2*fadeout);
    if size(out,1)==1
      han=han';
    end
    out(len-fadeout+1:end)=out(len-fadeout+1:end).*han(fadeout+1:end);
  end
else
    % fade out and extend
  if fadeout ~= 0
    if fadeout > ss
      error('FADEOUT is greater than signal length');
    end
    han = hanning(2*fadeout);
    inp2(ss-fadeout+1:end)=inp2(ss-fadeout+1:end).*han(fadeout+1:end);
  end
  out = [zeros(offset,1); inp2; zeros(len-ss,1)];
end
end

function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

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

function s = gain2slope(g)
% s = gain2slope(g)

s = rad2deg(acos(1./sqrt(g.^2+1)));
end

function [r,p] = corrcoeff(x,y)
% internal function to evaluate correlation coefficients in order to avoid 
% dependence on MATLAB statistics toolbox
%
%   Usage: r = corrcoeff(x,y)

if not(size(x) == size(y))
  error('corrcoeff: x and y must have same size!')
end

try
  [r,p] = corrcoef(x,y);
  r = r(2);
  p = p(2);
catch
  r = cov(x(:),y(:))/std(x(:))/std(y(:));
  r = r(2);
  n = length(y(:));
  df = n-2;
  t = r/sqrt((1-r^2)/df);
  p = 2*(1-tcdf(t,df));
end

end
