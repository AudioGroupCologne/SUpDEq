function varargout = exp_baumgartner2014(varargin)
%EXP_BAUMGARTNER2014 Results from Baumgartner et al. (2014)
%   Usage: data = exp_baumgartner2014(flag) 
%
%   EXP_BAUMGARTNER2014(flag) reproduces figures of the study from 
%   Baumgartner et al. (2014).
%
%   The following flags can be specified
%
%     'fig2'    Reproduce Fig.2:
%               Binaural weighting function best fitting results from 
%               Morimoto (2001) labeled as [1] and Macpherson and Sabin (2007) 
%               labeled as [2] in a least squared error sense.
%
%     'fig3'    Reproduce Fig.3:
%               Prediction examples. Actual responses and response predictions 
%               for three exemplary listeners when listening to median-plane 
%               targets in the baseline condition. 
%               Actual response angles are shown as open circles. 
%               Probabilistic response predictions are encoded by brightness 
%               according to the color bar to the right. Actual (A:) and 
%               predicted (P:) quadrant error rates (QE) and local polar 
%               RMS errors (PE) are listed above each panel.
%               
%     'fig4'    Reproduce Fig.4:
%               Model parametrization. Partial and joint prediction residues 
%               as functions of the degree of selectivity and the motoric 
%               response scatter. Residuum functions are normalized to the 
%               minimum residuum obtained for the optimal parameter value.
%
%     'fig5'    Reproduce Fig.5:
%               Effect of band limitation and spectral warping. Actual 
%               responses and response predictions for listener NH12 when 
%               listening to broadband (BB), low-pass filtered (LP), or 
%               spectrally warped (W) DTFs of the median plane. Data were 
%               pooled within pm15^circ of lateral angle.
%               All other conventions are as in Fig.3.
%
%     'fig6'    Reproduce Fig.6:
%               Effect of band limitation and spectral warping. 
%               Listeners were tested with broadband (BB), low-pass 
%               filtered (LP), and spectrally warped (W) DTFs. 
%               Actual: experimental results from Majdak et al. (2013). 
%               Part.: Model predictions for the actual eight participants 
%               based on the actually tested target positions. Pool: Model 
%               predictions for our listener pool based on all possible 
%               target positions. Symbols and whiskers show median values 
%               and inter-quartile ranges, respectively. Symbols were 
%               horizontally shifted to avoid overlaps. Dotted horizontal 
%               lines represent chance rate. Correlation coefficients, r, 
%               and prediction residues, e, specify the correspondence 
%               between actual and predicted listener-specific performances.
%
%     'fig7'    Reproduce Fig.7:
%               Effect of spectral resolution in terms of varying the number 
%               of spectral channels of a channel vocoder. Actual responses 
%               and response predictions for exemplary listener NH12. 
%               Results for 24, 9, and 3 channels are shown. All other 
%               conventions are as in Fig.3.
%
%     'fig8'    Reproduce Fig.8:
%               Effect of spectral resolution in terms of varying the number 
%               of spectral channels of a channel vocoder. Actual experimental 
%               results are from Goupell et al. (2010). Stimulation with broadband 
%               click trains (CL) represents an unlimited number of channels. 
%               All other conventions are as in Fig.6.
%
%     'fig9'    Reproduce Fig.9:
%               Effect of non-individualized HRTFs in terms of untrained 
%               localization with others' instead of own ears. Statistics 
%               summaries with open symbols represent actual experimental 
%               results replotted from Fig.,13 of Middlebrooks (1999), 
%               statistics with filled symbols represent predicted results.
%               Horizontal lines represent 25th, 50th, and 75th percentiles, 
%               the whiskers represent 5th and 95th percentiles, and crosses 
%               represent minima and maxima. Circles and squares denote mean values.
%               Dimensions of 
%
%     'fig10'   Reproduce Fig.10:
%               Effect of spectral ripples. Actual experimental results (circles) 
%               are from Macpherson and Middlebrooks (2003). Predicted results (filled circles) 
%               were modeled for our listener pool (squares). Either the ripple 
%               depth of 40,dB (top) or the ripple density of one ripple/octave 
%               (bottom) was kept constant. Ordinates show the listener-specific 
%               difference in error rate between a test and the baseline condition. 
%               Baseline performances are shown in the bottom right panel.
%               Symbols and whiskers show median values and inter-quartile ranges, 
%               respectively. Symbols were horizontally shifted to avoid overlaps. 
%               Diamonds with dashed lines show predictions (P) of the model 
%               without positive spectral gradient extraction (PSGE).
%
%     'fig11'   Reproduce Fig.11:
%               Effect of high-frequency attenuation in speech localization. 
%               Actual experimental results are from Best et al. (2005). 
%               Absolute polar angle errors (top) and QE (bottom) were averaged 
%               across listeners. Circles and squares show actual and predicted 
%               results, respectively. Diamonds with dashed lines show predictions 
%               of the model without positive spectral gradient extraction.
%
%     'fig12'   Reproduce Fig.12:
%               Listener-specific likelihood statistics used to evaluate 
%               target-specific predictions for baseline condition. Bars 
%               show actual likelihoods, dots show mean expected likelihoods, 
%               and whiskers show tolerance intervals with 99% confidence 
%               level of expected likelihoods.
%
%     'fig13'   Reproduce Fig.13:
%               Exemplary baseline predictions. Same as Fig.3 but for listeners 
%               where actual likelihoods were outside the tolerance intervals.    
%
%     'fig14'   Reproduce Fig.14:
%               Baseline performance as a function of the magnitude of the 
%               lateral response angle. Symbols and whiskers show median 
%               values and inter-quartile ranges, respectively. Open symbols 
%               represent actual and closed symbols predicted results. Symbols 
%               were horizontally shifted to avoid overlaps. Triangles with 
%               dashed lines show predictions (P) of the model without the 
%               sensomotoric mapping (SMM) stage.  
%
%     'tab1'    Reproduce Tab.1:
%               Listener-specific sensitivity calibrated on the basis of N 
%               baseline targets in proximity of the median plane (+-30deg). 
%               Listeners are labeled as NHl. Actual and predicted quadrant 
%               errors (QE) and local polar RMS errors (PE) are shown pairwise 
%               (Actual | Predicted).
%
%     'tab2'    Reproduce Tab.2:
%               The effects of model configurations on the prediction residues. 
%               PSGE: model with or without positive spectral gradient extraction. 
%               MBA: model with or without manual bandwidth adjustment to the 
%               stimulus bandwidth. Prediction residues between actual and 
%               predicted PE and QE are listed for acute performance with 
%               the broadband (BB), low-passed (LP) and warping (W) conditions 
%               of the experiments from Majdak et al. (2013).
%
%     'tab3'    Reproduce Tab.3:
%               Performance predictions for binaural, ipsilateral, and 
%               contralateral listening conditions. The binaural weighting 
%               coefficient was varied in order to represent the three 
%               conditions: binaural: Phi = 13^circ; 
%               ipsilateral: Phi rightarrow +0^circ; 
%               contralateral: Phi rightarrow -0^circ. 
%               Prediction residues and correlation coefficients between 
%               actual and predicted results are shown together with predicted 
%               average performances. 
%
%
%   Further, cache flags (see amt_cache) and plot flags can be specified:
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'no_plot'  Don't plot, only return data.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v0.4.3 or higher from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in auxdata/baumgartner2014
%
%   3) Statistics Toolbox for Matlab (for some of the figures)
%
%   Examples:
%   ---------
%
%   To display Fig.2 use :
%
%     exp_baumgartner2014('fig2');
%
%   To display Fig.3 use :
%
%     exp_baumgartner2014('fig3');
%
%   To display Fig.4 use :
%
%     exp_baumgartner2014('fig4');
%
%   To display Fig.5 use :
%
%     exp_baumgartner2014('fig5');
%
%   To display Fig.6 use :
%
%     exp_baumgartner2014('fig6');
%
%   To display Fig.7 use :
%
%     exp_baumgartner2014('fig7');
%
%   To display Fig.8 use :
%
%     exp_baumgartner2014('fig8');
%
%   To display Fig.9 use :
%
%     exp_baumgartner2014('fig9');
%
%   To display Fig.10 use :
%
%     exp_baumgartner2014('fig10');
%
%   To display Fig.11 use :
%
%     exp_baumgartner2014('fig11');
%
%   To display Fig.12 use :
%
%     exp_baumgartner2014('fig12');
%
%   To display Fig.13 use :
%
%     exp_baumgartner2014('fig13');
%
%   To display Fig.14 use :
%
%     exp_baumgartner2014('fig14');
%
%
%   See also: baumgartner2014 data_baumgartner2014
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%     M. Morimoto. The contribution of two ears to the perception of vertical
%     angle in sagittal planes. J. Acoust. Soc. Am., 109:1596--1603, 2001.
%     
%     P. Majdak, T. Walder, and B. Laback. Effect of long-term training on
%     sound localization performance with spectrally warped and band-limited
%     head-related transfer functions. J. Acoust. Soc. Am., 134:2148--2159,
%     2013.
%     
%     M. J. Goupell, P. Majdak, and B. Laback. Median-plane sound
%     localization as a function of the number of spectral channels using a
%     channel vocoder. J. Acoust. Soc. Am., 127:990--1001, 2010.
%     
%     J. C. Middlebrooks. Virtual localization improved by scaling
%     nonindividualized external-ear transfer functions in frequency. J.
%     Acoust. Soc. Am., 106:1493--1510, 1999.
%     
%     E. A. Macpherson and J. C. Middlebrooks. Vertical-plane sound
%     localization probed with ripple-spectrum noise. J. Acoust. Soc. Am.,
%     114:430--445, 2003.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_baumgartner2014.php


% #Author: Robert Baumgartner (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% ------ Check input options --------------------------------------------

definput.import={'amt_cache'};
definput.keyvals.FontSize = 12;
definput.keyvals.MarkerSize = 6;
definput.flags.type = {'missingflag','fig2','fig3','fig4','fig5','fig6',...
                       'fig7','fig8','fig9','fig10','fig11','fig12','fig13',...
                       'fig14','tab1','tab2','tab3'};%,...s
definput.flags.plot = {'plot','no_plot'};


[flags,kv]  = ltfatarghelper({'FontSize','MarkerSize'},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%% General Plot Settings
TickLength = [0.02,0.04];

%% ------ FIG 2 -----------------------------------------------------------
if flags.do_fig2
  
  % Original Data
  Morimoto_left = [1,0.5,0];
  Morimoto_right = 1-Morimoto_left;
  Morimoto_ang = [60,0,-60];

  Macpherson_left = [3.75,0.5,-4]/10+0.5;
  Macpherson_right = [-3.75,-0.5,4.75]/10+0.5;
  Macpherson_ang = [30,0,-30];

  data =  [Morimoto_left, -Macpherson_right, Macpherson_left];
  lat =   [Morimoto_ang ,  Macpherson_ang  , Macpherson_ang];

  % Fit Slope
  bwslope = 1:0.1:90;
  resid = zeros(size(bwslope));
  for ii = 1:length(bwslope)
    binw_left = 1./(1+exp(-lat/bwslope(ii))); 
    resid(ii) = rms(binw_left-data);
  end
  [~,idmin] = min(resid);
  contralateralGain = bwslope(idmin);
  amt_disp(['Phi: ' num2str(contralateralGain,'%2.0f') ' deg'],'documentation');
  
  % Calculate specific weights to plot
  lat = -90:5:90;
  binw_left = 1./(1+exp(-lat/contralateralGain)); % weight of left ear signal with 0 <= binw <= 1
  binw_right = 1-binw_left;
  
  if flags.do_plot
    figure;
    plot(lat,binw_left)
    hold on
    plot(lat,binw_right,'r--')
    plot(Morimoto_ang,Morimoto_left,'vk','MarkerSize',kv.MarkerSize,'MarkerFaceColor','w');
    plot(Macpherson_ang,Macpherson_left,'ok','MarkerSize',kv.MarkerSize,'MarkerFaceColor','w')

    plot(Morimoto_ang,Morimoto_left,'vb','MarkerSize',kv.MarkerSize,'MarkerFaceColor','w')
    plot(Morimoto_ang,Morimoto_right,'vr','MarkerSize',kv.MarkerSize,'MarkerFaceColor','w')
    plot(Macpherson_ang,Macpherson_left,'ob','MarkerSize',kv.MarkerSize,'MarkerFaceColor','w')
    plot(Macpherson_ang,Macpherson_right,'or','MarkerSize',kv.MarkerSize,'MarkerFaceColor','w')

    l = legend('\itL','\itR','[1]','[2]');
    set(l,'Location','East','FontSize',kv.FontSize-1)
    set(gca,'XLim',[lat(1) lat(end)],'YLim',[-0.05 1.05],'XTick',-60:30:60,'FontSize',kv.FontSize)
    xlabel('\phi_k (deg)','FontSize',kv.FontSize)
    ylabel('w_{\zeta}(\phi_k)','FontSize',kv.FontSize)
  end
  
  % Output
  clear data
  data.contralateralGain = contralateralGain;
  
end

%% ------ FIG 3&13 --------------------------------------------------------
if flags.do_fig3 || flags.do_fig13
  
  latseg = [-20,0,20]; ii = 2; % centers of lateral segments
%   dlat =  10;         % lateral range (+-) of each segment

  s = data_baumgartner2014('baseline',flags.cachemode);
  
  if flags.do_fig3
    idselect = ismember({s.id},{'NH15','NH22','NH62'});
  else
    idselect = ismember({s.id},{'NH12','NH39','NH18'});
  end
  s = s(idselect);

  %% LocaMo
  qe = zeros(length(s),length(latseg));
  pe = zeros(length(s),length(latseg));
  for ll = 1:length(s)

      s(ll).sphrtfs{ii} = 0;     % init
      s(ll).p{ii} = 0;        % init

      [s(ll).sphrtfs{ii},polang] = extractsp( latseg(ii),s(ll).Obj );
      [s(ll).p{ii},respangs] = baumgartner2014(...
          s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).fs,...
          'S',s(ll).S,'lat',latseg(ii),'polsamp',polang); 

      [ qe(ll,ii),pe(ll,ii) ] = baumgartner2014_pmv2ppp( ...
          s(ll).p{ii} , polang , respangs , s(ll).target{ii});

      if flags.do_plot
        if ll ==1; figure; end
        subplot(1,3,ll)
        Nmax = min(150,s(ll).Ntargets{ii});
        idplot = round(1:s(ll).Ntargets{ii}/Nmax:s(ll).Ntargets{ii});
        plot_baumgartner2014(s(ll).p{ii},polang,respangs,...
                  s(ll).target{ii}(idplot),s(ll).response{ii}(idplot),...
                  'MarkerSize',kv.MarkerSize,'cmax',0.05,'no_colorbar');
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

%% ------ FIG 4 -----------------------------------------------------------
if flags.do_fig4
  
  [perr,perr_exp,qerr,qerr_exp,gamma,mrs,Ntargets] = amt_cache('get','parametrization',flags.cachemode);
  
  if isempty(perr)
    amt_disp('Note that this procedure may last one or two hours!');
    
    gamma = [1,3,3,3,4,4,4,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,8,8,8,9,9,9,...
      10,10,10,12,12,12,16,16,16,30,30,30];%,100,100,100];
    mrs = [0,17,18,19,19,20,21,19,20,21, 0,10,100,16,17,18,19,20,21,23,30, 5,...
      19,20,21,19,20,21,19,20,21,19,20,21,19,20,21,19,20,21,19,20,21];%,19,20,21];

    latseg = -60:20:60; % centers of lateral segments
    dlat =  10;  % lateral range (+-) of each segment
        
    for g = 1:length(gamma)
      
      cname = ['result_baseline_g' num2str(gamma(g),'%u') '_mrs' num2str(mrs(g),'%u')];
      [s,qe, pe] = amt_cache('get',cname,flags.cachemode);
      if isempty(s)

        s = data_baumgartner2014('baseline','gamma',gamma(g),'mrsmsp',mrs(g),flags.cachemode);

        qe_exp = zeros(length(s),length(latseg));
        pe_exp = zeros(length(s),length(latseg));
        for ll = 1:length(s)

          s(ll).target = [];
          s(ll).response = [];
          s(ll).Nt = [];
          for ii = 1:length(latseg)

            latresp = s(ll).itemlist(:,7);
            idlat = latresp <= latseg(ii)+dlat & latresp > latseg(ii)-dlat;
            s(ll).mm2 = s(ll).itemlist(idlat,:);

            s(ll).mm2(:,7) = 0; % set lateral angle to 0deg such that localizationerror works also outside +-30deg

            pe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'rmsPmedianlocal'));
            qe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'querrMiddlebrooks'));

            s(ll).target{ii} = real(s(ll).mm2(:,6)); % polar angle of target
            s(ll).response{ii} = real(s(ll).mm2(:,8)); % polar angle of response
            s(ll).Nt{ii} = length(s(ll).target{ii});

          end
        end


        %% LocaMo
        qe = zeros(length(s),length(latseg));
        pe = zeros(length(s),length(latseg));
        for ll = 1:length(s)

          for ii = 1:length(latseg)

            s(ll).sphrtfs{ii} = 0;     % init
            s(ll).p{ii} = 0;        % init

            [s(ll).sphrtfs{ii},polang] = extractsp( latseg(ii),s(ll).Obj );
            [s(ll).p{ii},respangs] = baumgartner2014(...
                s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).fs,...
                'S',s(ll).S,'lat',latseg(ii),'polsamp',polang,...
                'gamma',gamma(g),'mrsmsp',mrs(g)); 

            if s(ll).Nt{ii} > 0
              [ qe(ll,ii),pe(ll,ii) ] = baumgartner2014_pmv2ppp( ...
                  s(ll).p{ii} , polang , respangs , s(ll).target{ii});
            else
              qe(ll,ii) = NaN; 
              pe(ll,ii) = NaN;
            end

          end

        end
        s = rmfield(s,{'Obj','itemlist','mm2','sphrtfs'}); % reduce file size
        amt_cache('set',cname,s,qe, pe, qe_exp, pe_exp);
      end
    end

    %% Combine results to single mat file
    perr = zeros(length(s),length(latseg),length(gamma));
    qerr = perr;
    for g = 1:length(gamma)
      fn = ['result_baseline_g' num2str(gamma(g),'%u') '_mrs' num2str(mrs(g),'%u')];
      [s,qerr(:,:,g), perr(:,:,g),qerr_exp,perr_exp] = amt_cache('get',fn,'cached');
    end
    
    % Number of targets for each listener and lateral segment
    Ntargets = zeros(length(s),7);
    for jj = 1:length(s)
      Ntargets(jj,:) = [s(jj).Nt{:}];
    end
    
    amt_cache('set','parametrization',perr,perr_exp,qerr,qerr_exp,gamma,mrs,Ntargets);
   
  end
  
  [qerr0,perr0] = baumgartner2014_pmv2ppp(ones(72,44)); % chance performances

  % extract all different gammas
  g = gamma;
  gamma = unique(gamma);

  Nset = size(perr,3);
  idnum = Ntargets ~= 0;
  relNt = Ntargets/sum(Ntargets(:));

  % Compute all residues
  resid.perr = zeros(length(gamma),1);
  resid.qerr = resid.perr;
  for ii = 1:Nset

    dperr = perr_exp - perr(:,:,ii);
    dqerr = qerr_exp - qerr(:,:,ii);
    resid.perr(ii) = sqrt( relNt(idnum)' * (dperr(idnum)).^2 );
    resid.qerr(ii) = sqrt( relNt(idnum)' * (dqerr(idnum)).^2 );

  end
  resid.total = resid.perr/perr0 + resid.qerr/qerr0;

  % Select optimal residues for various gamma
  id_g = zeros(length(gamma),1);
  etotal_g = zeros(length(gamma),1);
  for ii = 1:length(gamma)
    idgamma = find(g == gamma(ii));
    [etotal_g(ii),id] = min(resid.total(idgamma));
    id_g(ii) = idgamma(id);
  end
  eperr_g = resid.perr(id_g);
  eqerr_g = resid.qerr(id_g);

  [tmp,idopt] = min(etotal_g);
  etotal_g = etotal_g / etotal_g(idopt);
  eperr_g = eperr_g / eperr_g(idopt);
  eqerr_g = eqerr_g / eqerr_g(idopt);
  amt_disp(['Optimal Gamma: ' num2str(gamma(idopt),'%u') ' dB^-1'],'documentation');

  % Select residues for optimal gamma and various mrs
  idgammaopt = find(g == gamma(idopt));
  mrs_gopt = mrs(idgammaopt);
  [mrssort,idsort] = sort(mrs_gopt);
  idmrs = idgammaopt(idsort);
  [tmp,idopt_mrs] = min(resid.total(idmrs));
  idnorm = idmrs(idopt_mrs);
  etotal_gopt = resid.total(idmrs) / resid.total(idnorm);
  eperr_gopt = resid.perr(idmrs) / resid.perr(idnorm);
  eqerr_gopt = resid.qerr(idmrs) / resid.qerr(idnorm);
  amt_disp(['Optimal MRS: ' num2str(mrssort(idopt_mrs),'%u') ' deg'],'documentation');
  
  if flags.do_plot
    
    %% Plot residues for various gamma

    % Interpolate data
    gamma_int = logspace(0,2.1,1000);
    inttype = 'pchip';
    dperr_int = interp1(log10(gamma),eperr_g,log10(gamma_int),inttype);
    dqerr_int = interp1(log10(gamma),eqerr_g,log10(gamma_int),inttype);
    dtot_int = interp1(log10(gamma),etotal_g,log10(gamma_int),inttype);

    % Plot
    figure;
    subplot(1,2,1)
    semilogx(gamma_int,dperr_int,'k: ')
    hold on
    semilogx(gamma_int,dqerr_int,'k--')
    semilogx(gamma_int,dtot_int,'k-')
    semilogx(gamma(idopt),0.95,'vk','MarkerFaceColor','k','MarkerSize',kv.MarkerSize+1)

    leg = legend('PE','QE','PE&QE','\{\epsilon,\Gamma\}_{opt}');
    set(leg,'Location','northeast','FontSize',kv.FontSize-1)

    ylabel('e(\Gamma) / e(\Gamma_{opt})','FontSize',kv.FontSize)
    xlabel('\Gamma (dB^{-1})','FontSize',kv.FontSize)

    set(gca,'XLim',[gamma(1)-0.1 gamma(end)+20],'YLim',[0.91 1.7],'XMinorTick','on',...
      'FontSize',kv.FontSize-1)
    set(gca,'XTick',[1:10,20:10:100],'XTickLabel',{1,2,3,'',5,'','','','',10,20,30,'',50,'','','','',100})
    set(gca,'TickLength',TickLength)

    %% Plot residues for optimal gamma and various mrs

    % Interpolation
    mrs_int = 0:0.1:45;
    inttype = 'pchip';
    dperr_int = interp1(mrssort,eperr_gopt,mrs_int,inttype);
    dqerr_int = interp1(mrssort,eqerr_gopt,mrs_int,inttype);
    dtot_int = interp1(mrssort,etotal_gopt,mrs_int,inttype);

    % Plot
    subplot(1,2,2)
    plot(mrs_int,dperr_int,'k:')
    hold on
    plot(mrs_int,dqerr_int,'k--')
    plot(mrs_int,dtot_int,'k-')
    plot(mrssort(idopt_mrs),0.95,'vk','MarkerFaceColor','k','MarkerSize',kv.MarkerSize+1)

    ylabel('e(\epsilon) / e(\epsilon_{opt})','FontSize',kv.FontSize)
    xlabel('\epsilon (deg)','FontSize',kv.FontSize)

    set(gca,'XLim',[mrssort(1) 32],'YLim',[0.91 1.7],'XMinorTick','on',...
      'FontSize',kv.FontSize-1)
    set(gca,'TickLength',TickLength)
    
  end
end

%% ------ FIG 12 & TAB 1 -------------------------------------------------
if flags.do_fig12 || flags.do_tab1
  
  [s,qe,pe,qe_exp,pe_exp,latseg] = amt_cache('get','baseline',flags.cachemode);
  if isempty(s)
    
    latseg = -60:20:60; % centers of lateral segments
    dlat =  10;  % lateral range (+-) of each segment

    s = data_baumgartner2014('baseline',flags.cachemode);

    qe_exp = zeros(length(s),length(latseg));
    pe_exp = zeros(length(s),length(latseg));
    for ll = 1:length(s)

      s(ll).target = [];
      s(ll).response = [];
      s(ll).Nt = [];
      for ii = 1:length(latseg)
        
        latresp = s(ll).itemlist(:,7);
        idlat = latresp <= latseg(ii)+dlat & latresp > latseg(ii)-dlat;
        s(ll).mm2 = s(ll).itemlist(idlat,:);

        s(ll).mm2(:,7) = 0; % set lateral angle to 0deg such that localizationerror works outside +-30deg

        pe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'rmsPmedianlocal'));
        qe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'querrMiddlebrooks'));

        s(ll).target{ii} = real(s(ll).mm2(:,6)); % polar angle of target
        s(ll).response{ii} = real(s(ll).mm2(:,8)); % polar angle of response
        s(ll).Nt{ii} = length(s(ll).target{ii});

      end
    end


    %% LocaMo
    qe = zeros(length(s),length(latseg));
    pe = zeros(length(s),length(latseg));
    for ll = 1:length(s)

      for ii = 1:length(latseg)

        s(ll).sphrtfs{ii} = 0;     % init
        s(ll).p{ii} = 0;        % init

        [s(ll).sphrtfs{ii},polang{ii}] = extractsp( latseg(ii),s(ll).Obj );
        [s(ll).p{ii},respangs{ii}] = baumgartner2014(...
            s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).fs,...
            'S',s(ll).S,'lat',latseg(ii),'polsamp',polang{ii}); 

        if s(ll).Nt{ii} > 0
          [ qe(ll,ii),pe(ll,ii) ] = baumgartner2014_pmv2ppp( ...
              s(ll).p{ii} , polang{ii} , respangs{ii} , s(ll).target{ii});
        else
          qe(ll,ii) = NaN; 
          pe(ll,ii) = NaN;
        end

      end
      
      [s(ll).la,s(ll).le,s(ll).ci,s(ll).lr] = baumgartner2014_likelistat(s(ll).p,polang,respangs,s(ll).target,s(ll).response);
      
      
    end
%     sum( ci(:,1)-la' <=0 & ci(:,2)-la' >=0 )

    s = rmfield(s,{'Obj','itemlist','mm2','sphrtfs'}); % reduce file size 
    amt_cache('set','baseline',s,qe,pe,qe_exp,pe_exp,latseg)
  end
  
  varargout{1} = struct('s',s, 'qe',qe, 'pe',pe, 'qe_exp',qe_exp, 'pe_exp',pe_exp,...
      'latseg',latseg);
  
  flags.do_pm30deglat = true; % consider lateral range of +-30 deg

  
  Ns = length(s);
  relfreq = zeros(Ns,length(latseg));
  Ntall = nan(1,Ns);
  for jj = 1:Ns
    Ntlat = [s(jj).Nt{:}];
    Ntall(jj) = sum(Ntlat);
    relfreq(jj,:) = Ntlat/Ntall(jj);
  end

  if flags.do_pm30deglat
    idlat = find(latseg <= 30 & latseg >= -30);
  else % consider only median plane (+-10 deg)
    idlat = latseg == 0;
  end
  relfreqPerSubject = relfreq(:,idlat)./repmat(sum(relfreq(:,idlat),2),1,3);
  pe = sum(relfreqPerSubject .* pe(:,idlat) , 2);
  qe = sum(relfreqPerSubject .* qe(:,idlat) , 2);

  
  if flags.do_fig12

    % IDs for xlabel
    NHs = nan(Ns,4);
    for ll = 1:Ns
      NHs(ll,:) = s(ll).id;
    end

    if flags.do_plot

      [tmp,idsort] = sort([s.la]);
      la = [s.la];
      le = [s.le];
      ci = [s.ci]';
      lr = [s.lr]';
      
      fig=figure;

      plot_baumgartner2014_likelistat(la(idsort),le(idsort),ci(idsort,:),lr(idsort,:))
      ylabel({'Likelihood'},'FontSize',kv.FontSize)
      xlabel('Listener (NH)','FontSize',kv.FontSize)
      set(gca,'XLim',[0 Ns+1],'XTickLabel',NHs(idsort,3:4),...
          'YMinorTick','on','FontSize',kv.FontSize)
      
      % Bottom Line
      hold on
      ylim = [min(lr(:,1)) max(lr(:,4))];
      plot([0,Ns+1],[ylim(1),ylim(1)]+0.0015*diff(ylim),'k')
     
      set(fig,'PaperPosition',[1,1,9,3.5])

    end
    
  end
  
  if flags.do_tab1

    Labels = {'ID','N','S','actual QE','predicted QE','actual PE','predicted PE'};    

    mtx = zeros(length(Labels),length(s));
    for ll = 1:length(s)
      mtx(1,ll) = str2double(s(ll).id(3:end));
      mtx(2,ll) = sum([s(ll).Ntargets{:}]);
      mtx(3,ll) = s(ll).S;
      mtx(4,ll) = s(ll).qe_exp;
      mtx(5,ll) = qe(ll);
      mtx(6,ll) = s(ll).pe_exp;
      mtx(7,ll) = pe(ll);
    end
    
    [tmp,idsort] = sort(mtx(1,:)); % sort acc. to ID
    mtx = mtx(:,idsort);
    
    varargout{1} = mtx;
    varargout{2} = Labels;
    
  end

end

%% ------ FIG 14 ----------------------------------------------------------
if flags.do_fig14
  
  [s,qe,pe,qe_exp,pe_exp,latseg] = amt_cache('get','baseline',flags.cachemode);
  if isempty(s)
    exp_baumgartner2014('fig5',flags.cachemode);
    [s,qe,pe,qe_exp,pe_exp,latseg] = amt_cache('get','baseline',flags.cachemode);
  end
  
  [paradata.perr,~,paradata.qerr,~,paradata.g,paradata.mrs] = ...
    amt_cache('get','parametrization',flags.cachemode);

  idmrs0 = paradata.g == 6 & paradata.mrs == 0;
  mrs0.pe = paradata.perr(:,:,idmrs0);
  mrs0.qe = paradata.qerr(:,:,idmrs0);

  %% # of targets
  Ns = length(s);
  Nlat = length(latseg);
  Ntlat = zeros(Ns,Nlat);
  relfreq = zeros(Ns,Nlat);
  Ntall = zeros(Ns,1);
  for jj = 1:Ns
    Ntlat(jj,:) = [s(jj).Nt{:}];
    Ntall(jj) = sum(Ntlat(jj,:));
    relfreq(jj,:) = Ntlat(jj,:)/Ntall(jj);
  end
  relfreq = relfreq.*repmat(Ntall,1,Nlat)/sum(Ntall);

  %% Pooling to lateralization
  idlat0 = round(Nlat/2);
  idleft = idlat0-1:-1:1;
  idright = idlat0+1:Nlat;
  latseg = latseg(idlat0:end);
  relfreqLR = Ntlat(:,idleft) ./ (Ntlat(:,idleft) + Ntlat(:,idright) + eps);

  pe = [pe(:,idlat0) , relfreqLR.*pe(:,idleft) + (1-relfreqLR).*pe(:,idright)];
  pe_exp = [pe_exp(:,idlat0) , relfreqLR.*pe_exp(:,idleft) + (1-relfreqLR).*pe_exp(:,idright)];
  qe = [qe(:,idlat0) , relfreqLR.*qe(:,idleft) + (1-relfreqLR).*qe(:,idright)];
  qe_exp = [qe_exp(:,idlat0) , relfreqLR.*qe_exp(:,idleft) + (1-relfreqLR).*qe_exp(:,idright)];
  relfreq = [relfreq(:,idlat0) , relfreq(:,1:idlat0-1) + relfreq(:,Nlat:-1:idlat0+1)];

  mrs0.pe = [mrs0.pe(:,idlat0) , relfreqLR.*mrs0.pe(:,idleft) + (1-relfreqLR).*mrs0.pe(:,idright)];
  mrs0.qe = [mrs0.qe(:,idlat0) , relfreqLR.*mrs0.qe(:,idleft) + (1-relfreqLR).*mrs0.qe(:,idright)];

  %% Evaluation Metrics
  idnum = not(isnan(pe_exp) | isnan(pe));
  dpe = sqrt( relfreq(idnum)' * (pe_exp(idnum) - pe(idnum)).^2 );
  dqe = sqrt( relfreq(idnum)' * (qe_exp(idnum) - qe(idnum)).^2 );
  r_pe = corrcoef(pe_exp(idnum),pe(idnum));
  r_qe = corrcoef(qe_exp(idnum),qe(idnum));

  mrs0.dpe = sqrt( relfreq(idnum)' * (pe_exp(idnum) - mrs0.pe(idnum)).^2 );
  mrs0.dqe = sqrt( relfreq(idnum)' * (qe_exp(idnum) - mrs0.qe(idnum)).^2 );
  [mrs0.r_pe,mrs0.p_pe] = corrcoef(pe_exp(idnum),mrs0.pe(idnum));
  [mrs0.r_qe,mrs0.p_qe] = corrcoef(qe_exp(idnum),mrs0.qe(idnum));


  %% Quartiles
  quart_pe = zeros(3,length(latseg),2); % 1st dim: 25/50/75 quantiles; 2nd dim: lat; 3rd dim: model/experiment/mrs0
  quart_qe = zeros(3,length(latseg),2);
  qlow = 0.25;
  qhigh = 0.75;
  for ii = 1:length(latseg)

    id = not(isnan(pe(:,ii)));
    quart_pe(:,ii,1) = quantile(pe(id,ii),[qlow .50 qhigh]);
    quart_pe(:,ii,3) = quantile(mrs0.pe(id,ii),[qlow .50 qhigh]);
    id = not(isnan(pe_exp(:,ii)));
    quart_pe(:,ii,2) = quantile(pe_exp(id,ii),[qlow .50 qhigh]);

    id = not(isnan(qe(:,ii)));
    quart_qe(:,ii,1) = quantile(qe(id,ii),[qlow .50 qhigh]);
    quart_qe(:,ii,3) = quantile(mrs0.qe(id,ii),[qlow .50 qhigh]);
    id = not(isnan(qe_exp(:,ii)));
    quart_qe(:,ii,2) = quantile(qe_exp(id,ii),[qlow .50 qhigh]);

  end


  if flags.do_plot
     
    dx = 3;
    
    %% PE

    fig = figure;
    subplot(1,2,1)
    errorbar(latseg-dx,quart_pe(2,:,1),...
      quart_pe(2,:,1)-quart_pe(1,:,1),...
      quart_pe(3,:,1)-quart_pe(2,:,1),...
      'ok-','MarkerSize',kv.MarkerSize,'MarkerFaceColor','k');
    hold on
    errorbar(latseg+dx,quart_pe(2,:,3),...
      quart_pe(2,:,3)-quart_pe(1,:,3),...
      quart_pe(3,:,3)-quart_pe(2,:,3),...
      'vk--','MarkerSize',kv.MarkerSize,'MarkerFaceColor','k');
    errorbar(latseg,quart_pe(2,:,2),...
      quart_pe(2,:,2)-quart_pe(1,:,2),...
      quart_pe(3,:,2)-quart_pe(2,:,2),...
      'ok-','MarkerSize',kv.MarkerSize,'MarkerFaceColor','w');
    
    titstr = {['w/ SMM:  e_{PE} = ' num2str(dpe,'%0.1f') '\circ , r_{PE} = ' num2str(r_pe(2),'%0.2f')];...
      ['w/o SMM: e_{PE} = ' num2str(mrs0.dpe,'%0.1f') '\circ , r_{PE} = ' num2str(mrs0.r_pe(2),'%0.2f')]};
    amt_disp(titstr);
    title(titstr,'FontSize',kv.FontSize)
    set(gca,'XLim',[min(latseg)-2*dx,max(latseg)+2*dx],'YLim',[21.1,45.9],...
      'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))
    ylabel('Local Polar RMS Error (deg)','FontSize',kv.FontSize)
    xlabel('Magnitude of Lateral Angle (deg)','FontSize',kv.FontSize)


    %% QE

    subplot(1,2,2)
    errorbar(latseg-dx,quart_qe(2,:,1),...
      quart_qe(2,:,1)-quart_qe(1,:,1),...
      quart_qe(3,:,1)-quart_qe(2,:,1),...
      'ok-','MarkerSize',kv.MarkerSize,'MarkerFaceColor','k');
    hold on
    errorbar(latseg+dx,quart_qe(2,:,3),...
      quart_qe(2,:,3)-quart_qe(1,:,3),...
      quart_qe(3,:,3)-quart_qe(2,:,3),...
      'vk--','MarkerSize',kv.MarkerSize,'MarkerFaceColor','k');
    errorbar(latseg,quart_qe(2,:,2),...
      quart_qe(2,:,2)-quart_qe(1,:,2),...
      quart_qe(3,:,2)-quart_qe(2,:,2),...
      'ok-','MarkerSize',kv.MarkerSize,'MarkerFaceColor','w');
    titstr = {['w/ SMM:  e_{QE} = ' num2str(dqe,'%0.1f') '% , r_{QE} = ' num2str(r_qe(2),'%0.2f')];...
      ['w/o SMM: e_{QE} = ' num2str(mrs0.dqe,'%0.1f') '% , r_{QE} = ' num2str(mrs0.r_qe(2),'%0.2f')]};
    amt_disp(titstr);
    title(titstr,'FontSize',kv.FontSize)
    set(gca,'XLim',[min(latseg)-2*dx,max(latseg)+2*dx],'YLim',[2.1,26.9],...
      'XTick',latseg,'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))
    ylabel('Quadrant Error (%)','FontSize',kv.FontSize)
    xlabel('Magnitude of Lateral Angle (deg)','FontSize',kv.FontSize)

    l = legend('P with SMM','P w/o SMM','Actual');
    set(l,'FontSize',kv.FontSize-1,'Location','northwest')

    set(fig,'PaperPosition',[1,1,10,3.5])
    
  end
end

%% ------ FIG 5 -----------------------------------------------------------
if flags.do_fig5
  
  latdivision = 0;  % lateral angle
  dlat = 15;

  % Experimental Settings
  Conditions = {'BB','LP','W'};


  %% Computations
  s = data_baumgartner2014('pool',flags.cachemode);  
  s = s(ismember({s.id},'NH12'));
  amt_disp(['Listener: ' s.id],'documentation');
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
                s(ll).spdtfs_c{ii} = local_warphrtf(s(ll).spdtfs{ii},s(ll).fs);
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

        [s(ll).p{ii},rang] = baumgartner2014(...
              s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},s(ll).fs,...
              'S',s(ll).S,'lat',latdivision(ii),...
              'polsamp',s(ll).polang{ii});
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
                targets,responses,'MarkerSize',kv.MarkerSize,'cmax',0.05,'no_colorbar')
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
        
%         [la(C,ll),le(C,ll),ci(C,ll,:)] = baumgartner2014_likelistat(s(ll).p,s(ll).polang,respangs,s(ll).target,s(ll).response);
        
      end

    end

  end
  
  varargout{1} = s;

end


%% ------ FIG 6 -----------------------------------------------------------
if flags.do_fig6
  
  [s,cc] = amt_cache('get','spatstrat',flags.cachemode);
  if isempty(s)
    
    latdivision = [-20,0,20];            % lateral angle
    dlat = 10;

    % Experimental Settings
    Conditions = {'BB','LP','W'};

    %% Computations
    s = data_baumgartner2014('pool',flags.cachemode);
%     chance = [];
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
%             chance = [chance;mm2];
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
                  s(ll).spdtfs_c{ii} = local_warphrtf(s(ll).spdtfs{ii},s(ll).fs);
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

          [s(ll).p{ii},rang] = baumgartner2014(...
                s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},s(ll).fs,...
                'S',s(ll).S,'lat',latdivision(ii),...
                'polsamp',s(ll).polang{ii});
          respangs{ii} = rang;

          [ qe(ii),pe(ii) ] = baumgartner2014_pmv2ppp(s(ll).p{ii} , s(ll).polang{ii} , rang);

          if sum(ismember({data.id},s(ll).id)) % if actual participant actual targets
            [ qe_t(ii),pe_t(ii) ] = baumgartner2014_pmv2ppp( ...
                s(ll).p{ii} , s(ll).polang{ii} , rang , s(ll).target{ii} );

          end

        end

        % Model results of pool
        wlat = cos(deg2rad(latdivision)); % weighting compensating lateral compression
        wlat = wlat/sum(wlat);
        s(ll).qe_pool(C,1) = wlat * qe(:); 
        s(ll).pe_pool(C,1) = wlat * pe(:);

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

          [la(C,ll),le(C,ll),ci(C,ll,:)] = baumgartner2014_likelistat(s(ll).p,s(ll).polang,respangs,s(ll).target,s(ll).response);
        end

      end
    end
    s = rmfield(s,{'Obj','spdtfs_c','spdtfs'});% reduce file size

    %% Compute Chance Performance
%     chance = repmat(chance,10,1);
%     id_chance = randi(size(chance,1),size(chance,1),1);
%     chance(:,8) = chance(id_chance,6);
%     pe_chance = localizationerror(chance,'rmsPmedianlocal');
%     qe_chance = localizationerror(chance,'querrMiddlebrooks');

    [r,p] =  corrcoef([s.qe_exp],[s.qe_part]);
    cc.qe.r = r(2);
    cc.qe.p = p(2);
    amt_disp(['QE: r = ' num2str(r(2),'%0.2f') ', p = ' num2str(p(2),'%0.3f')],'documentation');

    [r,p] =  corrcoef([s.pe_exp],[s.pe_part]);
    cc.pe.r = r(2);
    cc.pe.p = p(2);
    amt_disp(['PE: r = ' num2str(r(2),'%0.2f') ', p = ' num2str(p(2),'%0.3f')],'documentation');

    amt_cache('set','spatstrat',s,cc)
  end
  varargout{1} = s;
  varargout{2} = cc;
  
  %% Measures

  % Quartiles
  quart_pe_part = quantile([s.pe_part]',[.25 .50 .75]);
  quart_qe_part = quantile([s.qe_part]',[.25 .50 .75]);

  quart_pe_pool = quantile([s.pe_pool]',[.25 .50 .75]);
  quart_qe_pool = quantile([s.qe_pool]',[.25 .50 .75]);

  quart_pe_exp = quantile([s.pe_exp]',[.25 .50 .75]);
  quart_qe_exp = quantile([s.qe_exp]',[.25 .50 .75]);

  % RMS Differences
  % individual:
  Ntargets = [s.Nt]'; % # of targets
  relfreq = Ntargets/sum(Ntargets(:));
  sd_pe = ([s.pe_part]'-[s.pe_exp]').^2; % squared differences
  dpe = sqrt(relfreq(:)' * sd_pe(:));    % weighted RMS diff.
  sd_qe = ([s.qe_part]'-[s.qe_exp]').^2;
  dqe = sqrt(relfreq(:)' * sd_qe(:));

  % Chance performance
%   qe0 = qe_chance;
%   pe0 = pe_chance;
  [qe0,pe0] = baumgartner2014_pmv2ppp('chance');

  if flags.do_plot
    
    dx = 0.15;
    
    figure 

    subplot(121)
    errorbar((1:3)+dx,quart_pe_part(2,:),...
        quart_pe_part(2,:) - quart_pe_part(1,:),...
        quart_pe_part(3,:) - quart_pe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar((1:3)-dx,quart_pe_pool(2,:),...
        quart_pe_pool(2,:) - quart_pe_pool(1,:),...
        quart_pe_pool(3,:) - quart_pe_pool(2,:),...
        'ks-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    errorbar((1:3),quart_pe_exp(2,:),...
        quart_pe_exp(2,:) - quart_pe_exp(1,:),...
        quart_pe_exp(3,:) - quart_pe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');

    plot([0,4],[pe0,pe0],'k:')

    title(['e_{PE} = ' num2str(dpe,'%0.1f') '\circ , r_{PE} = ' num2str(cc.pe.r,'%0.2f')],...
      'FontSize',kv.FontSize)
    ylabel('Local Polar RMS Error (deg)','FontSize',kv.FontSize)
    set(gca,...
        'XLim',[0.5 3.5],...
        'XTick',1:3,...
        'YLim',[27 54.9],...
        'XTickLabel',{'BB';'LP';'W'},...
        'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))

    subplot(122)
    errorbar((1:3)+dx,quart_qe_part(2,:),...
        quart_qe_part(2,:) - quart_qe_part(1,:),...
        quart_qe_part(3,:) - quart_qe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar((1:3)-dx,quart_qe_pool(2,:),...
        quart_qe_pool(2,:) - quart_qe_pool(1,:),...
        quart_qe_pool(3,:) - quart_qe_pool(2,:),...
        'ks-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    errorbar((1:3),quart_qe_exp(2,:),...
        quart_qe_exp(2,:) - quart_qe_exp(1,:),...
        quart_qe_exp(3,:) - quart_qe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');

    l = legend('Part.','Pool','Actual');
    set(l,'Location','northwest','FontSize',kv.FontSize-1)
    
    plot([0,4],[qe0 qe0],'k:')

    title(['e_{QE} = ' num2str(dqe,'%0.1f') '% , r_{QE} = ' num2str(cc.qe.r,'%0.2f')],...
      'FontSize',kv.FontSize)
    ylabel('Quadrant Error (%)','FontSize',kv.FontSize)
    set(gca,...
        'XLim',[0.5 3.5],...
        'XTick',1:3,...
        'YLim',[0.1 54],...
        'XTickLabel',{'BB';'LP';'W'},...
        'YAxisLocation','left',...
        'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))

    set(gcf,'PaperPosition',[1,1,10,3.5])
    
  end
end

%% ------ FIG 7 -----------------------------------------------------------
if flags.do_fig7
  
  % Model Settings
  latdivision = 0;            % lateral angle
  dlat = 10;

  % Experimental Settings
  Conditions = {'N24','N9','N3'};

  % Vocoder Settings 
  flow = 300;     % lowest corner frequency
  fhigh = 16000;  % highest corner frequency
  N = [24,9,3];

  %% Computations
  s = data_baumgartner2014('pool',flags.cachemode);
  s = s(ismember({s.id},'NH12')); 
  amt_disp(['Listener: ' s.id],'documentation');
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

        if N(C)==Inf
            s(ll).spdtfs_c{ii} = s(ll).spdtfs{ii};

        else
          n = N(C);

          [syncrnfreq, GETtrain] = local_getVocoder('',imp,n,flow,fhigh,0,100,stimPar);
          corners = [syncrnfreq(1);syncrnfreq(:,2)];

          ref = s(ll).spdtfs{ii};

          cond = zeros(length(imp),size(ref,2),2);

          for ch = 1:size(ref,3)
              for ang = 1:size(ref,2)
                  cond(:,ang,ch) = local_channelize('', 0.5*ref(:,ang,ch), ref(:,1), imp, n, corners, [], ...
                                  GETtrain, stimPar, 1, 0.01*s(ll).fs, 0.01*s(ll).fs);
              end

          end

          s(ll).spdtfs_c{ii} = cond;

        end
      end
    end


    %% Run Model

    for ll = 1:length(s)
      clear qe pe qe_t pe_t
      for ii = 1:length(latdivision)

        [p,rang] = baumgartner2014(...
              s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},s(ll).fs,...
              'S',s(ll).S,'lat',latdivision(ii),...
              'polsamp',s(ll).polang{ii});

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

          if C==1; figure; end
          subplot(1,3,C)
          plot_baumgartner2014(p,s(ll).polang{ii},rang,...
                s(ll).target{ii},s(ll).response{ii},...
                    'MarkerSize',kv.MarkerSize,'cmax',0.05,'no_colorbar');
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

%% ------ FIG 8 ----------------------------------------------------------
if flags.do_fig8
  
  [s,cc,N] = amt_cache('get','numchan',flags.cachemode);
  if isempty(s)
    
    % Model Settings
    latdivision = 0; % lateral angle
    dlat = 10;

    % Experimental Settings
    Conditions = {'CL','N24','N18','N12','N9','N6','N3'};

    % Vocoder Settings 
    N = fliplr([3,6,9,12,18,24,30]);	% # of vocoder channels
    flow = 300;     % lowest corner frequency
    fhigh = 16000;  % highest corner frequency


    %% Computations
    s = data_baumgartner2014('pool',flags.cachemode);
%     chance = [];
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

          if C==1
            s(ll).spdtfs_c{ii} = s(ll).spdtfs{ii};

          else
            n = N(C);

            [syncrnfreq, GETtrain] = local_getVocoder('',imp,n,flow,fhigh,0,100,stimPar);
            corners = [syncrnfreq(1);syncrnfreq(:,2)];

            ref = s(ll).spdtfs{ii};

            cond = zeros(length(imp),size(ref,2),2);

            for ch = 1:size(ref,3)
              for ang = 1:size(ref,2)
                  cond(:,ang,ch) = local_channelize('', 0.5*ref(:,ang,ch), ref(:,1), imp, n, corners, [], ...
                                  GETtrain, stimPar, 1, 0.01*s(ll).fs, 0.01*s(ll).fs);
              end
            end
            s(ll).spdtfs_c{ii} = cond;

          end
        end
      end


      %% Run Model

      for ll = 1:length(s)
        clear qe pe qe_t pe_t
        for ii = 1:length(latdivision)

          [p,rang] = baumgartner2014(...
                s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},s(ll).fs,...
                'S',s(ll).S,'lat',latdivision(ii),...
                'polsamp',s(ll).polang{ii});
              
          s(ll).p{ii} = p;    
          respangs{ii} = rang;

          [ qe(ii),pe(ii) ] = baumgartner2014_pmv2ppp(p , s(ll).polang{ii} , rang);

          if sum(ismember({data.id},s(ll).id)) % if actual participant actual targets
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

          [la(C,ll),le(C,ll),ci(C,ll,:)] = baumgartner2014_likelistat(s(ll).p,s(ll).polang,respangs,s(ll).target,s(ll).response);
          
        end

      end
      amt_disp(['Condition ' Cond ' completed.']);
    end
    
%     Crange = max(chance(:,8))-min(chance(:,8));
%     chance(:,8) = Crange*rand(size(chance,1),1)-min(chance(:,8));
%     pe_chance = localizationerror(chance,'rmsPmedianlocal');
%     qe_chance = localizationerror(chance,'querrMiddlebrooks');

    [r,p] =  corrcoef([s.qe_exp],[s.qe_part]);
    cc.qe.r = r(2);
    cc.qe.p = p(2);
    
    amt_disp(['QE: r = ' num2str(r(2),'%0.2f') ', p = ' num2str(p(2),'%0.3f')],'documentation');

    [r,p] =  corrcoef([s.pe_exp],[s.pe_part]);
    cc.pe.r = r(2);
    cc.pe.p = p(2);
    amt_disp(['PE: r = ' num2str(r(2),'%0.2f') ', p = ' num2str(p(2),'%0.3f')],'documentation');
    
    s = rmfield(s,{'spdtfs','spdtfs_c','Obj','itemlist'});
    
    amt_cache('set','numchan',s,cc,N)
  end
  varargout{1} = s;
  varargout{2} = cc;
  varargout{3} = N;
  
  %% Measures

  % Quartiles
  quart_pe_part = fliplr(quantile([s.pe_part]',[.25 .50 .75]));
  quart_qe_part = fliplr(quantile([s.qe_part]',[.25 .50 .75]));

  quart_pe_pool = fliplr(quantile([s.pe_pool]',[.25 .50 .75]));
  quart_qe_pool = fliplr(quantile([s.qe_pool]',[.25 .50 .75]));

  quart_pe_exp = fliplr(quantile([s.pe_exp]',[.25 .50 .75]));
  quart_qe_exp = fliplr(quantile([s.qe_exp]',[.25 .50 .75]));

  % RMS Differences
  % individual:
  Ntargets = [s.Nt]'; % # of targets
  relfreq = Ntargets/sum(Ntargets(:));
  sd_pe = ([s.pe_part]'-[s.pe_exp]').^2; % squared differences
  dpe = sqrt(relfreq(:)' * sd_pe(:));    % weighted RMS diff.
  sd_qe = ([s.qe_part]'-[s.qe_exp]').^2;
  dqe = sqrt(relfreq(:)' * sd_qe(:));

  % Chance performance
%   qe0 = qe_chance;
%   pe0 = pe_chance;
  [qe0,pe0] = baumgartner2014_pmv2ppp('chance');

  
  if flags.do_plot
    
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
    errorbar(fliplr(N)-dx,quart_pe_pool(2,:),...
        quart_pe_pool(2,:) - quart_pe_pool(1,:),...
        quart_pe_pool(3,:) - quart_pe_pool(2,:),...
        'ks-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    errorbar(fliplr(N),quart_pe_exp(2,:),...
        quart_pe_exp(2,:) - quart_pe_exp(1,:),...
        quart_pe_exp(3,:) - quart_pe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');
    plot([0,2*max(N)],[pe0,pe0],'k:')
    xlabel('Num. of Channels','FontSize',kv.FontSize)
    ylabel('Local Polar RMS Error (deg)','FontSize',kv.FontSize)

    title(['e_{PE} = ' num2str(dpe,'%0.1f') '\circ , r_{PE} = ' num2str(cc.pe.r,'%0.2f')],...
      'FontSize',kv.FontSize)
    set(gca,'XLim',[1 32],'XTick',[3 6 9 12 18 24 30],...
        'XTickLabel',{3;6;9;12;18;24;'CL'},...
        'YLim',[27 54.9],...
        'YMinorTick','on','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))

    %% QE
    subplot(122)
    errorbar(fliplr(N)+dx,quart_qe_part(2,:),...
        quart_qe_part(2,:) - quart_qe_part(1,:),...
        quart_qe_part(3,:) - quart_qe_part(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    errorbar(fliplr(N)-dx,quart_qe_pool(2,:),...
        quart_qe_pool(2,:) - quart_qe_pool(1,:),...
        quart_qe_pool(3,:) - quart_qe_pool(2,:),...
        'ks-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','k');
    errorbar(fliplr(N),quart_qe_exp(2,:),...
        quart_qe_exp(2,:) - quart_qe_exp(1,:),...
        quart_qe_exp(3,:) - quart_qe_exp(2,:),...
        'ko-','MarkerSize',kv.MarkerSize,...
        'MarkerFaceColor','w');

    l = legend('Part.','Pool','Actual');
    set(l,'Location','northeast','FontSize',kv.FontSize-1)
      
    plot([0,2*max(N)],[qe0,qe0],'k:')
      
    title(['e_{QE} = ' num2str(dqe,'%0.1f') '% , r_{QE} = ' num2str(cc.qe.r,'%0.2f')],...
      'FontSize',kv.FontSize)
    xlabel('Num. of Channels','FontSize',kv.FontSize)
    ylabel('Quadrant Error (%)','FontSize',kv.FontSize)
    set(gca,'XLim',[1 32],'XTick',[3 6 9 12 18 24 30],...
        'XTickLabel',{3;6;9;12;18;24;'CL'},...
        'YLim',[0.1 54],...
        'YMinorTick','on',...
        'YAxisLocation','left','FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))
      

    set(gcf,'PaperPosition',[1,1,10,3.5])

  end
end

%% ------ FIG 9 ----------------------------------------------------------
if flags.do_fig9
  
  [qe_pool,pe_pool,pb_pool] = amt_cache('get','nonindividual',flags.cachemode);
  if isempty(qe_pool)
    
    % Settings
    latdivision = [-20,0,20];  % lateral center angles of SPs
    flow = 4e3;

    s = data_baumgartner2014('pool',flags.cachemode);
    ns = length(s);

    % DTFs of the SPs
    for ll = 1:ns
      for ii = 1:length(latdivision)

        s(ll).latang{ii} = latdivision(ii);
        s(ll).polangs{ii} = [];
        s(ll).spdtfs{ii} = [];
        [s(ll).spdtfs{ii},s(ll).polangs{ii}] = extractsp(...
            s(ll).latang{ii},s(ll).Obj);

      end
    end

    amt_disp('Please wait a moment!');
    qe = zeros(ns,ns,length(latdivision)); % init QEs
    pe = qe;           % init PEs
    pb = qe;           % init Polar Biases
    for ll = 1:ns    % listener
        for jj = 1:ns    % ears
            for ii = 1:length(latdivision) % SPs

              s(ll).p{jj,ii} = [];
              s(ll).respangs{ii} = [];
              [s(ll).p{jj,ii},s(ll).respangs{ii}] = baumgartner2014(...
                  s(jj).spdtfs{ii},s(ll).spdtfs{ii},s(ll).fs,...
                  'S',s(ll).S,'lat',s(ll).latang{ii},...
                  'polsamp',s(ll).polangs{ii},'flow',flow);

              [ qe(ll,jj,ii),pe(ll,jj,ii),pb(ll,jj,ii) ] = baumgartner2014_pmv2ppp( ...
                  s(ll).p{jj,ii} , s(jj).polangs{ii} , s(ll).respangs{ii});

            end
        end
        amt_disp([' Subject ' num2str(ll,'%2u') ' of ' num2str(ns,'%2u')]);
    end

    lat_weight = cos(pi*latdivision/180);     %lateral weight compensating compression of polar dimension
    lat_weight = lat_weight/sum(lat_weight);  % normalize
    lat_weight = repmat(reshape(lat_weight,[1,1,length(latdivision)]),[ns,ns,1]);
    qe_pool = sum(qe.*lat_weight,3);
    pe_pool = sum(pe.*lat_weight,3);
    pb_pool = sum(pb.*lat_weight,3);

    amt_cache('set','nonindividual',qe_pool,pe_pool,pb_pool);
  end
  varargout{1} = struct('qe',qe_pool,'pe',pe_pool,'pb',pb_pool,'dimensions',{'listener (template)','ears (target)'});
  
  data = data_middlebrooks1999;
  
  %% Model outcomes
  ns = size(pe_pool,1);
  own = eye(ns) == 1;
  other = not(own);
  pb_pool = abs(pb_pool);
  qe_own.quantiles = quantile(qe_pool(own),[0,0.05,0.25,0.5,0.75,0.95,1]);
  pe_own.quantiles = quantile(pe_pool(own),[0,0.05,0.25,0.5,0.75,0.95,1]);
  pb_own.quantiles = quantile(pb_pool(own),[0,0.05,0.25,0.5,0.75,0.95,1]);
  qe_own.mean = mean(qe_pool(own));
  pe_own.mean = mean(pe_pool(own));
  pb_own.mean = mean(pb_pool(own));
  
  qe_other.quantiles = quantile(qe_pool(other),[0,0.05,0.25,0.5,0.75,0.95,1]);
  pe_other.quantiles = quantile(pe_pool(other),[0,0.05,0.25,0.5,0.75,0.95,1]);
  pb_other.quantiles = quantile(pb_pool(other),[0,0.05,0.25,0.5,0.75,0.95,1]);
  qe_other.mean = mean(qe_pool(other));
  pe_other.mean = mean(pe_pool(other));
  pb_other.mean = mean(pb_pool(other));
  
  if flags.do_plot
    dx = -0.2;
    Marker = 'ks';
    data.Marker = 'ko';
    MFC = 'k'; % Marker Face Color
    data.MFC = 'w';
    
    figure;
    subplot(131)
    local_middlebroxplot(1-dx,qe_own.quantiles,kv.MarkerSize)
    plot(1-dx,qe_own.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',MFC)
    local_middlebroxplot(1+dx,data.qe_own.quantiles,kv.MarkerSize)
    plot(1+dx,data.qe_own.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.MFC)
    local_middlebroxplot(2-dx,qe_other.quantiles,kv.MarkerSize)
    plot(2-dx,qe_other.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',MFC)
    local_middlebroxplot(2+dx,data.qe_other.quantiles,kv.MarkerSize)
    plot(2+dx,data.qe_other.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.MFC)
    ylabel('Quadrant Errors (%)','FontSize',kv.FontSize)
    set(gca,'YLim',[-2 43],'XLim',[0.5 2.5],...
      'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))

    subplot(132)
    plot(1-dx,pe_own.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',MFC)
    hold on
    plot(1+dx,data.pe_own.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.MFC)

    local_middlebroxplot(1-dx,pe_own.quantiles,kv.MarkerSize)
    plot(1-dx,pe_own.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',MFC)
    local_middlebroxplot(1+dx,data.pe_own.quantiles,kv.MarkerSize)
    plot(1+dx,data.pe_own.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.MFC)
    local_middlebroxplot(2-dx,pe_other.quantiles,kv.MarkerSize)
    plot(2-dx,pe_other.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',MFC)
    local_middlebroxplot(2+dx,data.pe_other.quantiles,kv.MarkerSize)
    plot(2+dx,data.pe_other.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.MFC)
    ylabel('Local Polar RMS Error (deg)','FontSize',kv.FontSize)
    set(gca,'YLim',[-2 62],'XLim',[0.5 2.5],...
      'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))

    subplot(133)
    local_middlebroxplot(1-dx,pb_own.quantiles,kv.MarkerSize)
    plot(1-dx,pb_own.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',MFC)
    local_middlebroxplot(1+dx,data.pb_own.quantiles,kv.MarkerSize)
    plot(1+dx,data.pb_own.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.MFC)
    local_middlebroxplot(2-dx,pb_other.quantiles,kv.MarkerSize)
    plot(2-dx,pb_other.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',MFC)
    local_middlebroxplot(2+dx,data.pb_other.quantiles,kv.MarkerSize)
    plot(2+dx,data.pb_other.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.MFC)

    ylabel('Magnitude of Elevation Bias (deg)','FontSize',kv.FontSize)
    set(gca,'YLim',[-2 55],'XLim',[0.5 2.5],...
      'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
        'TickLength',2*get(gca,'TickLength'))
  end
end

%% ------ FIG 10 ----------------------------------------------------------
if flags.do_fig10
  
  [pe_exp1,pe_exp2,pe_flat,noDCN] = amt_cache('get','ripples',flags.cachemode);
  if isempty(pe_exp1)
    
    do_exp1 = true;
    do_exp2 = true;
    plotpmv = false;

    density = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8]; % ripples/oct
    depth =   10:10:40;        % ripple depth (peak-to-trough) in dB

    %% Stimulus: 
    % 250-ms bursts, 20-ms raised-cosine fade in/out, flat from 0.6-16kHz

    fs = 48e3;    % sampling rate
    flow = 1e3;   % lower corner frequency of ripple modification in Hz
    fhigh = 16e3; % upper corner frequency of ripple modification in Hz
    Nf = 2^10;    % # Frequency bins

    f = 0:fs/2/Nf:fs/2;	% frequency bins
    id600 = find(f<=600,1,'last'); % index of 600 Hz (lower corner frequency of stimulus energy)
    idlow = find(f<=flow,1,'last'); % index of flow (ripples)
    idhigh = find(f>=fhigh,1,'first');  % index of fhigh (ripples)
    N600low = idlow - id600 +1;   % # bins without ripple modification
    Nlowhigh = idhigh - idlow +1; % # bins with ripple modification     % 
    O = log2(f(idlow:idhigh)/1e3);   % freq. trafo. to achieve equal ripple density in log. freq. scale

    % Raised-cosine "(i.e., cos^2)" ramp 1/8 octave wide
    fup = f(idlow)*2^(1/8);       % upper corner frequency of ramp upwards 
    idup = find(f<=fup,1,'last');
    Nup = idup-idlow+1;
    rampup = cos(-pi/2:pi/2/(Nup-1):0).^2;
    fdown = f(idhigh)*2^(-1/8);  % lower corner frequency of ramp downwards
    iddown = find(f>=fdown,1,'first');
    Ndown = idhigh-iddown+1;
    rampdown = cos(0:pi/2/(Ndown-1):pi/2).^2;
    ramp = [rampup ones(1,Nlowhigh-Nup-Ndown) rampdown];
    ramp = [-inf*ones(1,id600-1) zeros(1,N600low) ramp -inf*ones(1,Nf - idhigh)];

    % Ripples of Experiment I
    Sexp1 = zeros(Nf+1,length(density),2);  % 3rd dim: 1:0-phase 2:pi-phase
    Sexp1(idlow:idhigh,:,1) = (40/2* sin(2*pi*density'*O+ 0))';  % depth: 40dB, 0-phase
    Sexp1(idlow:idhigh,:,2) = (40/2* sin(2*pi*density'*O+pi))';  % depth: 40dB, pi-phase
    Sexp1 = repmat(ramp',[1,length(density),2]) .* Sexp1;
    Sexp1 = [Sexp1;Sexp1(Nf:-1:2,:,:)];
    Sexp1(isnan(Sexp1)) = -100;
%     sexp1 = ifftreal(10.^(Sexp1/20),2*Nf);
    sexp1 = real(ifft(10.^(Sexp1/20),2*Nf));
    sexp1 = circshift(sexp1,Nf);  % IR corresponding to ripple modification
    % Ripples of Experiment II
    Sexp2 = zeros(Nf+1,length(depth),2);  % 3rd dim: 1:0-phase 2:pi-phase
    Sexp2(idlow:idhigh,:,1) = (depth(:)/2*sin(2*pi*1*O+ 0))';  % density: 1 ripple/oct, 0-phase
    Sexp2(idlow:idhigh,:,2) = (depth(:)/2*sin(2*pi*1*O+pi))';  % density: 1 ripple/oct, pi-phase
    Sexp2 = repmat(ramp',[1,length(depth),2]) .* Sexp2;
    Sexp2 = [Sexp2;Sexp2(Nf-1:-1:2,:,:)];
    Sexp2(isnan(Sexp2)) = -100;
%     sexp2 = ifftreal(10.^(Sexp2/20),2*Nf);
    sexp2 = real(ifft(10.^(Sexp2/20),2*Nf));
    sexp2 = circshift(sexp2,Nf);  % IR corresponding to ripple modification


    %% Modeling
    for psge = 0:1

      if psge == 1
        s = data_baumgartner2014('pool',flags.cachemode);
      else % recalib
        s = data_baumgartner2014('pool','do',psge,flags.cachemode);
      end

    latseg = 0;   % centers of lateral segments
    runs = 5;     % # runs of virtual experiments

    pe_exp1 = zeros(length(latseg),length(s),length(density),2);
    pe_exp2 = zeros(length(latseg),length(s),length(depth),2);
    pe_flat = zeros(length(latseg),length(s));
    for ss = 1:length(s)
      for ll = 1:length(latseg)

        [spdtfs,polang] = extractsp(latseg(ll),s(ss).Obj);

        % target elevation range of +-60 deg
        idt = find( polang<=60 | polang>=120 );
        targets = spdtfs(:,idt,:);
        tang = polang(idt);

        [pflat,rang] = baumgartner2014(targets,spdtfs,...
            'S',s(ss).S,'polsamp',polang,...
            'lat',latseg(ll),'stim',[1;0],'do',psge); % Impulse
        mflat = baumgartner2014_virtualexp(pflat,tang,rang,'runs',runs);
        [f,r] = localizationerror(mflat,'sirpMacpherson2000');
        pe_flat(ll,ss) = localizationerror(mflat,f,r,'perMacpherson2003');

        if plotpmv, 
          figure; 
          plot_baumgartner2014(pflat,tang,rang,mflat(:,6),mflat(:,8));title(num2str(pe_flat(ll,ss),2));pause(0.5); 
        end 

        if do_exp1  % Exp. I
        for ii = 1:2*length(density)

          [p,rang] = baumgartner2014(targets,spdtfs,...
            'S',s(ss).S,'polsamp',polang,...
            'lat',latseg(ll),'stim',sexp1(:,ii),'do',psge);
          m = baumgartner2014_virtualexp(p,tang,rang,'runs',runs);
          pe_exp1(ll,ss,ii) = localizationerror(m,f,r,'perMacpherson2003');% - pe_flat(ll,ss);

          if plotpmv; figure; plot_baumgartner2014(p,tang,rang,m(:,6),m(:,8));title([num2str(density(mod(ii-1,10)+1)) 'ripples/oct; PE:' num2str(pe_exp1(ll,ss,ii),2) '%']);pause(0.5); end

        end
        end

        if do_exp2 % Exp. II
        for ii = 1:2*length(depth)

          [p,rang] = baumgartner2014(targets,spdtfs,...
            'S',s(ss).S,'polsamp',polang,...
            'lat',latseg(ll),'stim',sexp2(:,ii),'do',psge);
          m = baumgartner2014_virtualexp(p,tang,rang,'runs',runs);
          pe_exp2(ll,ss,ii) = localizationerror(m,f,r,'perMacpherson2003');% - pe_flat(ll,ss);

          if plotpmv; plot_baumgartner2014(p,tang,rang,m(:,6),m(:,8));title([num2str(depth(mod(ii-1,4)+1)) 'dB; PE:' num2str(pe_exp2(ll,ss,ii),2) '%']);pause(0.5); end

        end
        end

      end
      amt_disp([num2str(ss,'%2u') ' of ' num2str(length(s),'%2u') ' subjects completed']);

    end
    amt_disp();

    if length(latseg) > 1
      pe_exp1 = squeeze(mean(pe_exp1));
      pe_exp2 = squeeze(mean(pe_exp2));
      pe_flat = squeeze(mean(pe_flat));
    else 
      pe_exp1 = squeeze(pe_exp1);
      pe_exp2 = squeeze(pe_exp2);
      pe_flat = squeeze(pe_flat);
    end

    %% Save
      if psge==0
        noDCN.pe_exp1 = pe_exp1;
        noDCN.pe_exp2 = pe_exp2;
        noDCN.pe_flat = pe_flat;
        delete(which('baumgartner2014calibration.mat'))
      end
    end

    amt_cache('set','ripples',pe_exp1,pe_exp2,pe_flat,noDCN)
  end
  varargout{1} = {pe_exp1,pe_exp2,pe_flat,noDCN};
  
  amt_disp(...
    {'Note: Predicted results slightly deviate from the original publication ';...
    'because of a mismatch in the selection of responses considered for the ';...
    'evaluation of polar error rates. Macpherson and Middlebrooks (2003) restricted ';...
    'the target angles whereas Baumgartner et al. (2014) restricted the response angles.';...
    'Now results are calculated for restricted target angles in line with ';...
    'Macpherson and Middlebrooks (2003).'})
  
  dcn_flag = true;
  
  % Original data:
  data = data_macpherson2003;


  %% Phase condition handling
  pe_exp1 = mean(pe_exp1,3);
  data.pe_exp1 = mean(data.pe_exp1,3);
  pe_exp2 = mean(pe_exp2,3);
  data.pe_exp2 = mean(data.pe_exp2,3);
  if dcn_flag
      noDCN.pe_exp1 = mean(noDCN.pe_exp1,3);
      noDCN.pe_exp2 = mean(noDCN.pe_exp2,3);
  end
  idphase = 1;
  
  
  %% Increase
  pe_exp1 = pe_exp1 - repmat(pe_flat(:),1,size(pe_exp1,2));
  pe_exp2 = pe_exp2 - repmat(pe_flat(:),1,size(pe_exp2,2));
  if dcn_flag
      noDCN.pe_exp1 = noDCN.pe_exp1 - repmat(noDCN.pe_flat(:),1,size(noDCN.pe_exp1,2));
      noDCN.pe_exp2 = noDCN.pe_exp2 - repmat(noDCN.pe_flat(:),1,size(noDCN.pe_exp2,2));
  end


  %% Statistics
  quart_pe_flat = quantile(pe_flat,[.25 .50 .75]);
  quart_pe_data_flat = quantile(data.pe_flat,[.25 .50 .75]);

  quart_pe_exp1 = quantile(pe_exp1,[.25 .50 .75]);
  quart_pe_data_exp1 = quantile(data.pe_exp1,[.25 .50 .75]);

  quart_pe_exp2 = quantile(pe_exp2,[.25 .50 .75]);
  quart_pe_data_exp2 = quantile(data.pe_exp2,[.25 .50 .75]);

  if dcn_flag
      noDCN.quart_pe_flat = quantile(noDCN.pe_flat,[.25 .50 .75]);
      noDCN.quart_pe_exp1 = quantile(noDCN.pe_exp1,[.25 .50 .75]);
      noDCN.quart_pe_exp2 = quantile(noDCN.pe_exp2,[.25 .50 .75]);
  end

  
  if flags.do_plot
    
    dx = 1.05;
    FontSize = kv.FontSize;
    MarkerSize = kv.MarkerSize;
    
    % Exp1
    figure;
    
    subplot(2,8,1:8)
    errorbar(data.density/dx,quart_pe_exp1(2,:,idphase),...
        quart_pe_exp1(2,:,idphase) - quart_pe_exp1(1,:,idphase),...
        quart_pe_exp1(3,:,idphase) - quart_pe_exp1(2,:,idphase),...
        'ks-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    if dcn_flag
        errorbar(data.density*dx,noDCN.quart_pe_exp1(2,:,idphase),...
        noDCN.quart_pe_exp1(2,:,idphase) - noDCN.quart_pe_exp1(1,:,idphase),...
        noDCN.quart_pe_exp1(3,:,idphase) - noDCN.quart_pe_exp1(2,:,idphase),...
        'kd--','MarkerSize',MarkerSize-1,...
        'MarkerFaceColor','k');
    end
    errorbar(data.density,quart_pe_data_exp1(2,:,idphase),...
        quart_pe_data_exp1(2,:,idphase) - quart_pe_data_exp1(1,:,idphase),...
        quart_pe_data_exp1(3,:,idphase) - quart_pe_data_exp1(2,:,idphase),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','w');
    set(gca,'XScale','log','YMinorTick','on')
    set(gca,'XLim',[0.25/1.2 8*1.2],'XTick',data.density,'YLim',[-16 59],'FontSize',FontSize)
    xlabel('Ripple Density (ripples/octave)','FontSize',FontSize)
    ylabel({'Increase in';'Polar Error Rate (%)'},'FontSize',FontSize)

    if dcn_flag
        leg = legend('P with PSGE','P w/o PSGE','Actual');
    else
        leg = legend('Predicted','Actual');
    end
    set(leg,'FontSize',FontSize-1,'Location','southwest')
    legpos = get(leg,'Position');
    legpos(1) = legpos(1)+0.05;
    set(leg,'Position',legpos)

    %% Exp2

    subplot(2,8,9:13)
    errorbar(data.depth-1,quart_pe_exp2(2,:,idphase),...
        quart_pe_exp2(2,:,idphase) - quart_pe_exp2(1,:,idphase),...
        quart_pe_exp2(3,:,idphase) - quart_pe_exp2(2,:,idphase),...
        'ks-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    if dcn_flag
        errorbar(data.depth+1,noDCN.quart_pe_exp2(2,:,idphase),...
        noDCN.quart_pe_exp2(2,:,idphase) - noDCN.quart_pe_exp2(1,:,idphase),...
        noDCN.quart_pe_exp2(3,:,idphase) - noDCN.quart_pe_exp2(2,:,idphase),...
        'kd--','MarkerSize',MarkerSize-1,...
        'MarkerFaceColor','k');
    end
    errorbar(data.depth,quart_pe_data_exp2(2,:,idphase),...
        quart_pe_data_exp2(2,:,idphase) - quart_pe_data_exp2(1,:,idphase),...
        quart_pe_data_exp2(3,:,idphase) - quart_pe_data_exp2(2,:,idphase),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','w');
    set(gca,'XLim',[data.depth(1)-5 data.depth(end)+5],'XTick',data.depth,...
      'YLim',[-16 59],'YMinorTick','on','FontSize',FontSize)
    xlabel('Ripple Depth (dB)','FontSize',FontSize)
    ylabel({'Increase in';'Polar Error Rate (%)'},'FontSize',FontSize)
    ytick = get(gca,'YTick');
    ticklength = get(gca,'TickLength');

    %% Baseline
    subplot(2,8,14:15)
    errorbar(-0.5,quart_pe_flat(2),...
        quart_pe_flat(2) - quart_pe_flat(1),...
        quart_pe_flat(3) - quart_pe_flat(2),...
        'ks-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','k');
    hold on
    if dcn_flag
        errorbar(0.5,noDCN.quart_pe_flat(2),...
        noDCN.quart_pe_flat(2) - noDCN.quart_pe_flat(1),...
        noDCN.quart_pe_flat(3) - noDCN.quart_pe_flat(2),...
        'kd-','MarkerSize',MarkerSize-1,...
        'MarkerFaceColor','k');
    end
    errorbar(0,quart_pe_data_flat(2),...
        quart_pe_data_flat(2) - quart_pe_data_flat(1),...
        quart_pe_data_flat(3) - quart_pe_data_flat(2),...
        'ko-','MarkerSize',MarkerSize,...
        'MarkerFaceColor','w');
    set(gca,'XLim',[-3 3],'XTick',0,'XTickLabel',{'Baseline'},...
      'YLim',[-15 59],'YTick',ytick,'TickLength',3*ticklength,...
      'FontSize',FontSize,'YAxisLocation','right')
    xlabel(' ','FontSize',FontSize)
    ylabel({'Polar Error Rate (%)'},'FontSize',FontSize)

    %% Overall correlation between actual and predicted median values
    if dcn_flag
      m_pe_pred = [quart_pe_exp1(2,:,idphase) quart_pe_exp2(2,:,idphase)];
      m_pe_pred_noDCN = [noDCN.quart_pe_exp1(2,:,idphase) noDCN.quart_pe_exp2(2,:,idphase)];
      m_pe_actual = [quart_pe_data_exp1(2,:,idphase) quart_pe_data_exp2(2,:,idphase)];
      r = corrcoef(m_pe_pred,m_pe_actual);
      rDCN = r(2);
      r = corrcoef(m_pe_pred_noDCN,m_pe_actual);
      rnoDCN = r(2);
      
%       r = corrcoef(m_pe_pred,m_pe_pred_noDCN);
%       rInter = r(2);
%       [t,p] = corrdifftest(rDCN,rnoDCN,rInter,14,'hotelling')
%       z = corrdifftest(rDCN,rnoDCN,rInter,14,'steiger')

      amt_disp('Correlation between actual and predicted median values (15 conditions):','documentation');
      amt_disp(['w/  PSGE: r = ' num2str(rDCN,'%0.2f')],'documentation');
      amt_disp(['w/o PSGE: r = ' num2str(rnoDCN,'%0.2f')],'documentation');
    end
    
  end
end

%% ------ FIG 11 ----------------------------------------------------------
if flags.do_fig11
  
  [ape_all,qe_all,ape_BBnoise,qe_BBnoise] = amt_cache('get','highfreqatten_do1',flags.cachemode);
  noDCN = amt_cache('get','highfreqatten_do0',flags.cachemode);
  if isempty(ape_all) || isempty(noDCN.ape_all)
    
    fnHarvard = fullfile(amt_basepath,'auxdata','baumgartner2014','HarvardWords');
      
    amt_disp('Note that this computation may take several hours!')
    
    %% Settings
    latseg = 0;%[-20,0,20];   % centers of lateral segments
    NsampModel = 260; % # of modeled speech samples (takes 30min/sample); max: 260
    startSamp = 1; 

    plotpmv = false;
    plotspec = false;


    %% Load Data

    % Speech Samples from Harvard Word list
    speechsample = amt_cache('get','best2005speechSamples');
    if isempty(speechsample)
      fs_orig = 80e3; % Hz
      fs = 48e3;   % Hz
      p_resamp = fs/fs_orig;
      kk = 1; 
      if NsampModel <= 51
        Nsamp = NsampModel;
        Nlists = 1;
      else
        Nsamp = 260;
        Nlists = 5;
      end
      lsamp = 120000*p_resamp;
      speechsample = cell(Nsamp,1);
      for ii = 1:Nlists
        tmp.list = ['list' num2str(ii,'%1.0u')];
        tmp.path = fullfile(fnHarvard,tmp.list);
        tmp.dir = dir(fullfile(tmp.path,'*.mat'));
        for jj = 1:length(tmp.dir)
          if jj > Nsamp; break; end
%           load(fullfile(tmp.path,tmp.dir(jj).name))
          sig = amt_load('baumgartner2014',fullfile('HarvardWords',tmp.list,tmp.dir(jj).name));
          signal = resample(sig.word,p_resamp*10,10);
          gcurve = exp(-0.5 * (0:0.001:10).^2) ./ (sqrt(2*pi));
          env = filter(gcurve,1,signal.^2);
          idon = max(find(env > 5e7,1,'first')-1e3,1);
          idoff = min(find(env > 5e7,1,'last')+1e3,lsamp);
          lwin = idoff-idon+1;
          speechsample{kk} = signal(idon:idoff) .* tukeywin(lwin,0.01)';
          kk = kk + 1;
        end
      end
      amt_cache('set','best2005speechSamples',speechsample);
    end

    % FIR Low-pass filters at 8kHz
    % Brick-wall (aka sinc-filter): fir1(200,1/3) -> -60 dB
    x=amt_load('baumgartner2014','highfreqatten_filters.mat');
    lp{1} = [1 zeros(1,100)];
    lp{2} = x.fir20db;
    lp{3} = x.fir40db;
    lp{4} = x.fir60db;

    %% Model Data
    for psge = 0:1

      s = data_baumgartner2014('pool','do',psge,flags.cachemode);
     
      cname = ['result_best2005noise_do' num2str(psge,'%u')];
      [ape_BBnoise,qe_BBnoise] = amt_cache('get',cname,flags.cachemode);
      if isempty(ape_BBnoise)
        ape_BBnoise = zeros(1,length(s),length(latseg));
        qe_BBnoise = ape_BBnoise;
        for ss = 1:length(s)
          for ll = 1:length(latseg)
            [spdtfs,polang] = extractsp(latseg(ll),s(ss).Obj);
            [p,rang] = baumgartner2014(spdtfs,spdtfs,'do',psge,...
                  'S',s(ss).S,'polsamp',polang,'lat',latseg(ll),'notprint');
            ape_BBnoise(1,ss,ll) = baumgartner2014_pmv2ppp(p,polang,rang,'absPE');
            qe_BBnoise(1,ss,ll) = baumgartner2014_pmv2ppp(p,polang,rang);

            if plotpmv; figure; plot_baumgartner2014(p,polang,rang); title(num2str(ape_BBnoise(1,ss,ll),2)); end

          end
        end
        % Pool Lateral Segments
        if length(latseg) > 1
          ape_BBnoise = mean(ape_BBnoise,3);
          qe_BBnoise = mean(qe_BBnoise,3);
        end
        amt_cache('set',cname,ape_BBnoise,qe_BBnoise);
      end
    end
    
    ape_all = zeros(length(lp),length(s),NsampModel-startSamp+1,2);
    qe_all = ape_all;
    for kk = startSamp:NsampModel 
      for psge = 0:1
        cname = ['result_best2005speech_samp' num2str(kk) '_do' num2str(psge,'%u')];
        [ape_lat,qe_lat] = amt_cache('get',cname,flags.cachemode);
        if isempty(ape_lat)

          s = data_baumgartner2014('pool','do',psge,flags.cachemode);

          ape_lat = zeros(length(lp),length(s),length(latseg));
          qe_lat = ape_lat;
          for ss = 1:length(s)
            for ll = 1:length(latseg)
              for ilp = 1:length(lp)

                stim = filter(lp{ilp},1,speechsample{kk});

                if plotspec; figure; audspecgram(stim(:),fs,'dynrange',150); end

                [spdtfs,polang] = extractsp(latseg(ll),s(ss).Obj);
                [p,rang] = baumgartner2014(spdtfs,spdtfs,'do',psge,...
                  'S',s(ss).S,'polsamp',polang,...
                  'lat',latseg(ll),'stim',stim,'notprint');
                ape_lat(ilp,ss,ll) = baumgartner2014_pmv2ppp(p,polang,rang,'absPE');
                qe_lat(ilp,ss,ll) = baumgartner2014_pmv2ppp(p,polang,rang);

                if plotpmv; figure; plot_baumgartner2014(p,polang,rang); title(num2str(ape_lat(ilp,ss,ll),2)); end

              end
            end
          end
          % Pool Lateral Segments
          if length(latseg) > 1
            ape_lat = mean(ape_lat,3);
            qe_lat = mean(qe_lat,3);
          end
          amt_cache('set',cname,ape_lat,qe_lat);
          amt_disp([num2str(kk,'%1.0u') ' of ' num2str(NsampModel,'%2.0u') ' samples completed']);
        end
        ape_all(:,:,kk,psge+1) = ape_lat;
        qe_all(:,:,kk,psge+1) = qe_lat;

      end
    end

    noDCN.ape_all = ape_all(:,:,:,1);
    noDCN.qe_all = qe_all(:,:,:,1);
    [noDCN.ape_BBnoise,noDCN.qe_BBnoise] = amt_cache('get','result_best2005noise_do1');
    amt_cache('set','highfreqatten_do0',noDCN)
    
    ape_all = ape_all(:,:,:,2);
    qe_all = qe_all(:,:,:,2);
    [ape_BBnoise,qe_BBnoise] = amt_cache('get','result_best2005noise_do1');
    amt_cache('set','highfreqatten_do1',ape_all,qe_all,ape_BBnoise,qe_BBnoise)
    
  end
  
  varargout{1} = {ape_all,qe_all,ape_BBnoise,qe_BBnoise};
  varargout{2} = noDCN;
  
  data = data_best2005;
  
  % Pool Samples
  ape_pooled = mean(ape_all,3);
  qe_pooled = mean(qe_all,3);
  noDCN.ape_pooled = mean(noDCN.ape_all,3);
  noDCN.qe_pooled = mean(noDCN.qe_all,3);

  % Confidence Intervals or standard errors
  df_speech = size(ape_all,2)-1;%*size(ape_all,3)-1;
  tquant_speech = 1;%icdf('t',.975,df_speech);
  seape_speech = std(ape_pooled,0,2)*tquant_speech/(df_speech+1);
  df_noise = size(ape_BBnoise,2)-1;
  tquant_noise = 1;%icdf('t',.975,df_noise);
  seape_noise = std(ape_BBnoise,0,2)*tquant_noise/(df_noise+1);
  seape = [seape_noise;seape_speech];
  % DCN
  df_speech = size(noDCN.ape_all,2)-1;%*size(ape_all,3)-1;
  tquant_speech = 1;%icdf('t',.975,df_speech);
  seape_speech = std(noDCN.ape_pooled,0,2)*tquant_speech/(df_speech+1);
  df_noise = size(noDCN.ape_BBnoise,2)-1;
  tquant_noise = 1;%icdf('t',.975,df_noise);
  seape_noise = std(noDCN.ape_BBnoise,0,2)*tquant_noise/(df_noise+1);
  noDCN.seape = [seape_noise;seape_speech];

  % Means
  ape = mean([ape_BBnoise ; ape_pooled],2);
  qe = mean([qe_BBnoise ; qe_pooled],2);
  noDCN.ape = mean([noDCN.ape_BBnoise ; noDCN.ape_pooled],2);
  noDCN.qe = mean([noDCN.qe_BBnoise ; noDCN.qe_pooled],2);


  if flags.do_plot
    
    dx = 0;
    MarkerSize = kv.MarkerSize;
    FontSize = kv.FontSize;
    
    xticks = 0:size(ape_all,1);
    ape0 = baumgartner2014_pmv2ppp('absPE','chance');
    
    figure;
    subplot(211)
    h(1) = errorbar(xticks-dx,ape,seape,'ks');
    set(h(1),'MarkerFaceColor','k','MarkerSize',MarkerSize,'LineStyle','-')
    hold on
    h(3) = errorbar(xticks+dx,noDCN.ape,noDCN.seape,'kd');
    set(h(3),'MarkerFaceColor','k','MarkerSize',MarkerSize-1,'LineStyle','--')
    h(2) = errorbar(xticks,data.ape,data.seape,'ko');
    set(h(2),'MarkerFaceColor','w','MarkerSize',MarkerSize,'LineStyle','-')
    plot([-0.5 4.5],[ape0 ape0],'k:') % chance performance
%     ylabel('| \theta - \vartheta | (deg)','FontSize',FontSize)
    ylabel('Polar Error (deg)','FontSize',FontSize)
    set(gca,'XTick',xticks,'XTickLabel',[],'FontSize',FontSize)
    set(gca,'XLim',[-0.5 4.5],'YLim',[12 95],'YMinorTick','on')

    pos = get(gca,'Position');
    pos(2) = pos(2)-0.11;
    set(gca,'Position',pos)

%     leg = legend('Predicted with edge extraction','Predicted without edge extraction','Actual');
%     set(leg,'FontSize',FontSize-2,'Location','northoutside')
%     pos = get(leg,'Position');
%     pos(2) = pos(2)+0.14;
%     set(leg,'Position',pos)

    qe0 = baumgartner2014_pmv2ppp('QE','chance');
    
    subplot(212)
    h(1) = plot(xticks-dx,qe,'ks');
    set(h(1),'MarkerFaceColor','k','MarkerSize',MarkerSize,'LineStyle','-')
    hold on
    h(3) = plot(xticks+dx,noDCN.qe,'kd');
    set(h(3),'MarkerFaceColor','k','MarkerSize',MarkerSize-1,'LineStyle','--')
    h(2) = plot(xticks([1 2 5]),data.qe([1 2 5]),'ko'); 
    % In Baumgartner et al. (2014), we accidentially missed the actual data 
    % for the -20dB and -40dB conditions. Uncomment the following line to   
    % also show these data points.
%     h(2) = plot(xticks,data.qe,'ko');
    set(h(2),'MarkerFaceColor','w','MarkerSize',MarkerSize,'LineStyle','-')
    plot([-0.5 4.5],[qe0 qe0],'k:') % chance performance
    ylabel('Quadrant Err. (%)','FontSize',FontSize)
    set(gca,'XTick',xticks,'XTickLabel',data.meta,'FontSize',FontSize,...
      'XLim',[-0.5 4.5],'YLim',[-3 54],'YMinorTick','on')

  end
end


%% ------ TAB 2 ---------------------------------------------------------- 
if flags.do_tab2
  
  [qe_exp,pe_exp,qe_part,pe_part] = amt_cache('get','spatstrat_do0',flags.cachemode);
  if isempty(qe_exp)
    
    latdivision = [-20,0,20];            % lateral angle
    dlat = 10;

    % Experimental Settings
    Conditions = {'BB','LP','W'};

    %% Computations
    for ido = 0:1

      if ido == 1
        s = data_baumgartner2014('pool',flags.cachemode);
      else % recalib
        s = data_baumgartner2014('pool','do',0,flags.cachemode);
      end
    
      for C = 1:length(Conditions)

        Cond = Conditions{C};

        %% Data

        % Experimental data
        data = data_majdak2013(Cond);
        
        % Consider only actual participants
        idpart = [];
        for ii = 1:length(data)
          idpart = [idpart,find(ismember({s.id},data(ii).id))];
        end
        s = s(idpart);
        
        for ll = 1:length(s)
          s(ll).itemlist=data(ismember({data.id},s(ll).id)).mtx; 
          for ii = 1:length(latdivision)
            latresp = s(ll).itemlist(:,7);
            idlat = latresp <= latdivision(ii)+dlat & latresp > latdivision(ii)-dlat;
            mm2 = s(ll).itemlist(idlat,:);
            s(ll).target{ii} = mm2(:,6); % polar angle of target
            s(ll).response{ii} = mm2(:,8); % polar angle of response
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
                s(ll).spdtfs_c{ii} = local_warphrtf(s(ll).spdtfs{ii},s(ll).fs);
            end

          end
        end


      %% Run Model

        for ll = 1:length(s)
            
          fh = [8500,18000]; % Hz
          for ff=1:length(fh)
            qe = zeros(1,length(latdivision));
            pe = zeros(1,length(latdivision));
            qe_t = zeros(1,length(latdivision));
            pe_t = zeros(1,length(latdivision));
            for ii = 1:length(latdivision)

              if C == 1
                fhigh = 18000;
              else
                fhigh = fh(ff);
              end
              
              [s(ll).p{ii},rang] = baumgartner2014(...
                    s(ll).spdtfs_c{ii},s(ll).spdtfs{ii},s(ll).fs,...
                    'S',s(ll).S,'lat',latdivision(ii),...
                    'polsamp',s(ll).polang{ii},'do',ido,'fhigh',fhigh);
              respangs{ii} = rang;

              [ qe_t(ii),pe_t(ii) ] = baumgartner2014_pmv2ppp( ...
                  s(ll).p{ii} , s(ll).polang{ii} , rang , s(ll).target{ii} );

            end

            % Model results of participants
            if length(latdivision) == 3
              qe_part(ll,C,2*ido+ff) = (qe_t(1)*length(s(ll).target{1}) + ...
                  qe_t(2)*length(s(ll).target{2}) + ...
                  qe_t(3)*length(s(ll).target{3}))/...
                  (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
              pe_part(ll,C,2*ido+ff) = (pe_t(1)*length(s(ll).target{1}) + ...
                  pe_t(2)*length(s(ll).target{2}) + ...
                  pe_t(3)*length(s(ll).target{3}))/...
                  (length(s(ll).target{1})+length(s(ll).target{2})+length(s(ll).target{3}));
            else 
              s(ll).qe_part(C,2*ido+ff) = mean(qe_t);
              s(ll).pe_part(C,2*ido+ff) = mean(pe_t);
            end
          end

          % Actual experimental results
          qe_exp(ll,C) = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
          pe_exp(ll,C,1) = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
%           s(ll).Nt(C,1) = size(s(ll).itemlist,1);
            
        end
      end
    end
    s = rmfield(s,{'Obj','spdtfs_c','spdtfs'});% reduce file size
    amt_cache('set','spatstrat_do0',qe_exp,pe_exp,qe_part,pe_part)
  end
  
  result = struct('qe_exp',qe_exp,'pe_exp',pe_exp,'qe_part',qe_part,'pe_part',pe_part);
  
  meta = {'DCN no,  BWA yes';...
          'DCN no,  BWA no ';...
          'DCN yes, BWA yes';...
          'DCN yes, BWA no ';...
          };
  
  [qe0,pe0] = baumgartner2014_pmv2ppp('chance');

  qe_exp = permute(result.qe_exp,[2,1]);
  pe_exp = permute(result.pe_exp,[2,1]);
  qe_all = permute(result.qe_part,[2,1,3]);
  pe_all = permute(result.pe_part,[2,1,3]);

  % Statistics
  for cond = 1:3

    if cond == 1
      amt_disp('BB:','documentation');
    elseif cond == 2
      amt_disp('LP:','documentation');
    else
      amt_disp('W:','documentation');
    end

    Ns = size(pe_exp,2);

    group{1} = repmat(meta{1},Ns,1);
    for im = 2:length(meta)
      group{1} = [group{1} ; repmat(meta{im},Ns,1)];
    end
    group{2} = repmat(1:Ns,1,length(meta));

      data = sqrt(((qe_all(cond,:,:) - repmat(qe_exp(cond,:),[1,1,length(meta)]))/qe0).^2) + ...
             sqrt(((pe_all(cond,:,:) - repmat(pe_exp(cond,:),[1,1,length(meta)]))/pe0).^2);

      [p,t,stat] = friedman(squeeze(data),1);
%       [p,t,stat] = anovan(data(:),group,'display','off');

      amt_disp(['  Chi-sq = ' num2str(t{2,5},'%5.2f') ', p = ' num2str(p(1),'%3.3f')],'documentation');
      
      if p(1) < 0.05
        figure
        [c,m,h,nms] = multcompare(stat,'display','on');
        set(gca,'YTickLabel',flipud(meta))
      end

  end


  e_qe = zeros(length(meta),3);  % model, BB/LP/W
  e_pe = e_qe;
  r_pe = e_qe;
  r_qe = e_qe;
  for m = 1:length(meta)
    for c = 1:3

      e_qe(m,c) = rms(qe_all(c,:,m) - qe_exp(c,:));
      e_pe(m,c) = rms(pe_all(c,:,m) - pe_exp(c,:));
      tmp.r = corrcoef(qe_all(c,:,m) , qe_exp(c,:));
      r_qe(m,c) = tmp.r(2);
      tmp.r = corrcoef(pe_all(c,:,m) , pe_exp(c,:));
      r_pe(m,c) = tmp.r(2);

    end
  end

  %% Write Table
  mtx = [e_pe e_qe];
  mtx(:,1:2:end) = e_pe;
  mtx(:,2:2:end) = e_qe;
  mtx = round(mtx*10)/10;
  columnLabels = {'Spect. Proc.',...
    '$e_\mathrm{PE}$','$e_\mathrm{QE}$',...
    '$e_\mathrm{PE}$','$e_\mathrm{QE}$',...
    '$e_\mathrm{PE}$','$e_\mathrm{QE}$'};
  rowLabels = meta;
  
  varargout{1} = mtx;
  varargout{2} = rowLabels; 
  varargout{3} = columnLabels; 
  
end
 

%% ------ TAB 3 ----------------------------------------------------------    
if flags.do_tab3
  
  [s,qe,pe,qe_exp,pe_exp,latseg,bwcoef] = amt_cache('get','binWeighting',flags.cachemode);
  if isempty(s)
    
    bwcoef = [13 eps -eps Inf];
    latseg = -60:20:60; % centers of lateral segments
    dlat =  10;  % lateral range (+-) of each segment

    s = data_baumgartner2014('baseline',flags.cachemode);

    qe_exp = zeros(length(s),length(latseg));
    pe_exp = zeros(length(s),length(latseg));
    for ll = 1:length(s)

      s(ll).target = [];
      s(ll).response = [];
      s(ll).Nt = [];
      for ii = 1:length(latseg)
        
        latresp = s(ll).itemlist(:,7);
        idlat = latresp <= latseg(ii)+dlat & latresp > latseg(ii)-dlat;
        s(ll).mm2 = s(ll).itemlist(idlat,:);

        s(ll).mm2(:,7) = 0; % set lateral angle to 0deg such that localizationerror works outside +-30deg

        pe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'rmsPmedianlocal'));
        qe_exp(ll,ii) = real(localizationerror(s(ll).mm2,'querrMiddlebrooks'));

        s(ll).target{ii} = real(s(ll).mm2(:,6)); % polar angle of target
        s(ll).response{ii} = real(s(ll).mm2(:,8)); % polar angle of response
        s(ll).Nt{ii} = length(s(ll).target{ii});

      end
    end


    %% LocaMo
    
    qe = zeros(length(s),length(latseg),length(bwcoef));
    pe = zeros(length(s),length(latseg),length(bwcoef));
    for b = 1:length(bwcoef)
      for ll = 1:length(s)

        for ii = 1:length(latseg)

          s(ll).sphrtfs{ii} = 0;     % init
          s(ll).p{ii} = 0;        % init

          [s(ll).sphrtfs{ii},polang{ii}] = extractsp( latseg(ii),s(ll).Obj );
          [s(ll).p{ii},respangs{ii}] = baumgartner2014(...
              s(ll).sphrtfs{ii},s(ll).sphrtfs{ii},s(ll).fs,...
              'S',s(ll).S,'lat',latseg(ii),'polsamp',polang{ii},...
              'bwcoef',bwcoef(b)); 

          if s(ll).Nt{ii} > 0
            [ qe(ll,ii,b),pe(ll,ii,b) ] = baumgartner2014_pmv2ppp( ...
                s(ll).p{ii} , polang{ii} , respangs{ii} , s(ll).target{ii});
          else
            qe(ll,ii,b) = NaN; 
            pe(ll,ii,b) = NaN;
          end

        end

      end
    
    amt_disp([num2str(b,'%1.0f') ' of ' num2str(length(bwcoef),'%1.0f') ' completed']);
    end
    
    s = rmfield(s,{'Obj','itemlist','mm2','sphrtfs'}); % reduce file size 
    amt_cache('set','binWeighting',s,qe,pe,qe_exp,pe_exp,latseg,bwcoef)
    
  end
  
%   load(fn);
  varargout{1} = {s,qe,pe,qe_exp,pe_exp,latseg,bwcoef};
  
  %% # of targets
  Ns = length(s);
  Nlat = length(latseg);
  Ntlat = zeros(Ns,Nlat);
  relfreq = zeros(Ns,Nlat);
  Ntall = zeros(Ns,1);
  for jj = 1:Ns
    Ntlat(jj,:) = [s(jj).Nt{:}];
    Ntall(jj) = sum(Ntlat(jj,:));
    relfreq(jj,:) = Ntlat(jj,:)/Ntall(jj);
  end
  relfreq = relfreq.*repmat(Ntall,1,Nlat)/sum(Ntall);

  %% Pooling to lateralization
  idlat0 = round(Nlat/2);
  idleft = idlat0-1:-1:1;
  idright = idlat0+1:Nlat;
  latseg = latseg(idlat0:end);
  relfreqLR = Ntlat(:,idleft) ./ (Ntlat(:,idleft) + Ntlat(:,idright) + eps);

  relfreq = [relfreq(:,idlat0) , relfreq(:,1:idlat0-1) + relfreq(:,Nlat:-1:idlat0+1)];
  pe_exp = [pe_exp(:,idlat0) , relfreqLR.*pe_exp(:,idleft) + (1-relfreqLR).*pe_exp(:,idright)];
  qe_exp = [qe_exp(:,idlat0) , relfreqLR.*qe_exp(:,idleft) + (1-relfreqLR).*qe_exp(:,idright)];

  %% Evaluation Metrics  
  for b=1:length(bwcoef)
    
    % pooling lat.
    pe_b = [pe(:,idlat0,b) , relfreqLR.*pe(:,idleft,b) + (1-relfreqLR).*pe(:,idright,b)];
    qe_b = [qe(:,idlat0,b) , relfreqLR.*qe(:,idleft,b) + (1-relfreqLR).*qe(:,idright,b)];
    
    idnum = not(isnan(pe_exp) | isnan(pe_b));
    dpe(b) = sqrt( relfreq(idnum)' * (pe_exp(idnum) - pe_b(idnum)).^2 );
    dqe(b) = sqrt( relfreq(idnum)' * (qe_exp(idnum) - qe_b(idnum)).^2 );
    Mpe(b) = relfreq(idnum)' * pe_b(idnum);
    Mqe(b) = relfreq(idnum)' * qe_b(idnum);
    r = corrcoef(pe_exp(idnum),pe_b(idnum));
    r_pe(b) = r(2);
    r = corrcoef(qe_exp(idnum),qe_b(idnum));
    r_qe(b) = r(2);
    
  end

  %% Table

  % Write Table
  mtx = [round([dpe' dqe']*10)/10 , round([r_pe' r_qe']*100)/100 , round([Mpe' Mqe']*10)/10];
  columnLabels = {'','$e_\mathrm{PE}$','$e_\mathrm{QE}$','$r_\mathrm{PE}$',...
    '$r_\mathrm{QE}$','$\overline{\mathrm{PE}}$','$\overline{\mathrm{QE}}$'};
  rowLabels = {'$\Phi = 13^\circ$','$\Phi \rightarrow +0^\circ$','$\Phi \rightarrow -0^\circ$','$\Phi \rightarrow \infty^\circ$'};

  varargout{1} = mtx;
  varargout{2} = rowLabels;
  varargout{3} = columnLabels;
  
end
  


end



%% ------------------------------------------------------------------------
%  ---- INTERNAL FUNCTIONS ------------------------------------------------
%  ------------------------------------------------------------------------
function hM_warped = local_warphrtf(hM,fs)
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
        hwinfade = local_FWfade(hwin,512,24,96,192);
        hM_warped(1:end,el,canal)=hwinfade;
            
    end
end
end

function [syncrnfreq, GETtrain] = local_getVocoder(filename,in,channum,lower,upper,alpha,GaussRate,stimpar)
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
N = duration*GaussRate; % number of pulses
fN = floor(N);
Genv=zeros(fN,nsamples);
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
        Genv=zeros(fN,nsamples);
        Genv(1,:) = sqrt(Gamma(i)) * exp(-pi*(Gamma(i)*(t-T)).^2) .* sin(2*pi*cf(i)*t - T + pi/4); %!!! (t-T)
        Genv=repmat(Genv(1,:),[fN 1]);
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

function out=local_channelize(fwavout, h, h0, in, channum, corners, syncrnfreq, ...
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

ii=max(max(abs(out)));
if ii>=1
  error(['Maximum amplitude value is ' num2str(20*log10(ii)) 'dB. Set the HRTF scaling factor lower to avoid clipping']);
end
out=local_FWfade(out,0,fadein,fadeout);
% wavwrite(out,srate,stimpar.Resolution,fwavout);
end

function out = local_FWfade(inp, len, fadein, fadeout, offset)
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

function local_middlebroxplot(x,quantiles,MarkerSize)

lilen = 0.14; % length of horizontal lines

% Symbols
plot(x,quantiles(1),'kx','MarkerSize',MarkerSize) % min
hold on
plot(x,quantiles(7),'kx','MarkerSize',MarkerSize) % max

% Horizontal lines
line(x+0.5*[-lilen,lilen],repmat(quantiles(2),2),'Color','k') % lower whisker
line(x+[-lilen,lilen],repmat(quantiles(3),2),'Color','k') % 25% Quartile
line(x+[-lilen,lilen],repmat(quantiles(4),2),'Color','k') % Median
line(x+[-lilen,lilen],repmat(quantiles(5),2),'Color','k') % 75% Quartile
line(x+0.5*[-lilen,lilen],repmat(quantiles(6),2),'Color','k') % upper whisker

% Vertical lines
line([x,x],quantiles(2:3),'Color','k') % connector lower whisker
line([x,x],quantiles(5:6),'Color','k') % connector upper whisker
line([x,x]-lilen,quantiles([3,5]),'Color','k') % left box edge
line([x,x]+lilen,quantiles([3,5]),'Color','k') % left box edge

end



