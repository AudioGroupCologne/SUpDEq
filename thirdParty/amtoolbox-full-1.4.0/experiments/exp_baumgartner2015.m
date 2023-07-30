function varargout = exp_baumgartner2015(varargin)
%EXP_BAUMGARTNER2015 Results from  Baumgartner and Majdak (2015)
%   Usage: data = exp_baumgartner2015(flag) 
%
%   EXP_BAUMGARTNER2015(flag) reproduces figures of the studies from 
%   Baumgartner et al. (2015).
%
%
%   The following flags can be specified
%
%
%     'fig2'                       Reproduce Fig.2 of Baumgartner and Majdak (2015):
%                                  Example showing the spectral discrepancies obtained by VBAP. 
%                                  The targeted spectrum is the HRTF for 20 deg polar angle. 
%                                  The spectrum obtained by VBAP is the superposition of two 
%                                  HRTFs from directions 40 deg polar angle apart of each 
%                                  other with the tar- geted source direction centered in between.
%
%     'fig4'                       Reproduce Fig.4 of Baumgartner and Majdak (2015):
%                                  Response predictions to sounds created by VBAP with two 
%                                  loudspeakers in the median plane positioned at polar 
%                                  angles of -15 and 30 deg, respectively. Predictions for 
%                                  two exemplary listeners and pooled across all listeners. 
%                                  Each column of a panel shows the predicted PMV of 
%                                  polar-angle responses to a certain sound. Note the 
%                                  inter-individual differences and the generally small 
%                                  probabilities at response angles not occupied by the loudspeakers.
%
%     'fig5'                       Reproduce Fig.5 of Baumgartner and Majdak (2015):
%                                  Listener-specific increases in polar error as a function of 
%                                  the panning angle. Increase in polar error defined as 
%                                  the difference between the polar error obtained by the 
%                                  VBAP source and the polar error obtained by the real 
%                                  source at the corresponding panning angle. Same loudspeaker 
%                                  arrangement as for Fig. 4. Note the large inter-individual 
%                                  differences and the increase in polar error being largest 
%                                  at panning angles centered between the loudspeakers, i.e., 
%                                  at panning ratios around R = 0 dB.
%
%     'fig6'                       Reproduce Fig.6 of Baumgartner and Majdak (2015):
%                                  Panning angles for the loudspeaker arrangement of Fig. 4 
%                                  judged best for reference sources at polar angles of 
%                                  0 or 15 deg in the median plane. Comparison between 
%                                  experimental results from [2] and simulated results 
%                                  based on various response strategies: PM, CM, and both 
%                                  mixed. Dotted horizontal line: polar angle of the reference 
%                                  source. Hor- izontal line within box: median; 
%                                  box: inter-quartile range (IQR); 
%                                  whisker: within quartile +-1.5 IQR; 
%                                  star: outlier. 
%                                  Note that the simulations predicted a bias similar to 
%                                  the results from Pulkki (2001) for the reference source at 0 deg.
%
%     'tab1'                       Reproduce Tab.1 of Baumgartner and Majdak (2015):
%                                  Means and standard deviations of responded panning angles for the 
%                                  two reference sources (Ref.) together with corresponding GOFs 
%                                  evaluated for the actual results from Pulkki (2001) and 
%                                  predicted results based on various response strategies. 
%                                  Note the relatively large GOFs for the simulations based on 
%                                  mixed response strategies indicating a good correspondence 
%                                  between actual and predicted results.
%
%     'fig7'                       Reproduce Fig.7 of Baumgartner and Majdak (2015):
%                                  Increase in polar error (defined as in Fig. 5) as a function 
%                                  of loudspeaker span in the median plane with panning ratio 
%                                  R = 0 dB. Black line with gray area indicates mean 
%                                  +-1 standard deviation across listeners. Note that the 
%                                  increase in polar error monotonically increases with 
%                                  loudspeaker span.
%
%     'fig8'                       Reproduce Fig.8 of Baumgartner and Majdak (2015):
%                                  Effect of loudspeaker span in the median plane on coefficient 
%                                  of determination, r^2, for virtual source directions 
%                                  created by VBAP. Separate analysis for frontal, rear, 
%                                  and overall (frontal and rear) targets. Data pooled 
%                                  across listeners. Note the correspondence with the 
%                                  results obtained by Bremen et al. (2010).
%
%     'tab3'                       Reproduce Tab.3 of Baumgartner and Majdak (2015):
%                                  Predicted across-listener average of increase in polar 
%                                  errors as referred to a reference system containing 
%                                  loudspeakers at all considered directions. Distinction 
%                                  between mean and maximum degradation across directions. 
%                                  N: Number of loudspeakers. Ele.: Elevation of second layer. 
%                                  Notice that this elevation has a larger effect on mean 
%                                  and maximum degradation than N.
%
%     'fig9'   Reproduce Fig.9 of Baumgartner and Majdak (2015):
%                                  Predicted polar error as a function of the lateral and 
%                                  polar angle of a virtual source created by VBAP in 
%                                  various multichannel systems. Open circles indicate 
%                                  loudspeaker directions. Reference shows polar error 
%                                  predicted for a real source placed at the virtual 
%                                  source directions investigated for systems A, ..., F.
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
%   2) Data in hrtf/baumgartner2014
%
%   3) Statistics Toolbox for Matlab (for some of the figures)
%
%   Examples:
%   ---------
%
%
%   To display Fig.2 of Baumgartner and Majdak (2015) use :
%
%     exp_baumgartner2015('fig2');
%
%   To display Fig.4 of Baumgartner and Majdak (2015) use :
%
%     exp_baumgartner2015('fig4');
%
%   To display Fig.5 of Baumgartner and Majdak (2015) use :
%
%     exp_baumgartner2015('fig5');
%
%   To display Fig.6 of Baumgartner and Majdak (2015) use :
%
%     exp_baumgartner2015('fig6');
%
%   To display Fig.7 of Baumgartner and Majdak (2015) use :
%
%     exp_baumgartner2015('fig7');
%
%   To display Fig.8 of Baumgartner and Majdak (2015) use :
%
%     exp_baumgartner2015('fig8');
%
%   To display Fig.9 of Baumgartner and Majdak (2015) use :
%
%     exp_baumgartner2015('fig9');
%
%   To display Tab.1 of Baumgartner and Majdak (2015) use :
%
%     exp_baumgartner2015('tab1');
%
%   To display Tab.3 of Baumgartner and Majdak (2015) use :
%
%     exp_baumgartner2015('tab3');
%
%   See also: baumgartner2014 data_baumgartner2014
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%     R. Baumgartner and P. Majdak. Modeling Localization of Amplitude-Panned
%     Virtual Sources in Sagittal Planes. J. Audio Eng. Soc.,
%     63(7/8):562--569, Aug. 2015. [1]http ]
%     
%     P. Bremen, M. M. van Wanrooij, and A. J. van Opstal. Pinna cues
%     determine orientation response modes to synchronous sounds in
%     elevation. J. Neurosci., 30 (1):194--204, 2010.
%     
%     V. Pulkki. Localization of Amplitude-Panned Virtual Sources II: Two-
%     and Three-Dimensional Panning. J. Audio Eng. Soc., 49(4):753--767,
%     2001.
%     
%     References
%     
%     1. http://www.aes.org/e-lib/browse.cfm?elib=17842
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_baumgartner2015.php


%   #Author: Robert Baumgartner (2015)

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
% definput.flags.type = {'missingflag',...% 'fig5_baumgartner2015aro',...
%                        'fig2_baumgartner2015jaes','fig4_baumgartner2015jaes',...
%                        'fig5_baumgartner2015jaes','fig6_baumgartner2015jaes',...
%                        'fig7_baumgartner2015jaes','fig8_baumgartner2015jaes',...
%                        'fig9_baumgartner2015jaes','tab1_baumgartner2015jaes',...
%                        'tab3_baumgartner2015jaes',...
%                        };
definput.flags.type = {'missingflag',...
                       'fig2','fig4',...
                       'fig5','fig6',...
                       'fig7','fig8',...
                       'fig9','tab1',...
                       'tab3',...
                       };
definput.flags.plot = {'plot','no_plot'};


[flags,kv]  = ltfatarghelper({'FontSize','MarkerSize'},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%% General Plot Settings
TickLength = [0.02,0.04];



%% ------------------------------------------------------------------------
%  ---- baumgartner2015jaes -----------------------------------------------
%  ------------------------------------------------------------------------
%% ------ FIG 2 of baumgartner2015jaes ------------------------------------
if flags.do_fig2
  pol1 = 0;
  pol2 = 40;

  s = data_baumgartner2014('pool');

  fs = s(1).Obj.Data.SamplingRate;

  [dtfs,pol] = extractsp(0,s(1).Obj);

  polphant = pol1 + (pol2-pol1)/2;

  dtf1 = 10*dtfs(:,pol==pol1,1);
  dtf2 = 10*dtfs(:,pol==pol2,1);
  dtfreal = 10*dtfs(:,pol==polphant,1);
  dtfphant = (dtf1+dtf2)/2;

  %%
  if flags.do_plot
    
    figure

    set(gca,'LineWidth',1)
    plotfft(fft(dtfphant),fs,'posfreq')
    hold on
    plotfft(fft(dtfreal),fs,'posfreq')
    plotfft(fft(dtf1),fs,'posfreq')
    plotfft(fft(dtf2),fs,'posfreq');
    
    % Set line styles
    Color =  {[0.4660    0.6740    0.1880];...
              [0.3010    0.7450    0.9330];...
              [0.8500    0.3250    0.0980];...
              [     0    0.4470    0.7410]};
    LineWidth = [.5,.5,1,1];
    LineStyle = {':',':','--','-'};
    ch = get(gca,'Children');
    for ii = 1:length(ch)
      set(ch(ii),'Color',Color{ii},'LineStyle',LineStyle{ii},'LineWidth',LineWidth(ii));
    end
    
    leg = legend([num2str(polphant) '\circ VBAP'],...
      [num2str(polphant) '\circ source'],...
      ['  ' num2str(pol1) '\circ source'],...
      [num2str(pol2) '\circ source']);
    set(leg,'Location','southeast','FontSize',kv.FontSize)

    set(gca,'XLim',[500,17500],'YLim',[-39.9,9.9],'FontSize',kv.FontSize)

    xticks = get(gca,'XTick');
    set(gca,'XTickLabel',xticks/1000)
    xlabel({'Frequency (kHz)';' '},'FontSize',kv.FontSize)
    
  end
  
end

%% ------ FIG 4 of baumgartner2015jaes ------------------------------------
if flags.do_fig4
  
  [peI,s,pol0,DL,respang] = amt_cache('get','panningangle',flags.cachemode);
  
  if isempty(peI)
  
    MRS = 0;

    pol1 = -15; % polar angle of lower Lsp.
    pol2 = 30; % polar angle of higher Lsp.

    lat = 0; % must be 0 otherwise VBAP wrong!

    s = data_baumgartner2014('pool');

    dPol = pol2-pol1;
    s(1).spdtfs = [];
    [s(1).spdtfs,polang] = extractsp(lat,s(1).Obj); 
    idtest = find(polang >= pol1 & polang <= pol2);
    qeI = zeros(length(s),length(idtest));
    peI = qeI;

    for ii = 1:length(idtest)

      id1 = idtest(1); % ID lower pos.
      id2 = idtest(end); % ID higher pos.
      id0 = idtest(ii);  % ID pos. of phantom source
      pol0(ii) = polang(id0);

      % VBAP
      [p(1,1),p(1,2),p(1,3)] = sph2cart(lat,deg2rad(pol0(ii)),1);
      [L(1,1),L(1,2),L(1,3)] = sph2cart(lat,deg2rad(pol1),1);
      [L(2,1),L(2,2),L(2,3)] = sph2cart(lat,deg2rad(pol2),1);
      g = p/L;
      g = g/norm(g);

      DL(ii) = db(g(1)) - db(g(2));

      for ll = 1:length(s)

          s(ll).spdtfs = extractsp(lat,s(ll).Obj);

          % superposition
          s(ll).dtfs2{ii} = g(1)*s(ll).spdtfs(:,id1,:) + g(2)*s(ll).spdtfs(:,id2,:);

          [s(ll).p1(:,ii),respang] = baumgartner2014(...
            s(ll).spdtfs(:,id0,:),s(ll).spdtfs,s(ll).fs,'S',s(ll).S,...
            'mrsmsp',MRS,'polsamp',polang);
          s(ll).p2(:,ii) = baumgartner2014(...
            s(ll).dtfs2{ii},s(ll).spdtfs,s(ll).fs,'S',s(ll).S,...
            'mrsmsp',MRS,'polsamp',polang);

          [s(ll).qe1(ii),s(ll).pe1(ii)] = baumgartner2014_pmv2ppp(...
            s(ll).p1(:,ii),polang(id0),respang);
          [s(ll).qe2(ii),s(ll).pe2(ii)] = baumgartner2014_pmv2ppp(...
            s(ll).p2(:,ii),polang(id0),respang);

          % Increse of error
          qeI(ll,ii) = s(ll).qe2(ii) - s(ll).qe1(ii);
          peI(ll,ii) = s(ll).pe2(ii) - s(ll).pe1(ii);
      end
      amt_disp([num2str(ii) ' of ' num2str(length(idtest)) ' completed']);
    end

    amt_cache('set','panningangle',peI,s,pol0,DL,respang);
    
  end
  
  id1 = 'NH71'; % positive example
  id2 = 'NH62'; % negative example
  
  if flags.do_plot
    
    cmp= [0.2081    0.1663    0.5292;
          0.2052    0.2467    0.6931;
          0.0843    0.3472    0.8573;
          0.0157    0.4257    0.8789;
          0.0658    0.4776    0.8532;
          0.0777    0.5300    0.8279;
          0.0356    0.5946    0.8203;
          0.0230    0.6443    0.7883;
          0.0485    0.6793    0.7341;
          0.1401    0.7085    0.6680;
          0.2653    0.7327    0.5916;
          0.4176    0.7471    0.5142;
          0.5624    0.7487    0.4529;
          0.6872    0.7433    0.4029;
          0.7996    0.7344    0.3576;
          0.9057    0.7261    0.3105;
          0.9944    0.7464    0.2390;
          0.9847    0.8141    0.1734;
          0.9596    0.8869    0.1190;
          0.9763    0.9831    0.0538]; % parula colormap (defined for compatibility with older Matlab versions)
    
    figure

    p_pool = nan(size(s(1).p2,1),size(s(1).p2,2),length(s));
    for ll = 1:length(s)
      p_pool(:,:,ll) = s(ll).p2;
      if strcmp(s(ll).id,id1)
        subplot(1,4,1)
        plot_baumgartner2014(s(ll).p2,pol0,respang,'cmax',0.08)
        colormap(cmp)
        title(s(ll).id,'FontSize',kv.FontSize)
        xlabel(''); 
        ylabel('Response angle (deg)','FontSize',kv.FontSize);
        colorbar off
      elseif strcmp(s(ll).id,id2)
        subplot(1,4,2)
        plot_baumgartner2014(s(ll).p2,pol0,respang,'cmax',0.08)
        colormap(cmp)
        xlabel('Panning angle (deg)','FontSize',kv.FontSize);      
        ylabel(''); set(gca,'YTickLabel',[]);
        title(s(ll).id,'FontSize',kv.FontSize)
        colorbar off
      end
    end
    p_pool = mean(p_pool,3);
    subplot(1,4,3)
    plot_baumgartner2014(p_pool,pol0,respang,'cmax',0.08)
    colormap(cmp)
    title('Pool','FontSize',kv.FontSize)
    xlabel(''); 
    ylabel(''); set(gca,'YTickLabel',[]);
    colorbar off

    subplot(1,4,4)
    ymax = 8;
    y = -0.15:0.1:ymax+0.15;
    Ny = length(y);
    pcolor(1:2,y,repmat(y(:),1,2))
    shading flat
    axis tight
    colormap(cmp)
    title({' ';' '})
    set(gca,'XTick',[],'YTick',0:1:ymax,... 
      'YDir','normal','YAxisLocation','right','FontSize',kv.FontSize)
    ylabel({'Probability density (% per 5\circ)'},'FontSize',kv.FontSize)
    set(gca,'Position',get(gca,'Position').*[1.05,1,0.3,1])
  
  end
  
end

%% ------ FIG 5 of baumgartner2015jaes ------------------------------------
if flags.do_fig5
  
  [peI,s,pol0,DL,respang] = amt_cache('get','panningangle',flags.cachemode);
  
  if isempty(peI)
    exp_baumgartner2015('fig4','no_plot',flags.cachemode);
    [peI,s,pol0,DL,respang] = amt_cache('get','panningangle',flags.cachemode);
  end
  
  figure

  plot(peI')
  ylabel('Increase in polar error (deg)')

  % Panning angle axis
  set(gca,'XTick',1:10,'XTickLabel',round(pol0),...
      'YMinorTick','on','XLim',[1,10],'YLim',[-9,59])
  xlabel('Panning angle (deg)')
  
end

%% ------ FIG 6 of baumgartner2015jaes ------------------------------------
if flags.do_fig6
  
  results = amt_cache('get','replicatePulkki2001',flags.cachemode);
   
  if isempty(results)
  
    amt_disp('Results may slightly vary from simulation to simulation because noise stimulus is not fixed.');
    
    MRS = 0;

    pol1 = -15; % polar angle of lower Lsp.
    pol2 = 30; % polar angle of higher Lsp.

    polphant = [0,15]; % polar angles of phantom sources
    Pmax = nan(length(polphant),23); % max Probabilities
    panang_Pmax = nan(length(polphant),23); % panning angle selected by max P
    panang_Cen = nan(length(polphant),23); % panning angle selected by centroid
    for pp = 1:length(polphant)

    lat = 0; % must be 0 otherwise VBAP wrong!

    s = data_baumgartner2014('pool');

    s(1).spdtfs = [];
    [s(1).spdtfs,polang] = extractsp(lat,s(1).Obj); 

    % restrict response range
    idrang = find(polang >= pol1 & polang <= pol2); % to range between loudspeakers
    % idrang = find(polang <= 90); % to the front

    idtest = find(polang >= pol1 & polang <= pol2);
    id1 = find(polang >= pol1,1); % ID lower pos.
    id2 = find(polang <= pol2,1,'last'); % ID higher pos.
    tang = polang(id1:id2);

    for ll = 1:length(s)
      s(ll).spdtfs = extractsp(lat,s(ll).Obj);

      for ii = 1:length(idtest)

        id0 = idtest(ii);  % ID pos. of phantom source

        % VBAP
        [p(1,1),p(1,2),p(1,3)] = sph2cart(lat,deg2rad(polang(id0)),1);
        [L(1,1),L(1,2),L(1,3)] = sph2cart(lat,deg2rad(pol1),1);
        [L(2,1),L(2,2),L(2,3)] = sph2cart(lat,deg2rad(pol2),1);
        g = p/L;
        g = g/norm(g);

        DL(ii) = db(g(1)) - db(g(2));

        % superposition
        s(ll).dtfs2{ii} = g(1)*s(ll).spdtfs(:,id1,:) + g(2)*s(ll).spdtfs(:,id2,:);    

        [s(ll).p(:,ii),rang] = baumgartner2014(...
              s(ll).dtfs2{ii},s(ll).spdtfs(:,idrang,:),s(ll).fs,'S',s(ll).S,...
              'mrsmsp',MRS,'polsamp',polang(idrang),'rangsamp',5,...
              'stim',noise(10000,1,'pink')); % phantom source

      end

      id_rang = rang == polphant(pp);

      % interpolation between target angles
      tang_int = tang(1):1:tang(end);
      p_int = interp2(tang(:)',rang(:),s(ll).p,tang_int(:)',rang(:),'spline');
      p_int = p_int./repmat(sum(p_int,1),size(p_int,1),1); % normalize to PMVs

      % Variant 1: max P at source direction
      [Pmax(pp,ll),id_best_pan] = max(p_int(id_rang,:));
      panang_Pmax(pp,ll) = tang_int(id_best_pan);

      % Variant 2: centroid closest to source direction
      M = rang*p_int;
      [tmp,id_best_pan] = min(abs(M-polphant(pp)));
      panang_Cen(pp,ll) = tang_int(id_best_pan);

    end
      fprintf([num2str(pp) ' of ' num2str(length(polphant)) ' completed \n']);
    end

    results = struct('panang_Pmax',panang_Pmax,'Pmax',Pmax,...
      'panang_Cen',panang_Cen,'polphant',polphant,'DL',DL,'rang',rang);
    
    amt_cache('set','replicatePulkki2001',results);
    
  end
  
  [panang_varStrat,nCM_varStrat,p_varStrat,muhat,sigmahat] = amt_cache('get','replicatePulkki2001_varStrat',flags.cachemode);
  if isempty(panang_varStrat)
    
    pulkki01 = data_pulkki2001;
    [muhat(1),sigmahat(1)] = normfit(pulkki01(1,:));
    [muhat(2),sigmahat(2)] = normfit(pulkki01(2,:));
    
    Nsub = size(results.panang_Cen,2);
    panang_all = [];
    p = [];
    nCM = [];
    tmp = [];
    ii = 1;
    for inCM = 1:Nsub+1
      c = nchoosek(1:Nsub,inCM-1); % listeners with CM strategy
      lenC = size(c,1);
      nCM = [nCM;repmat(inCM,lenC,1)];
      panang_all = cat(3,panang_all , nan(2,Nsub,lenC));
      p = [p ; nan(lenC,2)];
      for ic = 1:lenC
        idCM = false(1,Nsub);
        idCM(c(ic,:)) = true;
        panang_all(:,:,ii) = [results.panang_Cen(:,idCM) results.panang_Pmax(:,not(idCM))];
        [tmp.h1,p(ii,1)] = kstest((panang_all(1,:,ii)-muhat(1))/sigmahat(1)); % center data acc. to target distribution and then test similarity to standard normal distribution
        [tmp.h2,p(ii,2)] = kstest((panang_all(2,:,ii)-muhat(2))/sigmahat(2));
        ii = ii+1;
      end
      disp([num2str(inCM) ' of ' num2str(Nsub+1) ' done'])
    end
    [tmp.min,idmax] = max(sum(p,2)); % best fit
    panang_varStrat = panang_all(:,:,idmax);
    nCM_varStrat = nCM(idmax);
    p_varStrat = p(idmax,:);
    
    amt_cache('set','replicatePulkki2001_varStrat',panang_varStrat,nCM_varStrat,p_varStrat,muhat,sigmahat)
    
  end
  
  if flags.do_plot
    
    pulkki01 = data_pulkki2001;
    Nsub = size(results.panang_Cen,2);
    
    figure
    for ii = 1:2

      subplot(1,2,ii)

      X = nan(size(pulkki01,2)*size(pulkki01,3),4);
      X(:,1) = pulkki01(ii,:);
      X(1:Nsub,2) = results.panang_Pmax(ii,:)';
      X(1:Nsub,3) = results.panang_Cen(ii,:)';

      X(1:Nsub,4) = panang_varStrat(ii,:)';

      plot([0,5],(ii-1)*15*[1,1],'k:') 
      hold on

      boxplot(X,'symbol','k*','outliersize',3)

      set(gca,'YLim',[-17,32], 'XTickLabel',{'[2]','PM','CM','Mixed'});
	  if ~verLessThan('matlab','8.4'), set(gca,'XTickLabelRotation',45); end
      if ii==1; 
        ylabel('Panning angle (deg)')
        text(0.7,27.5,'0\circ')
      else
        set(gca,'YTickLabel',[]); 
        text(0.7,27.5,'15\circ')
      end
    end
    
  end
  
end

%% ------ FIG 7 of baumgartner2015jaes ------------------------------------
if flags.do_fig7
  
  [peI,dPol] = amt_cache('get','loudspeakerspan',flags.cachemode);
  
  if isempty(peI)
    
    MRS = 0;

    flags.do_fig20 = false;
    flags.do_fig19 = false;

    s = data_baumgartner2014('pool');

    if flags.do_fig19
      dPol = [0 30,60];
      s = s(2); % NH12
    else
      dPol = 10:10:90; 
    end
    lat = 0;

    s(1).spdtfs = [];
    [s(1).spdtfs,polang] = extractsp(lat,s(1).Obj); 
    peI = zeros(length(s),length(dPol));
    peA = zeros(length(s),length(dPol)+1);
    ii = 0;
    while ii < length(dPol)
      ii = ii + 1;

      % find comparable angles
      id0 = [];
      id1 = [];
      id2 = [];
      for jj = 1: length(polang)
          t0 = find( round(polang) == round(polang(jj)+dPol(ii)/2) );
          t2 = find( round(polang) == round(polang(jj)+dPol(ii)) );
          if ~isempty(t0) && ~isempty(t2)
              id0 = [id0 t0];
              id1 = [id1 jj];
              id2 = [id2 t2];
          end
      end
      pol2{ii} = (polang(id1)+polang(id2)) /2;

      amt_disp([' Span: ' num2str(dPol(ii)) 'deg']);
      for ll = 1:length(s)

          s(ll).spdtfs = extractsp(lat,s(ll).Obj);

          % superposition
          s(ll).dtfs2{ii} = s(ll).spdtfs(:,id1,:) + s(ll).spdtfs(:,id2,:);

          [s(ll).p1{ii},respang] = baumgartner2014(...
            s(ll).spdtfs(:,id0,:),s(ll).spdtfs,s(ll).fs,'S',s(ll).S,...
            'mrsmsp',MRS,'polsamp',polang);
          s(ll).p2{ii} = baumgartner2014(...
            s(ll).dtfs2{ii},s(ll).spdtfs,s(ll).fs,'S',s(ll).S,...
            'mrsmsp',MRS,'polsamp',polang);

          % RMS Error
          [s(ll).qe1{ii},s(ll).pe1{ii}] = baumgartner2014_pmv2ppp(...
            s(ll).p1{ii},polang(id0),respang);
          [s(ll).qe2{ii},s(ll).pe2{ii}] = baumgartner2014_pmv2ppp(...
            s(ll).p2{ii},pol2{ii},respang);

          % Increse of error
          peI(ll,ii) = s(ll).pe2{ii} - s(ll).pe1{ii};

      end

    end

    amt_cache('set','loudspeakerspan',peI,dPol)
    
  end
    
  peIm = mean(peI,1);
  peIstd = std(peI,1,1);
    
  if flags.do_plot
    
    figure
    p = patch([dPol,fliplr(dPol)],[peIm+peIstd,fliplr(peIm-peIstd)],.7*[1,1,1]);
    set(p,'EdgeColor',.7*[1,1,1]);
    hold on
    h = plot(dPol,peIm,'k');
    set(h,'LineWidth',1)
    set(gca,'XTick',dPol,'XTickLabel',dPol,...
          'YMinorTick','on','Box','on','Layer','top')
    axis([9.8,90.2,-2,27])
    ylabel({'Increase in polar error (deg)'})
    xlabel({'Loudspeaker span (deg)';''})
    
  end
  
end

%% ------ FIG 8 of baumgartner2015jaes ------------------------------------
if flags.do_fig8
  
  [r2,dPol] = amt_cache('get','loudspeakerspan_r2',flags.cachemode);
  
  if isempty(r2)
    MRS = 0;

    s = data_baumgartner2014('pool');

    dPol = 0:10:105; 

    lat = 0;
    runs = 100;

    
    DL = -13:5:7; % panning ratios in dB
    
    s(1).spdtfs = [];
    [s(1).spdtfs,polang] = extractsp(lat,s(1).Obj); 
    r2 = zeros(length(s),length(dPol));
    r2 = [];
    ii = 0;
    while ii < length(dPol) % various spans
      ii = ii + 1;

      disp([' Span: ' num2str(dPol(ii)) 'deg']);

      r2total = nan(length(s),length(DL));
      r2front = nan(length(s),length(DL));
      r2rear = nan(length(s),length(DL));
      for idl = 1:length(DL)

      % find comparable angles
      id1 = []; % id of speaker with smaller polar angle
      id2 = []; % id of speaker with larger polar angle
      for jj = 1: length(polang)
%           t0 = find( round(polang) == round(polang(jj)+ipf*dPol(ii)/polFrac) );
          t2 = find( round(polang) == round(polang(jj)+dPol(ii)) );
          if ~isempty(t2)
%               id0 = [id0 t0];
              id1 = [id1 jj];
              id2 = [id2 t2];
          end
      end

      g = inv([1,1;1,-10^(DL(idl)/20)]) * [1;0]; % derived from db(g1/g2)=DL and g1+g2=1 (chosen arbitrarily since energy preservation is not relevant here)
      
      pol0 = nan(length(id1),1); % panning angles
      for jj = 1:length(id1)
        [L(1,1),L(1,2),L(1,3)] = sph2cart(0,deg2rad(polang(id1(jj))),1);
        [L(2,1),L(2,2),L(2,3)] = sph2cart(0,deg2rad(polang(id2(jj))),1);
        t = g'*L;
        [azi,ele,tmp.r] = cart2sph(t(1),t(2),t(3));
        [tmp.lat,pol0(jj)] = sph2hor(rad2deg(azi),rad2deg(ele));
      end

      for ll = 1:length(s) % various listeners

          s(ll).spdtfs = extractsp(lat,s(ll).Obj);

          % superposition
          s(ll).dtfs2{ii} = g(1)*s(ll).spdtfs(:,id1,:) + g(2)*s(ll).spdtfs(:,id2,:);

          [s(ll).p2{ii},rang] = baumgartner2014(...
            s(ll).dtfs2{ii},s(ll).spdtfs,s(ll).fs,'S',s(ll).S,...
            'mrsmsp',MRS,'polsamp',polang,'rangsamp',1); % phantom

          % total polar range 
          m2 = baumgartner2014_virtualexp(s(ll).p2{ii},pol0,rang,'runs',runs);

          % R2 via correlation coefficient
          r = corrcoef(m2(:,8),m2(:,6));
          r2total(ll,idl) = r(2);

          % restricted to frontal polar range 
          idfront = rang <= 90;
          respangfront = rang(idfront);
          idpol0front = pol0 <= 90;
          pol0front = pol0(idpol0front);
          p = s(ll).p2{ii}(idfront,idpol0front);

          m2front = baumgartner2014_virtualexp(p,pol0front,respangfront,'runs',runs);

          r = corrcoef(m2front(:,8),m2front(:,6));
          r2front(ll,idl) = r(2);

          % restricted to rear polar range 
          idrear = rang >= 90;
          respangrear = rang(idrear);
          idpol0rear = pol0 >= 90;
          pol0rear = pol0(idpol0rear);
          p = s(ll).p2{ii}(idrear,idpol0rear);

          m2rear = baumgartner2014_virtualexp(p,pol0rear,respangrear,'runs',runs);

          r = corrcoef(m2rear(:,8),m2rear(:,6));
          r2rear(ll,idl) = r(2);

      end
      end
      r2.total(:,ii) =  nanmean(r2total,2);
      r2.front(:,ii) =  nanmean(r2front,2); 
      r2.rear(:,ii) =  nanmean(r2rear,2);
    end
 
    amt_cache('set','loudspeakerspan_r2',r2,dPol)

  end
    
  % data extracted from Fig.6 (data:AVG) of Bremen et al. (2010, J Neurosci,
  % 30:194-204) via http://arohatgi.info/WebPlotDigitizer
  bremen2010.pol = 15:15:105;
  bremen2010.r2 = [.85 , .81 , .63 , .38 , .19 , .11 , .21];
  
  if flags.do_plot

    figure
    h(1) = plot(bremen2010.pol,bremen2010.r2,'bo-');
    hold on
    h(2) = plot(dPol,mean(r2.front),'bo-');

    h(3) = plot(dPol,mean(r2.rear),'rs-');
    h(4) = plot(dPol,mean(r2.total),'kd-');

    set(h(1),'MarkerSize',kv.MarkerSize,'MarkerFaceColor','b')
    set(h(2:4),'MarkerSize',kv.MarkerSize,'MarkerFaceColor','w')
    set(gca,'XLim',[-5,110],'YLim',[-.05,1.05])

    leg = legend('Frontal from [18]','Frontal','Rear','Overall','Location','best');
    set(leg,'FontSize',kv.FontSize)

    xlabel({'Loudspeaker span (deg)';''})
    ylabel('\it{r}^{ 2}')
  end
  
end

%% ------ FIG 9 of baumgartner2015jaes ------------------------------------
if flags.do_fig9
  
  SysName{1}  = 'NHK 22.2 (without bottom layer)';
  LSPsetup{1} = [ 0,0 ; 30,0 ; 60,0 ;  90,0 ; 135,0 ; ...
                180,0 ;-30,0 ;-60,0 ; -90,0 ;-135,0 ; ...
                  0,45; 45,45; 90,45; 135,45; 180,45; ...
               -135,45;-90,45;-45,45;   0,90];

  SysName{2}  = 'Samsung 11.2';
  LSPsetup{2} = [ 0,0 ; 60,0 ;  90,0 ; 135,0 ; ...
                       -60,0 ; -90,0 ;-135,0 ; ...
                 45,45;135,45; -45,45;-135,45];

  SysName{3}  = 'Samsung 10.2';
  LSPsetup{3} = [ 0,0 ; 60,0 ;  90,0 ; 135,0 ; ...
                       -60,0 ; -90,0 ;-135,0 ; ...
                 45,45;180,45; -45,45];

  SysName{4}  = 'USC 10.2';
  LSPsetup{4} = [ 0,0 ; 30,0 ; 60,0 ;  115,0 ; ...
                180,0 ;-30,0 ;-60,0 ; -115,0 ; ...
                 45,45;-45,45];

  SysName{5}  = 'Auro-3D 10.1';
  LSPsetup{5} = [ 30,0 ; 30,30 ; 135,30 ; 135,0;...
            0,0 ;-30,0 ;-30,30 ;-135,30 ;-135,0; 0,90];        

  SysName{6}  = 'Auro-3D 9.1';
  LSPsetup{6} = [ 30,0 ; 30,30 ; 135,30 ; 135,0;...
            0,0 ;-30,0 ;-30,30 ;-135,30 ;-135,0]; 
  
  latall = -45:5:45;
  polall = 0:10:180;
          
  pe = amt_cache('get','locaVBAP',flags.cachemode);
  
  if isempty(pe)
    
    MRS = 0;

    s = data_baumgartner2014('pool');
    
    pe = zeros(length(latall),length(polall),length(s),length(LSPsetup),2); % predicted local polar RMS errors
    for ll = 1:length(LSPsetup)
      for aa = 1:length(latall)
        lat = latall(aa);

        for pp = 1:length(polall)
          pol = polall(pp);

          % Select LSP-Triangle and compute VBAP gains
          [source_pos(1),source_pos(2)] = hor2sph(lat,pol);
          Nlsp = length(LSPsetup{ll});
          [g,IDsp] = local_vbap(LSPsetup{ll},source_pos);


          for jj = 1:length(s)

            % Compute binaural impulse response of loudspeaker triple
            dtfs = permute(double(s(jj).Obj.Data.IR),[3 1 2]);
            lsp_hrirs = zeros(length(IDsp),size(dtfs,1),2);
            for ii = 1:length(IDsp)
    %           [lat_sp,pol_sp] = sph2horpolar(LSPsetup{ll}(IDsp(ii),1),LSPsetup{ll}(IDsp(ii),2));
              idx = local_findNearestPoslocaVBAP(...
                LSPsetup{ll}(IDsp(ii),1:2),s(jj).Obj.SourcePosition(:,1:2));
              lsp_hrirs(ii,:,:) = squeeze(dtfs(:,idx,:));
            end
            target(:,1,1) = g*lsp_hrirs(:,:,1);
            target(:,1,2) = g*lsp_hrirs(:,:,2);

            % SP-template
            [spdtfs,polang] = extractsp(lat,s(jj).Obj);

            % Run loca model
            [p,rang] = baumgartner2014(...
                    target,spdtfs,s(jj).fs,'S',s(jj).S,...
                    'lat',lat,'polsamp',polang,'mrsmsp',MRS);

            m = baumgartner2014_virtualexp(p,pol,rang,'runs',1000);
            pe(aa,pp,jj,ll,1) = localizationerror(m,'precPnoquerr');
            [~,pe(aa,pp,jj,ll,2)] = baumgartner2014_pmv2ppp(p,pol,rang);

          end
        end
      end
      amt_disp([num2str(ll) ' of ' num2str(length(LSPsetup)) ' done']);
    end
    amt_cache('set','locaVBAP',pe)
  end
 
  MRS = 0;
  pe_ref = amt_cache('get','locaVBAP_ref',flags.cachemode);
  if isempty(pe_ref)
    
    s = data_baumgartner2014('pool');

    latall = -45:5:45;
    polall = 0:10:180;
    pe_ref = zeros(length(latall),length(polall),length(s),1,2); % predicted local polar RMS errors

    %% Computations

    for jj = 1:length(s)

      dtfs = permute(double(s(jj).Obj.Data.IR),[3 1 2]);

      for aa = 1:length(latall)
        lat = latall(aa);
        for pp = 1:length(polall)
          pol = polall(pp);
          [lat_sp,pol_sp,idx] = local_findNearestPoslocaVBAPref(lat,pol,s(jj).Obj.SourcePosition(:,1:2));
          target = dtfs(:,idx,:);

          % SP-template
          [spdtfs,polang] = extractsp(lat,s(jj).Obj);

          % Run loca model
          [p,rang] = baumgartner2014(...
                  target,spdtfs,s(jj).fs,'S',s(jj).S,...
                  'mrsmsp',MRS,'lat',lat,'polsamp',polang);

          m = baumgartner2014_virtualexp(p,pol,rang,'runs',1000);
          pe_ref(aa,pp,jj,1,1) = localizationerror(m,'precPnoquerr');
          [~,pe_ref(aa,pp,jj,1,2)] = baumgartner2014_pmv2ppp(p,pol,rang);

        end
      end
    end
    
    amt_cache('set','locaVBAP_ref',pe_ref)
  end

  N = length(LSPsetup);
  eRMS = pe_ref(:,:,:,:,2);
  eRMS(:,:,:,2:N+1) = pe(:,:,:,:,2);

  if flags.do_plot

    labels = {'Reference';'\it A';'\it B';'\it C';'\it D';'\it E';'\it F'};
    labels = labels(1:N+1,:);
    
    figure 
    for ll = 1:N+1

      pemean = squeeze(mean(eRMS(:,:,:,ll),3));
      subplot(1,N+2,ll)
      imagesc(latall,polall,pemean')
      set(gca,'YTick',0:30:180,'XTick',-30:30:30)
      axis equal tight
      colormap hot
      if MRS == 0
        ymin = 15.5;
        ymax = 49.5;
      else
        ymin = 15.5;
        ymax = 49.5;
      end
      caxis([ymin ymax])
      if ll==1; 
        ylabel('Polar angle (deg)','FontName','Helvetica'); 
      else 
        set(gca,'YTickLabel',[])
      end

      if ll==4; xlabel('Lateral angle (deg)','FontName','Helvetica'); end

      title(labels{ll},'FontName','Helvetica')

      % Loudspeaker positions
      if ll > 1
        [lat_lsp,pol_lsp] = sph2hor(LSPsetup{ll-1}(:,1),LSPsetup{ll-1}(:,2));
        idlat = abs(lat_lsp) <= max(abs(latall))+1;
        hold on
        h2 = plot(lat_lsp(idlat),pol_lsp(idlat),'wo');
        set(h2,'MarkerSize',2*3.5);
        h1 = plot(lat_lsp(idlat),pol_lsp(idlat),'ko');
        set(h1,'MarkerSize',2*3);
      end
    end

    subplot(1,N+2,N+2)
    y = ymin:1:ymax;
    Ny = length(y);
    pcolor(1:2,y,repmat(y(:),1,2))
    shading flat
    axis tight
    colormap hot
    title({' ';' '})
    set(gca,'XTick',[],'YTick',20:5:45,... %,'TickLength',[0.12,0.12]
      'YDir','normal','YAxisLocation','right','FontSize',kv.FontSize)
    ylabel({'Polar error (deg)'},'FontSize',kv.FontSize)
    set(gca,'Position',get(gca,'Position').*[1.03,1.8,0.3,0.8])
  end
  
end

%% ------ Tab 1 of baumgartner2015jaes ------------------------------------
if flags.do_tab1
  
  results = amt_cache('get','replicatePulkki2001',flags.cachemode);
  [panang_varStrat,nCM_varStrat,p_varStrat,muhat,sigmahat] = amt_cache('get','replicatePulkki2001_varStrat',flags.cachemode);
  if isempty(panang_varStrat)
    exp_baumgartner2015('fig6','no_plot',flags.cachemode)
  end
  
  pulkki01 = data_pulkki2001;
  
  [h1_pul,p1_pul] = kstest((pulkki01(1,:)-muhat(1))/sigmahat(1)); % center data acc. to target distribution and then test similarity to standard normal distribution
  [h2_pul,p2_pul] = kstest((pulkki01(2,:)-muhat(2))/sigmahat(2));
  [h1_pm,p1_pm] = kstest((results.panang_Pmax(1,:)-muhat(1))/sigmahat(1)); % center data acc. to target distribution and then test similarity to standard normal distribution
  [h2_pm,p2_pm] = kstest((results.panang_Pmax(2,:)-muhat(2))/sigmahat(2));
  [h1_cm,p1_cm] = kstest((results.panang_Cen(1,:)-muhat(1))/sigmahat(1)); % center data acc. to target distribution and then test similarity to standard normal distribution
  [h2_cm,p2_cm] = kstest((results.panang_Cen(2,:)-muhat(2))/sigmahat(2));

  amt_disp('p-values of K.S.-test:','documentation');
  amt_disp('Real source at 0 deg:','documentation');
  amt_disp(['  Results Pulkki (2001): p = ' num2str(p1_pul,'%3.2f')],'documentation');
  amt_disp(['  Probability Maximiz.:  p = ' num2str(p1_pm,'%3.2f')],'documentation');
  amt_disp(['  Centroid Match:        p = ' num2str(p1_cm,'%3.2f')],'documentation');
  amt_disp(['  Individual strategy:   p = ' num2str(p_varStrat(1),'%3.2f') ' (#CM = ' num2str(nCM_varStrat) ')'],'documentation');
  amt_disp('Real source at 15 deg:','documentation');
  amt_disp(['  Results Pulkki (2001): p = ' num2str(p2_pul,'%3.2f')],'documentation');
  amt_disp(['  Probability Maximi.:   p = ' num2str(p2_pm,'%3.2f')],'documentation');
  amt_disp(['  Centroid Match:        p = ' num2str(p2_cm,'%3.2f')],'documentation');
  amt_disp(['  Individual strategy:   p = ' num2str(p_varStrat(2),'%3.2f') ' (#CM = ' num2str(nCM_varStrat) ')'],'documentation');
  
end

%% ------ Tab 3 of baumgartner2015jaes ------------------------------------
if flags.do_tab3
  
  pe = amt_cache('get','locaVBAP',flags.cachemode);
  pe_ref = amt_cache('get','locaVBAP_ref',flags.cachemode);
  if isempty(pe)
    exp_baumgartner2014('fig9_baumgartner2015jaes','no_plot',flags.cachemode);
  end
  
  N = size(pe,4); % # loudspeakers
  eRMS = pe_ref(:,:,:,:,2);
  eRMS(:,:,:,2:N+1) = pe(:,:,:,:,2);
  
  labels = {'Reference';'\it A';'\it B';'\it C';'\it D';'\it E';'\it F'};
  labels = labels(1:N+1,:);

  amt_disp('RMS error difference from reference averaged across directions','documentation');
  amt_disp('System  min  mean  max','documentation');
  for ll = 2:N+1
    IeRMS = eRMS(:,:,:,ll) - eRMS(:,:,:,1);
    IeRMS = mean(IeRMS,3); % average across listeners
    amt_disp([labels{ll} ' ' num2str(min(IeRMS(:)),'%2.1f') ' ' num2str(mean(IeRMS(:)),'%2.1f') ' ' num2str(max(IeRMS(:)),'%2.1f')],'documentation');
  end
  
end

end



%% ------------------------------------------------------------------------
%  ---- INTERNAL FUNCTIONS ------------------------------------------------
%  ------------------------------------------------------------------------


function [idx,posN] = local_findNearestPoslocaVBAP(posdesired,posexist)
% FINDNEARESTPOS_LOCAVABAP finds nearest position. Data given in spherical
% coordinates
%
% Usage:    [idx,posN] = findNearestPos(posdesired,posexist)

ds = deg2rad(posdesired);
es = deg2rad(posexist);

[d(:,1),d(:,2),d(:,3)] = sph2cart(ds(:,1),ds(:,2),ones(size(ds,1),1));
[e(:,1),e(:,2),e(:,3)] = sph2cart(es(:,1),es(:,2),ones(size(es,1),1));

D = e-repmat(d,length(es),1);
[Dmin,idx] = min(sum(D.^2,2));

posN = posexist(idx,:);

end

function [latN,polN,idx] = local_findNearestPoslocaVBAPref(lat,pol,aziele)
% FINDNEARESTPOS_LOCAVBAP_REF finds nearest position according to lat. and pol. angle
%
% Usage:    [latN,polN,idx] = findNearestPos(lat,pol,aziele)

[positions(:,1),positions(:,2)] = sph2hor(aziele(:,1),aziele(:,2));
d_lat = abs(lat-positions(:,1));
d_pol = abs(pol-positions(:,2));
d = d_lat+d_pol;
[d_min,idx] = min(d);
latN = positions(idx,1);
polN = positions(idx,2);

end

function [g,IDspeaker] = local_vbap(lsp_coord,source_pos)
%VBAP Returns array of gains for VBAP triplet of speakers.         
%   Usage: [g,IDspeaker] = vbap(speaker_coord,source_pos)
%
%   Input Parameters:
%     lsp_coord  : spherical (azi,ele) coordinates of loudspeakers
%     source_pos : spherical (azi,ele) coordinates of phantom source
%
%   Output Parameters:
%     g          : VBAP gains of loudspeaker triplet
%     IDspeaker  : indices of selected loudspeaker triplet

Nlsp = size(lsp_coord,1);

% Convert to spherical coordinates winth angles in radians
speaker_sph = deg2rad(lsp_coord);
source_sph = deg2rad(source_pos);

% Convert to cartesian coordinates
[speaker_cart(:,1),speaker_cart(:,2),speaker_cart(:,3)] = sph2cart(...
  speaker_sph(:,1),speaker_sph(:,2),ones(Nlsp,1));
[source_cart(1),source_cart(2),source_cart(3)] = sph2cart(...
  source_sph(1),source_sph(2),1);

% Add imaginary speaker below
if min(lsp_coord(:,2)) >= 0
  speaker_cart = [speaker_cart;0,0,-10];
end

% Select lsp. triplet
dt = DelaunayTri(speaker_cart);
ch = convexHull(dt);
% figure; trisurf(ch, dt.X(:,1),dt.X(:,2),dt.X(:,3), 'FaceColor', 'cyan')
d = zeros(length(ch),1);
for ii = 1:length(ch)
  d(ii) = sum(local_dist(source_cart,speaker_cart(ch(ii,:),:)));
end
[~,IDch] = min(d);
IDspeaker = ch(IDch,:);

% Compute lsp. gains
g = source_cart / speaker_cart(IDspeaker,:);
g = g / norm(g);


end

function d = local_dist(x,Y)

d = sqrt(sum((repmat(x,size(Y,1),1)-Y).^2,2));

end


