function varargout=exp_baumgartner2013(varargin)
%EXP_BAUMGARTNER2013   Figures from Baumgartner et al. (2013)
%   Usage: data = exp_baumgartner2013(flag)
%
%   EXP_BAUMGARTNER2013(flags) reproduces figures of the book 
%   chapter from Baumgartner et al. (2013).
%
%   Optional fields of output data structure:
%
%   data.id
%      listener ID
%
%   data.u
%      listener-specific uncertainty
% 
%   data.dtfs
%      matrix containing DTFs. Dimensions: time, position, channel 
%      (more details see doc: HRTF format)
%
%   data.fs
%      sampling rate of impulse responses
%
%   data.pos
%      source-position matrix referring to 2nd dimension of 
%      hM and formated acc. to meta.pos (ARI format).
%      6th col is the lateral angle. 7th col is the polar angle.
% 
%   data.spdtfs
%      DTFs of specific SPs
%  
%   data.polangs
%      polar angles corresponding to data.spdtfs*
%   
%   data.pmv
%      predicted polar angular response PMVs
%
%   data.respangs
%      polar response angules corresponding to data.p*
%   
%   data.pe
%      predicted local polar RMS error in degrees
% 
%   data.qe
%      predicted quadrant error
%
%   The following flags can be specified:
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%     'fig12'   Reproduce Fig. 12:
%               DTFs of median SP (left ear). 
%               Left: NH12. Center: NH58. Right: NH33. 
%               Brightness: Spectral magnitude.
%
%     'fig13'   Reproduce Fig. 13:
%               Listener's localization performance for non-individualized 
%               versus individualized DTFs. 
%               Bars: Individualized DTFs. 
%               Circles: Non-individualized DTFs averaged over 16 DTF sets. 
%               Error bars: Â±1 standard deviation of the average. 
%               Horizontal line: Chance performance corresponding to 
%               guessing the direction of the sound.
%
%     'fig14'   Reproduce Fig. 14:
%               Bars: Increase in localization errors when listening to a 
%               binaural recording created with the DTFs from NH58. 
%               Asterisk: Difference between actual (experiment) and 
%               predicted performance when listening to individualized DTFs. 
%               Horizontal lines: chance performance corresponding to 
%               guessing the direction of the sound.
%
%     'fig15'   Reproduce Fig. 15:
%               Bars: Localization performance of the pool listening to 
%               selected DTFs. 
%               Circles: DTFs from NH12. 
%               Squares: DTFs from NH58. 
%               Triangles: DTFs from NH33.
%
%     'fig18'   Reproduce Fig. 18:
%               Predicted response probabilities (PMVs) as a function of 
%               the amplitude panning ratio between a rear and a front 
%               loudspeaker in the same SP (lateral angle: 30 deg.). 
%               Left: Results for NH12. 
%               Center: Results for NH15. 
%               Right: Results for our pool of listeners.
%               Circle: Maximum of a PMV. 
%               Panning ratio of 0: Only front loudspeaker active. 
%               Panning ratio of 1: Only rear loudspeaker active. 
%
%     'fig19'   Reproduce Fig. 19:
%               Prediction matrices (PMVs) for different loudspeaker spans
%               and NH12. 
%               Left: Span of 0 deg. (single-loudspeaker reproduction, baseline). 
%               Center: Span of 30 deg. 
%               Right: Span of 60 deg. 
%               P: Predicted performance.
%
%     'fig20'   Reproduce Fig. 20:
%               Increase in localization errors as a function of the
%               loudspeaker span. 
%               Circles: Averages over all listeners from the pool. 
%               Error bars: +/- 1 standard deviation of the averages.
%
%     'fig22'   Reproduce Fig. 22:
%               Predictions for the surround setups in the VBAP configuration. 
%               Left: 5.1 setup, panning between the loudspeakers L and LS. 
%               Center: 10.2 setup (DSX), panning from L (polar angle of 0 deg.) 
%               via LH2 (55 deg.) to LS (180 deg.)
%               Right: 9.1 setup (Auro-3D), panning from L (0 deg.) via 
%               LH1 (34 deg.) and LSH (121 deg.) to LS (180 deg.). 
%               Desired polar angle: Continuous scale representing the VBAP 
%               across pair-wise contributing loudspeakers. 
%               All other conventions as in Fig. 18.
%
%     'fig23'   Reproduce Fig. 23:
%               Predictions for two modifications to the 9.1 setup (Auro 3D). 
%               Left: Original setup, loudspeakers LS and LSH 
%               at the azimuth of 110 deg. 
%               Center: LSH at the azimuth of 130 deg. 
%               Right: LS and LSH at the azimuth of 130 deg. 
%               All other conventions as in Fig. 22.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v0.4.3 or higher from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2013
%
%   See also: baumgartner2013, data_baumgartner2013
%
%   Examples:
%   ---------
%
%   To display Fig. 12 use :
%
%     exp_baumgartner2013('fig12'); 
%
%   To display Fig. 13 use :
%
%     exp_baumgartner2013('fig13');
%
%   To display Fig. 14 use :
%
%     exp_baumgartner2013('fig14');
%
%   To display Fig. 15 use :
%
%     exp_baumgartner2013('fig15');
%
%   To display Fig. 18 use :
%
%     exp_baumgartner2013('fig18');
%
%   To display Fig. 19 use :
%
%     exp_baumgartner2013('fig19');
%
%   To display Fig. 20 use :
%
%     exp_baumgartner2013('fig20');
%
%   To display Fig. 22 use :
%
%     exp_baumgartner2013('fig22');
%
%   To display Fig. 23 use :
%
%     exp_baumgartner2013('fig23');
%
%   References:
%     R. Baumgartner. Modeling sagittal-plane sound localization with the
%     application to subband-encoded head related transfer functions.
%     Master's thesis, University of Music and Performing Arts, Graz, June
%     2012.
%     
%     R. Baumgartner, P. Majdak, and B. Laback. Assessment of Sagittal-Plane
%     Sound Localization Performance in Spatial-Audio Applications,
%     chapter 4, pages 93-119. Springer-Verlag GmbH, 2013.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/experiments/exp_baumgartner2013.php

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

  definput.flags.type = {'missingflag',...
    'fig12','fig13','fig14','fig15',... % Ch. 4.1 (binaural rec.)
    'fig18','fig19','fig20','fig22','fig23'}; % Ch. 4
  definput.flags.plot = {'plot','noplot'};
  definput.import={'amt_cache'};

  % Parse input options
  [flags,kv]  = ltfatarghelper({},definput,varargin);
        
if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%% Figure 12

if flags.do_fig12
  
  s = data_baumgartner2013('pool', flags.cachemode);
  
  s = s([2 1 14]);  % NH12, NH58, NH33
  
  nfft = 2^12;
  range = 50;   % dB
  lat = 0;     % median sagittal plane
  fhigh = 18e3; % Hz
  da = 3;       % white frame for ticks (angles)
  df = 0.1;     % white frame for ticks (freq)
  
  figure
  for ll = 1:length(s)
  
    f = 0:s(1).fs/nfft:s(ll).fs/2;
    f = f/1e3; % in kHz
    
    s(ll).spdtfs = [];
    s(ll).polangs = [];
    [s(ll).spdtfs,s(ll).polangs] = extractsp(lat,s(ll).Obj);
    mag = 20*log10(abs(fft(s(ll).spdtfs,nfft)));
    mag = mag(1:nfft/2+1,:,:);
    maxmag = max(max(max(mag)));
    mag = mag - maxmag;       % normalize to 0dB
    maxmag = 0;
    mag(mag < maxmag-range) = maxmag-range;  % limit dynamic range
    crange = [floor(maxmag-range), ceil(maxmag)];% range of color code

    idf = f<fhigh/1e3;
    magl = mag(idf,:,1);
    fp = f(idf);
    polp = s(ll).polangs;

    subplot(1,3,ll)
    h = pcolor(fp,polp,magl');
    caxis(crange)
    shading interp
    colormap('bone')
    xlabel('Frequency (kHz)')
    ylabel('Polar Angle (\circ)')
    set(gca,...
        'XLim',[fp(1)-df,fp(end)+df],...
        'YLim',[polp(1)-da,polp(end)+da],...    
        'XTick',1:2:17,...
        'YTick',-30:30:polp(end),...
        'XMinorTick','on',...
        'YMinorTick','on'...
        )
    
  end
  
  if nargout == 1
    varargout{1} = s;
  end
  
end


%% Figures 13-15 (non-individual binaural recordings)

if flags.do_fig13 || flags.do_fig14 || flags.do_fig15
  
  latdivision = [-20,0,20];  % lateral center angles of SPs

  [s,qe,pe]=amt_cache('get','fig13to15',flags.cachemode);
  if isempty(s)

    s = data_baumgartner2013('pool', flags.cachemode);
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

    %fprintf('\n Please wait a little! \n');
    qe = zeros(ns,ns,length(latdivision)); % init QEs
    pe = qe;                               % init PEs
    for ll = 1:ns    % listener
        for jj = 1:ns    % ears
            for ii = 1:length(latdivision) % SPs

              s(ll).pmv{jj,ii} = [];
              s(ll).respangs{ii} = [];
              [s(ll).pmv{jj,ii},s(ll).respangs{ii}] = baumgartner2013(...
                  s(jj).spdtfs{ii},s(ll).spdtfs{ii},s(ll).fs,...
                  'u',s(ll).u,'lat',s(ll).latang{ii},...
                  'polsamp',s(ll).polangs{ii});

              [ qe(ll,jj,ii),pe(ll,jj,ii) ] = baumgartner2013_pmv2ppp( ...
                  s(ll).pmv{jj,ii} , s(jj).polangs{ii} , s(ll).respangs{ii});

            end
        end
        amt_disp([num2str(ll,'%2u') ' of ' num2str(ns,'%2u') ' completed'],'progress');
    end

    qe = mean(qe,3);
    pe = mean(pe,3);
    for ll = 1:length(s)
      s(ll).pe = pe(ll,:);
      s(ll).qe = qe(ll,:);
    end 
    amt_cache('set','fig13to15',s,qe,pe);
  end  
  ns = length(s);
  
  if nargout == 1
    varargout{1} = s; % SAME for Fig. 13-15!!!!
  end
  
  if flags.do_fig13 % -----------------------------------------------------
  
    if flags.do_plot

      markersize = 4;

      [qechance,pechance] = baumgartner2013_pmv2ppp(ones(49,44));
      
      % IDs for XTick of plots
      for ll = 1:ns
          nhs(ll,:) = s(ll).id;
      end

      % Polar error

      figure
      subplot(121)
      peinc = zeros(ns,ns-1);
      for l = 1:ns
          peinc(l,:) = pe(l,[1:l-1,l+1:ns]);
      end

      mpeinc = mean(peinc,2);
      speinc = std(peinc,0,2);

      h = bar(diag(pe));
      set(h,'FaceColor','white')

      hold on

      errorbar(mpeinc,speinc,'ko',...
        'MarkerSize',markersize,'MarkerFaceColor','w');

      plot(0:ns+1,pechance*ones(ns+2,1),'k--')

      ylabel('Local Polar RMS Error (\circ)')
      xlabel({'NH'})

      set(gca,'YLim',[20.1 48],'XLim',[0 ns+1],...
        'XTick',1:ns,'XTickLabel',nhs(:,3:4),...
          'YMinorTick','on')

      % Quadrant error

      subplot(122)
      qeinc = zeros(ns,ns-1);
      for l = 1:ns
          qeinc(l,:) = qe(l,[1:l-1,l+1:ns]);
      end

      mqeinc = mean(qeinc,2);
      sqeinc = std(qeinc,0,2);

      h = bar(diag(qe));
      set(h,'FaceColor','white')

      hold on
      errorbar(mqeinc,sqeinc,'ko',...
        'MarkerSize',markersize,'MarkerFaceColor','w');

      plot(0:ns+1,qechance*ones(ns+2,1),'k--')

      set(gca,'YLim',[floor(min(diag(qe))-1) ceil(qechance+1)],...
        'XLim',[0 ns+1],'XTick',1:ns,'XTickLabel',nhs(:,3:4),...
        'YMinorTick','on')

      ylabel('Quadrant Error (%)')
      xlabel({'NH'})

    end
    
  elseif flags.do_fig14 % -------------------------------------------------
    
    ears = 1;  % NH58
    s(ears).peexp = 23.5;
    s(ears).qeexp = 0.775;
    
    if flags.do_plot
      
      [qechance,pechance] = baumgartner2013_pmv2ppp(ones(49,44));
      
      % IDs for XTick of plots
      for ll = 1:ns
          nhs(ll,:) = s(ll).id;
      end
      
      figure
      
      % Polar error
      subplot(121)
      h = bar(pe(:,ears)-diag(pe));
      hold on
      plot(ears,s(ears).peexp-pe(ears,ears),'k*')

      dx = 0.42;
      peincchance = pechance-diag(pe);
      for ll = 1:length(s)
          plot(ll+[-dx,+dx],ones(2,1)*peincchance(ll),'k--')
      end

      ylabel('Increase in Local Polar RMS Error (\circ)')
      xlabel('NH')
      set(h,'FaceColor','white')
      set(gca,'XLim',[0 length(s)+1],...
        'XTick',1:length(s),'XTickLabel',nhs(:,3:4),...
        'YLim',[-1 21.9],'YMinorTick','on')

      % Quadrant error
      subplot(122)
      h = bar(qe(:,ears)-diag(qe));
      hold on
      plot(ears,s(ears).qeexp-qe(ears,ears),'k*')

      dx = 0.4;
      qeincchance = qechance-diag(qe);
      for ll = 1:length(s)
          plot(ll+[-dx,+dx],ones(2,1)*qeincchance(ll),'k--')
      end

      set(h,'FaceColor','white')
      set(gca,'XLim',[0 length(s)+1],...
        'XTick',1:length(s),'XTickLabel',nhs(:,3:4),...
        'YLim',[-5 18],'YMinorTick','on')
      xlabel('NH')
      ylabel('Increase in Quadrant Error (%)')
      
    end
    
  elseif flags.do_fig15
    
    pepoor = pe(:,14); % NH33 (14)
    pemed = pe(:,1);   % NH58 (1)
    pegood = pe(:,2);  % NH12 (2)

    qepoor = qe(:,14);
    qemed = qe(:,1);
    qegood = qe(:,2);
    
    if flags.do_plot
      
      markersize = 4;

      [qechance,pechance] = baumgartner2013_pmv2ppp(ones(49,44));
      
      % IDs for XTick of plots
      for ll = 1:ns
          nhs(ll,:) = s(ll).id;
      end
      
      figure
      
      % Polar error
      subplot(121)
      h = bar(diag(pe));
      hold on
      plot(pegood,'ko','MarkerSize',markersize,'MarkerFaceColor','k')
      plot(pemed, 'ks','MarkerSize',markersize,'MarkerFaceColor',0.5*ones(3,1))
      plot(pepoor,'kv','MarkerSize',markersize,'MarkerFaceColor','w')
      plot(0:length(s)+1,pechance*ones(length(s)+2,1),'k--')

      set(gca,'XLim',[0 length(s)+1],'XTick',1:length(s),'XTickLabel',nhs(:,3:4),...
              'YLim',[20.1 48],'YMinorTick','on')
      set(h,'FaceColor','white')%,'BarWidth',0.5)

      ylabel('Local Polar RMS Error (\circ)')
      xlabel('NH')

      % Quadrant error
      subplot(122)
      h = bar(diag(qe));
      hold on
      plot(qegood,'ko','MarkerSize',markersize,'MarkerFaceColor','k')
      plot(qemed,'ks','MarkerSize',markersize,'MarkerFaceColor',0.5*ones(3,1))
      plot(qepoor,'kv','MarkerSize',markersize,'MarkerFaceColor','w')
      plot(0:length(s)+1,qechance*ones(length(s)+2,1),'k--')

      set(gca,'XLim',[0 length(s)+1],'XTick',1:length(s),'XTickLabel',nhs(:,3:4),...
          'YLim',[1 44],'YMinorTick','on')
      set(h,'FaceColor','white')%,'BarWidth',0.5)

      ylabel('Quadrant Error (%)')
      xlabel('NH')
      
    end
    
  end
    
end


%% Figures 18, 22, 23 (Amplitude panning between loudspeakers)

if flags.do_fig18 || flags.do_fig22 || flags.do_fig23

  % Coordinates (azi,ele) of loudspeakers
  if flags.do_fig18
    LSPsetup{1} = [ 30,0 ; 150,0 ];                  % lat: 30 deg.
    package='fig18';
  elseif flags.do_fig22
    LSPsetup{1} = [ 30,0 ; 110,0 ];                  % 5.1
    LSPsetup{2} = [ 30,0 ; 45,45 ; 110,0 ];          % 10.2
    LSPsetup{3} = [ 30,0 ; 30,30 ; 110,30 ; 110,0 ]; % 9.1
    package='fig22';
  elseif flags.do_fig23
    LSPsetup{1} = [ 30,0 ; 30,30 ; 110,30 ; 110,0 ]; % 9.1
    LSPsetup{2} = [ 30,0 ; 30,30 ; 140,30 ; 110,0 ]; % 9.1 (prop. LSH)
    LSPsetup{3} = [ 30,0 ; 30,30 ; 140,30 ; 140,0 ]; % 9.1 (prop. LSH&LS)
    package='fig23';
  end
  
  dpol = 5;    % step size in deg

  [s, D]=amt_cache('get',package,flags.cachemode);
  if isempty(s),
    s = data_baumgartner2013('pool', flags.cachemode);
    ns = length(s);

    for ss = 1:length(LSPsetup)

      LSPsph = LSPsetup{ss};
      nLSP = size(LSPsph,1);
      LSPhp = zeros(nLSP,2);
      [LSPhp(:,1),LSPhp(:,2)] = sph2horpolar(LSPsph(:,1),LSPsph(:,2));
      LSPcart = zeros(nLSP,3);
      [LSPcart(:,1),LSPcart(:,2),LSPcart(:,3)] = ...
            sph2cart(LSPsph(:,1),LSPsph(:,2),ones(nLSP,1));

      D = LSPhp(1,2):dpol:LSPhp(nLSP,2);     % Desired polar angle    

      %fprintf('\n Please wait a little! \n');
      for ll = 1:ns

        s(ll).pos(:,1)=bsxfun(@times,s(ll).Obj.SourcePosition(:,1),ones(s(ll).Obj.API.M,1));
        s(ll).pos(:,2)=bsxfun(@times,s(ll).Obj.SourcePosition(:,2),ones(s(ll).Obj.API.M,1));
        [poscart(:,1),poscart(:,2),poscart(:,3)] = ...
            sph2cart(s(ll).pos(:,1),s(ll).pos(:,2),ones(size(s(ll).pos,1),1));

        % DTF indices of LSPs
        idLSP = zeros(nLSP,1);
        for ii = 1:nLSP
          e = poscart-repmat(LSPcart(ii,:),size(poscart,1),1);
          [tmp,idLSP(ii)] = min( sqrt(sum(e.^2,2)) ); % minimum euclidean distance
        end
        clear poscart

        Lp = 1;  % LSP pair index
        for d = 1:length(D)

          if D(d) > LSPhp(Lp+1,2)  % Switch to next LSP pair?
            Lp = Lp + 1;
          end
          pol1 = LSPhp(Lp,2);
          pol2 = LSPhp(Lp+1,2);

          if abs(pol1-pol2) < 180

            pol0 = mean([pol1 pol2]);   % polar angle of center
            phi0 = abs(pol1-pol2)/2;    % basis angle
            M = (1-tan(deg2rad(pol0-D(d)))/tan(deg2rad(phi0)))/2;

          else

            M = (D(d)-pol1) / (pol2-pol1);

          end

          lat1 = LSPhp(Lp,1);
          lat2 = LSPhp(Lp+1,1);
          lat0 = mean([lat1 lat2]);   % lateral angle of center
          phi0 = abs(lat1-lat2)/2;    % basis angle
          phiM = rad2deg( atan( tan(deg2rad(phi0)) * (1-2*M) ) );
          latM = lat0 - phiM;

          s(ll).dtfs = permute(double(s(ll).Obj.Data.IR),[3 1 2]);

          dtfLSP = (1-M)*s(ll).dtfs(:,idLSP(Lp),:) + M*s(ll).dtfs(:,idLSP(Lp+1),:);
          [spdtfs,polangs] = extractsp(latM,s(ll).Obj);

          [pmv,respangs] = baumgartner2013(dtfLSP,spdtfs,s(ll).fs,...
              'lat',latM,'u',s(ll).u,...
              'polsamp',polangs);
          idP = respangs >= -30 & respangs <= 210;
          s(ll).pmv{ss}(:,d) = pmv(idP)/sum(pmv(idP));
          s(ll).respangs = respangs(idP);

        end

        amt_disp([num2str(ll,'%2u') ' of ' num2str(ns,'%2u') ' completed'],'progress')

      end

      % Pool
      ppool = zeros(size(s(1).pmv{ss}));
      for ll = 1:ns
          ppool = ppool + s(ll).pmv{ss};
      end
      s(ns+1).pmv{ss} = ppool/length(s);
      s(ns+1).respangs = s(1).respangs;
      s(ns+1).id = 'Pool';

    end
    amt_cache('set',package,s,D);
  end
  % Output
  if nargout == 1
    varargout{1} = s;
  end
  
  % Plots
  if flags.do_plot
    
    if flags.do_fig18
      subject = [9 5 18];
      setup = [1 1 1];
    elseif flags.do_fig22
      subject = [18 18 18];
      setup = 1:3;
    elseif flags.do_fig23
      subject = [18 18 18];
      setup = 1:3;
    end
    
    figure
    
    for ii = 1:3
      ll = subject(ii);
      ss = setup(ii);
      subplot(1,3,ii)
    
      h = pcolor(D,s(ll).respangs,s(ll).pmv{ss});
      if not(isoctave) % does not work with x11 (mac osx)
        axis equal
      end
      if flags.do_fig23 || (flags.do_fig22 && ss > 1)
        set(gca,'XLim',[D(1)-3 D(end)+3],'XMinorTick','on',...
          'XTick',0:30:180,'XTickLabel',{0:30:150,180})
        xlabel('Desired Polar Angle (\circ)')
      else
        set(gca,'XLim',[D(1)-3 D(end)+3],'XMinorTick','on',...
          'XTick',18:36:180,'XTickLabel',{0.1:0.2:0.9})
        xlabel('Panning Ratio')
      end
      set(gca,'YLim',[s(ll).respangs(1)-3 s(ll).respangs(end)+3],...
      	'YMinorTick','on','YTick',(0:30:180)+2.5,'YTickLabel',0:30:180)
      colormap bone
      shading flat
      caxis([0 0.1])
      ylabel('Response Angle (\circ)')

      [tmp,Imax] = max(s(ll).pmv{ss});
      hold on
      h1 = plot( D,s(ll).respangs(Imax)+2.5, 'wo');  % shadow
      set(h1,'MarkerSize',4,'MarkerFaceColor','none') 
      h2 = plot( D,s(ll).respangs(Imax)+2.5, 'ko'); 
      set(h2,'MarkerSize',3,'MarkerFaceColor','none') 
      
    end
    
  end
  
end


%% Figures 19, 20 (Effect of loudspeaker span)

if flags.do_fig19 || flags.do_fig20
  
  s = data_baumgartner2013('pool', flags.cachemode);
  
  if flags.do_fig19
    dPol = [0 30,60];
    s = s(2); % NH12
  elseif flags.do_fig20
    dPol = 10:10:90; 
  end
  lat = 0;
  
  s(1).spdtfs = [];
  [s(1).spdtfs,polang] = extractsp(lat,s(1).Obj); 
  qeI = zeros(length(s),length(dPol));
  peI = qeI;
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
    
    amt_disp([' Span: ' num2str(dPol(ii)) ' deg.'],'progress');
    for ll = 1:length(s)

        s(ll).spdtfs = extractsp(lat,s(ll).Obj);

        % superposition
        s(ll).dtfs2{ii} = s(ll).spdtfs(:,id1,:) + s(ll).spdtfs(:,id2,:);
        
        [s(ll).pmv1{ii},respang] = baumgartner2013(...
          s(ll).spdtfs(:,id0,:),s(ll).spdtfs,s(ll).fs,'u',s(ll).u,...
          'polsamp',polang);
        s(ll).pmv2{ii} = baumgartner2013(...
          s(ll).dtfs2{ii},s(ll).spdtfs,s(ll).fs,'u',s(ll).u,...
          'polsamp',polang);

        [s(ll).qe1{ii},s(ll).pe1{ii}] = baumgartner2013_pmv2ppp(...
          s(ll).pmv1{ii},polang(id0),respang);
        [s(ll).qe2{ii},s(ll).pe2{ii}] = baumgartner2013_pmv2ppp(...
          s(ll).pmv2{ii},pol2{ii},respang);
        
        % Increse of error
        qeI(ll,ii) = s(ll).qe2{ii} - s(ll).qe1{ii};
        peI(ll,ii) = s(ll).pe2{ii} - s(ll).pe1{ii};

    end

  end
  
  if nargout == 1
    varargout{1} = s;
  end

  if flags.do_plot
    
    if flags.do_fig19
      
      figure;
      subplot(1,3,1)
      plot_baumgartner2013(s(1).pmv1{1},polang,respang,'cmax',0.1,'nocolorbar');
      title(['P: PE = ' num2str(s(1).pe1{1},2) '\circ, QE = ' num2str(s(1).qe1{1},2) '%'])
      
      subplot(1,3,2)
      plot_baumgartner2013(s(1).pmv2{2},pol2{2},respang,'cmax',0.1,'nocolorbar');
      title(['P: PE = ' num2str(s(1).pe2{2},2) '\circ, QE = ' num2str(s(1).qe2{2},2) '%'])
            
      subplot(1,3,3)
      plot_baumgartner2013(s(1).pmv2{3},pol2{3},respang,'cmax',0.1,'nocolorbar');
      title(['P: PE = ' num2str(s(1).pe2{3},2) '\circ, QE = ' num2str(s(1).qe2{3},2) '%'])
      
    elseif flags.do_fig20
      
      qeIm = mean(qeI);
      peIm = mean(peI);
      qeIs = std(qeI);
      peIs = std(peI);

      MarkerSize = 5;

      figure
      subplot(121)
      errorbar(dPol,peIm,peIs,'ko-',...
          'MarkerSize',MarkerSize,...
          'MarkerFaceColor','k');
      ylabel('Increase in Local Polar RMS Error (\circ)')
      xlabel('Loudspeaker Span (\circ)')
      set(gca,...
          'XLim',[dPol(1)-10 dPol(end)+10],...
          'XTick',dPol,...
          'YLim',[-1,12.9],...
          'YMinorTick','on')
      axis square

      subplot(122)
      errorbar(dPol,qeIm,qeIs,'ko-',...
          'MarkerSize',MarkerSize,...
          'MarkerFaceColor','k');
      ylabel('Increase in Quadrant Error (%)')
      xlabel('Loudspeaker Span (\circ)')
      set(gca,...
          'XLim',[dPol(1)-10 dPol(end)+10],...
          'XTick',dPol,...
          'YLim',[-1,12.9],...
          'YMinorTick','on')
      axis square
      
    end
      
  end
  
end


end
