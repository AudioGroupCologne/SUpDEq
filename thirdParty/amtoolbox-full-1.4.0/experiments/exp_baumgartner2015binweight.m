function varargout = exp_baumgartner2015binweight(varargin)
%EXP_BAUMGARTNER2015BINWEIGHT Results from various publications of Baumgartner et al. (2015)
%   Usage: data = exp_baumgartner2015(flag) 
%
%   exp_baumgartner2015(flag) reproduces figures of the studies from 
%   Baumgartner et al. (2015).
%
%
%   The following flags can be specified
%
%
%     'fig5'    Reproduce Fig.5 of Baumgartner et al. (2015):
%                                  Effect of background noise on reliability of contralateral  
%                                  cues for various lateral eccentricities. Top row:
%                                  Across-listener averages of performance measures for
%                                  contralateral ear. Bottom row: Contralateral re ipsilateral
%                                  averages of performance measures.
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
%   To display Fig.5 of Baumgartner et al. (2015) use :
%
%     exp_baumgartner2015binweight('fig5');
%
%   See also: baumgartner2014 data_baumgartner2014
%
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. The reliability of
%     contralateral spectral cues for sound localization in sagittal planes.
%     In Midwinter Meeting of the Association for Research in Otolaryngology,
%     Baltimore, MD, Feb 2015.
%     
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_baumgartner2015binweight.php


% #Author: Robert Baumgartner
% #Author: Clara Hollomey (2021): adapted for AMT 1.1
% #Author: Piotr Majdak (2021): adapted for AMT 1.1

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
definput.flags.type = {'missingflag', 'fig5'};

definput.flags.plot = {'plot','no_plot'};


[flags,kv]  = ltfatarghelper({'FontSize','MarkerSize'},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%% General Plot Settings
TickLength = [0.02,0.04];
% ltfatsetdefaults('dbspl','dboffset',100); % PM: not required
%% ------ FIG 5 of baumgartner2015aro -------------------------------------
if flags.do_fig5
  
  [perr,qerr,snrFront,bwcoef,lat] = amt_cache('get','fig5_baumgartner2015aro',flags.cachemode);
  if isempty(perr)
    
    snrFront = -20:2:40; % in dB

    latecc = [10,30,50];% lateral eccentricities

    mrs = 17; % no sensorimotor mapping

    s = data_baumgartner2014('pool','cached'); % PM: Always cached. For recalculation call data_* directly.
    
    bwcoef = [13,+eps,-eps]; % configuration of binaural weighting stage (binaural, ipsilateral, contralateral

    lat = [-fliplr(latecc),latecc];

    maskerNoise = noise(0.05*s(1).Obj.Data.SamplingRate,1,'white');
    targetNoise = noise(0.05*s(1).Obj.Data.SamplingRate-255,1,'white');

    perr = nan(length(snrFront),length(bwcoef),length(lat),length(s));
    qerr = nan(length(snrFront),length(bwcoef),length(lat),length(s));
    for isub=1:length(s)
      idfrontal = find(s(isub).Obj.SourcePosition(:,1)==0 & s(isub).Obj.SourcePosition(:,2)==0);
      frontalDtfs = shiftdim(s(isub).Obj.Data.IR(idfrontal,:,:),2);
      frontalTarget = convolve(targetNoise,frontalDtfs);
      lvl = mean(dbspl(frontalTarget)); % level of frontal target stimulus in dB
      for ilat = 1:length(lat)
        [spdtfs,tang] = extractsp(lat(ilat),s(isub).Obj);
        targets = convolve(targetNoise,spdtfs);
        targets = reshape(targets,[length(targets),size(targets,2)/2,2]);
        for isnr = 1:length(snrFront)
          targetsPlusMasker = targets + ...
            repmat(scaletodbspl(maskerNoise,lvl-snrFront(isnr)),[1,size(targets,2),2]);
          for ibwc = 1:length(bwcoef)
            [p,rang] = baumgartner2014(targetsPlusMasker,s(isub).Obj,...
              'S',s(isub).S,'mrsmsp',mrs,...
              'lat',lat(ilat),'bwcoef',bwcoef(ibwc));
            [ qerr(isnr,ibwc,ilat,isub) , perr(isnr,ibwc,ilat,isub) ] = ...
              baumgartner2014_pmv2ppp(p,tang,rang);
          end
        end
      end
      amt_disp([num2str(isub) ' of ' num2str(length(s)) ' completed']);
    end

    amt_cache('set','fig5_baumgartner2015aro',perr,qerr,snrFront,bwcoef,lat)
    
  end
  
  r = struct('perr',perr,'qerr',qerr,'snrFront',snrFront,'bwcoef',bwcoef,'lat',lat);
  varargout{1} = r;
  
  if flags.do_plot
    
    % pool left/right
    perr = (r.perr(:,:,length(r.lat)/2:-1:1,:) + r.perr(:,:,1+length(r.lat)/2:length(r.lat),:))/2;
    qerr = (r.qerr(:,:,length(r.lat)/2:-1:1,:) + r.qerr(:,:,1+length(r.lat)/2:length(r.lat),:))/2;
    latecc = r.lat(1+length(r.lat)/2:length(r.lat)); % lateral eccentricity

    perr_ipsipro = squeeze(perr(:,3,:,:) - perr(:,2,:,:)); % contra minus ipsi
    qerr_ipsipro = squeeze(qerr(:,3,:,:) - qerr(:,2,:,:)); 

    snr_int = r.snrFront(1):r.snrFront(end);

    figure
    
    % display 0-error line
    for ii = 1:2
      subplot(2,2,2+ii)
      plot(snr_int,zeros(length(snr_int),1),'k:')
      hold on
    end
    
    color = [ 0.2081    0.1663    0.8292;...
              0.8292    0.1663    0.2081;...
              0.8081    0.6081    0.2081];
    for ii=1:length(latecc)
  
      subplot(2,2,1)
      abspecontra_int = interp1(r.snrFront,mean(perr(:,ii,3,:),4),snr_int,'spline');
      h(ii) = plot(snr_int,abspecontra_int,'Color',color(ii,:));  hold on
      xlabel('SNR (dB)','FontSize',kv.FontSize)
      ylabel('PE_{contra} (deg)','FontSize',kv.FontSize)
      axis([-20,40,31,54])
      set(gca,'FontSize',kv.FontSize)

      subplot(2,2,2)
      absqecontra_int = interp1(r.snrFront,mean(qerr(:,ii,3,:),4),snr_int,'spline');
      h(ii) = plot(snr_int,absqecontra_int,'Color',color(ii,:));  hold on
      xlabel('SNR (dB)','FontSize',kv.FontSize)
      ylabel('QE_{contra} (deg)','FontSize',kv.FontSize)
      axis([-20,40,6,49])
      set(gca,'FontSize',kv.FontSize)

      subplot(2,2,3)
      pe_int = interp1(r.snrFront,mean(perr_ipsipro(:,ii,:),3),snr_int,'spline');
      h(ii) = plot(snr_int,pe_int,'Color',color(ii,:));
      xlabel('SNR (dB)','FontSize',kv.FontSize)
      ylabel('PE_{contra} - PE_{ipsi} (deg)','FontSize',kv.FontSize)
      axis([-20,40,-4,29])
      set(gca,'FontSize',kv.FontSize)

      subplot(2,2,4)
      qe_int = interp1(r.snrFront,mean(qerr_ipsipro(:,ii,:),3),snr_int,'spline');
      plot(snr_int,qe_int,'Color',color(ii,:))
      xlabel('SNR (dB)','FontSize',kv.FontSize)
      ylabel('QE_{contra} - QE_{ipsi} (deg)','FontSize',kv.FontSize)
      axis([-20,40,-4,29])
      set(gca,'FontSize',kv.FontSize)
    end

    subplot(2,2,3)
    legendentries = [repmat('\phi = \pm',length(latecc),1) num2str(latecc(:)) repmat('\circ',length(latecc),1)];
    leg = legend(h,legendentries,'Location','north');
    set(leg,'FontSize',kv.FontSize)
    
  end
% ltfatsetdefaults('dbspl','dboffset',93.98);   % PM: not required.
end


