function data = exp_hassager2016(varargin)
%EXP_HASSAGER2016 Experiments of Hassager et al. (2016)
%   Usage: data = exp_hassager2016(flag) 
%
%   EXP_HASSAGER2016(flag) reproduces figures of the study from 
%   Hassager et al. (2016).
%
%   Optional fields of output data structure:
% 
%     'contralateralGain'   contralateral gain of binaural weighting function
%
%
%   The following flags can be specified
%
%     'fig6'    Reproduce Fig.6:
%               The mean of the seven listeners' perceived sound source 
%               location (open symbols) as a function of the bandwidth factor 
%               and the corresponding model predictions (filled symbols). 
%               The model predictions have been shifted slightly to the right 
%               for a better visual interpretation. The error bars are one 
%               standard error of the mean.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v0.4.3 or higher from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/...
%
%   3) Statistics Toolbox for Matlab (for some of the figures)
%
%   Examples:
%   ---------
%
%   To display Fig.6 use :
%
%     exp_hassager2016('fig6');
%
%   References:
%     H. G. Hassager, F. Gran, and T. Dau. The role of spectral detail in the
%     binaural transfer function on perceived externalization in a
%     reverberant environment. J. Acoust. Soc. Am., 139(5):2992--3000, 2016.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_hassager2016.php


%   #Author: Robert Baumgartner (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.import={'amt_cache'};
definput.flags.type = {'fig6'};
definput.flags.plot = {'plot','no_plot'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

% if flags.do_missingflag
%   flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
%              sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
%   error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
% end

if flags.do_fig6
  azi = [0,50];
  
  Pext_A = data_hassager2016;
  B = Pext_A.B;
  
  data = data_baumgartner2017;
  data = data(1:5);
  fs = data(1).Obj.Data.SamplingRate;
  fncache = 'hassager2016';
  Pext = amt_cache('get',fncache,flags.cachemode);
  if isempty(Pext)
    Pext = nan(length(B),length(data),length(azi));
    for isubj = 1:length(data)
      Obj = data(isubj).Obj;
      for iazi = 1:length(azi)
    %     templateSound = SOFAspat(in,Obj,azi(iazi),0);
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
    % plotfftreal(iB*db2mag(10)*fftreal(target(:,1)),fs,'flog'); hold on
          Pext(iB,isubj,iazi) = hassager2016(target,template,'fs',fs);
        end
      end
      amt_disp([num2str(isubj),' of ',num2str(length(data)),' subjects completed.']);
    end
    amt_cache('set',fncache,Pext);
  end
  
  %% Output
  for isub = 1:length(data)
    data(isub).Externalization = squeeze(Pext(:,isub,:));
    data(isub).BandwidthFactor = B;
    data(isub).Azimuth = azi;
  end
  
  %%
  if flags.do_plot
    BplotTicks = logspace(log10(0.25),log10(64),9);
    BplotTicks = round(BplotTicks*100)/100;
    BplotStr = ['Ref.';num2str(BplotTicks(:))];
    BplotTicks = [BplotTicks(1)/2,BplotTicks];
    
    B(1) = B(2)/2; % just for illustration
    figure 
    symb = {'-ko','-k^'};
    for iazi = 1:length(azi)
      subplot(1,2,iazi)
      h(1) = plot(B,Pext_A.rating(:,iazi),symb{iazi});
      set(h(1),'MarkerFaceColor','w')
      hold on
      h(2) = errorbar(B*1.1,mean(Pext(:,:,iazi),2),std(Pext(:,:,iazi),0,2)/sqrt(length(data)),symb{iazi});
      set(h(2),'MarkerFaceColor','k')
      set(gca,'XTick',BplotTicks,'XTickLabel',BplotStr,'XScale','log')
      axis([BplotTicks(1)/1.5,BplotTicks(end)*1.5,0.8,5.2])
      xlabel('Bandwidth Factor [ERB]')
      ylabel('Mean Externalization Rating')
      title([num2str(azi(iazi)),'\circ'])
      grid on
      leg = legend({'dir - data','dir - pred'},'Location','southwest');
    end
  end
  
end

end


