function exp_enzner2008(varargin) 
%EXP_ENZNER2008  Creates figures like [Enzner2008, Fig. 2], [Enzner2009, Fig. 4]
%   Usage: exp_enzner2008(flag)
%
%   Required data: hrtf/enzner2008
% 
%   The following flags can be specified:
%
%     'fig2'   Plot Fig. 2 from Enzner et al. (2008)
%
%     'fig4_enzner2009'   Plot Fig. 4 from Enzner et al. (2009)
%
%   Examples:
%   ---------
%
%   To display Figure 2 from the 2008 paper use :
%
%     exp_enzner2008('fig2');
%
%   To display Figure 4 from the 2009 paper use :
%
%     exp_enzner2008('fig4_enzner2009');
%
%   See also: enzner2008
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_enzner2008.php


%   #Author: Michael Weinert (2013)
%   #Author: Gerald Enzner (2013)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%addpath(fullfile(amt_basepath,'hrtf','continuous-azimuth HRIR'))

% enzner2008(1,1,varargin); % enzner2008(mu, delta_phi, varargin)

definput.flags.type = {'missingflag', 'fig2','fig4_enzner2009'};
[flags,kv]  = ltfatarghelper({},definput,varargin);
        
if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


rec_filename = 'example_1ch_white_noise_earsignals.wav';
ref_filename = 'example_1ch_white_noise_reference.wav';
P.adapt = 20000; % depends on the recording (overhead at the end and the beginnig)
P.sys_latency = 30;
P.mu=1;
P.delta_phi=1;
P.h_length = 256;

%% read signals
x = amt_load('enzner2008',ref_filename);
[y, fs] = amt_load('enzner2008',rec_filename);

%% Fig. 2 from Enzner (2008)
if flags.do_fig2
  
  hrir_data = enzner2008(x,y,P);

  h0=squeeze(hrir_data(:,1,271));
  h0(:,2)=squeeze(hrir_data(:,2,271));
  
  figure
  subplot(2,1,1)
  plot(20*log10(abs(h0(:,1)./max(max(abs(h0))))))
  xlim([1 length(h0(:,1))])
  ylim([-80 3])
  ylabel('|h_1(\kappa,\theta_k)| [dB]')
  title(['azimuth \theta_k = 270 deg, left ear'])
  grid on
  subplot(2,1,2)
  plot(20*log10(abs(h0(:,2)./max(max(abs(h0))))))
  xlim([1 length(h0(:,2))])
  ylim([-80 3])
  ylabel('|h_2(\kappa,\theta_k)| [dB]')
  xlabel('impulse response lag \kappa')
  title(['azimuth \theta_k = 270 deg, right ear'])
  grid on
  
end
%% Fig. 4 from Enzner (2009)
if flags.do_fig4_enzner2009

    [hrir_data, hrir_angles, errorsig] = enzner2008(x,y,P);
    
    SNR_l = 10*log10(var(errorsig(P.adapt+1:(length(y)-P.adapt),1))/var(y(P.adapt+1:(length(y)-P.adapt),1))/P.h_length);
    SNR_r = 10*log10(var(errorsig(P.adapt+1:(length(y)-P.adapt),2))/var(y(P.adapt+1:(length(y)-P.adapt),2))/P.h_length);
    t = linspace(0, length(y)/fs, length(y));
    figure
    subplot(2,1,1)
    txt1 = ['(\sigma_{e^l}^2/\sigma_{y^l}^2)/N = ',num2str(SNR_l),' dB'];
    title({txt1});
    hold on
    p1 = plot(t, y(:,1),'b');
    p2 = plot(t, errorsig(:,1),'m');
    xlim([(P.adapt+1)/fs (length(y)-P.adapt)/fs]);
    xlabel('time [s]');
    ylim([-0.5 0.5]);
    ylabel('amplitude');
    legend([p1 p2],{'left ear recording y^l(k)' 'error signal e^l(k)'}, 'location', 'northwest')
    grid on;
    subplot(2,1,2)
    txt2 = ['(\sigma_{e^r}^2/\sigma_{y^r}^2)/N = ',num2str(SNR_r),' dB'];
    title({txt2});
    hold on
    p1 = plot(t,y(:,2),'b');
    p2 = plot(t,errorsig(:,2),'m');
    xlim([(P.adapt+1)/fs (length(y)-P.adapt)/fs]);
    xlabel('time [s]');
    ylim([-0.5 0.5]);
    ylabel('amplitude');
    legend([p1 p2],{'right ear recording y^r(k)' 'error signal e^r(k)'}, 'location', 'northeast')
    grid on;
end;





