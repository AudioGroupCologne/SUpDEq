function example_output = exp_dietz2011(varargin)
%EXP_DIETZ2011 Experiments from Dietz et al. 2011
%
%   EXP_DIETZ2011(fig) reproduce Fig. no. fig from the Dietz et
%   al. 2011 paper.
%
%   *Note**: The input signals used in this routine are not identical to
%   the ones used in the original paper.
%
%   The following flags can be specified;
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%     'fig3'    Reproduce Fig. 3 panels a + b.
%
%     'fig4'    Reproduce Fig. 4.
%
%     'fig5'    Reproduce Fig. 5.
%
%     'fig6'    Reproduce Fig. 6.
%
%     'fig6_spille2013'    Reproduce Fig. 6 from Spille et al. (2013).
%
%   Examples:
%   ---------
%
%   To display Fig. 3 use :
%
%     exp_dietz2011('fig3');
%
%   To display Fig. 4 use :
%
%     exp_dietz2011('fig4');
%
%   To display Fig. 5 use :
%
%     exp_dietz2011('fig5');
%
%   To display Fig. 6 use :
%
%     exp_dietz2011('fig6');
%
%   To display Fig. 6 from Spille et al. (2013) use :
%
%     exp_dietz2011('fig6_spille2013');
%
%   References:
%     M. Dietz, S. D. Ewert, and V. Hohmann. Auditory model based direction
%     estimation of concurrent speakers from binaural signals. Speech
%     Communication, 53(5):592-605, 2011. [1]http ]
%     
%     C. Spille, B. Meyer, M. Dietz, and V. Hohmann. Binaural scene analysis
%     with multi-dimensional statistical filters, chapter 6. Springer-Verlag
%     GmbH, 2013.
%     
%     References
%     
%     1. http://www.sciencedirect.com/science/article/pii/S016763931000097X
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/experiments/exp_dietz2011.php

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

% AUTHOR: Mathias Dietz

definput.flags.type = {'missingflag','fig3','fig4','fig5','fig6','fig6_spille2013'};
definput.flags.plot = {'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
               sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

% Load polynomial lookup data for converting ITD to azimuth
lookup = amt_load('dietz2011','itd2angle_lookuptable.mat');

if flags.do_fig3

    signal=repmat(sig_competingtalkers('five_speakers'),2,[]);
    fs = 16000;
    signal=signal(1:5*fs,:);
    s_pos = [-80 -30 0 30 80];
    ic_threshold=0.98;
    panellabel = 'ab';

    % run IPD model on signal
    [hairc_fine, fc, hairc_ild]=dietz2011(signal,fs,'fhigh',1400);

    % convert interaural information into azimuth
    itd_unwrapped = ...
        dietz2011_unwrapitd(hairc_fine.itd_lp,hairc_ild,hairc_fine.f_inst,2.5);
    angl=itd2angle(itd_unwrapped,lookup);

    h_ic=zeros(91,12);
    h_all=histc(angl,-90:2:90);
    for n=1:12
        h_ic(:,n)=histc(angl(hairc_fine.ic(:,n)>ic_threshold&[diff(hairc_fine.ic(:,n))>0; 0],n),-90:2:90);
    end
    example_output.angle_fine = angl;
    example_output.IVS_fine = hairc_fine.ic;
    example_output.histogram_angle_label = -90:2:90;
    example_output.histograms_with_IVS = h_ic;
    example_output.histograms_without_IVS = h_all;

    if flags.do_plot
        figure;
        fontsize = 14;
        set(gcf,'Position',[100 100 1170 700]);
        
        for panel = 1:2
            subplot(1,2,panel)
            switch panel
                
              case 1
                bar(-90:2:90,sum(h_all,2),'r');
                title('Mean histogram of all fine-structure channels','Fontsize',fontsize);
								axis([-90 90 0 21900]);
								set(gca,'YTick',[5000 10000 15000 20000],'YTickLabel',{'5k','10k','15k','20k'});
                ymax = max(sum(h_all,2));
              case 2
								bar(-90:2:90,sum(h_ic,2));
								title('Mean histogram with VS filter','Fontsize',fontsize);
								axis([-90 90 0 5600]);
								set(gca,'YTick',[1000:1000:5000],'YTickLabel',{'1k','2k','3k','4k','5k'});
								ymax = max(sum(h_ic,2));
            end
            set(gca,'Fontsize',fontsize);
            set(gca,'XTick',s_pos);
%             xlim([-93 93])
%             ylim([0 ymax*1.1])
            xlabel('Azimuth [deg]','Fontsize',fontsize);
            ylabel('Frequency of occurence','Fontsize',fontsize);
            text (-80,ymax*.95,panellabel(panel),'Fontsize',fontsize+1,'FontWeight','bold');
        end
    end;
end;


if flags.do_fig4
    
    % This reproduces Figure 4 from Dietz et al Speech Comm. 2011

    signal=repmat(sig_competingtalkers('two_speakers'),2,[]);
    fs = 16000;
    s_pos = [-80 -30 0 30 80];
    ic_threshold=0.98;
    cn = [10 1]; % channel numbers for separate plots (1st entry also for time plot)
    panellabel = 'abc';

    % run IPD model on signal
    [hairc_fine, fc, hairc_ild, hairc_mod216]=dietz2011(signal,fs,'mod_center_frequency_hz',216);
    [~, ~, ~, hairc_mod135]=dietz2011(signal,fs,'mod_center_frequency_hz',135);
    % convert interaural information into azimuth
    itd_unwrapped = ...
        dietz2011_unwrapitd(hairc_fine.itd_lp,hairc_ild(:,1:12),hairc_fine.f_inst,2.5);
    angl=itd2angle(itd_unwrapped,lookup);
    angl_fmod216=hairc_mod216.itd_lp*140000; %linear approximation. paper version is better than this
    angl_fmod135=hairc_mod135.itd_lp*140000; %linear approximation. paper version is better than this

    h_ic=zeros(61,12);
    h_fmod216=zeros(61,11);
    h_fmod135=zeros(61,11);
    h_all=histc(angl,-60:2:60);
    for n=1:12
        h_ic(:,n)=histc(angl(hairc_fine.ic(:,n)>ic_threshold&[diff(hairc_fine.ic(:,n))>0; 0],n),-60:2:60);
    end
    for n=1:11
        h_fmod216(:,n)=histc(angl_fmod216(hairc_mod216.ic(:,n)>ic_threshold&[diff(hairc_mod216.ic(:,n))>0; 0],n),-60:2:60);
        h_fmod135(:,n)=histc(angl_fmod135(hairc_mod135.ic(:,n)>ic_threshold&[diff(hairc_mod135.ic(:,n))>0; 0],n),-60:2:60);
    end
    example_output.angle_fine = angl;
    example_output.IVS_fine = hairc_fine.ic;
    example_output.histogram_angle_label = -60:2:60;
    example_output.histogram_panel1 = h_fmod216;
    example_output.histogram_panel2 = h_fmod135;
    example_output.histogram_panel3 = sum(h_ic,2);

    if flags.do_plot
        figure;
        fontsize = 14;
        set(gcf,'Position',[100 100 1170 400])
        
        for panel = 1:3
            subplot(1,3,panel)
            switch panel
                
              case 1
                bar(-60:2:60,sum(h_fmod216,2))
                title('histogram of mod ITD channels 13-23','Fontsize',fontsize)
								axis([-50 50 0 7600]);
								set(gca,'YTick',[2500 5000 7500],'YTickLabel',{'2.5k','5k','7.5k'});
                ymax = max(sum(h_fmod216,2));
%                 ylim([0 ymax*1.15])
              case 2
                bar(-60:2:60,sum(h_fmod135,2))
								axis([-50 50 0 7600]);
								set(gca,'YTick',[2500 5000 7500],'YTickLabel',{'2.5k','5k','7.5k'});
                title('histogram of mod ITD channels 13-23','Fontsize',fontsize)
                ymax = max(sum(h_fmod135,2));
%                 ylim([0 ymax*1.15])
              case 3
                bar(-60:2:60,sum(h_ic,2))
                title('histogram of fine ITD channels 1-12','Fontsize',fontsize)
								axis([-50 50 0 76000]);
								set(gca,'YTick',[25000 50000 75000],'YTickLabel',{'25k','50k','75k'});
                ymax = max(sum(h_ic,2));
%                 ylim([0 ymax*1.15])
            end
            set(gca,'Fontsize',fontsize)
            set(gca,'XTick',s_pos)
%             xlim([-63 63])
            
            xlabel('Azimuth [deg]','Fontsize',fontsize)
            ylabel('Frequency of occurence','Fontsize',fontsize)
            text (-40,ymax*1.2,panellabel(panel),'Fontsize',fontsize+1,'FontWeight','bold')
            
        end
    end;
end;

if flags.do_fig5
    % mix signals
    signal1=sig_competingtalkers('one_of_three');
    signal2=sig_competingtalkers('two_of_three');
    signal3=sig_competingtalkers('three_of_three');
    noise  =sig_competingtalkers('bnoise');
    noise = noise(1:40000,:);
    fs = 16000;

    % derive histograms
    ic_threshold=0.98;
    k = zeros(14,38);
    for n = 1:7
        switch n
          case 1
            signal = signal1+noise;
            
          case 2
            signal = signal2+noise;
          case 3
            signal = signal3+noise;
          case 4
            signal = signal1+2*noise;
          case 5
            signal = signal2+2*noise;
          case 6
            signal = signal3+2*noise;
          case 7
            signal = noise;
        end
        % run IPD model on signal
        [hairc_fine, fc, hairc_ild]=dietz2011(signal,fs,'fhigh',1400);
        % convert interaural information into azimuth
        itd_unwrapped = ...
            dietz2011_unwrapitd(hairc_fine.itd_lp,hairc_ild,hairc_fine.f_inst,2.5);
        angl=itd2angle(itd_unwrapped,lookup);

        h_ic=zeros(38,12);
        k(2*n,:)=sum(histc(angl,-92.5:5:92.5),2);
        for erb=1:12
            h_ic(:,erb)=histc(angl(hairc_fine.ic(:,erb)>ic_threshold&[diff(hairc_fine.ic(:,erb))>0; 0],erb),-92.5:5:92.5);
        end
        k(2*n-1,:)=sum(h_ic,2);
        example_output.angle_fine(:,:,n)=angl;
        example_output.IVS_fine = hairc_fine.ic;
        example_output.histogram_angle_label = -92.5:5:92.5;
        example_output.histograms_with_IVS(:,n)=sum(h_ic,2);
        example_output.histograms_without_IVS(:,n)=sum(histc(angl,-92.5:5:92.5),2);
    end

    if flags.do_plot
        % plot
        figure;
        set(gcf,'Position',[100 100 990 500])
        y=[0.14 0.57];
        x=[0.1 0.22 0.34 0.49 0.61 0.73 0.88];
        cols = 2;
        condi={'1S 0dB','','2S 0dB','',...
               '3S 0dB','','1S -6dB','',...
               '2S -6dB','','3S -6dB','','noise',''};
        
        for n = 1:14
            axes('position',[x(ceil(n/cols)) y(mod(n-1,cols)+1) 0.11 0.42],...
                 'box','on','fontSize',12)
            
            if mod(n,2)==0
                set(gca,'YTick',[10 20 30],'YTickLabel',{'','',''},'ylim',[0 32])
                set(gca,'Visible','on','XTick',[-30 0 30],...
                        'xlim',[-50 50],'XTickLabel',{'','','',''})
            else
                set(gca,'YTick',[1 2 3],'YTickLabel',{'','',''},'ylim',[0 3.33])
                set(gca,'Visible','on','XTick',[-30 0 30],...
                        'xlim',[-50 50],'XTickLabel',{'-30','0','+30'})
                xlabel(condi(n))
            end
            
            if n==2
                set(gca,'YTick',[10 20 30],'YTickLabel',{'10k','20k','30k'})
                ylabel('no IVS mask')
            elseif n==1        
                set(gca,'YTick',[1 2 3],'YTickLabel',{'1k','2k','3k'})
                ylabel('with IVS mask')
            end
            hold on
            bar(-92.5:5:92.5,k(n,:)/1000)
            colormap gray 
            
        end
        hold off
    end;
end;

if flags.do_fig6

    signal=sig_competingtalkers('one_speaker_reverb');
    fs = 16000;
    s_pos = [-45 0 45];
    ic_threshold=0.98;
    cn = [10 1]; % channel numbers for separate plots (1st entry also for time plot)
    panellabel = 'abcd';

    % run IPD model on signal
    [hairc_fine, fc, hairc_ild, hairc_mod]=dietz2011(signal,fs);
    % convert interaural information into azimuth
    itd_unwrapped = ...
        dietz2011_unwrapitd(hairc_fine.itd_lp,hairc_ild(:,1:12),hairc_fine.f_inst,2.5);
    angl=itd2angle(itd_unwrapped,lookup);
    angl_fmod=hairc_mod.itd_lp*140000; %linear approximation. paper version is better than this

    h_ic=zeros(71,12);
    h_fmod=zeros(71,11);
    h_all=histc(angl,-70:2:70);
    for n=1:12
        h_ic(:,n)=histc(angl(hairc_fine.ic(:,n)>ic_threshold&[diff(hairc_fine.ic(:,n))>0; 0],n),-70:2:70);
    end
    for n=1:11
        h_fmod(:,n)=histc(angl(hairc_mod.ic(:,n)>ic_threshold&[diff(hairc_mod.ic(:,n))>0; 0],n),-70:2:70);
    end
    example_output.angle_fine = angl;
    example_output.IVS_fine = hairc_fine.ic;
    example_output.histogram_angle_label = -70:2:70;
    example_output.histograms_with_IVS = h_ic;
    example_output.histograms_without_IVS = h_all;
    example_output.histogram_panel1 = sum(h_ic,2);
    example_output.histogram_panel2 = sum(h_ic(:,6:12),2);
    example_output.histogram_panel3 = sum(h_ic(:,1:3),2);
    example_output.histogram_panel4 = sum(h_fmod,2);

    if flags.do_plot
        figure;
        fontsize = 14;
        set(gcf,'Position',[60 100 1370 350])
        
        for panel = 1:4
            subplot(1,4,panel)
            switch panel
                
              case 1
                bar(-70:2:70,sum(h_ic,2))
                title('fine ITD channels 200-1400 Hz','Fontsize',fontsize)
                ymax = max(sum(h_ic,2));
              case 2
                bar(-70:2:70,sum(h_ic(:,6:12),2))
                title('fine ITD channels 500-1400 Hz','Fontsize',fontsize)
                ymax = max(sum(h_ic,2));
              case 3
                bar(-70:2:70,sum(h_ic(:,1:3),2))
                title('fine ITD channels 200-400 Hz','Fontsize',fontsize)
                ymax = max(sum(h_ic,2));
              case 4
                bar(-70:2:70,sum(h_fmod,2))
                title('mod ITD channels 1400-5000 Hz','Fontsize',fontsize)
                ymax = max(sum(h_ic,2));
            end
            set(gca,'Fontsize',fontsize)
            set(gca,'XTick',s_pos)
%             xlim([-73 73])
%             ylim([0 ymax*1.1])
						axis([-73 73 0 5500]);
						set(gca,'YTick',[1000:1000:5000],'YTickLabel',{'1k','2k','3k','4k','5k'});

            xlabel('Azimuth [deg]','Fontsize',fontsize)
            ylabel('Frequency of occurence','Fontsize',fontsize)
            text (-50,ymax*.97,panellabel(panel),'Fontsize',fontsize+1,'FontWeight','bold')
            
        end        
    end;    
end;

if flags.do_fig6_spille2013
  
  signal=amt_load('dietz2011','s123456.wav');
  fs=44100;
  s_pos =[-75 -40 -10 10 40 75];

  ic_threshold=0.98;
  cn = [10 1]; % channel numbers for separate plots (1st entry also for time plot)
  panellabel = 'acbd';

  % run IPD model on signal
  [hairc_fine,fc,hairc_ild]=dietz2011(signal,fs,'fhigh',1400);

  % convert interaural information into azimuth
  itd_unwrapped = ...
      dietz2011_unwrapitd(hairc_fine.itd_lp,hairc_ild,hairc_fine.f_inst_lp,2.5);
  lookup = amt_load('dietz2011','itd2angle_lookuptable.mat');
  angl=itd2angle(itd_unwrapped,lookup);

  h_ic=zeros(91,12);
  h_all=histc(angl,-90:2:90);
  for n=1:12
      h_ic(:,n)=histc(angl(hairc_fine.ic(:,n)>ic_threshold&[diff(hairc_fine.ic(:,n))>0; 0],n),-90:2:90);
  end

  % plot
  figure;
  fontsize = 14;
  set(gcf,'Position',[100 100 1170 700])

  for panel = 1:4
      subplot(2,2,panel)
      switch panel
          case 1
              bar(-90:2:90,h_all(:,cn(1)))
              hold on
              bar(-90:2:90,h_ic(:,cn(1)),'r')
              title(['Azimuth histogram of ch. ' num2str(cn(1)) ' at cf = ' ...
                  num2str(round(fc(cn(1)))) ' Hz'],'Fontsize',fontsize)
              ymax = max(h_all(:,cn(1)));
          case 2
              bar(-90:2:90,h_all(:,cn(2)))
              hold on
              bar(-90:2:90,h_ic(:,cn(2)),'r')
              title(['Azimuth histogram of ch. ' num2str(cn(2)) ' at cf = ' ...
                  num2str(round(fc(cn(2)))) ' Hz'],'Fontsize',fontsize)
              ymax = max(h_all(:,cn(2)));
          case 3
              t=(1:length(signal))*1/fs;
              plot(t,angl(:,cn(1)),'b.');
              hold on;
              plot(t(hairc_fine.ic(:,cn(1))>ic_threshold),...
                  angl(hairc_fine.ic(:,cn(1))>ic_threshold,cn(1)),'r.');
              title(['Azimuth over time in ch. ' num2str(cn(1)) ' at cf = ' ...
                  num2str(round(fc(cn(1)))) ' Hz'],'Fontsize',fontsize)
          case 4
              bar(-90:2:90,mean(h_all,2))
              hold on
              bar(-90:2:90,mean(h_ic,2),'r')
              title('Mean histogram of all fine-structure channels','Fontsize',fontsize)
              ymax = max(mean(h_all,2));
      end
      set(gca,'Fontsize',fontsize)
      if panel ~= 3
          set(gca,'XTick',s_pos)
          xlim([-93 93])
          ylim([0 ymax*1.1])
          xlabel('Azimuth [deg]','Fontsize',fontsize)
          ylabel('Frequency of occurence','Fontsize',fontsize)
          text (-80,ymax*.95,panellabel(panel),'Fontsize',fontsize+1,'FontWeight','bold')
      else
          set(gca,'YTick',s_pos)
          xlim([0 5.5])
          ylim([-95 95])
          ylabel('Azimuth [deg]','Fontsize',fontsize)
          xlabel('Time [s]','Fontsize',fontsize)
          rectangle('Position',[0.32,60,0.3,21],'FaceColor','white')
          text (0.38,69,panellabel(panel),'Fontsize',fontsize+1,'FontWeight','bold')
      end
  end

  example_output.angle_fine = angl;
  example_output.IVS_fine = hairc_fine.ic;
  example_output.histogram_angle_label = -90:2:90;
  example_output.histograms_with_IVS = h_ic;
  example_output.histograms_without_IVS = h_all;
  
end
