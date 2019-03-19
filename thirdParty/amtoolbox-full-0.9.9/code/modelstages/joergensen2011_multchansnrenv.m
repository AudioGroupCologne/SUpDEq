function result =  multChanSNRenv(Mix,Noise,fs,T)
%
%
% Usage: result =  multChanSNRenv(signal,test,noise,fs,N,T)
%
%   signal  : clean input signal
%   test    : distorted signal (by noise, reverb etc.)
%   noise   : noise signal, processed in the same way as the test signal
%   fs      : sampling frequency
%   N       : Number of samples in the input signals
%   T       : Duration in seconds of the input signals
%
%  Output:
%   result.mod_fcs               : Center frequencies of the modulation filterbank
%   result.outSNRenvs            : Matrix with an SNRenv value for each modulation filter in each gammatone filter;  
%   result.sEPSM_ExcPtns         : (Optional) Modulation excitation patterns
%
%  multChanSNRenv: calculates the SNRenv in 7 modulation filters in 22
%  Gammatone filters with 1/3-octave spacing.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/joergensen2011_multchansnrenv.php

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



debug = 0; % 0 =  no debug; 1 = debug with figures


% Center frequencies of the gammatone filters:
midfreq=[63    80   100  125  160  200    250    315  400   500  630 800  1000  1250  1600  2000 2500 3150 4000 5000 6300  8000];

%  diffuse field hearing threshold in quiet: ISO 387-7:2005
HT_diffuse = [37.5 31.5  26.5  22.1 17.9 14.4 11.4 8.4 5.8  3.8  2.1  1.0  0.8  1.9  0.5 -1.5 -3.1 -4.0 -3.8 -1.8 2.5   6.8 ];

stim(1,:) = Mix ;
stim(2,:) = Noise;


%% calculating filterbank time domain output
for k = 1:2
    g = gammatonefir(midfreq,fs,'complex');
    tmp  = 2*real(ufilterbank(stim(k,:),g,1));
    GT_output(:,:,k)  = tmp; % time domain outputs
end

%
%% ------------------ determining which frequency bands of the noise alone that are above
%% the hearing threshold. They are used for further processing.
  
% Filtering the mix with a rectangular 1/3-octave filterbank
[fcs1 mix_rms_out] = thirdoctrmsanalysis24(Mix,fs);
mixRMS_dB =  20*log10(mix_rms_out)+100;

mix_spec_level_corr = mixRMS_dB - 10*log10(fcs1*0.231); % converted to spectrum level according to ANSI 1997.
%     Find the bands with level above diffuse field Hearing Threshold
Nbands_to_process = find(mix_spec_level_corr(1:22)>HT_diffuse);

%% ------------------ calculating envelopes of temporal outputs

for p = Nbands_to_process
    
    clear stim
    stim(1,:) = GT_output(:,p,1) ;
    stim(2,:) = GT_output(:,p,2);
    
    % lowpass filtering at 150 Hz
    [bb, aa] = butter(1, 150*2/fs);
    
    for k = 1:2
        tmp = abs(hilbert(stim(k,:))); %  envelope
        tmp = filter(bb,aa,tmp);
        
        env(k,:) = tmp;
       
    end
    
    %% ------------estimation of SNRenv in this audio channel
    [mod_fcs outSNRenvs(:,p) ] = joergensen2011snrenv(env(1,:),env(2,:),fs);
    
end

% Saving output data:
result.mod_fcs = mod_fcs;
result.outSNRenvs = outSNRenvs;

if debug == 1
    cut = 1;
    t = 0:1/fsNew:T;
    %     t = t(cut:end-1);
    t2 = 0:1/fs:T;
    t2 = t2(1:end-1);
    linwidth = 2;
    fontsize = 14;
    D2 = 1;
    xmax = 100;
    xmin = 0.4;
    ymax = 1.2;
    ymin = -1.2;
    %     x2_normed = x2/max(x2);
    %     test2_normed = test2/max(test2);
    Bandidx = 14;
    figure
    subplot(3,1,1)
    plot(t2,GT_output(:,p,1)/max(GT_output(:,p,1)),'k','linewidth',linwidth),hold on
    plot(t2,GT_output(:,p,3)/max(GT_output(:,p,3)),'r','linewidth',linwidth)
    plot(t2,GT_output(:,p,2)/max(GT_output(:,p,2)),'color',[.5 .5 .5],'linewidth',linwidth)
    title(['Output from GT-band centered at ', num2str(midfreq(Bandidx)), ' Hz'])
    xlabel('Time [s]','FontSize',fontsize)
    ylabel('Amplitude ','FontSize',fontsize)
    % xlim([xmin xmax])
    ylim([ymin ymax])
    set(gca,'fontsize',fontsize );
    subplot(3,1,2)
    %     plot(t(1:D2:end),x2_normed(1:D2:end),'color',[.5 .5 .5],'linewidth',linwidth),hold on
    plot(t, env_normed(:,1,Bandidx)','k','linewidth',linwidth), hold on
    plot(t, env_normed(:,2,Bandidx)','color',[.5 .5 .5],'linewidth',linwidth)
    %     title(['speech' ' SNR = ' num2str(SNRs(end)) ])
    xlabel('Time [s]','FontSize',fontsize)
    ylabel('Amplitude ','FontSize',fontsize)
    % xlim([xmin xmax])
    ylim([ymin ymax])
    set(gca,'fontsize',fontsize );
    subplot(3,1,3)
    %     plot(t(1:D2:end),test2_normed(1:D2:end),'color',[.5 .5 .5],'linewidth',linwidth),hold on
    plot(t, env_normed(:,3,Bandidx)','k','linewidth',linwidth)
    % plot(t,env_int(3,:),'g','linewidth',linwidth)
    title('speech + noise')
    xlabel('Time [s]','FontSize',fontsize)
    ylabel('Amplitude ','FontSize',fontsize)
    % xlim([xmin xmax])
    ylim([ymin ymax])
    set(gca,'fontsize',fontsize );
    
    Nnew = length(env_normed(:,1,Bandidx));
    for k = 1:3
        tmp =  abs(fft(env_normed(:,k,Bandidx))/(Nnew/2)).^2 ; % only norm with N/2, why?
        outSxx(:,k) = tmp(1:fix(Nnew/2)+1); %(1:fix(N/2)+1)
        
        pos_freqs= linspace(0,fsNew/2,length(outSxx(:,k)));
        [fcs1 mod_spec(:,k)] = octave_third_low(pos_freqs,outSxx(:,k));
        [fcs_EPSM, outSxxEPSM(:,k)] =  EPSM3(pos_freqs,outSxx(:,k)',0);
    end
    
    % tmp = find(isnan(mod_spec));
    % mod_spec(tmp) = 0.00101;
    octbands = 1:3:25;
    xmin = .5;
    xmax = 60;
    ymin = -60;
    ymax = 1;
    figure
    subplot(1,2,1)
    % plot(pos_freqs, 10*log10(specMean(3,:)),'- k','linewidth',linwidth),hold on
    % plot(pos_freqs, 10*log10(specMean(1,:)),'- k','linewidth',linwidth),hold on
    % plot(fcs1, 10*log10(mod_spec(:,1)),'- s k','linewidth',linwidth), hold on
    % plot(fcs1, 10*log10(mod_spec(:,3)),': * k','linewidth',linwidth),
    % plot(fcs1, 10*log10(mod_spec(:,2)),'color',[.5 .5 .5],'linewidth',linwidth),
    plot(pos_freqs, 10*log10(outSxx(:,1)),'- s k','linewidth',linwidth), hold on
    plot(pos_freqs, 10*log10(outSxx(:,3)),': * k','linewidth',linwidth),
    plot(pos_freqs, 10*log10(outSxx(:,2)),'color',[.5 .5 .5],'linewidth',linwidth),
    xlim([xmin xmax])
    ylim([ymin ymax])
    xlabel('Modulation frequency [Hz]','FontSize',fontsize)
    ylabel('Amplitude ','FontSize',fontsize)
    set(gca,'XTick',fcs1(octbands),'XScale','log','XMinorTick','on','ytick',ymin:10:ymax,'fontsize',12 );
    legend('clean','noisy','noise')
    ymax = 10;
    ymin = -25;
    subplot(1,2,2)
    plot([1 fcs_EPSM], 10*log10(outSxxEPSM(:,1)),'-s k','linewidth',linwidth),hold on
    plot([1 fcs_EPSM], 10*log10(outSxxEPSM(:,3)),': * k','linewidth',linwidth)
    plot([1 fcs_EPSM], 10*log10(outSxxEPSM(:,2)),'color',[.5 .5 .5],'linewidth',linwidth)
    xlim([xmin xmax])
    ylim([ymin ymax])
    xlabel('Modulation filter f_c [Hz]','FontSize',fontsize)
    ylabel('Envelope power','FontSize',fontsize)
    legend('clean','noisy','noise')
    set(gca,'XTick',[1 fcs_EPSM],'XScale','log','XMinorTick','on','ytick',-25:2:10,'fontsize',12);
end


