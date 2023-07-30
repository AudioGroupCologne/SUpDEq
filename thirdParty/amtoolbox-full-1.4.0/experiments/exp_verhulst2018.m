function data = exp_verhulst2018(varargin)
%EXP_VERHULST2018 Figures from Verhulst et al. (2018)
%
%   Usage: output = exp_verhulst2018(flag)
%
%   This function can be used to obtain Figures 3C, 5A, or 7A from the paper
%   by Verhulst, Altoe, and Vasilkov (2018). These simulations might take
%   some time depending on the computing power of your computer. On an average
%   laptop the time required to generate Figs 3C, 5A, is of about 5 minutes
%   and of about 8-9 minutes for Fig 7A.
%
%   - Fig. 3C compares the simulated auditory nerve rates using the verhulst2018
%     and verhulst2015 models to an AM tone with fc=4 kHz, fmod=100 Hz,
%     and level of 60 dB SPL.
%   - Fig. 5A Computes on-frequency average spiking rates of high-, medium-, 
%     and low-spontaneous rate neurons for pure tones of different level
%     with either carriers of 1 kHz or 4 kHz.
%   - Fig. 7A shows the envelope-following response amplitudes to 4 kHz AM tone
%     fmod = 98 Hz, and levels between 45 and 80 dB, when the cochlear
%     profile associated to a normal audiogram (Flat00) is used compared
%     to a cochlear profile with a normal audiogram up to 1 kHz and
%     with a sloping hearing loss that produces a hearing threshold
%     of 35 dB HL at 8 kHz (Flat00_Slope35).
%     The EFR amplitude with the Flat00 profile is also compared between
%     a no-synaptopathy condition (all neurons: 13-3-3) and a synaptopathy
%     profile where the low and middle spontaneous rate neurons have
%     been removed (13-0-0).
%
%   Examples:
%   ---------
%
%   To display Figure 3c from Verhulst et al. (2018) use:
%
%     exp_verhulst2018('fig3c');
%
%   To display Figure 5a from Verhulst et al. (2018) use:
%
%     exp_verhulst2018('fig5a');
%
%   To display Figure 7a from Verhulst et al. (2018) use:
%
%     exp_verhulst2018('fig7a');
%
%   License:
%   --------
%
%   This model is licensed under the UGent Academic License. Further usage details are provided 
%   in the UGent Academic License which can be found in the AMT directory "licences" and at 
%   <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.
%
%   References:
%     S. Verhulst, H. Bharadwaj, G. Mehraei, C. Shera, and
%     B. Shinn-Cunningham. Functional modeling of the human auditory
%     brainstem response to broadband stimulation. jasa, 138(3):1637--1659,
%     2015.
%     
%     S. Verhulst, A. Altoè, and V. Vasilkov. Functional modeling of the
%     human auditory brainstem response to broadband stimulation.
%     hearingresearch, 360:55--75, 2018.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_verhulst2018.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   #Author: Alejandro Osses (2020): primary implementation for the AMT
%   #Author: Piotr Majdak (2021): adapted to the AMT 1.0
%   #License: ugent

data = [];

definput.import={'amt_cache'}; % from arg_amt_cache
definput.flags.disp = {'no_debug','debug'}; % flag to provide debugging information when model called with 'debug', see amt_disp
definput.flags.plot={'plot','no_plot'};

definput.flags.type={'missingflag','fig3c','fig5a','fig7a'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%%% Common parameters:
fs = 44100; % Hz, sampling frequency in Hz. Can be any sampling frequency
            % but in the model the input signals are always resampled to 100 kHz
cf_flag = 'abr';

dBFS = 94; % dB full scale. Amplitude 1 = 1 Pascal
if dBFS == 94
    p0 = 2e-5; % Pa, reference pressure
end

%% ------ FIG 3c Verhulst, Altoe, Vasilkov (2018) -------------------------
if flags.do_fig3c
    %%% Stimulus parameters
    fc       = 4000; % Hz, frequency of the carrier
    fmod     =  100; % Hz, frequency of the modulator
    mdepth   = 1; % value between 0 and 1
    dur      = 500e-3; % stimulus duration in seconds
    dur_ramp = 2.5e-3; % s, duration of the ramp

    lvl = 60; % level, dB
    
    N_samples = round(dur*fs);
    dur_ramp_samples = round((dur_ramp)*fs);

    % Creating a cosine ramp:
    ramp = ones(N_samples,1);
    ramp(1:dur_ramp_samples)         = rampup(dur_ramp_samples);
    ramp(end-dur_ramp_samples+1:end) = rampdown(dur_ramp_samples);

    % AM stimulus and calibration:
    t = (0:N_samples-1)/fs; t=t(:);
    carrier = sin(2*pi*fc*t); % starts at phase = 0
    env = (1 + mdepth * sin(2*pi*fmod*t-pi/2) ); % modulator starts at minimum (phase=-pi/2)
    insig = env .* carrier; % Amplitude-modulated signal
    insig = scaletodbspl(insig,lvl,dBFS);
    insig = ramp.*insig;

    %%% Model parameters (only one hearing profile, Flat00 and 13-3-3):
    hear_profile = 'Flat00';
    numH = 13; % Number of neurons, HSR
    numM =  3;
    numL =  3;
    
    cf_flag_here = fc; % on-frequency simulation only (1 channel)
    outputs = verhulst2018(insig,fs,cf_flag_here,'hearing_profile',hear_profile,...
				'numL',numL,'numM',numM,'numH',numH, ... % number of AN neurons
				'anfH',... % provide detailed information about the high-SR fibres
				'no_ihc', ... % detailed information about the IHC not required here
				'no_cn','no_ic',... % simulations of CN and IC not required here
				flags.disp);
				
  	out2015 = verhulst2015(insig,fs,cf_flag_here,'hearing_profile',hear_profile,...
        'numL',numL,'numM',numM,'numH',numH,...
        'anfH',... % provide detailed information about the high-SR fibres
				'no_ihc', ... % detailed information about the IHC not required here
        'no_cn','no_ic',flags.disp);                        
    
    fs_abr = outputs.fs_abr;
    t_anf = (1:length(outputs.anfH))'/outputs.fs_abr;
    
    figure;
    plot(t_anf*1000,outputs.anfH); hold on, grid on
    plot(t_anf*1000,out2015.anfH,'k');
    xlabel('Time [ms]')
    ylabel('Firing rate [spikes/s]');
    
    num_tot = numH+numM+numL;
    
    %%%
    L_bin = 15; % samples
    percent = 90; % percentage overlap
    psthbinwidth = L_bin/fs_abr;
    L_overlap = round((percent/100)*L_bin); 
    %%%
    
    t_psth = buffer(t_anf, L_bin, L_overlap,'nodelay');
    t_psth = t_psth(1,:);
    
    anfH     = buffer(outputs.anfH, L_bin, L_overlap,'nodelay'); anfH     = mean(anfH);
    anfH2015 = buffer(out2015.anfH, L_bin, L_overlap,'nodelay'); anfH2015 = mean(anfH2015);
    
    anf     = buffer(outputs.an_summed/num_tot, L_bin, L_overlap,'nodelay'); anf     = mean(anf);
    anf2015 = buffer(out2015.an_summed/num_tot, L_bin, L_overlap,'nodelay'); anf2015 = mean(anf2015);
    
    if flags.do_plot
        XL = [350 380];
        YL = [-30 330];
    
        figure;
        plot(t_psth*1000,anf,'b-','LineWidth',2); hold on, grid on
        plot(t_psth*1000,anf2015,'k--','LineWidth',2);
        xlim(XL);
        ylim(YL);
        xlabel('Time [ms]');
        ylabel('AN firing rate [spikes/s]');
        legend('Verhulst et al. (2018)','Verhulst et al. (2015)');
        title(sprintf('Auditory Nerve model output (%.0f-%.0f-%.0f neurons; bin size=%.2f ms, %.1f percent overlap)',numH,numM,numL,psthbinwidth*1000,percent));

        figure;
        plot(t_psth*1000,anfH,'b-','LineWidth',2); hold on, grid on
        plot(t_psth*1000,anfH2015,'k--','LineWidth',2);

        xlim(XL);
        ylim(YL);
        xlabel('Time [ms]');
        ylabel('Firing rate [spikes/s]');
        legend('Verhulst et al. (2018)','Verhulst et al. (2015)');
        title(sprintf('High-spontaneous rate neuron output (bin size=%.2f ms, %.1f percent overlap)',psthbinwidth*1000,percent));
    end
    
    data.anfH = anfH; % output one HSR neuron
    data.anf  = anf;  % average output using 19 neurons (13-3-3)
    data.anfH2015 = anfH2015; % output one HSR neuron using Verhulst et al. (2015)
    data.anf2015  = anf2015;  % average output using 19 neurons (13-3-3) using Verhulst et al. (2015)
    data.anf_unit = 'spikes/s';
    data.fc   = fc; 
end

%% ------ FIG 5a Verhulst, Altoe, Vasilkov (2018) -------------------------
if flags.do_fig5a
 
    %%% Model parameters (only one hearing profile, Flat00 and 13-3-3):
    hear_profile = 'Flat00';
    numH = 13; % Number of neurons, HSR
    numM =  3;
    numL =  3;
    
    %%% Stimulus parameters:
    fc   = [1000 4000]; % 8e3;  % CF in Hz;
    N_fc = length(fc);
    LW = [1 2]; % LineWidth to be used in the plots (one for ech CF)
    
    lvls = 0:10:100; % test levels
    L = length(lvls);

    dur   = 50e-3;  % stimulus duration in seconds
    dur_ramp = 2.5e-3; % s, duration of the up/down ramp
    
    % Memory allocation for output variable:
    spike_rates = zeros(N_fc,3,L); % 3: 1 for HSR, 1 for MSR, 1 for LSR
    
    t = 0:1/fs:dur-1/fs; % time vector
    N_samples = length(t); % length (duration) of the input signal in samples
    dur_ramp_samples = round(dur_ramp*fs);

    % Linear ramps:
    ramp_up = (0:(dur_ramp_samples-1))'/dur_ramp_samples;
    ramp_dn = (dur_ramp_samples:-1:0)'/dur_ramp_samples;

    for i = 1:N_fc
        CF = fc(i); % stimulus frequency in Hz
        insig_orig(:,i) = sin(2*pi*CF*t'); % column array

        if i == 1
            ramp4insig = ones(size(insig_orig(:,i)));
            idxs = 1:dur_ramp_samples;
            ramp4insig(idxs) = ramp4insig(idxs) .* ramp_up;
            idxs = (N_samples-dur_ramp_samples):N_samples;
            ramp4insig(idxs) = ramp4insig(idxs) .* ramp_dn;
        end

        insig_orig(:,i) = insig_orig(:,i) .* ramp4insig; % applying the ramp
    end

    if flags.do_plot
        % New (empty) figure where the results for each CF will be appended:
        figure
    end

    for i = 1:N_fc

        CF = fc(i); % stimulus frequency in Hz
        
        for k = 1:L
            insigs(:,k) = sqrt(2)*p0*10^(lvls(k)/20)*insig_orig(:,i); % calibrated stimulus
        end
        % insigs = repmat(insigs,nrep,num_stims);
        outputs = verhulst2018(insigs,fs,cf_flag,'hearing_profile',hear_profile,...
					'numL',numL,'numM',numM,'numH',numH, ... % number of AN neurons
					'anfL','anfM','anfH', ... % provide detailed information about all types of AN neurons
					'no_ihc', ... % detailed information about IHC not required here
					'no_cn', 'no_ic', ... % simulations of CN and IC not required here
					flags.disp);
					
        idx_cf = find(outputs(1).cf>CF,1,'last');

        fs_an = outputs(1).fs_an;

        idxi = round(15e-3*fs_an)+1; % start after 15 ms to skip the strong onset
        idxf = round(dur*fs_an);     % end 

        for k = 1:L
            % Reads the corresponding bin (idx_cf):
            psthL =outputs(k).anfL(:,idx_cf); 
            psthM =outputs(k).anfM(:,idx_cf);
            psthH =outputs(k).anfH(:,idx_cf);
            
            % Averaging
            spike_rates(i,1,k)= mean(psthH(idxi:idxf));
            spike_rates(i,2,k)= mean(psthM(idxi:idxf));
            spike_rates(i,3,k)= mean(psthL(idxi:idxf));
        end

        if flags.do_plot
            plot(lvls,squeeze(spike_rates(i,1,:)),'ro-' ,'LineWidth',LW(i)); hold on, grid on
            plot(lvls,squeeze(spike_rates(i,2,:)),'bs--','LineWidth',LW(i));
            plot(lvls,squeeze(spike_rates(i,3,:)),'md-' ,'LineWidth',LW(i));
        end
        %%%
    end
    
    if flags.do_plot
        xlim([min(lvls) max(lvls)])
        title(sprintf('Average firing rates at CF=%.1f Hz',CF))
        xlabel('Stimulus Level (dB SPL)')
        ylabel('Firing Rate (/s)')
        % legend(sprintf('LSR (%.0f units)',numsponts(1)),sprintf('MSR (%.0f units)',numsponts(2)),sprintf('HSR (%.0f units)',numsponts(3)))
    end
    
    data.figure_flag = 'do_fig5a';
    data.spikes_rates = spike_rates;
    data.spikes_rates_unit = 'spikes/s';
    data.spikes_rates_description = 'avg. spiking rates (onset excluded) for i=carrier frequency (x2); j=1,2,3 for HSR,MSR,LSR neurons; k=test levels';
    data.lvls = lvls;
end

%% ------ FIG 7a Verhulst, Altoe, Vasilkov (2018) -------------------------
if flags.do_fig7a
    
    %%% Model parameters:
    % For calculations after the simulations:
    pars = arg_verhulst2018; % load all verhulst2018 default parameters;
    M1 = pars.keyvals.M1; 
    M3 = pars.keyvals.M3;
    M5 = pars.keyvals.M5;
   
    %%% Stimulus parameters
    fc       = 4000; % Hz, frequency of the carrier
    fmod     =   98; % Hz, frequency of the modulator
    mdepth   = 0.85; % value between 0 and 1
    dur      = 100e-3; % stimulus duration in seconds
    dur_ramp = 2.5e-3; % s, duration of the ramp

    lvls = 45:5:80; % level, dB
    L  = length(lvls);

    N_samples = round(dur*fs);
    dur_ramp_samples = round((dur_ramp)*fs);

    % Creating a cosine ramp:
    ramp = ones(N_samples,1);
    ramp(1:dur_ramp_samples)         = rampup(dur_ramp_samples);
    ramp(end-dur_ramp_samples+1:end) = rampdown(dur_ramp_samples);

    t = (0:N_samples-1)/fs; t=t(:);
    carrier = sin(2*pi*fc*t); % starts at phase = 0
    env = (1 + mdepth * sin(2*pi*fmod*t-pi/2) ); % modulator starts at minimum (phase=-pi/2)
    insig_orig = env .* carrier; % Amplitude-modulated signal
    %%%

    insig = nan(round(dur*fs),L);
    for i = 1:L
        insig(:,i) = scaletodbspl(insig_orig,lvls(i),dBFS);
        insig(:,i) = ramp.*insig(:,i);
    end

    %%% Simulation for 3 cochlear+synaptopathy profiles:
    for k = 1:3

        % tic
        out = [];
        switch k
            case 1
                numH = 13;
                numM =  3;
                numL =  3;
                out = verhulst2018(insig,fs,cf_flag,'hearing_profile','Flat00','no_ihc','numH',numH,'numM',numM,'numL',numL,flags.disp);
            case 2
                numH = 13;
                numM =  3;
                numL =  3;
                out = verhulst2018(insig,fs,cf_flag,'hearing_profile','Flat00_Slope35','no_ihc','numH',numH,'numM',numM,'numL',numL,flags.disp);
            case 3
                numH = 13;
                numM =  0;
                numL =  0;
                out = verhulst2018(insig,fs,cf_flag,'hearing_profile','Flat00','no_ihc','numH',numH,'numM',numM,'numL',numL,flags.disp);
        end
        % toc

        fs_abr = out.fs_abr;
        K = fs_abr/2;
        N = length(out(1).w5);
        for i = 1:L
            AN = out(i).an_summed; 
            CN = out(i).cn;
            IC = out(i).ic;        

            EFR_sim(i,:,k) = M1*sum(AN,2)+M3*sum(CN,2)+M5*sum(IC,2); % EFRs in time domain

            [hAm(:,i,k),f] = freqz(EFR_sim(i,:,k),1,K,fs_abr);
            hAm(:,i,k) = hAm(:,i,k)/N; % Parseval's theorem: [rmsdb(EFR_sim(i,:,k)) rmsdb(hAm(:,i,k))]

            idx1 = fmod-20; % makes sure that the DC is not accounted for
            idx2 = K;
            [EFR_amp(i,k) ,idx_max(i,k)]  = max(abs( hAm(idx1:idx2,i,k))); 
        end
    end
    EFR_amp  = 20*log10(EFR_amp /1e-6); % converting to dB re 1uV
    
    if flags.do_plot
        figure;
        plot(lvls,EFR_amp(:,1),'s-' ,'LineWidth',2,'Color','b'), hold on, grid on;
        plot(lvls,EFR_amp(:,2),'o-' ,'LineWidth',2,'Color','r');
        plot(lvls,EFR_amp(:,3),'d--','LineWidth',2,'Color','m');
        hl = legend('NH', 'HI', 'HSR','Location','NorthWest');
        set(hl,'FontSize',10);
        legend('boxoff');
        xlabel('Stimulus Level [dB SPL]');
        ylabel('EFR Amplitude [dB re. 1\mu V]')
        ylim([-62 -18]);
        xlim([40 85]);
        set(gca,'XTick',45:5:80);
        set(gca,'YTick',-60:5:-20);
    end
    
    data.figure_flag  = 'do_fig7a';
    data.EFR_amp      = EFR_amp;
    data.EFR_amp_unit = 'dB re. 1uV';
    data.spikes_rates_description = 'amplitude of the envelope-following response (EFR) for AM tones of i=test levels; k=hearing profiles';
    data.lvls = lvls;
end


