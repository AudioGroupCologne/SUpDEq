function data = exp_osses2021(varargin)
%EXP_OSSES2021 monaural perceptual similarity
%
%   Usage: data = exp_osses2021(flags)
%
%   exp_osses2021 reproduces Figs. 4, 11 and 14 from Osses and Kohlrausch (2021), 
%   where a modified version of the dau1997 model is used. 
%   The figures are similar to Figs. 4.14, C.9B, and C.11B from Osses (2018).
%
%
%   The following flags can be specified:
%
%     'redo'    Recomputes data for specified figure
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'no_plot'  Don't plot, only return data.
%
%     'fig4_osses2020'    Reproduce Fig. 4 of Osses et al. (2020).
%
%     'fig14_osses2020'    Reproduce Fig. 14 of  Osses et al. (2020).
%     'fig4'     Reproduce Fig. 4 of Osses and Kohlrausch (2021). This is the
%                same figure as Fig. 4 of Osses and Kohlrausch (2020, preprint)
%
%     'fig11'    Reproduce Fig. 11 of Osses and Kohlrausch (2021). This is the
%                same figure as Fig. 14 of Osses and Kohlrausch (2020, preprint)
%
%     'fig14'    Reproduce Fig. 14 of Osses and Kohlrausch (2021).
%
%   Fig. 4 - Two internal representations of a piano sound ('P1') using the  
%   PEMO model with two configurations of the adaptation loops are shown:
%   Overshoot limitation with a factor of 5, as suggested in Osses and 
%   Kohlrausch (2021), and with a factor of 10 (see Dau et al., 1997).
%   To display Fig. 4 of Osses and Kohlrausch (2021) use :
%
%     out = exp_osses2021('fig4'); % same as: out = exp_dau1997('fig4_osses2020');
%
%   Fig. 11 - The effect of the overshoot limitation with factors of 5 and 10
%   are shown for a 4-kHz pure tone of 70 dB SPL that includes 2.5-ms up/down 
%   ramps. For these plots the outer and middle ear stages are skipped. One
%   gammatone filter at 4 kHz is used, followed by the IHC stage (ihc_breebaart2001),
%   and the adaptation loops (adt_osses2021 for lim=5, adt_dau1997 for lim=10).
%
%   To display Fig. 11 of Osses and (2020) use :
%
%     out = exp_osses2021('fig11'); % same as: out = exp_dau1997('fig14_osses2020');
%
%   Fig. 14 - Modulation transfer functions for the 12 modulation filters in
%   modfilterbank.m. This figure is obtained for a click of unit amplitude
%   while all modules in osses2021 are by-passed except for the modulation
%   filter bank. In the modulation filter bank, the phase insensitivity for
%   filters with mfc>10 is disabled (see App. C of Osses and Kohlrausch, 2021)
%   in a way that the outputs of the modulation filter are the complex valued
%   filtered signals for each mfc. 
%   To display Fig. 14 of Osses and Kohlrausch (2021) use :
%
%     out = exp_osses2021('fig14'); 
%
%
%   References:
%     A. Osses Vecchi and A. Kohlrausch. Perceptual similarity between piano
%     notes: Simulations with a template-based perception model. J. Acoust.
%     Soc. Am., 141(4), 2020.
%     
%     M. Jepsen, S. Ewert, and T. Dau. A computational model of human
%     auditory signal processing and perception. J. Acoust. Soc. Am., 124(1),
%     2008.
%     
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band
%     carriers. J. Acoust. Soc. Am., 102:2892--2905, 1997a.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_osses2021.php


%   #Author: Alejandro Osses (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 




data = [];

definput.import={'amt_cache'};
definput.flags.type={'missingflag','fig4','fig11','fig14'}; % Osses and Kohlrausch 2021, JASA

definput.flags.plot={'plot','no_plot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

%% ------ FIG 4 Osses and Kohlrausch 2021 ---------------------------------
if flags.do_fig4

    % Goes to 'http://sofacoustics.org/data/amt-0.10.0/auxdata/' and loads
    %   the data under the folder 'dau1997', with name 'P1-GH05-Cd5_1-dur-1300-ms.wav'
    [insig, fs] = amt_load('dau1997','P1-GH05-Cd5_1-dur-1300-ms.wav');
    
    tobs = 0.25; % s, only the first 0.25 s of the waveform
    insig = insig(1:fs*tobs);
    
    subfs = 16000; % sampling frequency for the internal representation

    
    %%% Using osses2021:
    flags_common = {'afb_osses2021','ihc_breebaart2001','mfb_jepsen2008'};
    [outsig05,fc,mfc] = osses2021(insig,fs,'subfs',subfs,flags_common{:},'adt_osses2021');
    outsig10          = osses2021(insig,fs,'subfs',subfs,flags_common{:},'adt_dau1997'); 
    
    N = length(fc); % number of audio frequencies
    K = length(mfc);
    %%% Memory allocation:
    Ik_05 = zeros(1,K); % one value for each modulation frequency
    Ik_10 = zeros(1,K); % one value for each modulation frequency
    
    Im_05 = zeros(1,N);  % one value for each audio frequency
    Im_10 = zeros(1,N);  % one value for each audio frequency
    %%%
    
    for j = 1:N
        Im_05(j) = Im_05(j)+1/subfs*sum(outsig05{j}(:).^2); % all mod. filters together
        Im_10(j) = Im_10(j)+1/subfs*sum(outsig10{j}(:).^2);
    end
    
    for i = 1:K
        for j = 1:N
            Nr_mod_filters = size(outsig05{j},2);
            if i <= Nr_mod_filters 
                % summed up only if the mod filter i is present.
                Ik_05(i) = Ik_05(i)+1/subfs*sum(outsig05{j}(:,i).^2);  
                Ik_10(i) = Ik_10(i)+1/subfs*sum(outsig10{j}(:,i).^2);
            else
                % Nothing to do, in this case mfc(i) > 1/4*fc(j)
            end  
        end
    end
    
    Itot_05 = sum(Ik_05);
    Itot_10 = sum(Ik_10);
    
    fc_erb = freqtoaud(fc);
    mfc_nr = 1:K;
    
    if flags.do_plot
        
        % Panel A
        figure;
        plot(fc_erb,100*Im_05/Itot_05,'ro-'); hold on; grid on;
        plot(fc_erb,100*Im_10/Itot_10,'bs--');

        set(gca,'XTick',3:33);
        set(gca,'YTick',0:1:10);
        ylim([-1 11])
        xlim([2 34])
        legend('lim=5 (Osses2021)','lim=10');
        Pos = get(gcf,'Position');
        Pos(3) = 800;
        set(gcf,'Position',Pos); % setting width of the figure

        xlabel('Audio centre frequency f_c (ERB_N)');
        ylabel('Percentage (%)');

        title('A. Information-based audio-frequency analysis')
        
        % Panel B
        figure;
        plot(mfc_nr,100*Ik_05/Itot_05,'ro-'); hold on; grid on;
        plot(mfc_nr,100*Ik_10/Itot_10,'bs--'); 

        set(gca,'XTick',1:12);
        set(gca,'YTick',0:3:24);
        xlim([0 13])
        ylim([-1 27])
        legend('lim=5 (Osses2021)','lim=10');

        xlabel('Modulation centre frequency mf_c (Nr.)');
        ylabel('Percentage (%)');
    
        title('B. Information-based modulation-frequency analysis')
    end
    
    data.figure_flag = 'do_fig4';
    data.fc = fc;
    data.mfc = mfc;
    data.Ik_05 = Ik_05;
    data.Ik_10 = Ik_10;
    data.Ik = 'Model Units (MU)';
    data.Im_05 = Im_05;
    data.Im_10 = Im_10;
    data.Im_unit = 'Model Units (MU)';
    data.Itot_05 = Itot_05;
    data.Itot_10 = Itot_10;
    data.Itot_unit = 'Model Units (MU)';
    data.description = 'Energy content (MU) for each audio frequency band n, and each modulation frequency band k';
end
  
%% ------ FIG 11 Osses and Kohlrausch 2021 --------------------------------
if flags.do_fig11
    
    % 1. Stimulus creation:
    fs = 44100;
    dur = 300e-3; % 300 ms
    lvl = 70;
    dBFS = 100; % AMT default
    
    t = (1:dur*fs)/fs; t = t(:); % creates 't' as a column array
    fc = 4000;
    dur_ramp_ms = 2.5;
    dur_ramp = round((dur_ramp_ms*1e-3)*fs); % duration ramp in samples

    insig = sin(2*pi*fc.*t);
    insig = scaletodbspl(insig,lvl,dBFS); % calibration before applying the ramp
    
    rp    = ones(size(insig)); 
    rp(1:dur_ramp) = rampup(dur_ramp);
    rp(end-dur_ramp+1:end) = rampdown(dur_ramp);
    insig = rp.*insig;

    insig = [zeros(50e-3*fs,1); insig; zeros(200e-3*fs,1)]; % 50 and 200 ms 
          % of silence before and after the sine tone
    t = (1:length(insig))/fs; t = t(:); % creates 't' as a column array
    
    kv_here    = {'basef',fc,'flow',fc,'fhigh',fc}; 
    flags_here = {'no_outerear','no_middleear','ihc','adt','no_mfb'};
        
    outsig05 = osses2021(insig,fs,kv_here{:},flags_here{:}); 
    outsig10 = osses2021(insig,fs,kv_here{:},'adt_dau1997',flags_here{:});
    
    %%%
    onset05 = max(outsig05);
    onset10 = max(outsig10);
    
    id_steady = find(t>=.330-dur_ramp_ms*1e-3 & t<=.350-dur_ramp_ms*1e-3);

    steady05 = mean(outsig05(id_steady));
    steady10 = mean(outsig10(id_steady));
    
    ra05 = onset05/steady05;
    ra10 = onset10/steady10;
    
    fprintf('Lim  5: Onset = %.1f MU, steady = %.1f MU\n',onset05,steady05);
    fprintf('Lim 10: Onset = %.1f MU, steady = %.1f MU\n',onset10,steady10);
    
    if flags.do_plot
        leg4plot{1} = sprintf('lim =  5: ratio onset/steady=%.1f',ra05);
        leg4plot{2} = sprintf('lim = 10: ratio onset/steady=%.1f',ra10);
    
        figure; 
        plot(t,outsig05,'r-','LineWidth',2); hold on, grid on;
        plot(t,outsig10,'b-')
       
        xlabel('Time (s)');
        ylabel('Amplitude \Psi (Model Units)')
        legend(leg4plot);
        
        set(gca,'XTick',[50:50:450]*1e-3);
        xlim([0.025 0.475]);
        ylim([-300 1480])
        set(gca,'YTick',-200:100:1400)
    end
    
    data.figure_flag = 'do_fig11';
    data.ra05 = ra05;
    data.onset05 = onset05;
    data.steady05 = steady05;
    data.ra10 = ra10;
    data.onset10 = onset10;
    data.steady10 = steady10;    
    data.onset_unit = 'Model Units (MU)';
    data.steady_unit = 'Model Units (MU)';
end

%% ------ FIG 14 Osses and Kohlrausch 2021 --------------------------------
if flags.do_fig14
    % From local file: g20210330_characterising_MFB_review2021.m    
    
    N = 2^16; % arbitrary number of samples: more samples = more FFT resolution,
              % (no zero-padding or whatsoever)
    K = N/2;
    fs = 44100;
    insig = [zeros(N/2-1,1); 1; zeros(N/2,1)]; % temporally centred dirac
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1 Calculation with impulse responses:
    kv = {'silent_mode',0}; % to activate the screen display
    flags_common = {'no_outerear','no_middleear','no_afb','no_ihc','no_adt','mfb','no_phase_insens'};
    [out1,~,~,~]        = osses2021(insig,fs,kv{:},'mfb_dau1997'           ,flags_common{:});
    [out2,~,mfc,params] = osses2021(insig,fs,kv{:},'mfb_jepsen2008'         ,flags_common{:});
    [out3,~,~,~]        = osses2021(insig,fs,kv{:},'mfb_osses2021_att_gain',flags_common{:});
     
    out1 = out1{1};
    out2 = out2{1};
    out3 = out3{1};
    
    for k = 1:size(out1,2)
        HdB(:,k)     = 20*log10(abs(freqz(out1(:,k),1,K))); % mfb_dau1997 % HdB     = To_dB(abs(hresp));
        HdB_tot(:,k) = 20*log10(abs(freqz(out2(:,k),1,K))); % mfb_jepsen2008
        HdB_new(:,k) = 20*log10(abs(freqz(out3(:,k),1,K))); % mfb_osses2021_att_gain
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.2 Calculation using the coefficients:
    f = (0:K-1)/K*(fs/2);
    for k = 1:size(out1,2)
        
        h_here = filter(params.mfb_b(k,:),params.mfb_a(k,:),insig);
        HdB_coeff(:,k) = 20*log10(abs(freqz(h_here,1,K))); % mfb_dau1997 % HdB     = To_dB(abs(hresp));
        
        insig_filt = filter(params.b_lp_150_Hz,params.a_lp_150_Hz,insig);
        h_here = filter(params.mfb_b(k,:),params.mfb_a(k,:),insig_filt);
        HdB_tot_coeff(:,k) = 20*log10(abs(freqz(h_here,1,K))); % mfb_jepsen2008
        
        hlpf = freqz(insig_filt,1,K);
        hlpf_dB = 20*log10(abs(hlpf));
        
        gain_mfcs = interp1(f,hlpf_dB,mfc); % gains according to the 150-Hz filter bank
        gain_lin = 10.^(gain_mfcs/20);
        h_here = filter(gain_lin(k)*params.mfb_b(k,:),params.mfb_a(k,:),insig);
        HdB_new_coeff(:,k) = 20*log10(abs(freqz(h_here,1,K))); % mfb_osses2021_att_gain
    end
    

    if flags.do_plot
        %%% 2. Plotting
        % Trick to add XTick labels
        XLT = [];
        XL = [0.1 1 2.5 5 10 20 50 125 250 500 1000 2000];
        for i = 1:length(XL)
            if XL(i) < 1000
                XLT{i} = num2str(XL(i));
            else
                XLT{i} = [num2str(XL(i)/1000) 'k'];
            end
        end

        h = [];

        N_mfb_types = 3;
        figure;
        for j = 1:N_mfb_types   
            % figure(j);
            subplot(3,1,j);
            for i = 1:12;
                switch j
                    case 1 % No LPF
                        HdB_here  = HdB(:,i);
                        HdB_here2 = HdB_coeff(:,i);
                    case 2 % With LPF
                        HdB_here  = HdB_tot(:,i);
                        HdB_here2 = HdB_tot_coeff(:,i);
                    case 3 % With LPF as suggested
                        HdB_here  = HdB_new(:,i);
                        HdB_here2 = HdB_new_coeff(:,i);
                end
                if i >= 10
                    LW = 2;
                else
                    LW = 1;
                end

                semilogx(f,HdB_here ,'k','LineWidth',LW); hold on, grid on
                plot(f,HdB_here2,'r--','LineWidth',2);

                switch j
                    case 1
                        title('MTFs: dau1997')
                    case 2
                        title('MTFs: osses2021 (as used by Osses & Kohlrausch, 2021)')
                    case 3
                        title('MTFs: osses2021-att-gain (as suggested for further development)')
                end
                if j == 2 || j == 3
                    plot(f,hlpf_dB,':','LineWidth',2,'Color',0.5*[1 1 1]); hold on
                end
            end
            h(end+1) = gcf;
            ylim([-53 3]);
            xlim([.9 2200])
            set(gca,'XTick',XL)
            set(gca,'XTickLabel',XLT)
            xlabel('Frequency (Hz)')
            ylabel({'','Amplitude (dB)'})


            ylim([-20 6])
            set(gca,'YTick',-18:3:3);
        end


        set(gcf,'Position',[0 0 570 700]);

        legend('from insig','from filter coeff.','150-Hz LPF');
    end
    
    data.figure_flag = 'do_fig14';
    data.figure_description = 'Modulation transfer function of the modulation filters using three configurations';

end


