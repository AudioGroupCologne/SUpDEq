function data = exp_osses2022(varargin)
%EXP_OSSES2022 comparison of eight monaural models
%
%   Usage: data = exp_osses2022(flags)
%          data = exp_osses2022(..., 'models', models)
%
%   exp_osses2022 reproduces the figures from Osses et al. (2022), where
%   the following models are compared:
%
%     'dau1997'           : Model #1 from Dau et al., (1997)
%
%     'zilany2014'        : Model #2 from Zilany et al., (2014)
%
%     'verhulst2015'      : Model #3 from Verhulst et al., (2015)
%
%     'verhulst2018'      : Model #4 from Verhulst et al., (2018)
%
%     'bruce2018'         : Model #5 from Bruce et al., (2018)
%
%     'king2019'          : Model #6 from King et al., (2019)
%
%     'relanoiborra2019'  : Model #7 from Relano-Iborra et al., (2019)
%
%     'osses2021'         : Model #8 from Osses and Kohlrausch (2021)
%
%   The following flags can be specified;
%
%     'redo'    Recomputes data for specified figure
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'no_plot' Doesn't plot, only return data.
%
%     'models'  Vector selecting the models, if not all models to process.
%               For example 'models',[1 3 6] shows the data for Dau et al.,
%               (1997), Verhulst et al., (2015), and King et al., (2019) only.
%
%
%   To display Fig. 3 of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig3');
%
%   To display Fig. 4 of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig4');
%
%   To display Fig. 5 of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig5');
%
%   To display Fig. 6 of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig6');
%
%   To display Figs. 7 of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig7'); % this flag will obtain both figures, 6 and 7
%
%   To display Figs. 8 of Osses et al., (2022) use :
%     out = exp_osses2022('fig8'); % this flag will obtain both figures, 6 and 7
%
%   To display Fig. 9 of Osses et al., (2022) use : 
%
%     out = exp_osses2022('fig9');
%
%   To display Fig. 10 of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig10');
%
%   To display Fig. 11 of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig11');
%
%   To display Fig. 12a or 12b of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig12a');
%     out = exp_osses2022('fig12b');
%
%   To display Fig. 13 of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig13');
%
%   To display Fig. 14a or 14b of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig14a');
%     out = exp_osses2022('fig14b');
%
%   To display Fig. 15 of Osses et al., (2022) use :
%
%     out = exp_osses2022('fig15');
%     out = exp_osses2022('fig15','redo'); % example where the data are recomputed
%
%   References:
%     H. Relaño-Iborra, J. Zaar, and T. Dau. A speech-based computational
%     auditory signal processing and perception model. J. Acoust. Soc. Am.,
%     146(5), 2019.
%     
%     A. Osses Vecchi and A. Kohlrausch. Perceptual similarity between piano
%     notes: Simulations with a template-based perception model. J. Acoust.
%     Soc. Am., 141(4), 2020.
%     
%     A. King, L. Varnet, and C. Lorenzi. Accounting for masking of frequency
%     modulation by amplitude modulation with the modulation filter-bank
%     concept. J. Acoust. Soc. Am., 145(2277), 2019.
%     
%     I. C. Bruce, Y. Erfani, and M. S. R. Zilany. A phenomenological model
%     of the synapse between the inner hair cell and auditory nerve:
%     Implications of limited neurotransmitter release sites. Hearing
%     Research, 360:40--54, 2018.
%     
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band
%     carriers. J. Acoust. Soc. Am., 102:2892--2905, 1997a.
%     
%     M. S. A. Zilany, I. C. Bruce, and L. H. Carney. Updated parameters and
%     expanded simulation options for a model of the auditory periphery. The
%     Journal of the Acoustical Society of America, 135(1):283--286, Jan.
%     2014.
%     
%     S. Verhulst, H. Bharadwaj, G. Mehraei, C. Shera, and
%     B. Shinn-Cunningham. Functional modeling of the human auditory
%     brainstem response to broadband stimulation. jasa, 138(3):1637--1659,
%     2015.
%     
%     S. Verhulst, A. Altoè, and V. Vasilkov. Functional modeling of the
%     human auditory brainstem response to broadband stimulation.
%     hearingresearch, 360:55--75, 2018.
%     
%     A. Osses, S. Verhulst, L. Varnet, T. Dau, L. H. Carney, I. C. Bruce,
%     and P. Majdak. A comparative study of eight human auditory models
%     including brainstem processing. Acta Acustica, 2022.
%     
%
%   See also: dau1997 verhulst2015 verhulst2018 zilany2014 bruce2018 osses2021 king2019 relanoiborra2019
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_osses2022.php


%   #Author: Alejandro Osses (2020-2021): Integration in the AMT
%   #Author: Piotr Majdak (2021): Adaptations for the AMT 1.0
%   #Author: Alejandro Osses (2023): Updating display options, as of AMT 1.3.1

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

data = [];

definput.import={'amt_cache'};
definput.flags.type={'missingflag','fig3','fig4','fig5','fig6','fig7','fig8', ...
                 'fig9','fig10','fig11','fig12a','fig12b','fig13','fig14a','fig14b','fig15'};

definput.flags.plot={'plot','no_plot'};
definput.keyvals.models=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

dispflag = amt_flags;

%%%
FontSize = 13; % FontSize for all figure axes: TODO - apply to all the figures
MediumFontSize = 11;
SmallFontSize = 9;
fs = 100e3; % 100 kHz of sampling rate, this is arbitrary
dBFS = 94; % i.e., amplitude 1 = 1 Pa = 94 dB SPL re 2x10^{-5} Pa

%%% for Figs 3, 4, 5, 6 (confirm for the first figures)
models = {'dau1997','zilany2014','verhulst2015','verhulst2018','bruce2018','king2019','relanoiborra2019','osses2021'};
Colours = {'b',local_rgb('Green'),'r',local_rgb('LightSkyBlue'),local_rgb('Maroon'),'m','k',local_rgb('mediumorchid')}; 
Markers = {'o','s','<','>','v','^','p','d'};
MarkersSize = [10 10 10 10 10 10 10 10];
LineStyle = {'-','-','-.','-','--','-','-','-'};
LineWidth = [2 2 2 2 2 2 2 2];

if ~isempty(keyvals.models)
	warning('Not all models selected');
	idxs=keyvals.models;
	Markers = Markers(idxs);
	MarkersSize = MarkersSize(idxs);
	Colours = Colours(idxs);
	LineStyle = LineStyle(idxs);
	LineWidth = LineWidth(idxs);
	models = models(idxs);
end

N_models = length(models);

figure_handle = []; % multiple figures will be generated
figure_name   = [];

%% ------ FIG 3 Osses et al. 2021 -----------------------------------------
if flags.do_fig3

    % Params for FFT:
    K = fs/2; % for FFT
    % Memory allocation
    h_dB    = nan(K,N_models);
    hmax_dB = nan(1,N_models);
    h0_dB   = nan(K,N_models);
    k2remove = [];
    label2use = [];

    for k = 1:N_models
        switch models{k}
            case {'dau1997','king2019'}
                % No middle ear filter
                me_type = [];
            case {'relanoiborra2019'}
                me_type = 'lopezpoveda2001';
                label2use{end+1} = 'relanoiborra2019, osses2021';
            case 'osses2021'
                % uses the same as relanoiborra (no need to calculate it again)
                me_type = [];
            case 'zilany2014'
                me_type = 'zilany2009';
                label2use{end+1} = 'zilany2014, bruce2018';
            case 'bruce2018'
                me_type = [];
            case 'verhulst2015'
                me_type = 'verhulst2015';
                label2use{end+1} = 'verhulst2015';
            case 'verhulst2018'
                me_type = 'verhulst2018';
                label2use{end+1} = 'verhulst2018';
        end

        if ~isempty(me_type)

            [b,a] = middleearfilter(fs,me_type);

            Nr_cascaded = size(b,1);
            if Nr_cascaded == 1 % no cascaded
                [h,w]=freqz(b,a,K);
            else
                h = ones(K,1);
                for j = 1:Nr_cascaded
                    [htmp,w]=freqz(b(j,:),a(j,:),K);
                    h = h.*htmp;
                end
            end
            f = w/pi*(fs/2); % freq. of the FFT bins in Hz
            h_dB(:,k)   = 20*log10(abs(h));
            hmax_dB(k)  = max(h_dB(:,k));
            h0_dB(:,k)  = h_dB(:,k) - hmax_dB(k);

        else
            k2remove = [k2remove k];
        end
    end

    % Preparing the plots with the four middle ear filters:
    models(k2remove) = [];
    Colours(k2remove) = [];
    Markers(k2remove) = [];
    MarkersSize(k2remove) = [];
    LineStyle(k2remove) = [];
    LineWidth(k2remove) = [];

    h_dB(:,k2remove) = [];
    hmax_dB(k2remove) = [];
    h0_dB(:,k2remove) = [];

    %%%
    b_oear = headphonefilter(fs); % Pralong filter
    htmp=freqz(b_oear,1,K);
    hOuter_ear = 20*log10(abs(htmp));

    %%%
    if flags.do_plot
        pl = [];
        figure;
    end
    for k = 1:length(models)
        if flags.do_plot
            opts_plot = {'Color',Colours{k},'LineWidth',LineWidth(k),'LineStyle',LineStyle{k}};
            pl(end+1) = semilogx(f,h0_dB(:,k),opts_plot{:}); hold on, grid on;

            if strcmp(models{k},'relanoiborra2019')
                hCombined = h0_dB(:,k) + hOuter_ear;
                semilogx(f,hCombined,'--','Color',[0.5 0.5 0.5],'LineWidth',4); hold on, grid on;

                %%% Plotting again
                semilogx(f,h0_dB(:,k),opts_plot{:}); hold on, grid on;
            end
        end
        L1 = [min(f) max(f); -3 -3];
        L2 = [f'; h0_dB(:,k)'];
        P = local_InterX(L1,L2);

        if ~isempty(P)
            f_min03_low(k)  = P(1,1);
            f_min03_high(k) = P(1,end);
        else
            f_min03_low(k)  = nan;
            f_min03_high(k) = nan;
        end

        L1 = [min(f) max(f); -15 -15];
        P = local_InterX(L1,L2);

        if ~isempty(P)
            f_min15_low(k)  = P(1,1);
            f_min15_high(k) = P(1,end);
        else
            f_min15_low(k)  = nan;
            f_min15_high(k) = nan;
        end

    end

    if flags.do_plot
        xlim([80 30000]);
        hl = legend([pl(2:end) pl(1)],label2use);
        set(hl ,'FontSize',SmallFontSize);
        set(gca,'FontSize',FontSize);

        f2tick = [63 125 250 500 1000 2000 4000 8000 16000];
        set(gca,'XTick',f2tick);
        set(gca,'YTick',-21:3:3);
        ylim([-19.5 4.5])
        xlabel('Frequency (Hz)');
        ylabel('Amplitude (dB re. max)');

        Pos = get(gcf,'Position');
        Pos(3:4) = [640 330]; % fix width
        set(gcf,'Position',Pos);

        data.figure_flag   = 'do_fig3';
        data.figure_handle = gcf;
        data.figure_name   = 'fig03-me-resp'; % middle ear response
    end
    
    data.description = 'Frequency response of the middle ears used in the selected models';
    data.hmax_dB = hmax_dB;
    data.f_min03_low  = f_min03_low;
    data.f_min03_high = f_min03_high;
    data.f_min15_low  = f_min15_low;
    data.f_min15_high = f_min15_high;
    data.evaluated_models = models;
end

%% ------ FIG 4 Excitation patterns (Fig. 5) / FFT resp. (Fig. 6)----------
if flags.do_fig4
    
    % 1. Top panels, based on data from Fig. 6: ---------------------------
    data5 = exp_osses2022('fig5','no_plot','models',keyvals.models); % data from Fig. 5
    vref = data5.v_ref;
    
    idx_models = 1:length(models);
    idx2remove = find(strcmp(models,'bruce2018'));
    idx_models(idx2remove) = [];
    
    lvls = data5.lvls;
    lvls2look = [40 70 100];
    for i = 1:length(lvls2look)
        idxs2look(i) = find(lvls==lvls2look(i));
    end
    models = data5.models;
    
    K = length(models);
    
    if flags.do_plot
        Markers = {'s','s','s','s','s','s','s','s'}; 
        if length(idx_models) ~= length(Markers) % it will always be the case (bruce2018 is not processed)
            % if true, then the user requested a subset of the models
            Markers = Markers(idx_models);
            Colours = Colours(idx_models);
        end

        num_rows = 4;
        num_cols = 8;
        YL = [-95 12];

        ha1 = [];
        ha2 = [];

        str_extra = {'a)','b)','c)','d)','e)','f)','g)'};
        for k = 1:K
            models_str{k} = [str_extra{k} ' ' models{k}];
        end
        format_str = {'','-','-'};
        
        Colour_low1  = [0.8 0.8 0.8];
        Colour_high1 = [1 1 0]; % rgb('Yellow')
        Colour_low2  = [0.6 0.6 0.6];
        Colour_high2 = [1 0.6445 0]; % rgb('Orange');
        
        Ylabel = 'Output level (dB)';
        Xlabel = 'Frequency (Hz)';
    end
    
    flow1  = audtofreq(freqtoaud(500)-1);
    fhigh1 = audtofreq(freqtoaud(500)+1);
    flow2  = audtofreq(freqtoaud(4000)-1);
    fhigh2 = audtofreq(freqtoaud(4000)+1);
        
    fm1ERB = flow1; % 425;
    fp1ERB = fhigh1; % 583;
    
    fm1ERB = freqtoaud(fm1ERB); % converting from Hz to ERB_N
    fp1ERB = freqtoaud(fp1ERB); % converting from Hz to ERB_N
    
    i = 2; % i = 2 is the index of the tested 500 Hz tone
    if flags.do_plot
        figure;
    end
    
    for idx_lvl = idxs2look 
        for k = 1:K
            fc_resp = data5.data_fc_resp{k,i};
            fc_resp = freqtoaud(fc_resp); % converting from Hz to ERB_N
            
            freq_resp = data5.data_out_freq_resp_all{k,i}(:,idx_lvl);
            curve = freq_resp - vref(k);
            
            if flags.do_plot
                subplot(num_rows,num_cols,k)
                plot(fc_resp,curve,format_str{i},'Color',Colours{k},'LineWidth',2); hold on; grid on;
            end
                        
            L1 = [fc_resp(:) curve(:)]';
            L2 = [fm1ERB*[1 1]; min(curve)-3 max(curve)+3];
            P = local_InterX(L1,L2);
            if dispflag.do_verbose
                amt_disp(sprintf('%s: freq=%.0f, amp=%.1f dB re max',models{k},fm1ERB,P(2)));
            end
            
            if flags.do_plot
                plot(P(1),P(2),'Color',Colours{k},'Marker',Markers{k},'MarkerFaceColor',Colour_low1,'LineWidth',1,'MarkerSize',6);
            end
            
            L2 = [fp1ERB*[1 1]; min(curve)-3 max(curve)+3];
            P = local_InterX(L1,L2);
            if dispflag.do_verbose
              amt_disp(sprintf('%s: freq=%.0f, amp=%.1f dB re max',models{k},fp1ERB,P(2)));
            end
            
            if flags.do_plot
                plot(P(1),P(2),'Color',Colours{k},'Marker',Markers{k},'MarkerFaceColor',Colour_high1,'LineWidth',1,'MarkerSize',6);
                            
                title(models{k});
                ylabel(Ylabel);
            
                if idx_lvl == 5
                    ha1(end+1) = gca;
                    % plot(flow1*[1 1],YL,'k-');
                    % plot(fhigh1*[1 1],YL,'k-');
                end
            end
        end
    end
    if flags.do_plot
        subplot(num_rows,num_cols,8)
        semilogx(fc_resp,nan(size(freq_resp))); 
        ha1(end+1) = gca;
        ylabel(Ylabel);
    end
    
    i = 3; % i = 3 is the index of the tested 4000 Hz tone
    fm1ERB = flow2; % 3567;
    fp1ERB = fhigh2; % 4482;
    
    fm1ERB = freqtoaud(fm1ERB); % converting from Hz to ERB_N
    fp1ERB = freqtoaud(fp1ERB); % converting from Hz to ERB_N
    
    for idx_lvl = idxs2look
        for k = 1:K
            fc_resp = data5.data_fc_resp{k,i};
            fc_resp = freqtoaud(fc_resp);
            
            freq_resp = data5.data_out_freq_resp_all{k,i}(:,idx_lvl);
            curve = freq_resp - vref(k);
            
            if flags.do_plot
                subplot(num_rows,num_cols,k + num_cols)
                plot(fc_resp,curve,format_str{i},'Color',Colours{k},'LineWidth',2); hold on; grid on;
            end
            
            L1 = [fc_resp(:) curve(:)]';
            L2 = [fm1ERB*[1 1]; min(curve)-3 max(curve)+3];
            P = local_InterX(L1,L2);
            if dispflag.do_verbose
                amt_disp(sprintf('%s: freq=%.0f, amp=%.1f dB re max',models{k},fm1ERB,P(2)));
            end
            
            if flags.do_plot
                plot(P(1),P(2),'Color',Colours{k},'Marker',Markers{k},'MarkerFaceColor',Colour_low2,'LineWidth',1,'MarkerSize',6);
            end
            
            L2 = [fp1ERB*[1 1]; min(curve)-3 max(curve)+3];
            P = local_InterX(L1,L2);
            if dispflag.do_verbose
                amt_disp(sprintf('%s: freq=%.0f, amp=%.1f dB re max',models{k},fp1ERB,P(2)));
            end
            
            if flags.do_plot
                plot(P(1),P(2),'Color',Colours{k},'Marker',Markers{k},'MarkerFaceColor',Colour_high2,'LineWidth',1,'MarkerSize',6);
            
                if k == 4
                    xlabel(Xlabel);
                end
                ylabel(Ylabel)
            
                if idx_lvl == 5
                    ha1(end+1) = gca;
                    % plot(flow2*[1 1],YL,'k-');
                    % plot(fhigh2*[1 1],YL,'k-');
                end
            end
        end
    end
    if flags.do_plot
        subplot(num_rows,num_cols,16)
        semilogx(fc_resp,nan(size(freq_resp))); 
        ha1(end+1) = gca;
        ylabel(Ylabel);
    end
    
    data6 = exp_osses2022('fig6','no_plot','models',keyvals.models); % data from Fig. 6
    
    %%% Remove 'data_fc_resp'
    %           data_out_freq_resp_basef
    lvls = [40 70 100]; % These were the tested levels (a priori knowledge)
    fc2look = [500 4000];
    
    idx_f0 = 1; % maybe remove
    for ii = 1:length(fc2look)
        for k = 1:K
            P1 = [];
            P2 = []; 
            val = [];
            idx = [];
            for idx_lvl = 1:3
                fc_resp = data6.data_fig4.ff{k};
                switch fc2look(ii)
                    case 500
                        freq_resp_MA = data6.data_fig4.freq_resp_MA0500{k}(:,idx_lvl);
                        freq_resp    = data6.data_fig4.freq_resp0500{k}(:,idx_lvl);
                    case 4000
                        freq_resp_MA = data6.data_fig4.freq_resp_MA4000{k}(:,idx_lvl);
                        freq_resp    = data6.data_fig4.freq_resp4000{k}(:,idx_lvl);
                    otherwise
                        error('Only tested using f=500 and 4000 Hz');
                end
                
                fc_resp = freqtoaud(fc_resp);

                % freq_resp = data5.data_out_freq_resp_all{k,idx_lvl};
                curve_broad = freq_resp - vref(k);
                if flags.do_plot
                    subplot(num_rows,num_cols,k + num_cols*(ii+1)) % ii+1 = 3 for 500 Hz and ii+1 = 4 for 4000 Hz 
                                                                   % i.e., third and fourth row of plots
                    plot(fc_resp,curve_broad,'-','Color',.9*[1 1 1],'LineWidth',3.5); hold on; grid on;
                end
                
                fc_resp_basef = fc_resp;
                freq_resp_basef = freq_resp_MA; % data5.data_out_freq_resp_basef{k,idx_lvl};

                curve = freq_resp_basef(:,idx_f0) - vref(k);
                if flags.do_plot
                    plot(fc_resp_basef,curve,'-','Color',Colours{k},'LineWidth',2); hold on;
                end

                L1 = [fc_resp_basef(:) curve(:)]';
                L2 = [min(fc_resp) max(fc_resp); max(curve)-3 max(curve)-3];
                P = local_InterX(L1,L2);

                pp = audtofreq(P(1,:));
                if dispflag.do_verbose
                    fc_here = (pp(2)+pp(1))/2;
                    bw_here = pp(2)-pp(1);
                    amt_disp(sprintf('%s: fc=%.0f Hz, BW=%.1f Hz',models{k},fc_here,bw_here));
                end
                
                if fc2look(ii) == 4000
                    switch models{k}
                        case {'verhulst2015','verhulst2018','relanoiborra2019'}
                        [val(idx_lvl),idx(idx_lvl)] = max(curve);
                    end
                end
                
                %%%
                P1(end+1,:) = P(1,1:2);
                P2(end+1,:) = P(2,1:2);

                if flags.do_plot
                    plot(P(1,:),P(2,:),'k-','LineWidth',1.5); hold on;
                
                    if fc2look(ii) == 4000
                        if k == 4
                            xlabel(Xlabel);
                        end
                    end

                    ylabel(Ylabel);
                    if idx_lvl == 1
                        ha2(end+1) = gca;
                    end
                end
            end

            if flags.do_plot
                if fc2look(ii) == 4000
                    if ~isempty(val)
                        yoff = +5; % off
                        for i = 1:3
                            plot(fc_resp_basef(idx(i)),val(i)+yoff,'kv','MarkerFaceColor','g','MarkerSize',6); hold on
                        end
                    end
                end
                for i = 1:2
                    plot(P1(1,i)*[1 1],[P2(1,i) -120],'k--','LineWidth',1.5);
                    plot(P1(2,i)*[1 1],[P2(2,i) P2(1,i)],'k--','LineWidth',1.5);
                    plot(P1(3,i)*[1 1],[P2(3,i) P2(2,i)],'k--','LineWidth',1.5);
                end
            end
        end
        
        if flags.do_plot
            subplot(num_rows,num_cols,8+num_cols*(ii+1)) % Extra empty panel
            semilogx(fc_resp,nan(size(freq_resp))); 
            ha2(end+1) = gca;
            ylabel(Ylabel);
        end
    end
    
    %%%
    if flags.do_plot
        XL = [250 8000];
        XL = freqtoaud(XL);
        linkaxes([ha1 ha2],'x');
        xlim(XL)

        set([ha1 ha2],'YTick',-120:20:0)
        linkaxes(ha1,'xy');
        set(ha1,'YLim',YL);

        linkaxes(ha2,'xy');
        set(ha2,'YLim',YL-20);

        XT_Hz = [audtofreq(freqtoaud(500)-1)   500 audtofreq(freqtoaud( 500)+1) ...
                 audtofreq(freqtoaud(4000)-1) 4000 audtofreq(freqtoaud(4000)+1)];
        XT = freqtoaud(XT_Hz);
        XTL = {'' '500' '' '' '4000' ''};

        set([ha1 ha2],'XTick',XT)
        set([ha1 ha2],'XTickLabel',XTL)

        Pos = get(gcf,'Position');
        Pos(3:4) = [1325 700];
        set(gcf,'Position',Pos);

        figure_handle = gcf;
        figure_name = 'fig04-freq-responses';

        data.figure_flag   = 'do_fig5';
        data.figure_handle = figure_handle;
        data.figure_name   = figure_name;
    end
    
end

%% ------ FIG 5 Osses et al. 2021 -----------------------------------------
if flags.do_fig5

    k2remove = [];
    N_models_here = N_models; % Bruce2018, is removed later...

    %%% 1. Generating the input signals: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Signal parameters:
    dur=100e-3; % 50 ms
    lvls=0:10:100; % Level of the test signals
    N_lvls = length(lvls);

    f0s = [1000 500 4000];

    dt = 1/fs; % \Delta t in s
    t=0:dt:dur-dt;

    cfs = local_m2hz(401); % default characteristic frequencies from Verhulst2012

    Pos34 = [500 350];

    for i = 1:length(f0s)
        f0_target = f0s(i);

        if f0_target == 1000
            idx_lvls = find(lvls==100); % 0-dB reference
        else
            idx_lvls = 1:N_lvls;
        end

        idx = find(cfs>f0_target,1,'last');
        f0 = cfs(idx);

        f0_lo = audtofreq(freqtoaud(f0_target)-1); % min 1 ERB
        idx_off(1) = find(cfs>f0_lo,1,'last');

        f0_hi = audtofreq(freqtoaud(f0_target)+1); % plus 1 ERB
        idx_off(2) = find(cfs>f0_hi,1,'last');

        if dispflag.do_verbose
            amt_disp(sprintf('actual f0=%.1f Hz using Verhulst''s mapping (bin=%.0f)',f0,idx));
            amt_disp(sprintf('f0 off=%.1f Hz (bin=%.0f)',f0_lo,idx_off(1)));
            amt_disp(sprintf('f0 off=%.1f Hz (bin=%.0f)',f0_hi,idx_off(2)));
        end

        % Basis input signal
        insig = zeros(length(t),N_lvls); % Memory allocation
        insig_orig = sin(2*pi*f0*t(:)); % the same sinusoid will be scaled at different levels

        % Up/down cosine ramp (fixed)
        dur_ramp_ms = 10;
        dur_ramp = round((dur_ramp_ms*1e-3)*fs); % duration ramp in samples

        rp    = ones(size(insig_orig));
        rp(1:dur_ramp) = rampup(dur_ramp);
        rp(end-dur_ramp+1:end) = rampdown(dur_ramp);
        %%%

        % Calibration of the input signals:
        for j = idx_lvls
            insig(:,j) = scaletodbspl(insig_orig,lvls(j),dBFS); % setdbspl
            insig(:,j) = rp.*insig(:,j);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N = length(lvls(idx_lvls));
        vrms        = nan(N,N_models); % for the on-frequency values
        vrms_off_lo = nan(N,N_models); % for the off-frequency values
        vrms_off_hi = nan(N,N_models); % for the off-frequency values

        k_test = 1;

        insig_here = insig(:,idx_lvls);

        for k = 1:N_models_here

            %%% Loading flags and keyvals
            [fg,kv] = local_get_flags(models{k});
            fc_flags   = fg.fc_flags;
            afb_flags  = fg.afb_flags;
            afb_kv     = kv.afb_keyvals;
            % ihc_flags  = fg.ihc_flags;
            noihc_flags  = fg.noihc_flags;
            noan_flags = fg.noan_flags; % No auditory nerve module
            % nomfb_flags = fg.nomfb_flags;

            fname = ['fig05_' models{k} '-f0-' num2str(f0_target) '-Hz'];
            c = amt_cache('get',fname,flags.cachemode);

            if ~isempty(c)
                bRun = 0;
                out_all    = c.out_all;
                out_lo_all = c.out_lo_all;
                out_hi_all = c.out_hi_all;
                
                if isfield(c,'out_freq_resp_all')
                    fc_resp = c.fc_resp;
                    out_freq_resp_all = c.out_freq_resp_all;
                end
                if k == 1
                    fc_resp_ref = fc_resp;
                end
                data_fc_resp{k,i} = fc_resp;
                data_out_freq_resp_all{k,i} = out_freq_resp_all;
            else
                bRun = 1;
                out_all    = [];
                out_lo_all = [];
                out_hi_all = [];
                
                out_freq_resp_all = [];
                
                switch models{k}
                    case {'verhulst2015','verhulst2018'}
                        % Nothing to do
                    case 'dau1997'
                        % Nothing to do
                    otherwise
                        fc_resp = fc_resp_ref;
                end
            end

            if bRun
                flow_here = 80;
                fhigh_here = 8000;
                bwmul_here = 0.125; % extra resolution for dau1997, osses2021
                fc_kv_extra = {'bwmul',bwmul_here};

                amt_disp(['Calculating ' models{k} '...']);
                switch models{k}
                    case 'dau1997'
                        for j = 1:N
                            fc_kv = {'flow',f0,'fhigh',f0,'basef',f0,'dboffset',dBFS};
                            out = dau1997(insig_here(:,j),fs,fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_all(:,end+1) = out;

                            fc_kv = {'flow',f0_lo,'fhigh',f0_lo,'basef',f0_lo,'dboffset',dBFS};
                            out = dau1997(insig_here(:,j),fs,fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_lo_all(:,end+1) = out;

                            fc_kv = {'flow',f0_hi,'fhigh',f0_hi,'basef',f0_hi,'dboffset',dBFS};
                            out_hi = dau1997(insig_here(:,j),fs,fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_hi_all(:,end+1) = out;
                            
                            fc_kv = {'flow',flow_here,'fhigh',fhigh_here,'basef',f0,'dboffset',dBFS};
                            [out,fc_resp] = dau1997(insig_here(:,j),fs,fc_kv{:},fc_kv_extra{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_freq_resp_all(:,end+1) = local_rmsdb(out);
                            fc_resp_ref = fc_resp;
                        end

                    case 'zilany2014'

                        kv={'numH',0,'numM',0,'numL',0,'fiberType',4,'nrep',1};
                        for j = 1:N
                            [~,~,~,c1,c2] = zilany2014(insig_here(:,j),fs,f0,kv{:},afb_flags{:});
                            idxs = 1:length(insig_here); % to remove the zero padding
                            out = c1(:)+c2(:);
                            out = out(idxs);
                            out_all(:,end+1) = out;

                            [~,~,~,c1,c2] = zilany2014(insig_here(:,j),fs,f0_lo,kv{:},afb_flags{:});
                            out = c1(:)+c2(:);
                            out = out(idxs);
                            out_lo_all(:,end+1) = out;

                            [~,~,~,c1,c2] = zilany2014(insig_here(:,j),fs,f0_hi,kv{:},afb_flags{:});
                            out = c1(:)+c2(:);
                            out = out(idxs);
                            out_hi_all(:,end+1) = out;
                            
                            out_tmp = [];
                            for kk = 1:length(fc_resp)
                                fc_here = fc_resp(kk);
                                if fc_here >= 125
                                    [~,~,~,c1,c2] = zilany2014(insig_here(:,j),fs,fc_here, kv{:},afb_flags{:});
                                    out = c1(:)+c2(:);
                                    out = out(idxs);
                                    out_tmp(:,end+1) = out;    
                                else
                                    last_idx_nan = kk;
                                end
                            end
                            out_tmp = local_rmsdb(out_tmp);
                            out_freq_resp_all(:,end+1) = out_tmp(:);
                        end
                        nan_array = nan(last_idx_nan, size(out_freq_resp_all,2));
                        out_freq_resp_all = [nan_array; out_freq_resp_all];
                        
                    case 'verhulst2015'
                        %idx_here = find(strcmp(afb_flags,'no_v'));
                        %if ~isempty(idx_here)
                        %    afb_flags{idx_here} = 'v'; % overwritting 'no_v' by 'v':
                        %end
                        out_model = verhulst2015(insig_here,fs,fc_flags{:},afb_flags{:},'v',noihc_flags{:},noan_flags{:});
                        for j = 1:N
                            out = out_model(j).v(:,idx);
                            out_all(:,end+1) = out;

                            out = out_model(j).v(:,idx_off(1));
                            out_lo_all(:,end+1) = out;

                            out = out_model(j).v(:,idx_off(2));
                            out_hi_all(:,end+1) = out;
                            
                            if j == 1
                                fc_resp = out_model(j).cf;
                            end
                            out_freq_resp_all(:,end+1) = local_rmsdb(out_model(j).v);
                        end

                    case 'verhulst2018'
                        %idx_here = find(strcmp(afb_flags,'no_v'));
                        %if ~isempty(idx_here)
                        %    afb_flags{idx_here} = 'v'; % overwritting 'no_v' by 'v':
                        %end
                        out_model = verhulst2018(insig_here,fs,fc_flags{:},afb_flags{:},'v',noihc_flags{:},noan_flags{:});
                        for j = 1:N
                            out = out_model(j).v(:,idx);
                            out_all(:,end+1) = out;

                            out = out_model(j).v(:,idx_off(1));
                            out_lo_all(:,end+1) = out;

                            out = out_model(j).v(:,idx_off(2));
                            out_hi_all(:,end+1) = out;
                            
                            if j == 1
                                fc_resp = out_model(j).cf;
                            end
                            out_freq_resp_all(:,end+1) = local_rmsdb(out_model(j).v);
                        end

                    case 'bruce2018'
                        % bruce2018 has the same preprocessing as zilany2014
                        if i == 1
                            k2remove = [k2remove k];
                        end

                    case 'king2019'
                        for j = 1:N
                            fc_kv = {'flow',f0,'fhigh',f0,'basef',f0,'dboffset',dBFS};
                            out = king2019(insig_here(:,j),fs,afb_kv{:},fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_all(:,end+1) = out;

                            fc_kv = {'flow',f0_lo,'fhigh',f0_lo,'basef',f0_lo,'dboffset',dBFS};
                            out = king2019(insig_here(:,j),fs,afb_kv{:},fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_lo_all(:,end+1) = out;

                            fc_kv = {'flow',f0_hi,'fhigh',f0_hi,'basef',f0_hi,'dboffset',dBFS};
                            out = king2019(insig_here(:,j),fs,afb_kv{:},fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_hi_all(:,end+1) = out;
                            
                            fc_kv = {'flow',flow_here,'fhigh',fhigh_here,'basef',f0,'dboffset',dBFS};
                            [out,fc_resp] = king2019(insig_here(:,j),fs,afb_kv{:},fc_kv{:},fc_kv_extra{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_freq_resp_all(:,end+1) = local_rmsdb(out);
                        end

                    case 'relanoiborra2019'
                        for j = 1:N
                            fc_kv = {'flow',f0,'fhigh',f0,'basef',f0,'no_ihc','no_an'};
                            [~,~,out] = relanoiborra2019_featureextraction(insig_here(:,j),fs,fc_kv{:},afb_flags{:});
                            out_all(:,end+1) = out;

                            fc_kv = {'flow',f0_lo,'fhigh',f0_lo,'basef',f0_lo,'no_ihc','no_an'};
                            [~,~,out] = relanoiborra2019_featureextraction(insig_here(:,j),fs,fc_kv{:},afb_flags{:});
                            out_lo_all(:,end+1) = out;

                            fc_kv = {'flow',f0_hi,'fhigh',f0_hi,'basef',f0_hi,'no_ihc','no_an'};
                            [~,~,out] = relanoiborra2019_featureextraction(insig_here(:,j),fs,fc_kv{:},afb_flags{:});
                            out_hi_all(:,end+1) = out;
                            
                            fc_kv = {'flow',flow_here,'fhigh',fhigh_here,'basef',f0,'no_ihc','no_an'};
                            [~,~,out,fc_resp] = relanoiborra2019_featureextraction(insig_here(:,j),fs,fc_kv{:},fc_kv_extra{:},afb_flags{:});
                            out_freq_resp_all(:,end+1) = local_rmsdb(out);
                        end

                    case 'osses2021'
                        for j = 1:N
                            fc_kv = {'flow',f0,'fhigh',f0,'basef',f0,'dboffset',dBFS};
                            out = osses2021(insig_here(:,j),fs,fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_all(:,end+1) = out;

                            fc_kv = {'flow',f0_lo,'fhigh',f0_lo,'basef',f0_lo,'dboffset',dBFS};
                            out = osses2021(insig_here(:,j),fs,fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_lo_all(:,end+1) = out;

                            fc_kv = {'flow',f0_hi,'fhigh',f0_hi,'basef',f0_hi,'dboffset',dBFS};
                            out = osses2021(insig_here(:,j),fs,fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            out_hi_all(:,end+1) = out;
                            
                            fc_kv = {'flow',flow_here,'fhigh',fhigh_here,'basef',f0,'dboffset',dBFS};
                            [out,fc_resp] = osses2021(insig_here(:,j),fs,fc_kv{:},fc_kv_extra{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                            
                            out_freq_resp_all(:,end+1) = local_rmsdb(out);
                        end
                end

                if ~isempty(out_all)
                    c = [];
                    c.out_all = out_all;
                    c.out_lo_all = out_lo_all;
                    c.out_hi_all = out_hi_all;
                    
                    c.fc_resp = fc_resp;
                    c.out_freq_resp_all = out_freq_resp_all;
                    amt_cache('set',fname,c);
                    
                    data_fc_resp{k,i} = fc_resp;
                    data_out_freq_resp_all{k,i} = out_freq_resp_all;
                end
            end

            if ~isempty(out_all)
                % Now calculation of RMS values:
                for j = 1:N
                    vrms(j,k)        = local_rmsdb(out_all(:,j));
                    vrms_off_lo(j,k) = local_rmsdb(out_lo_all(:,j));
                    vrms_off_hi(j,k) = local_rmsdb(out_hi_all(:,j));
                end
            end
        end

        if i == 1
            % Preparing the plots with the four middle ear filters:
            models(k2remove) = [];
            Colours(k2remove) = [];
            Markers(k2remove) = [];
            MarkersSize(k2remove) = [];
            LineStyle(k2remove) = [];
            LineWidth(k2remove) = [];

            N_models_here = length(models); % number of models after removing

            vrms(:,k2remove) = [];
            vrms_off_lo(:,k2remove) = [];
            vrms_off_hi(:,k2remove) = [];
            
            data_fc_resp(k2remove) = [];
            data_out_freq_resp_all(k2remove) = [];
        end

        if f0_target == 1000
            v_ref_0_dB = vrms;
            % v_ref_0_dB = [   -7.3561  -44.0567  -77.9298 -101.8625  -44.0567]; % calculated at 1 kHz
            data.v_ref = v_ref_0_dB;
        else
            if f0_target == 500
                prefix4title = {'a) ','b) ','c) '};
            elseif f0_target == 4000
                prefix4title = {'d) ','e) ','f) '};
            end
            v_ref = v_ref_0_dB;

            if flags.do_plot
                Format = [];
                for k = 1:N_models_here
                     Format{k} = {'Color',Colours{k},'LineStyle',LineStyle{k}, ...
                        'LineWidth',LineWidth(k),'Marker',Markers{k},'MarkerFaceColor','w', ...
                        'MarkerSize',MarkersSize(k)};
                end
                %%% End formatting options

                figure;
                for k = 1:N_models_here
                    pl(k)=plot(lvls,vrms(:,k)-v_ref(k),Format{k}{:});
                    grid on, hold on;
                end

                title(sprintf('\n%sCF_n at %.0f Hz, on-freq. (n=%.0f)',prefix4title{1},f0,idx));
                xlabel('Input level (dB SPL)');
                ylabel('Output level (dB)');
                ylim([-103 13])
                xlim([min(lvls)-5 max(lvls)+5])
                set(gca,'XTick', 0:10:100)
                set(gca,'YTick',-100:10:10)
                figure_handle(end+1) = gcf;
                figure_name{end+1}   = sprintf('fig05-IO-at-%.0f-Hz',f0);
                Pos = get(gcf,'Position');
                Pos(3:4) = Pos34;
                set(gcf,'Position',Pos);

                %%%%
                figure;
                for k = 1:N_models_here
                    plot(lvls,vrms_off_lo(:,k)-v_ref(k),Format{k}{:});
                    grid on, hold on
                end
                title(sprintf('\n%sCF_n at %.0f Hz, off-freq. (n=%.0f)',prefix4title{2},f0_lo,idx_off(1)));
                xlabel('Input level (dB SPL)');
                ylabel('Output level (dB)');
                ylim([-103 13])
                xlim([min(lvls)-5 max(lvls)+5])
                set(gca,'XTick', 0:10:100)
                set(gca,'YTick',-100:10:10)
                figure_handle(end+1) = gcf;
                figure_name{end+1}   = sprintf('fig05-IO-at-%.0f-Hz-off-lo',f0);
                Pos = get(gcf,'Position');
                Pos(3:4) = Pos34;
                set(gcf,'Position',Pos);

                figure;
                for k = 1:N_models_here
                    plot(lvls,vrms_off_hi(:,k)-v_ref(k),Format{k}{:});
                    grid on, hold on
                end
                title(sprintf('\n%sCF_n at %.0f Hz, off-freq. (n=%.0f)',prefix4title{3},f0_hi,idx_off(2)));
                xlabel('Input level (dB SPL)');
                ylabel('Output level (dB)');
                ylim([-103 13])
                xlim([min(lvls)-5 max(lvls)+5])
                set(gca,'XTick', 0:10:100)
                set(gca,'YTick',-100:10:10)
                Pos = get(gcf,'Position');
                Pos(3:4) = Pos34;
                set(gcf,'Position',Pos);

                if i == length(f0s)
                    text4leg = models;
                    if strcmp(models{2},'zilany2014')
                        text4leg{2} = [text4leg{2} ',bruce2018'];
                    end
                    hl = legend(text4leg,'Location','SouthEast');
                    set(hl,'FontSize',8);
                end

                figure_handle(end+1) = gcf;
                figure_name{end+1}   = sprintf('fig05-IO-at-%.0f-Hz-off-hi',f0);
            end
            data(k_test).vrms = vrms;
            data(k_test).vrms_off_lo = vrms_off_lo;
            data(k_test).vrms_off_hi = vrms_off_hi;
            data(k_test).f0 = f0;

            k_test = k_test+1;
        end

    end
    if flags.do_plot
        data.figure_flag   = 'do_fig5';
        data.figure_handle = figure_handle;
        data.figure_name   = figure_name;
    end
    data.models = models;

    data.v_ref_0_dB = v_ref_0_dB;
    data.insig = insig;
    data.insig_orig = insig_orig;
    data.fs    = fs;
    data.ramp  = rp;
    data.dBFS  = dBFS;
    data.lvls  = lvls;
    
    data.data_fc_resp = data_fc_resp;
    data.data_out_freq_resp_all = data_out_freq_resp_all;
end

%% ------ FIG 6 Frequency selectivity ------------------------------------
if flags.do_fig6 

    k2remove = [];

    fsig = 4000; % Hz
    basef = fsig;

    dt = 1/fs;
    lvls=[40 70 100];
	N_lvls = length(lvls);

    flow  = 20;
    fhigh = 10000;
    BW = fhigh-flow;
    fc = BW/2;
    dur_stim_total = 3; % s, Total duration of the stimulus

    bCreate_insig = 1;
    fname = 'fig06_dau1997-QERB-40';
    c = amt_cache('get',fname,flags.cachemode);
    if ~isempty(c)
        if isfield(c,'insig_not_scaled')
            insig_not_scaled = c.insig_not_scaled;
            bCreate_insig = 0;
        end
    end

    if bCreate_insig == 1
        insig_not_scaled = sig_bandpassnoise(fc,fs,dur_stim_total,50,BW);
    end

    % Bandpass noise (Stimulus bType = 2 in previous versions of this script)
    dur_each_section = 0.5; % s
    N_sections = dur_stim_total/dur_each_section;

    % Up/down cosine ramp (fixed)
    dur_ramp_ms = 10; % in ms
    dur_ramp = round((dur_ramp_ms*1e-3)*fs); % duration ramp in samples

    N_samples = round(dur_stim_total*fs);
    rp    = ones(N_samples,1);
    rp(1:dur_ramp) = rampup(dur_ramp);
    rp(end-dur_ramp+1:end) = rampdown(dur_ramp);

    for i=1:N_lvls
        lvl_tot = lvls(i); %  + 10*log10(BW);
        insigs(:,i) = rp.*scaletodbspl(insig_not_scaled,lvl_tot,dBFS);
    end
    % End: generating bandpass noise

    cf = local_m2hz(401); % Get characteristic frequencies from Verhulst models
    cf_max = 10000;
    cf_min =  125;
    %%%
    N_sigs = 31; % Approximate number of bands
    bin_cf_max = find(cf>cf_max,1,'last');
    bin_cf_min = find(cf>cf_min,1,'last');
    df_bin = floor(abs(bin_cf_max-bin_cf_min)/N_sigs);

    idx_cf = [bin_cf_min:-df_bin:bin_cf_max 305 112]; % 305 - bin for 502 Hz, 112 bin for 4013 Hz
    idx_cf = sort(idx_cf,'descend');
    fc_ref = cf(idx_cf);
    
    N_sigs = length(fc_ref); % Exact number of bands

    for k = 1:N_models
        %%% Loading flags and keyvals
        [fg,kv] = local_get_flags(models{k});
        fc_flags   = fg.fc_flags;
        afb_flags  = fg.afb_flags;
        afb_kv     = kv.afb_keyvals;
        noihc_flags  = fg.noihc_flags;
        noan_flags = fg.noan_flags; % No auditory nerve module
        nomfb_flags = fg.nomfb_flags;

        for j = 1:N_lvls
            fname = ['fig06_' models{k} '-QERB-' num2str(lvls(j))];
            c = amt_cache('get',fname,flags.cachemode);

            if ~isempty(c)
                bRun = 0;

                outsig = c.outsig;
                fc     = c.fc;
                if isfield(c,'insig_not_scaled')
                    % Only for dau1997:
                    insig_not_scaled = c.insig_not_scaled;
                end
            else
                bRun = 1;
                outsig = [];
            end

            if bRun
                amt_disp(['Calculating ' models{k} '...']);
                switch models{k}
                    case 'dau1997'
                        for ii = 1:N_sigs
                            fc_kv = {'flow',fc_ref(ii),'fhigh',fc_ref(ii),'basef',fc_ref(ii),'dboffset',dBFS};
                            [outsig(:,ii),fc(ii)] = dau1997(insigs(:,j),fs,fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                        end

                    case 'zilany2014'
                        outsig = [];
                        kv={'nrep',1}; % kv={'numH',0,'numM',0,'numL',0,'fiberType',4,'nrep',1,'reptime',1.2};
                        for ii = 1:N_sigs
                            [~,~,~,c1,c2] = zilany2014(insigs(:,j),fs,fc_ref(ii),kv{:},afb_flags{:});
                            if ii == 1
                                L = size(insigs,1);
                            end
                            outsig(:,ii) = c1(1:L)+c2(1:L);
                        end

                    case 'verhulst2015'
                        %idx = find(strcmp(afb_flags,'no_v'));
                        %if ~isempty(idx)
                        %    afb_flags{idx} = 'v'; % overwritting 'no_v' by 'v':
                        %end
                        outs = verhulst2015(insigs(:,j),fs,fc_flags{:},afb_flags{:},'v',noihc_flags{:},noan_flags{:});
                        if j == 1
                            CF_here = outs.cf;
                        end
                        outsig = [];
                        for ii = 1:N_sigs
                            outsig(:,ii) = outs.v(:,idx_cf(ii));
                        end

                    case 'verhulst2018'
                        %idx = find(strcmp(afb_flags,'no_v'));
                        %if ~isempty(idx)
                        %    afb_flags{idx} = 'v'; % overwritting 'no_v' by 'v':
                               % For most of the simulations in this script 'v' provides redundant (big-sized) information
                               %   for this reason the 'no_v' is only for this figure overwritten to 'v', because BM velocity
                               %   is what we actually need to obtain this figure.
                        %end
                        outs = verhulst2018(insigs(:,j),fs,fc_flags{:},afb_flags{:},'v',noihc_flags{:},noan_flags{:});
                        if j == 1
                            CF_here = outs.cf;
                        end
                        outsig = [];
                        for ii = 1:N_sigs
                            outsig(:,ii) = outs.v(:,idx_cf(ii));
                        end

                    case 'bruce2018'
                        % already in zilany2014
                        if j == 1
                            k2remove = [k2remove k];
                        end
                        outsig = [];
                        fc = [];

                    case 'king2019'
                        for ii = 1:N_sigs
                            fc_kv = {'flow',fc_ref(ii),'fhigh',fc_ref(ii),'basef',fc_ref(ii),'dboffset',dBFS};
                            [outsig(:,ii),fc(ii)] = king2019(insigs(:,j),fs,fc_kv{:},afb_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                        end

                    case 'relanoiborra2019'
                        for ii = 1:N_sigs
                            fc_kv = {'flow',fc_ref(ii),'fhigh',fc_ref(ii),'basef',fc_ref(ii),'erbspacebw','no_internalnoise','no_ihc','no_an'};
                            [~,~,outsig(:,ii),fc(ii)] = relanoiborra2019_featureextraction(insigs(:,j), fs,fc_kv{:});
                        end

                    case 'osses2021'
                        for ii = 1:N_sigs
                            fc_kv = {'flow',fc_ref(ii),'fhigh',fc_ref(ii),'basef',fc_ref(ii),'dboffset',dBFS};
                            [outsig(:,ii),fc(ii)] = osses2021(insigs(:,j),fs,fc_kv{:},afb_flags{:},noihc_flags{:},noan_flags{:});
                        end
                end
                if ~isempty(outsig)
                    c = [];
                    if strcmp(models{k},'dau1997')
                        c.insig_not_scaled = insig_not_scaled;
                    end
                    c.outsig = outsig;
                    c.fc     = fc;
                    amt_cache('set',fname,c);
                end
            end

            K = N_samples/2;

            if ~isempty(outsig)
                
              idx0500 = find(fc_ref >  500, 1, 'first'); % needed for Fig. 4
              idx4000 = find(fc_ref > 4000, 1, 'first'); % needed for Fig. 4
                
              switch j
                case 1 % Low intensity
                  if dispflag.do_verbose
                      amt_disp(sprintf('%s: Calculating low-level tuning',models{k}));
                  end
                  for ii = 1:length(fc_ref)
                    N_here = length(outsig(:,ii));
                    insig_here = reshape(outsig(:,ii),N_here/2,2);
                    [BW_low_each(ii,:,k), QERB_low_each(ii,:,k), ...
                     Q03_low_each(ii,:,k), Q10_low_each(ii,:,k), extra] ...
                     = local_Get_ERB_estimation_multi(outsig(:,ii), fc_ref(ii), fs, N_sections);
                 
                    if ii == idx0500
                        data_fig4.ff{k} = extra.freq_here;
                        data_fig4.freq_resp0500{k}(:,j) = extra.Amp_dB    + extra.max_dB;
                        data_fig4.freq_resp_MA0500{k}(:,j) = extra.Amp_dB_MA + extra.max_dB + extra.max_MA;
                    end
                    if ii == idx4000
                        data_fig4.freq_resp4000{k}(:,j) = extra.Amp_dB    + extra.max_dB;
                        data_fig4.freq_resp_MA4000{k}(:,j) = extra.Amp_dB_MA + extra.max_dB + extra.max_MA;
                    end 
                  end

                case 2 % Middle intensity
                  if dispflag.do_verbose
                      amt_disp(sprintf('%s: Calculating mid-level tuning',models{k}));
                  end
                  for ii = 1:length(fc_ref)
                    N_here = length(outsig(:,ii));
                    insig_here = reshape(outsig(:,ii),N_here/2,2);
                    [BW_mid_each(ii,1,k), QERB_mid_each(ii,1,k), ...
                    Q03_mid_each(ii,1,k), Q10_mid_each(ii,1,k), extra] ...
                    = local_Get_ERB_estimation_multi(outsig(:,ii), fc_ref(ii), fs, N_sections);
                    
                    if ii == idx0500
                        data_fig4.ff{k} = extra.freq_here;
                        data_fig4.freq_resp0500{k}(:,j) = extra.Amp_dB    + extra.max_dB;
                        data_fig4.freq_resp_MA0500{k}(:,j) = extra.Amp_dB_MA + extra.max_dB + extra.max_MA;
                    end
                    if ii == idx4000
                        data_fig4.freq_resp4000{k}(:,j) = extra.Amp_dB    + extra.max_dB;
                        data_fig4.freq_resp_MA4000{k}(:,j) = extra.Amp_dB_MA + extra.max_dB + extra.max_MA;
                    end 
                  end

                case 3 % Higher intensity
                  if dispflag.do_verbose
                      amt_disp(sprintf('%s: Calculating high-level tuning',models{k}));
                  end
                  for ii = 1:length(fc_ref)
                    N_here = length(outsig(:,ii));
                    insig_here = reshape(outsig(:,ii),N_here/2,2);
                    [BW_high_each(ii,1,k), QERB_high_each(ii,1,k), ...
                     Q03_high_each(ii,1,k), Q10_high_each(ii,1,k), extra] ...
                     = local_Get_ERB_estimation_multi(outsig(:,ii), fc_ref(ii), fs, N_sections);
                 
                    if ii == idx0500
                        data_fig4.ff{k} = extra.freq_here;
                        data_fig4.freq_resp0500{k}(:,j) = extra.Amp_dB    + extra.max_dB;
                        data_fig4.freq_resp_MA0500{k}(:,j) = extra.Amp_dB_MA + extra.max_dB + extra.max_MA;
                    end
                    if ii == idx4000
                        data_fig4.freq_resp4000{k}(:,j) = extra.Amp_dB    + extra.max_dB;
                        data_fig4.freq_resp_MA4000{k}(:,j) = extra.Amp_dB_MA + extra.max_dB + extra.max_MA;
                    end 
                  end % end ii
              end % end switch j
            end % end ~isempty
            
        end % end for j
    end % end for models

    Pos34 = [600 350];
    %%%
    % Preparing the plots
    models(k2remove) = [];
    Colours(k2remove) = [];
    Markers(k2remove) = [];
    MarkersSize(k2remove) = [];
    LineStyle(k2remove) = [];
    LineWidth(k2remove) = [];

    N_models_here = length(models); % number of models after removing

    BW_low_each(:,:,k2remove) = [];
    QERB_low_each(:,:,k2remove) = [];
    Q03_low_each(:,:,k2remove) = [];
    Q10_low_each(:,:,k2remove) = [];

    BW_mid_each(:,:,k2remove) = [];
    QERB_mid_each(:,:,k2remove) = [];
    Q03_mid_each(:,:,k2remove) = [];
    Q10_mid_each(:,:,k2remove) = [];

    BW_high_each(:,:,k2remove) = [];
    QERB_high_each(:,:,k2remove) = [];
    Q03_high_each(:,:,k2remove) = [];
    Q10_high_each(:,:,k2remove) = [];

    idx = find(fc_ref<= 8000);
    fc_ref = fc_ref(idx);

    QERB_low_each = QERB_low_each(idx,:,:);
    Q03_low_each = Q03_low_each(idx,:,:);
    Q10_low_each = Q10_low_each(idx,:,:);

    QERB_mid_each = QERB_mid_each(idx,:,:);
    Q03_mid_each = Q03_mid_each(idx,:,:);
    Q10_mid_each = Q10_mid_each(idx,:,:);

    QERB_high_each = QERB_high_each(idx,:,:);
    Q03_high_each = Q03_high_each(idx,:,:);
    Q10_high_each = Q10_high_each(idx,:,:);

    data_fig4.ff(k2remove) = [];
    data_fig4.freq_resp0500(k2remove) = [];
    data_fig4.freq_resp_MA0500(k2remove) = [];
    data_fig4.freq_resp4000(k2remove) = [];
    data_fig4.freq_resp_MA4000(k2remove) = [];
    %%% ---------------------------------------------------------------
    if flags.do_plot
        figure;

        % Analytical expressions:
        fc_ref_here = [125 fc_ref 10000];
        alpha = 0.3;  % Shera2002, Table I, col 'human'
        beta  = 12.7; % Shera2002, Table I, col 'human'
        Qerb_Shera_extended = beta*((fc_ref_here/1000).^alpha);
        Qerb_Shera = Qerb_Shera_extended(2:end-1);
        % Q10 = Qerb*0.505+0.2085; % Ibrahim2010, Eq. 40.2
        pl = [];
        pl(end+1) = semilogx(fc_ref_here,Qerb_Shera_extended,'-','Color',local_rgb('LightGray'),'LineWidth',5); hold on

        BW_Glasberg = 24.7*(4.37*(fc_ref_here/1000)+1);
        Qerb_Glasberg_extended = fc_ref_here./BW_Glasberg;
        Qerb_Glasberg = Qerb_Glasberg_extended(2:end-1);
        % Q10   = Qerb*0.505+0.2085; % Ibrahim2010, Eq. 40.2
        pl(end+1) = plot(fc_ref_here,Qerb_Glasberg_extended,'-','Color',local_rgb('Gray'),'LineWidth',5);

        for k = 1:N_models_here % little adjustments in format for the coming figures
            switch models{k}
                case 'dau1997'
                    Marker = [Markers{k} '-'];
                    MSize = 7;
                    LW = 1;
                case 'zilany2014'
                    Marker = [Markers{k} LineStyle{k}];
                    MSize = 7;
                    LW = LineWidth(k);
                otherwise
                    Marker = LineStyle{k};
                    MSize = MarkersSize(k);
                    LW = LineWidth(k);
            end
            format_kv = {Marker,'Color',Colours{k},'MarkerSize',MSize,'LineWidth',LW};

            if k == 1
                disp('Using Q03')
            end
            Q_low_here(:,k) = mean(Q03_low_each(:,:,k),2); % average of the positive and negative click
            pl(end+1) = semilogx(fc_ref,Q_low_here(:,k),format_kv{:},'MarkerFaceColor',Colours{k}); hold on
            % plot(fc_ref,QERB_higher(:,k),'s--',format_kv{:},'MarkerFaceColor','w');
        end

        grid on;
        XT = [20 125 250 500 1000 2000 4000 8000 16000];
        set(gca,'XTick',XT);
        YT = 2:2:26;
        set(gca,'YTick',YT);
        ylim([1 27])

        xlabel('Frequency (Hz)');
        ylabel('Q factor');

        Pos = get(gcf,'Position');
        Pos(3:4) = Pos34; % [600 350];
        set(gcf,'Position',Pos);

        figure_handle(end+1) = gcf;
        figure_name{end+1}   = sprintf('fig06-QERB-low-level');

        % ---------------------------------------------------------------------
        for k = 1:N_models_here

            X_ref_empirical = mean(QERB_low_each(:,:,k),2);
            Y_10 = Q10_low_each(:,1,k);
            Y_03 = Q03_low_each(:,1,k);

            % Fit a polynomial p of degree 1 to the (x,y) data:
            p10(k,:) = polyfit(X_ref_empirical,Y_10,1);
            p03(k,:) = polyfit(X_ref_empirical,Y_03,1);

            switch models{k}
                case {'dau1997','relanoiborra2019','king2019','osses2021'}
                    X_ref(:,k) = Qerb_Glasberg(:);
                case {'zilany2014','verhulst2015','verhulst2018'}
                    X_ref(:,k) = Qerb_Shera(:);
            end
            p10_with_formula(k,:) = polyfit(X_ref(:,k),Y_10,1);
            p03_with_formula(k,:) = polyfit(X_ref(:,k),Y_03,1); % figure; plot(X_ref,Y_10); hold on; plot(X_ref,Y_03,'r')

            Q03_recons(:,k) = X_ref(:,k)*p03_with_formula(k,1) + p03_with_formula(k,2);
        end

        for k = 1:N_models_here
            f_current = 125;
            f_high_lim = 8000;

            count(k) = 0;
            while f_current < f_high_lim
                count(k) = count(k)+1;
                f_low(count(k),k) = f_current;
                if count(k)>10000; warning('Break'); break; end % break the loop if got stuck
                Q_here = interp1(fc_ref, mean(Q03_low_each(:,:,k),2), f_low(count(k),k),'linear','extrap');

                BW_here = f_low(count(k),k)/Q_here;
                f_high(count(k),k) = f_low(count(k),k) + BW_here;
                f_current = f_high(count(k),k);
            end

            fc_max = max(f_high);
            fc_min = f_low(1,:);
        end
        % Values reported in Table II:
        counts   = count;
        erb_step = (freqtoaud(fc_max)-freqtoaud(fc_min))./(count-1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure; % Now middle levels only

        % Analytical expressions:
        pl = [];
        pl(end+1) = semilogx(fc_ref_here,Qerb_Shera_extended   ,'-','Color',local_rgb('LightGray'),'LineWidth',5); hold on
        pl(end+1) = plot(    fc_ref_here,Qerb_Glasberg_extended,'-','Color',local_rgb('Gray')     ,'LineWidth',5);

        for k = 1:N_models_here
            switch models{k}
                case 'dau1997'
                    Marker = [Markers{k} '-'];
                    MSize = 7;
                    LW = 1;
                case 'zilany2014'
                    Marker = [Markers{k} LineStyle{k}];
                    MSize = 7;
                    LW = LineWidth(k);
                otherwise
                    Marker = LineStyle{k};
                    MSize = MarkersSize(k);
                    LW = LineWidth(k);
            end
            format_kv = {Marker,'Color',Colours{k},'MarkerSize',MSize, ...
                'LineWidth',LW};
            if k == 1 % bUse_Q03
                disp('Using Q03')
            end
            Q_mid_here(:,k) = mean(Q03_mid_each(:,:,k),2); % average of the positive and negative click
            pl(end+1) = semilogx(fc_ref,Q_mid_here(:,k),format_kv{:},'MarkerFaceColor',Colours{k}); hold on
            % plot(fc_ref,QERB_higher(:,k),'s--',format_kv{:},'MarkerFaceColor','w');
        end
        grid on;
        set(gca,'XTick',XT);
        YT = 2:2:26;
        set(gca,'YTick',YT);

        ylim([1 27])

        leg = models;
        k = 2;
        if strcmp(leg{k},'zilany2014')
            leg{k} = 'zilany2014,bruce2018';
        end
        leg{end+1} = 'Shera (2001)';
        leg{end+1} = 'G & M (1990)';

        legend([pl(3:end) pl(1:2)],leg,'Location','NorthWest');

        xlabel('Frequency (Hz)');
        ylabel('Q factor');

        Pos = get(gcf,'Position');
        Pos(3:4) = Pos34;
        set(gcf,'Position',Pos);

        figure_handle(end+1) = gcf;
        figure_name{end+1}   = sprintf('fig06-QERB-mid-level');

        f_low = [];
        f_high = [];
        count = [];
        for k = 1:N_models_here
            f_current = 125;
            f_high_lim = 8000;

            count(k) = 0;
            while f_current < f_high_lim
                count(k) = count(k)+1;
                f_low(count(k),k) = f_current;

                Q_here = interp1(fc_ref, mean(Q03_mid_each(:,:,k),2), f_low(count(k),k),'linear','extrap');
                if isnan(Q_here)
                    % Q_here = Q03_recons_hi(1,k);
                    warning('Warning what happens here')
                end

                BW_here = f_low(count(k),k)/Q_here;
                f_high(count(k),k) = f_low(count(k),k) + BW_here;

                f_current = f_high(count(k),k);
            end

            fc_max = max(f_high);
            fc_min = f_low(1,:);
        end
        counts(2,:)   = count;
        erb_step(2,:) = (freqtoaud(fc_max)-freqtoaud(fc_min))./(count-1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure; % Now high levels only

        % Analytical expressions:
        pl = [];
        pl(end+1) = semilogx(fc_ref_here,Qerb_Shera_extended   ,'-','Color',local_rgb('LightGray'),'LineWidth',5); hold on
        pl(end+1) = plot(    fc_ref_here,Qerb_Glasberg_extended,'-','Color',local_rgb('Gray')     ,'LineWidth',5);

        for k = 1:N_models_here
            switch models{k}
                case 'dau1997'
                    Marker = [Markers{k} '-'];
                    MSize = 7;
                    LW = 1;
                case 'zilany2014'
                    Marker = [Markers{k} LineStyle{k}];
                    MSize = 7;
                    LW = LineWidth(k);
                otherwise
                    Marker = LineStyle{k};
                    MSize = MarkersSize(k);
                    LW = LineWidth(k);
            end
            format_kv = {Marker,'Color',Colours{k},'MarkerSize',MSize, ...
                'LineWidth',LW};
            if k == 1
                disp('Using Q03')
            end
            Q_high_here(:,k) = mean(Q03_high_each(:,:,k),2); % average of the positive and negative click
            pl(end+1) = semilogx(fc_ref,Q_high_here(:,k),format_kv{:},'MarkerFaceColor',Colours{k}); hold on
        end
        grid on;
        set(gca,'XTick',XT);
        YT = 2:2:26;
        set(gca,'YTick',YT);

        ylim([1 27])

        xlabel('Frequency (Hz)');
        ylabel('Q factor');

        Pos = get(gcf,'Position');
        Pos(3:4) = Pos34;
        set(gcf,'Position',Pos);

        figure_handle(end+1) = gcf;
        figure_name{end+1}   = sprintf('fig06-QERB-high-level');

        % -----------------------------------------------------------------
        figure; % Difference Low - high

        % Analytical expressions:
        pl = [];
        for k = 1:N_models_here
            switch models{k}
                case 'dau1997'
                    Marker = [Markers{k} '-'];
                    MSize = 7;
                    LW = 1;
                case 'zilany2014'
                    Marker = [Markers{k} LineStyle{k}];
                    MSize = 7;
                    LW = LineWidth(k);
                otherwise
                    Marker = LineStyle{k};
                    MSize = MarkersSize(k);
                    LW = LineWidth(k);
            end
            format_kv = {Marker,'Color',Colours{k},'MarkerSize',MSize,'LineWidth',LW};
            pl(end+1) = semilogx(fc_ref,Q_low_here(:,k)-Q_high_here(:,k),format_kv{:},'MarkerFaceColor',Colours{k}); hold on
        end
        grid on;

        Pos = get(gcf,'Position');
        Pos(3:4) = [600 250];
        set(gcf,'Position',Pos);

        set(gca,'XTick',XT);
        ylim([-4.5 18.5])
        YT = -4:2:18;
        set(gca,'YTick',YT);

        xlabel('Frequency (Hz)');
        ylabel('Q factor difference');

        figure_handle(end+1) = gcf;
        figure_name{end+1}   = sprintf('fig06-QERB-difference-high-low');

        X_ref_empirical = [];
        for k = 1:N_models_here
            X_ref_empirical(:,k) = Q_high_here(:,k);
            Y_10 = Q10_high_each(:,1,k);
            Y_03 = Q03_high_each(:,1,k);
        end

        f_low = [];
        f_high = [];
        count = [];
        for k = 1:N_models_here
            f_current = 125;
            f_high_lim = 8000;

            count(k) = 0;
            while f_current < f_high_lim
                count(k) = count(k)+1;
                f_low(count(k),k) = f_current;
                if count(k)>10000; warning('break - this should not happen'); break; end % break the loop if got stuck
                % Q_here = interp1(fc_ref,X_ref(:,k),f_low(count(k),k),'linear','extrap');
                Q_here = interp1(fc_ref, mean(Q03_high_each(:,:,k),2), f_low(count(k),k),'linear','extrap');
                if isnan(Q_here)
                    warning('Warning what happens here') % Q_here = Q03_recons_hi(1,k);
                end

                BW_here = f_low(count(k),k)/Q_here;
                f_high(count(k),k) = f_low(count(k),k) + BW_here;

                f_current = f_high(count(k),k);
            end

            fc_max = max(f_high);
            fc_min = f_low(1,:);
        end
        counts(3,:)   = count;
        erb_step(3,:) = (freqtoaud(fc_max)-freqtoaud(fc_min))./(count-1);
        data.figure_handle = figure_handle;
        data.figure_name   = figure_name;
        data.counts        = counts;
        data.counts_description = 'Number of bands required to cover a frequency range from 125 to 8000 Hz';
        data.erb_step      = erb_step;
        data.erb_step_description = 'ERB spacing that would ensure filters overlapped at their -3 dB points';
    end
    
    data.data_fig4     = data_fig4;
    
    data.figure_flag   = 'do_fig6';
    data.models        = models;
end

%% ------ FIG 7 or 8 ------------------------------------------------------
if flags.do_fig7 || flags.do_fig8

    k2remove = [];
    N_models_here = N_models; % Bruce2018, is removed later...

    t_ms_all = [];
    vihc_all = [];
    ti_all = [];

    %%% 1. Generating the input signals: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Signal parameters:
    dur = 0.08; % 80 ms
    lvl = 80; % 80 dB SPL
    dur_i = 0.02; % starting time for AC assessment
    dur_f = 0.07; % ending time for AC assessment

    % Arbitrary zero padding before and after the signal:
    dur_sil = 50e-3;
    sil_bef = zeros(round(dur_sil*fs),1);

    idxi = length(sil_bef)+1; % index corresponding to the first signal sample after the zero padding

    % AC will be calculated between dur_i and dur_f, excluding the silence:
    idxi_ac = idxi+round(dur_i*fs)-1;
    idxf_ac = idxi+round(dur_f*fs)-1;

    [cf,x] = local_m2hz(401); % Get characteristic frequencies from Verhulst models

    cf_max = 4000;
    cf_min =  125;
    %%%
    N_sigs = 12;

    vac  = nan(N_models, N_sigs);
    vdc  = nan(size(vac));

    %%% Getting f0 ('f0_method'):
    df_bin = 25;
    idx_cf(N_sigs) = find(cf>cf_max,1,'last');

    for i = N_sigs-1:-1:1
        idx_cf(i) = idx_cf(i+1)+df_bin;
    end
    f0 = cf(idx_cf);
    %%%

    dt = 1/fs; % \Delta t in s
    t=0:dt:dur-dt;

    for i = 1:N_sigs
        insig(:,i) = sin(2*pi*f0(i)*t); % the same sinusoid will be scaled at different levels
    end
    % Up/down cosine ramp (fixed)
    dur_ramp_ms = 5; % in ms
    dur_ramp = round((dur_ramp_ms*1e-3)*fs); % duration ramp in samples

    rp    = ones(size(insig(:,1)));
    rp(1:dur_ramp) = rampup(dur_ramp);
    rp(end-dur_ramp+1:end) = rampdown(dur_ramp);

    % Adjusting the signal level before the ramp and silences are added:
    insig = scaletodbspl(insig,lvl,dBFS);
    % Adding the silence and ramps to the signals:
    insig = [repmat(sil_bef,1,N_sigs); repmat(rp,1,N_sigs).*insig];
    %%%

    for k = 1:N_models_here

        %%% Loading flags and keyvals
        [fg,kv] = local_get_flags(models{k});
        fc_flags   = fg.fc_flags;
        afb_flags  = fg.afb_flags;
        afb_kv     = kv.afb_keyvals;
        ihc_flags  = fg.ihc_flags;
        noan_flags = fg.noan_flags; % No auditory nerve module
        nomfb_flags = fg.nomfb_flags;
        %%%

        fname = ['fig07_and_08_ihc-' models{k}];
        c = amt_cache('get',fname,flags.cachemode);

        if ~isempty(c)
            bRun = 0;

            vihc = c.vihc;
            fs   = c.fs;
            fc   = c.fc;
        else
            bRun = 1;
            vihc = [];
        end

        amt_disp(['Calculating ' models{k} '...']);
        switch models{k}
            case 'dau1997'
                zoom_higher(k) = 1;
                if bRun
                    for j = 1:N_sigs
                        fc_here = f0(j);
                        fc_kv = {'flow',fc_here,'fhigh',fc_here,'basef',fc_here,'dboffset',dBFS};
                        [outsig,fc] = dau1997(insig(:,j),fs,fc_kv{:},afb_flags{:}, ...
                            ihc_flags{:},noan_flags{:},nomfb_flags{:});
                        vihc(:,j) = outsig;
                    end

                    c = [];
                    c.vihc = vihc;
                    c.fs   = fs;
                    c.fc   = fc;
                    amt_cache('set',fname,c);
                end

            case 'zilany2014'
                zoom_higher(k) = 3;
                if bRun
                    kv={'nrep',1}; % kv={'numH',0,'numM',0,'numL',0,'fiberType',4,'nrep',1};
                    for j = 1:N_sigs
                        [~,~,vihc(:,j)] = zilany2014(insig(:,j),fs,f0(j),kv{:},afb_flags{:});
                    end
                    c = [];
                    c.vihc = vihc;
                    c.fs   = fs;
                    c.fc   = fc;
                    amt_cache('set',fname,c);
                end

            case 'verhulst2015'
                zoom_higher(k) = 3;
                if bRun
                    outsig = verhulst2015(insig,fs,fc_flags{:},afb_flags{:}, ...
                        ihc_flags{:},noan_flags{:},nomfb_flags{:});
                    for j = 1:N_sigs
                        vihc(:,j) = outsig(j).ihc(:,idx_cf(j));
                    end
                    c = [];
                    c.vihc = vihc;
                    c.fs   = fs;
                    c.fc   = fc;
                    amt_cache('set',fname,c);
                end

            case 'verhulst2018'
                zoom_higher(k) = 3;
                if bRun
                    outsig = verhulst2018(insig,fs,fc_flags{:},afb_flags{:}, ...
                        ihc_flags{:},noan_flags{:},nomfb_flags{:});

                    for j = 1:N_sigs
                        vihc(:,j) = outsig(j).ihc(:,idx_cf(j));
                    end
                    c = [];
                    c.vihc = vihc;
                    c.fs   = fs;
                    c.fc   = fc;
                    amt_cache('set',fname,c);
                end

            case 'bruce2018'
                % already in zilany2014
                k2remove = [k2remove k];

            case 'king2019'
                zoom_higher(k) = 1;
                if bRun
                    for j = 1:N_sigs
                        fc_here = f0(j);
                        fc_kv = {'flow',fc_here,'fhigh',fc_here,'basef',fc_here,'dboffset',dBFS};
                        [outsig,fc] = king2019(insig(:,j),fs,fc_kv{:},afb_kv{:}, ...
                            afb_flags{:},ihc_flags{:},noan_flags{:},nomfb_flags{:});
                        vihc(:,j) = outsig;
                    end

                    c = [];
                    c.vihc = vihc;
                    c.fs   = fs;
                    c.fc   = fc;
                    amt_cache('set',fname,c);
                end

            case 'relanoiborra2019'
                zoom_higher(k) = 3;
                if bRun
                    for j = 1:N_sigs
                        fc_here = f0(j);
                        fc_kv = {'flow',fc_here,'fhigh',fc_here,'basef',fc_here};
                        [~,~,outsig,fc] = relanoiborra2019_featureextraction(insig(:,j), fs, ...
                            fc_kv{:},afb_flags{:},noan_flags{:});
                        vihc(:,j) = outsig;
                    end
                    c = [];
                    c.vihc = vihc;
                    c.fs   = fs;
                    c.fc   = fc;
                    amt_cache('set',fname,c);
                end

            case 'osses2021'
                zoom_higher(k) = 3;
                if bRun
                    for j = 1:N_sigs
                        fc_here = f0(j);
                        fc_kv = {'flow',fc_here,'fhigh',fc_here,'basef',fc_here,'dboffset',dBFS};
                        [outsig,fc] = osses2021(insig(:,j),fs,fc_kv{:},afb_flags{:}, ...
                            ihc_flags{:},noan_flags{:},nomfb_flags{:});
                        vihc(:,j) = outsig;
                    end
                    c = [];
                    c.vihc = vihc;
                    c.fs   = fs;
                    c.fc   = fc;
                    amt_cache('set',fname,c);
                end
        end

        dc_reg = 1:length(sil_bef);
        ac_reg = idxi_ac:idxf_ac;

        if k~=5,
          [ac_target,dc_target,vr] = local_get_ac_dc_osses2022(vihc,ac_reg,dc_reg);

          vac(k,:) = ac_target;
          vdc(k,:) = dc_target;
          vrest(k,1) = vr(1);

          % Plotting: -------------------------------------------------------
          % % Normalisation to 1.5 the maximum value:
          v_norm_factor(k) = 1.5*max(max(abs(vihc)));

          t_ms = 1000*(1:size(vihc,1))/fs;
          t2plot = [35 105];
          idxi_plot = find(t_ms<=t2plot(1),1,'last');
          idxf_plot = find(t_ms<=t2plot(2),1,'last');
          idxs = idxi_plot:idxf_plot-1;

          v_dec_norm = (vihc-vrest(k))/v_norm_factor(k); % normalisation to maximum value

          idxs_N = 9:12;
          v_dec_norm(:,idxs_N) = zoom_higher(k)*v_dec_norm(:,idxs_N);

          idxs2store = idxs(1:5000);
          t_ms_all = [t_ms_all t_ms(idxs2store)];
          ti_all(end+1) = size(vihc_all,1)+1;
          vihc_all = [vihc_all; v_dec_norm(idxs2store,:)];

          vmin(k,:) = min(vihc(idxs2store,:));
        end

    end

    models(k2remove) = [];
    Colours(k2remove) = [];
    Markers(k2remove) = [];
    MarkersSize(k2remove) = [];
    LineStyle(k2remove) = [];
    LineWidth(k2remove) = [];

    zoom_higher(k2remove)  = [];
    vac(k2remove,:)  = []; vdc(k2remove,:)  = []; vmin(k2remove,:) = [];
    v_norm_factor(k2remove) = [];

    N_models_here = length(models); % number of models after removing

    if flags.do_plot
        %%% All formatting options in one variable. The formatting options
        %     for the bruce2018 are overwritten for visibility reasons,
        %     given that it provides the same outputs of zilany2014
        Format = [];
        for k = 1:N_models_here
            Format{k} = {'Color',Colours{k},'LineStyle',LineStyle{k}, ...
            'LineWidth',LineWidth(k),'Marker',Markers{k},'MarkerFaceColor','w', ...
            'MarkerSize',MarkersSize(k)};
        end
        %%% End formatting options
    end

    ti_all(end+1) = length(t_ms_all);
    toffset = t_ms(idxs(1));
    tf = 0;
    idxi = 1;
    if flags.do_plot
        YL = [-.5 N_sigs];
        plt = [];
        idx_split = 8;

        panel_labels = {'a) ','b) ','c) ','d) ','e) ','f) ','g) '};
        for k = 1:N_models_here
            if strcmp(models{k},'zilany2014')
                leg{k} = [panel_labels{k} 'zilany2014,bruce2018'];
            else
                leg{k} = [panel_labels{k} models{k}];
            end
        end
    end

    for k = 1:N_models_here

        if flags.do_plot
            if k == 1
                figure;
                figure_handle(end+1) = gcf;
                figure_name{end+1} = 'fig07-ihc-waveforms';
            end
        end

        idxi = ti_all(k);
        idxf = ti_all(k+1)-1;
        t_here    = t_ms_all(idxi:idxf);
        vihc_here = vihc_all(idxi:idxf,:);
        
        if flags.do_plot
            for j = 1:N_sigs
                offy = j-1; % (j-1) lower to higher freqs from bottom to top. Use (N_sigs-j) otherwise
                if j == 1
                    plt(end+1) = plot(t_here-toffset+tf,vihc_here(:,j)+offy,'Color',Colours{k});
                else
                    plot(t_here-toffset+tf,vihc_here(:,j)+offy,'Color',Colours{k});
                end
                hold on
                if j == 1 && (k == 1  || k == idx_split)
                    ylim(YL);
                end
                plot(tf*[1 1],YL,'k-');

                if (k == 1 || k == idx_split)
                    text(2,offy+.3,sprintf('%.0f Hz',f0(j)),'Color','k','FontSize',8);
                    grid on

                    dur_plot = 1000*length(idxs2store)/fs;
                    XL = [0 (idx_split-1)*dur_plot];
                    xlim(XL);
                    Pos = get(gcf,'Position');
                    Pos(3:4) = [1300 530];
                    set(gcf,'Position',Pos);

                    XT  = 0:10:XL(end);
                    XTL = [0 repmat([10:10:dur_plot],1,(idx_split-1))];

                    XT(dur_plot/10+1:dur_plot/10:end) = [];
                    XTL(dur_plot/10+1:dur_plot/10:end) = [];
                    set(gca,'XTick',XT);
                    set(gca,'XTickLabel',XTL);
                    YT = 0:N_sigs;
                    set(gca,'YTick',YT);
                    set(gca,'YTickLabel','');
                    set(gca,'FontSize',MediumFontSize);
                end

                if j >= 9
                    % For the last three waveforms...
                    text((k*dur_plot-8),offy-.07,sprintf('x %.0f',zoom_higher(k)),'Color','k','FontSize',SmallFontSize);
                end
            end % end for j

            x_coor = ((mod(k,idx_split))-1)*1/(idx_split-1)+.01;
            y_coor = .98;
            text(x_coor,y_coor,leg{k},'Units','Normalize','FontSize',SmallFontSize);

            if k == idx_split-1
                tf = 0;
            else
                tf = tf+(t_here(end)+dt*1000-toffset);
            end
            ylabel('Amplitude (a.u.)');
            xlabel('Time (ms)');
            
            plot(tf*[1 1],YL,'k-'); % last time for the second figure
        end
    end
    
    ra  = vac ./vdc;

    if flags.do_plot
        % Fig. 8 ----------------------------------------------------------
        figure;
        for k = 1:N_models_here
            pl(k) = loglog(f0,ra(k,:),Format{k}{:}); grid on, hold on
            xlabel('Frequency (Hz)');
            ylabel('IHC AC/DC ratio');
        end
        set(gca,'XTick',[80 125 250 500 1000 2000 4000 8000]);
        YT = get(gca,'YTick');
        set(gca,'YTickLabel',YT);
        xlim([100 6000]);
        ylim([0.008 180]);

        Pos = get(gcf,'Position');
        Pos(3:4) = [600 300];
        set(gcf,'Position',Pos);

        leg = models;
        if strcmp(leg{2},'zilany2014')
            leg{2} = 'zilany2014,bruce2018';
        end
        legend(leg,'Location','SouthWest');

        figure_handle(end+1) = gcf;
        figure_name{end+1} = 'fig08-IHC-AC-DC';

        data.figure_flag   = 'do_fig7 or do_fig8';
        data.figure_handle = figure_handle;
        data.figure_name   = figure_name;
    end
    
    data.insig = insig;
    data.fs    = fs;
    data.ramp  = rp;

    data.models = models;
    data.ACDCratio = ra;
    data.AC    = vac;
    data.DC    = vdc;

    data.dBFS  = dBFS;
    data.lvl   = lvl;
    data.f0_approx = round(f0');
    data.f0_exact  = round(cf(idx_cf)');

    data.vrest = vrest;
    data.vmin  = vmin;
    data.v_norm_factor = v_norm_factor;
end

%% ------ FIG. 9 ----------------------------------------------------------
if flags.do_fig9
    outsig_all = [];
    t_an_all = [];
    
    % 1. Stimulus creation:
    dur = 300e-3; % 300 ms
    lvl = 70;

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

	ti_steady = 300*1e-3; % ms
	tf_steady = 340*1e-3; % ms

    for k = 1:N_models
        %%% Loading flags and keyvals
        [fg,kv] = local_get_flags(models{k});
        afb_flags = fg.afb_flags;
        afb_kv    = kv.afb_keyvals;
        ihc_flags = fg.ihc_flags;
        an_flags  = fg.an_flags;
        an_kv     = kv.an_keyvals;
        nomfb_flags = fg.nomfb_flags;
        
        fname = ['fig09_' models{k} '-an-wave'];
        c = amt_cache('get',fname,flags.cachemode);

        if ~isempty(c)
            bRun  = 0;

            outsig = c.outsig;
            fs_an  = c.fs_an;
            fc     = c.fc;
            if isfield(c,'out_psth')
                out_psth = c.out_psth;
            end
        else
            bRun = 1;
            outsig = [];
        end

        YTL = []; % empty tick label
        
        if bRun && ~isempty(an_kv)
            % Loading common parameters
            switch models{k}
                case {'verhulst2015','verhulst2018','zilany2014','bruce2018'}
                    if strcmp(an_kv{1},'numH')
                        numH = an_kv{2};
                    else
                        error('an_kv should contain numH')
                    end
             
                    if strcmp(an_kv{3},'numM')
                        numM = an_kv{4};
                    else
                        error('an_kv should contain numM')
                    end
             
                    if strcmp(an_kv{5},'numL')
                        numL = an_kv{6};
                    else
                        error('an_kv should contain numM')
                    end
            end 
            switch models{k}
                case {'zilany2014','bruce2018'}
                    idx = find(strcmp(an_kv,'psth_binwidth'));
                    if ~isempty(idx)
                        % Only true for bruce2018:
                        psth_binwidth = an_kv{idx+1};
                    end
                    nrep = 100;
            end
        end
        
        if bRun
            amt_disp(['Calculating ' models{k} '...']);
            switch models{k}
                case 'dau1997'
                    fc_kv = {'flow',fc,'fhigh',fc,'basef',fc,'dboffset',dBFS};

                    outsig = dau1997(insig,fs,fc_kv{:},afb_flags{:}, ...
                        ihc_flags{:},an_flags{:},nomfb_flags{:});
                    fs_an = fs;

                    c = [];
                    c.outsig = outsig;
                    c.fs_an  = fs_an;
                    c.fc     = fc;
                    amt_cache('set',fname,c);

                case 'zilany2014'
                    kv={'fiberType',4,'numH',numH,'numM',numM,'numL',numL,'nrep',nrep,'psth_binwidth',psth_binwidth};
                    [outsig,psth,~,~,~,~,out] = zilany2014(insig,fs,fc,kv{:},afb_flags{:});
                    fs_an = fs;
                    out_psth.psth = psth;
                    out_psth.psth_binwidth = psth_binwidth;

                    if isfield(out,'psth_LSR')
                        out_psth.psth_LSR = out.psth_LSR;
                    end
                    if isfield(out,'psth_MSR')
                        out_psth.psth_MSR = out.psth_MSR;
                    end
                    if isfield(out,'psth_HSR')
                        out_psth.psth_HSR = out.psth_HSR;
                    end

                    c = [];
                    c.outsig = outsig;
                    if isfield(out,'meanrate_LSR')
                        c.outsig_LSR = out.meanrate_LSR;
                    end
                    if isfield(out,'meanrate_MSR')
                        c.outsig_MSR = out.meanrate_MSR;
                    end
                    if isfield(out,'meanrate_HSR')
                        c.outsig_HSR = out.meanrate_HSR;
                    end
                    c.fs_an  = fs_an;
                    c.fc     = fc;
                    c.out_psth = out_psth;

                    amt_cache('set',fname,c);

                case 'verhulst2015'
                    fc_flag = fc; % one frequency will be simulated only
                    out = verhulst2015(insig,fs,fc_flag,afb_flags{:},an_flags{:},'anfL','anfM','anfH',an_kv{:},nomfb_flags{:});
                    outsig = out(1).an_summed/(numL + numM + numH);
                    fs_an = out(1).fs_an;

                    c = [];
                    c.outsig = outsig;
                    c.fs_an  = fs_an;
                    c.fc     = fc;
                    amt_cache('set',fname,c);

                case 'verhulst2018'
                    fc_flag = fc; % one frequency will be simulated only
                    out = verhulst2018(insig,fs,fc_flag,afb_flags{:},an_flags{:},an_kv{:},nomfb_flags{:});
                    outsig = out(1).an_summed/(numL + numM + numH);
                    fs_an = out(1).fs_an;

                    c = [];
                    c.outsig = outsig;
                    c.fs_an  = fs_an;
                    c.fc     = fc;
                    amt_cache('set',fname,c);

                case 'bruce2018'
                    kv = {'numH',numH,'numM',numM,'numL',numL,'psthbinwidth_mr',psth_binwidth, ...
                        'nrep',nrep,'specificSRautoTiming'};
                    out = bruce2018(insig,fs,fc,kv{:},afb_flags{:});

                    outsig = out.meanrate;
                    out_psth.psth = out.psth;
                    out_psth.psth_binwidth = psth_binwidth;

                    if isfield(out,'psth_LSR')
                        out_psth.psth_LSR = out.psth_LSR;
                    end
                    if isfield(out,'psth_MSR')
                        out_psth.psth_MSR = out.psth_MSR;
                    end
                    if isfield(out,'psth_HSR')
                        out_psth.psth_HSR = out.psth_HSR;
                    end

                    fs_an = fs;

                    c = [];
                    c.outsig = outsig;
                    if isfield(out,'meanrate_LSR')
                        c.outsig_LSR = out.meanrate_LSR;
                    end
                    if isfield(out,'meanrate_MSR')
                        c.outsig_MSR = out.meanrate_MSR;
                    end
                    if isfield(out,'meanrate_HSR')
                        c.outsig_HSR = out.meanrate_HSR;
                    end
                    c.fs_an  = fs_an;
                    c.fc     = fc;
                    c.out_psth = out_psth;

                    amt_cache('set',fname,c);

                case 'king2019'
                    fc_kv = {'flow',fc,'fhigh',fc,'basef',fc,'dboffset',dBFS};

                    outsig = king2019(insig,fs,fc_kv{:},afb_flags{:}, afb_kv{:}, ...
                        ihc_flags{:},an_flags{:},nomfb_flags{:});
                    fs_an = fs;

                    c = [];
                    c.outsig = outsig;
                    c.fs_an  = fs_an;
                    c.fc     = fc;
                    amt_cache('set',fname,c);

                case 'relanoiborra2019'
                    fc_kv = {'flow',fc,'fhigh',fc,'basef',fc};
                    [~,~,outsig] = relanoiborra2019_featureextraction(insig, fs,fc_kv{:},afb_flags{:});
                    fs_an = fs;

                    c = [];
                    c.outsig = outsig;
                    c.fs_an  = fs_an;
                    c.fc     = fc;
                    amt_cache('set',fname,c);

                case 'osses2021'
                    fc_kv = {'flow',fc,'fhigh',fc,'basef',fc,'dboffset',dBFS};

                    outsig = osses2021(insig,fs,fc_kv{:},afb_flags{:}, ...
                        ihc_flags{:},an_flags{:},nomfb_flags{:});
                    fs_an = fs;

                    c = [];
                    c.outsig = outsig;
                    c.fs_an  = fs_an;
                    c.fc     = fc;
                    amt_cache('set',fname,c);
            end
        end

        switch models{k}
            case {'zilany2014','bruce2018'}
                pst_hist = out_psth.psth;
                dt = out_psth.psth_binwidth;
                t_psth = (1:length(pst_hist))*dt;
        end

        t_an = (1:length(outsig))/fs_an;
        t_an = t_an(:); % creates 't_an' as a column array

        outsig_all{k} = outsig;
        t_an_all{k}   = t_an;

        id_steady = find(t_an>=ti_steady & t_an<=tf_steady);

        onset  = max(outsig);
        steady = mean(outsig(id_steady));

        ra(k) = onset/steady;

        data_sim(k).onset = onset;
        data_sim(k).onset = steady;
        fs_an_all(k) = fs_an;
    end

    if flags.do_plot
        for k = 1:N_models

            psth_offx = 20;
            switch models{k}
                case {'dau1997','relanoiborra2019','osses2021'}
                    figure(1);
                    if strcmp(models{k},'osses2021') % osses2021 will be shifted
                        offx = psth_offx;
                    else
                        offx = 0;
                    end
                    handle_data(k) = plot(t_an_all{k}*1000+offx,outsig_all{k},'Color',Colours{k},'LineWidth',2); hold on, grid on;

                    if k == 1
                        figure_handle(end+1) = gcf; % multiple figures will be generated
                        figure_name{end+1}   = ['fig09-tone-4-kHz-' models{k}];
                    end

                    units_amplitude = '(MU)';
                    YL = [-300 1500]; % a priori knowledge
                    stepY = 100;
                    YT = YL(1)+stepY:stepY:YL(2)-stepY;

                    if k == 1
                        YTL = [];
                        for ii = 1:length(YT)
                            if mod(YT(ii),200)==0
                                YTL{ii} = num2str(YT(ii));
                            else
                                YTL{ii} = '';
                            end
                        end
                    end

                case {'zilany2014','bruce2018'}
                    figure;
                    handle_data(k) = stairs(t_psth*1000 + psth_offx,pst_hist,'Color',0.5*Colours{k},'LineWidth',2); hold on, grid on;
                    plot(t_an_all{k}*1000,outsig_all{k},'Color',Colours{k},'LineWidth',2); hold on, grid on;

                    % Trick to add labels later:
                    plot(t_psth(1)*1000 + psth_offx,pst_hist(1),'Color',0.5*Colours{k},'LineWidth',2)
                    % End: Trick
                    figure_handle(end+1) = gcf; % multiple figures will be generated
                    figure_name{end+1}   = ['fig09-tone-4-kHz-' models{k}];

                    units_amplitude = '(spikes/s)';
                    YL = [0 1000]; % a priori knowledge
                    stepY = 100;
                    YT = YL(1)+stepY:stepY:YL(2)-stepY;

                case {'verhulst2015','verhulst2018'}
                    figure(3)

                    if strcmp(models{k},'verhulst2015')
                        offx = psth_offx;
                    else
                        offx = 0;
                    end
                    handle_data(k) = plot(t_an_all{k}*1000 + offx,outsig_all{k},'Color',Colours{k},'LineWidth',2); hold on, grid on;
                    if k == 3
                        figure_handle(end+1) = gcf; % multiple figures will be generated
                        figure_name{end+1}   = ['fig09-tone-4-kHz-' models{k}];
                    end

                    units_amplitude = '(spikes/s)';
                    YL = [0 1000]; % a priori knowledge
                    stepY = 100;
                    YT = YL(1)+stepY:stepY:YL(2)-stepY;

                case 'king2019'
                    factor = 1e3;
                    figure;
                    handle_data(k) = plot(t_an_all{k}*1000,factor*outsig_all{k},'Color',Colours{k},'LineWidth',2); hold on, grid on;
                    if k == 6
                        figure_handle(end+1) = gcf; % multiple figures will be generated
                        figure_name{end+1}   = ['fig09-tone-4-kHz-' models{k}];
                    end

                    units_amplitude = '(arbitrary units x 10^{-3})';
                    YL = factor*1.8e-3*[-1 1]; % a priori knowledge
                    stepY = factor*.3e-3;
                    YT = YL(1)+stepY:stepY:YL(2)-stepY;
            end
            xlabel('Time (ms)');
            ylabel(['Amplitude ' units_amplitude])

            set(gca,'XTick',50:50:450);
            xlim([25 475]);

            ylim(YL);
            set(gca,'YTick',YT);
            switch k
                case {1,7,8}
                    set(gca,'YTickLabel',YTL);
            end

            offy = 5;

            Xhere = 1000*[t_an(min(id_steady)) t_an(max(id_steady))];
            Yhere = (YL(1)+offy)*[1 1];
            plot(Xhere,Yhere,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);

            Yhere = (YL(2)-offy)*[1 1];
            plot(Xhere,Yhere,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);

            Yhere = YL;
            Xhere = 1000*t_an(min(id_steady))*[1 1];
            plot(Xhere,Yhere,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);

            Xhere = 1000*t_an(max(id_steady))*[1 1];
            plot(Xhere,Yhere,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);

            Pos = get(gcf,'Position');
            Pos(3:4) = [600 250];
            set(gcf,'Position',Pos);
        end
    end

    if flags.do_plot
        data.figure_flag = 'do_fig9';
    end
    data.ra = ra;

    data.insig = insig;
    data.fs    = fs;
    data.data_sim = data_sim;
    data.models = models;
    if flags.do_plot
        data.figure_handle = figure_handle;
        data.handle_data   = handle_data;
        data.figure_name   = figure_name;
    end
    data.outsig_all = outsig_all;
    data.fs_an_all  = fs_an_all;
end

%% ------ FIG 10 or FIG11 -------------------------------------------------
if flags.do_fig10 || flags.do_fig11
    % Code adapted from testPopRateLevel_BEZ2018_mine.m
    %%% 1. Generating the input signals: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Signal parameters:
    CF    = 4000; % CF in Hz;
    lvls = 0:10:100;
    N_lvls = length(lvls);

    dur     = 300e-3; % stimulus duration in seconds
    dur_sil =  50e-3; % duration of silence at the beginning

    dur_ramp_ms = 2.5; % rise/fall time in seconds
    dur_ramp = round((dur_ramp_ms*1e-3)*fs); % duration ramp in samples

    t = 0:1/fs:dur-1/fs; t = t(:); % time vector

    insig_orig = sin(2*pi*CF*t);
    insig = zeros(length(t),N_lvls); % memory allocation

    rp    = ones(size(insig_orig));
    rp(1:dur_ramp) = rampup(dur_ramp);
    rp(end-dur_ramp+1:end) = rampdown(dur_ramp);

    for kk = 1:N_lvls
        lvl = lvls(kk);
        insig(:,kk) = scaletodbspl(insig_orig,lvl,dBFS); % calibration before applying the ramp
    end
    insig = repmat(rp,1,N_lvls).*insig;
    %%% End generating the signal (the level adjustment is done later)

    insig = [zeros(dur_sil*fs,N_lvls); insig]; % 50 ms of silence before and after the sine tone

    ti_steady_ms = (dur+dur_sil)*1e3 - 50; % 300 ms (if dur=300 ms)
    tf_steady_ms = (dur+dur_sil)*1e3 - 10; % 340 ms (if dur=340 ms)

    numH = 12; % Automate this
    numM =  4;
    numL =  4;
    numTot = numH+numM+numL;
    numCum = cumsum([numL numM numH]);

    rates_steady = nan(numTot,N_lvls,N_models);
    rates        = nan(numTot,N_lvls,N_models);

    for k = 1:N_models

        meanAN_H = [];
        meanAN_M = [];
        meanAN_L = [];

        %%% Loading flags and keyvals
        [fg,kv] = local_get_flags(models{k});
        afb_flags = fg.afb_flags;
        afb_kv    = kv.afb_keyvals;
        ihc_flags = fg.ihc_flags;
        an_flags = fg.an_flags;
        nomfb_flags = fg.nomfb_flags;
        an_kv    = kv.an_keyvals;
        %%%

        fname = ['fig10_rate-level-' models{k}];
        c = amt_cache('get',fname,flags.cachemode);

        if ~isempty(c)
            bRun  = 0;

            meanAN_H = c.meanAN_H;
            meanAN_M = c.meanAN_M;
            meanAN_L = c.meanAN_L;
            rates  = c.rates;
            fs_an  = c.fs_an;
            idxi   = c.idxi;
            idxf   = c.idxf;
            outsig = c.outsig;

            if isfield(c,'psth_L')
                psth_L = c.psth_L;
            end
            if isfield(c,'psth_M')
                psth_M = c.psth_M;
            end
            if isfield(c,'psth_H')
                psth_H = c.psth_H;
            end
            if isfield(c,'rates_psth')
                rates_psth = c.rates_psth;
            end
        else
            bRun = 1;
        end

        if ~isempty(an_kv)
            switch models{k}
                case {'zilany2014','bruce2018'}
                    idx = find(strcmp(an_kv,'psth_binwidth'));
                    if ~isempty(idx)
                        % Only true for bruce2018:
                        psth_binwidth = an_kv{idx+1};
                    end
                    nrep = 100;
            end
        end
            
        if bRun
            amt_disp(['Calculating ' models{k} '...']);
            
            % Initialisation:
            meanAN_H = []; 
            meanAN_M = []; 
            meanAN_L = []; 
            rates = [];
            switch models{k}
                case {'zilany2014','bruce2018'}
                    psth_L = []; psth_M = []; psth_H = []; rates_psth = [];
            end
            switch models{k}
                case 'dau1997'
                    fc_kv = {'flow',CF,'fhigh',CF,'basef',CF,'dboffset',dBFS};
                    for j = 1:N_lvls
                        meanAN_H(:,j) = dau1997(insig(:,j),fs,fc_kv{:}, ...
                            afb_flags{:},ihc_flags{:},an_flags{:},nomfb_flags{:});
                    end
                    outsig = meanAN_H;
                    fs_an = fs;

                    units_amplitude = '(Model Units)';
                    idxi = round(ti_steady_ms*1e-3*fs_an)+1;
                    idxf = round(tf_steady_ms*1e-3*fs_an);

                case 'zilany2014'
                    L = size(insig,1);
                    for j = 1:N_lvls

                        % out_psTH.psth_binwidth = psth_binwidth;
                        kv={'fiberType',4,'numH',numH,'numM',numM,'numL',numL,'nrep',nrep,'psth_binwidth',psth_binwidth};
                        [outsig,psth,~,~,~,~,out] = zilany2014(insig(:,j),fs,CF,kv{:},afb_flags{:});

                        meanAN_L(:,j) = out.meanrate_LSR(1:L);
                        meanAN_M(:,j) = out.meanrate_MSR(1:L);
                        meanAN_H(:,j) = out.meanrate_HSR(1:L);
                        rates(:,j) = outsig(1:L);

                        psth_L(:,j) = out.psth_LSR;
                        psth_M(:,j) = out.psth_MSR;
                        psth_H(:,j) = out.psth_HSR;
                        rates_psth(:,j) = psth;

                        outsig = rates;
                    end
                    units_amplitude = '(spikes/s)';
                    fs_an = fs;

                    idxi = round(ti_steady_ms*1e-3*fs_an)+1; % round(15e-3*fs_an)+1; % start after 15 ms (strong onset);
                    idxf = round(tf_steady_ms*1e-3*fs_an); % end

                case 'verhulst2015'
                    fc_flag = CF;
                    out = verhulst2015(insig,fs,fc_flag,afb_flags{:},an_flags{:},'anfL','anfM','anfH',an_kv{:},nomfb_flags{:});
                    for j = 1:N_lvls
                        % psTH =out(kk).an_summed/numTot;
                        meanAN_L(:,j) = out(j).anfL;
                        meanAN_M(:,j) = out(j).anfM;
                        meanAN_H(:,j) = out(j).anfH;
                    end
                    units_amplitude = '(spikes/s)';
                    outsig = (numL*meanAN_L + numM*meanAN_M + numH*meanAN_H)/numTot;
                    fs_an = out(1).fs_an;

                    idxi = round(ti_steady_ms*1e-3*fs_an)+1;
                    idxf = round(tf_steady_ms*1e-3*fs_an);

                case 'verhulst2018'
                    fc_flag = CF;
                    out = verhulst2018(insig,fs,fc_flag,afb_flags{:},an_flags{:},'anfL','anfM','anfH',an_kv{:},nomfb_flags{:});
                    for j = 1:N_lvls
                        % psTH =out(kk).an_summed/numTot;
                        meanAN_L(:,j) = out(j).anfL;
                        meanAN_M(:,j) = out(j).anfM;
                        meanAN_H(:,j) = out(j).anfH;
                    end
                    units_amplitude = '(spikes/s)';

                    outsig = (numL*meanAN_L + numM*meanAN_M + numH*meanAN_H)/numTot;
                    fs_an = out(1).fs_an;

                    idxi = round(ti_steady_ms*1e-3*fs_an)+1;
                    idxf = round(tf_steady_ms*1e-3*fs_an);

                case 'bruce2018'
                    psth_L = []; psth_M = []; psth_H = [];
                    L = size(insig,1);
                    kv={'numH',numH,'numM',numM,'numL',numL,'nrep',nrep,'psthbinwidth_mr',psth_binwidth, ...
                        'specificSRautoTiming'};
                    for j = 1:N_lvls
                        out = bruce2018(insig(:,j),fs,CF,kv{:},afb_flags{:});

                        meanAN_L(:,j) = out.meanrate_LSR(1:L);
                        meanAN_M(:,j) = out.meanrate_MSR(1:L);
                        meanAN_H(:,j) = out.meanrate_HSR(1:L);
                        rates(:,j) = out.meanrate(1:L);

                        psth_L(:,j) = out.psth_LSR;
                        psth_M(:,j) = out.psth_MSR;
                        psth_H(:,j) = out.psth_HSR;
                        rates_psth(:,j) = out.psth;
                    end

                    % bruce2018 returns the 'weighted LSR, MSR, and HSR', unweighting (for visual scaling),
                    %  to obtain an average PSTH per neurone type...
                    disp('Applying scaling, assuming an_summed = numL*Raw_L + numM*Raw_M + numH*Raw_H')
                    psth_L = psth_L/(numL/numTot);
                    psth_M = psth_M/(numM/numTot);
                    psth_H = psth_H/(numH/numTot);

                    outsig = rates;
                    units_amplitude = '(spikes/s)';

                    fs_an = fs;

                    idxi = round(ti_steady_ms*1e-3*fs_an)+1;
                    idxf = round(tf_steady_ms*1e-3*fs_an);

                case {'king2019'}
                    fc_kv = {'flow',CF,'fhigh',CF,'basef',CF,'dboffset',dBFS};
                    % clear meanAN_H
                    for j = 1:N_lvls
                        meanAN_H(:,j) = king2019(insig(:,j),fs,fc_kv{:}, ...
                            afb_flags{:},afb_kv{:},ihc_flags{:},an_flags{:},nomfb_flags{:});
                    end
                    outsig = meanAN_H;
                    fs_an = fs;

                    factor = 1e3;
                    if factor == 1e3
                        units_amplitude = '(a.u. x 10^{-3})';
                    else
                        units_amplitude = '(a.u.)';
                    end

                    idxi = round(ti_steady_ms*1e-3*fs_an)+1;
                    idxf = round(tf_steady_ms*1e-3*fs_an);

                case {'relanoiborra2019'}
                    fc_kv = {'flow',CF,'fhigh',CF,'basef',CF};
                    for j = 1:N_lvls
                        [~,~,meanAN_H(:,j)] = relanoiborra2019_featureextraction(insig(:,j), fs,fc_kv{:},afb_flags{:});
                    end
                    outsig = meanAN_H;
                    fs_an = fs;

                    units_amplitude = '(Model Units)';
                    idxi = round(ti_steady_ms*1e-3*fs_an)+1;
                    idxf = round(tf_steady_ms*1e-3*fs_an);

                case 'osses2021'
                    fc_kv = {'flow',CF,'fhigh',CF,'basef',CF,'dboffset',dBFS};
                    for j = 1:N_lvls
                        meanAN_H(:,j) = osses2021(insig(:,j),fs,fc_kv{:}, ...
                            afb_flags{:},ihc_flags{:},an_flags{:},nomfb_flags{:});
                    end
                    outsig = meanAN_H;
                    fs_an = fs;

                    units_amplitude = '(Model Units)';

                    idxi = round(ti_steady_ms*1e-3*fs_an)+1;
                    idxf = round(tf_steady_ms*1e-3*fs_an);
            end

            c = [];
            c.meanAN_H = meanAN_H; 
            c.meanAN_M = meanAN_M; 
            c.meanAN_L = meanAN_L;
            c.rates = rates;
            c.outsig = outsig;

            switch models{k}
                case {'zilany2014','bruce2018'}
                    c.psth_binwidth = psth_binwidth;

                    c.psth_L = psth_L;
                    c.psth_M = psth_M;
                    c.psth_H = psth_H;
                    c.rates_psth  = rates_psth;
            end
            c.fs_an = fs_an;
            c.idxi = idxi;
            c.idxf = idxf;
            c.units_amplitude = units_amplitude;

            amt_cache('set',fname,c);
        end

        rates_max_tot(:,k) = max(outsig)';
        rates_mean_all(:,k) = mean(outsig(idxi:idxf,:))';

        for j = 1:N_lvls
            if ~isempty(meanAN_L)
                rates_max(          1:numCum(1),j,k)= max(meanAN_L(:,j));
                rates_steady(       1:numCum(1),j,k)= mean(meanAN_L(idxi:idxf,j));
            else
                rates_max(          1:numCum(1),j,k)= nan;
                rates_steady(       1:numCum(1),j,k)= nan;
            end
            if ~isempty(meanAN_M)
                rates_max(   numCum(1)+1:numCum(2),j,k)= max(meanAN_M(:,j));
                rates_steady(numCum(1)+1:numCum(2),j,k)= mean(meanAN_M(idxi:idxf,j));
            else
                rates_max(numCum(1)+1:numCum(2),j,k) = nan;
                rates_steady(numCum(1)+1:numCum(2),j,k) = nan;
            end
            rates_max(numCum(2)+1:numCum(3),j,k)= max(meanAN_H(:,j));
            rates_steady(numCum(2)+1:numCum(3),j,k)= mean(meanAN_H(idxi:idxf,j));

            switch models{k}
                case {'zilany2014','bruce2018'}

                    dt = psth_binwidth;
                    % t_psth = (1:size(rates_psth,1))*dt;

                    idxi_psth = round(ti_steady_ms*1e-3/dt)+1; % start after 15 ms (exludes the first 15 ms of strong onset)
                    idxf_psth = round(tf_steady_ms*1e-3/dt);

                    rates_steady_PSTH(1,j,k)= mean(psth_L(idxi_psth:idxf_psth,j));
                    rates_steady_PSTH(2,j,k)= mean(psth_M(idxi_psth:idxf_psth,j));
                    rates_steady_PSTH(3,j,k)= mean(psth_H(idxi_psth:idxf_psth,j));


                    rates_max_PSTH(1,j,k)= max(psth_L(:,j));
                    rates_max_PSTH(2,j,k)= max(psth_M(:,j));
                    rates_max_PSTH(3,j,k)= max(psth_H(:,j));

                    rates_max_PSTH(4,j,k)= max(rates_psth(:,j));
            end
        end

        Raw_L = squeeze(rates_steady(1:numCum(1),:,k));
        Raw_M = squeeze(rates_steady(numCum(1)+1:numCum(2),:,k));
        Raw_H = squeeze(rates_steady(numCum(2)+1:numCum(3),:,k));
        rates_avg_all{k} = squeeze(rates_steady(numCum(2)+1:numCum(3),:,k));
        if k == 1
            rates_avg_description = 'rate levels for ''HSR''';
        end
        if k == 1
            rates_max_description = 'onset rate levels for ''HSR''';
        end
        if flags.do_plot

            do_PSTH = 0;
            Pos34 = [320 250]; % width and height of the resulting plot
            factor = 1;

            if flags.do_fig10
                switch models{k}
                    case {'dau1997','relanoiborra2019','osses2021'}
                        figure(1);
                        if k == 1
                            figure_handle(end+1) = gcf;
                            figure_name{end+1} = ['fig10-AN-firing-rates-model-' models{k}];
                        else
                            hold on; % grid on
                        end
                        YL = 'Amplitude (MU)';
                        YLim = [-10 130];
                        YT = 0:20:130;

                    case {'zilany2014','bruce2018'}
                        do_PSTH = 1; % to do after...
                        figure;
                        YL = 'Firing rate (spikes/s)';
                        figure_handle(end+1) = gcf;
                        figure_name{end+1} = ['fig10-AN-firing-rates-model-' models{k}];
                        if strcmp(models{k},'zilany2014')
                            YLim = [-20 380];
                            YT = 0:30:360;
                        else
                            YLim = [-20 230];
                            YT = 0:30:300;
                        end

                    case {'verhulst2015','verhulst2018'}
                        figure

                        figure_handle(end+1) = gcf;
                        figure_name{end+1} = ['fig10-AN-firing-rates-model-' models{k}];

                        YL = 'Firing rate (spikes/s)';
                        YLim = [-20 300];
                        YT = 0:30:280;

                    case 'king2019'
                        factor = 1e5;

                        figure;

                        if factor == 1e5
                            YL = 'Amplitude (a.u. x 10^{-5})';
                        else
                            YL = 'Amplitude (a.u.)';
                        end
                        figure_handle(end+1) = gcf;
                        figure_name{end+1} = ['fig10-AN-firing-rates-model-' models{k}];
                        YLim = [-0.2 3.2];
                        YT = (0:.3:3);

                end

                switch models{k}
                    case {'zilany2014','bruce2018','verhulst2015','verhulst2018'}
                        % No change in Colour
                        Colour_here = local_rgb('Gray'); % Colours{k}
                        hl1 = plot(lvls,factor*mean(Raw_L),'o-','Color',Colour_here,'MarkerFaceColor','w','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w'); hold on, grid on
                        hl2 = plot(lvls,factor*mean(Raw_M),'^--','Color',Colour_here,'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',[0.8 0.8 0.8]);
                        hl3 = plot(lvls,factor*mean(Raw_H),'s-','Color',Colour_here,'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',Colour_here);
                    otherwise
                        hl1 = 0;
                        hl2 = 0;
                        hl3 = 0;
                end
                hl4 = plot(lvls,factor*rates_mean_all(:,k),'d-','Color',Colours{k},'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',Colours{k}); grid on

                handle_data(k,:) = [hl1(1) hl2(1) hl3(1) hl4(1)];

                xlim([min(lvls)-2 max(lvls)+2])
                ylim(YLim)

                xlabel('Stimulus Level (dB SPL)')
                ylabel(YL)

                set(gca,'XTick',0:10:100);
                set(gca,'YTick',YT);

                Pos = get(gcf,'Position');
                Pos(3:4) = Pos34;
                set(gcf,'Position',Pos);

                factors(k) = factor;
                %%%
                if do_PSTH

                    factor = 1;
                    figure;
                    YL = 'Firing rate (spikes/s)';
                    figure_handle(end+1) = gcf;
                    figure_name{end+1} = ['fig10-AN-firing-rates-model-' models{k} '-PSTH'];

                    Raw_L = squeeze(rates_steady_PSTH(1,:,k));
                    Raw_M = squeeze(rates_steady_PSTH(2,:,k));
                    Raw_H = squeeze(rates_steady_PSTH(3,:,k));

                    hl1 = plot(lvls,factor*Raw_L,'o-','Color' ,local_rgb('Gray'),'MarkerFaceColor','w','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w'); hold on, grid on
                    hl2 = plot(lvls,factor*Raw_M,'^--','Color',local_rgb('Gray'),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',[0.8 0.8 0.8]);
                    hl3 = plot(lvls,factor*Raw_H,'s-','Color' ,local_rgb('Gray'),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',local_rgb('Gray'));
                    psth_tot(k,:) = (numL*Raw_L+numM*Raw_M+numH*Raw_H)/(numL+numM+numH);
                    hl4 = plot(lvls,factor*psth_tot(k,:),'d-','Color',0.5*Colours{k},'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',0.5*Colours{k});
                    handle_data_psth(k,:) = [hl1(1) hl2(1) hl3(1) hl4(1)];

                    xlim([min(lvls)-2 max(lvls)+2])
                    ylim(YLim)

                    xlabel('Stimulus Level (dB SPL)')
                    ylabel(YL)

                    set(gca,'XTick',0:10:100);
                    set(gca,'YTick',YT);

                    Pos = get(gcf,'Position');
                    Pos(3:4) = Pos34;
                    set(gcf,'Position',Pos);
                end
            end

            %%%%
            if flags.do_fig11
                switch models{k}
                    case {'dau1997','relanoiborra2019','osses2021'}
                        figure(1);
                        if k == 1
                            figure_handle(end+1) = gcf;
                            figure_name{end+1} = sprintf('fig11-tone-4-kHz-IO-onset-%s',models{k});
                        else
                            hold on; % grid on
                        end
                        YL = 'Amplitude (MU)';
                        YLim = [-50 1650];
                        stepY = 100;
                        YT = 0:stepY:1600;

                        if k == 1
                            YTL = [];
                            for ii = 1:length(YT)
                                if mod(YT(ii),200)==0
                                    YTL{ii} = num2str(YT(ii));
                                else
                                    YTL{ii} = '';
                                end
                            end
                        end

                    case {'zilany2014','bruce2018'}

                        do_PSTH = 1; % to do after...
                        figure;
                        YL = 'Firing rate (Spikes/s)';
                        figure_handle(end+1) = gcf; % multiple figures will be generated
                        figure_name{end+1}   = sprintf('fig10-tone-4-kHz-IO-onset-%s',models{k});

                        if strcmp(models{k},'zilany2014')
                            YLim = [-50 1650];
                        else
                            YLim = [-50 1050];
                        end
                        stepY = 100;
                        YT = 0:stepY:1600;

                    case {'verhulst2015','verhulst2018'}
                        figure

                        figure_handle(end+1) = gcf; % multiple figures will be generated
                        figure_name{end+1}   = sprintf('fig10-tone-4-kHz-IO-onset-%s',models{k});

                        YL = 'Firing rate (Spikes/s)';
                        YLim = [-50 1650];
                        stepY = 100;
                        YT = 0:stepY:1600;

                    case 'king2019'
                        factor = 1e3;

                        figure;

                        if factor == 1e3
                            YL = '(a.u. x 10^{-3})';
                        else
                            YL = '(a.u.)';
                        end

                        figure_handle(end+1) = gcf; % multiple figures will be generated
                        figure_name{end+1}   = sprintf('fig10-tone-4-kHz-IO-onset-%s',models{k});
                        YLim = [-0.2 5.2];
                        YT = (0:.5:5);

                end
                Raw_L = squeeze(rates_max(1:numCum(1),:,k));
                Raw_M = squeeze(rates_max(numCum(1)+1:numCum(2),:,k));
                Raw_H = squeeze(rates_max(numCum(2)+1:numCum(3),:,k));

                switch models{k}
                    case {'zilany2014','bruce2018','verhulst2015','verhulst2018'}
                        % No change in Colour
                        Colour_here = local_rgb('Gray'); % Colours{k}
                        hl1 = plot(lvls,factor*mean(Raw_L),'o-','Color',Colour_here,'MarkerFaceColor','w','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w'); hold on, grid on
                        hl2 = plot(lvls,factor*mean(Raw_M),'^--','Color',Colour_here,'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',[0.8 0.8 0.8]);
                        hl3 = plot(lvls,factor*mean(Raw_H),'s-','Color',Colour_here,'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',Colour_here);
                    otherwise
                        hl1 = 0;
                        hl2 = 0;
                        hl3 = 0;
                end
                hl4 = plot(lvls,factor*rates_max_tot(:,k),'d-','Color',Colours{k},'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',Colours{k}); grid on

                handle_data(k,:) = [hl1(1) hl2(1) hl3(1) hl4(1)];

                xlim([min(lvls)-2 max(lvls)+2])
                ylim(YLim)

                xlabel('Stimulus Level (dB SPL)')
                ylabel(YL);
                set(gca,'YTick',YT);

                switch models{k}
                    case {'king2019'}

                    otherwise
                        set(gca,'YTickLabel',YTL);
                end

                set(gca,'XTick',0:10:100);

                Pos = get(gcf,'Position');
                Pos(3:4) = Pos34;
                set(gcf,'Position',Pos);

                factors(k) = factor;
                %%%
                if do_PSTH

                    factor = 1;

                    % case {'zilany2014','bruce2018'}
                    figure;
                    YL = 'Firing rate (Spikes/s)';
                    figure_handle(end+1) = gcf;
                    figure_name{end+1}   = sprintf('fig10-tone-4-kHz-IO-onset-%s-PSTH',models{k});

                    Raw_L = squeeze(rates_max_PSTH(1,:,k));
                    Raw_M = squeeze(rates_max_PSTH(2,:,k));
                    Raw_H = squeeze(rates_max_PSTH(3,:,k));
                    psth_tot(k,:) = squeeze(rates_max_PSTH(4,:,k))';
                    % psth_totTest = (numL*Raw_L+numM*Raw_M+numH*Raw_H)/(numL+numM+numH);

                    hl1 = plot(lvls,factor*Raw_L,'o-','Color' ,local_rgb('Gray'),'MarkerFaceColor','w','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w'); hold on, grid on
                    hl2 = plot(lvls,factor*Raw_M,'^--','Color',local_rgb('Gray'),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',[0.8 0.8 0.8]);
                    hl3 = plot(lvls,factor*Raw_H,'s-','Color' ,local_rgb('Gray'),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',local_rgb('Gray'));

                    hl4 = plot(lvls,factor*psth_tot(k,:),'d-','Color',0.5*Colours{k},'LineWidth',3,'MarkerSize',10,'MarkerFaceColor',0.5*Colours{k});
                    handle_data_psth(k,:) = [hl1(1) hl2(1) hl3(1) hl4(1)];

                    YLim = [-50 1650];

                    xlim([min(lvls)-2 max(lvls)+2])
                    ylim(YLim)

                    xlabel('Stimulus Level (dB SPL)')
                    ylabel(YL)

                    set(gca,'XTick',0:10:100);
                    set(gca,'YTick',YT);
                    set(gca,'YTickLabel',YTL);

                    Pos = get(gcf,'Position');
                    Pos(3:4) = Pos34;
                    set(gcf,'Position',Pos);
                end
            end
            %%%
        end
    end
    
    data.rate_avg_all  = rates_avg_all;
    data.rate_avg_description = rates_avg_description;
    if flags.do_fig10
        data.figure_flag   = 'do_fig10';
        % data.rate_avg_all  = rates_avg_all;
        % data.rate_avg_description = rates_avg_description;
    end
    data.rate_max_tot = rates_max_tot;
    data.rate_max_description = rates_max_description;
    if flags.do_fig11
        if flags.do_plot
            data.figure_flag   = 'do_fig11';
        end
    end
    
    if flags.do_plot
        data.factors = factors;
        data.figure_handle = figure_handle;
        data.figure_name   = figure_name;
        data.handle_data   = handle_data;
        if exist('handle_data_psth','var')
            data.handle_data_psth = handle_data_psth; 
        end
    end
    data.models = models;
    data.lvls = lvls;
end

%% ------ FIG 12 Osses et al. 2021 ----------------------------------------
if flags.do_fig12a || flags.do_fig12b

    figure_handle_ax = [];

    % This is similar to Verhulst et al. 2018, Fig. 3C:
    %%% 1. Generating the input signals: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Signal parameters:
    fc_target = 4000; % Hz, frequency of the carrier
    cf = local_m2hz(401); % Get characteristic frequencies from Verhulst models
    idx_cf = find(cf>fc_target,1,'last');
    fc = cf(idx_cf);

    fmod     = 100; % Hz, frequency of the modulator
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

    %%%
    psth_bin4plot = 0.5e-3; % 0.5 ms
    percent4plot = 90; % percentage overlap
    %%%
    outs_anf = [];

    for k = 1:N_models
        %%% Loading flags and keyvals:
        [fg,kv] = local_get_flags(models{k});
        fc_flags   = fg.fc_flags;
        afb_flags  = fg.afb_flags;
        afb_kv     = kv.afb_keyvals;
        ihc_flags  = fg.ihc_flags;
        an_flags   = fg.an_flags; % No auditory nerve module
        nomfb_flags= fg.nomfb_flags;
        an_kv = kv.an_keyvals;

        fname = ['fig12_' models{k} '-an'];
        c = amt_cache('get',fname,flags.cachemode);

        if ~isempty(c)
            bRun  = 0;

            an_summed = c.an_summed;
            fs_an     = c.fs_an;
            anfH      = c.anfH;
            anfM      = c.anfM;
            anfL      = c.anfL;
            if isfield(c,'out_psth')
                out_psth = c.out_psth;
            end
        else
            bRun = 1;
            an_summed = [];
            amt_disp(['Calculating ' models{k} '...']);
        end

        if bRun && ~isempty(an_kv)
            % Loading common parameters
            switch models{k}
                case {'verhulst2015','verhulst2018','zilany2014','bruce2018'}
                    if strcmp(an_kv{1},'numH')
                        numH = an_kv{2};
                    else
                        error('an_kv should contain numH')
                    end
             
                    if strcmp(an_kv{3},'numM')
                        numM = an_kv{4};
                    else
                        error('an_kv should contain numM')
                    end
             
                    if strcmp(an_kv{5},'numL')
                        numL = an_kv{6};
                    else
                        error('an_kv should contain numM')
                    end
            end 
            switch models{k}
                case {'zilany2014','bruce2018'}
                    idx = find(strcmp(an_kv,'psth_binwidth'));
                    if ~isempty(idx)
                        % Only true for bruce2018:
                        psth_binwidth = an_kv{idx+1};
                    end
                    nrep = 100;
            end
        end
        
        %%% Figure settings (overloaded later for dau1997):
        text4ylabel='(spikes/s)';
        YL = [-48 330];
        factor = 1;
        %%%
        switch models{k}
            case 'dau1997'
                if bRun
                    fc_kv = {'flow',fc,'fhigh',fc,'basef',fc,'dboffset',dBFS};
                    out = dau1997(insig,fs,fc_kv{:},afb_flags{:}, ...
                        ihc_flags{:},an_flags{:},nomfb_flags{:});
                    fs_an = fs;
                    an_summed = out;
                    anfH = out;
                    anfM = []; % 0*out;
                    anfL = []; % 0*out;

                    c = [];
                    c.an_summed = an_summed;
                    c.fs_an = fs_an;
                    c.anfH  = anfH;
                    c.anfM  = anfM;
                    c.anfL  = anfL;
                    amt_cache('set',fname,c);
                end

                text4ylabel='(MU)';
                YL = [-230 1450];

            case {'zilany2014'}
                if bRun
                    kv={'fiberType',4,'numH',numH,'numM',numM,'numL',numL,'nrep',nrep,'psth_binwidth',psth_binwidth};

                    [an_summed,out_psth.psth,~,~,~,~,out] = zilany2014(insig,fs,fc,kv{:},afb_flags{:});
                    fs_an = fs;
                    anfH = out.meanrate_HSR;
                    anfM = out.meanrate_MSR;
                    anfL = out.meanrate_LSR;
                    out_psth.psth_binwidth = psth_binwidth;

                    c = [];
                    c.an_summed = an_summed;
                    c.fs_an = fs_an;
                    c.anfH  = anfH;
                    c.anfM  = anfM;
                    c.anfL  = anfL;
                    c.out_psth = out_psth;
                    amt_cache('set',fname,c);
                end
                YL = [-20 780];
                
            case {'verhulst2015'}
                if bRun
                    hear_profile = 'Flat00'; % NH audiogram, default
                    fc_kv = {'hearing_profile',hear_profile};
                    out = verhulst2015(insig,fs,fc_flags{:},fc_kv{:},afb_flags{:},an_kv{:},an_flags{:},'anfL','anfM','anfH',nomfb_flags{:});
                    fs_an = out.fs_abr;

                    num_tot = numH+numM+numL;
                    an_summed = out.an_summed(:,idx_cf)/num_tot;
                    anfH = out.anfH;
                    anfM = out.anfM;
                    anfL = out.anfL;

                    c = [];
                    c.an_summed = an_summed;
                    c.fs_an = fs_an;
                    c.anfH  = anfH;
                    c.anfM  = anfM;
                    c.anfL  = anfL;
                    amt_cache('set',fname,c);
                end
                YL = [-20 780];
            case {'verhulst2018'}
                if bRun
                    %%% Model parameters (only one hearing profile, Flat00 and 13-3-3):
                    hear_profile = 'Flat00';
                    fc_kv = {'hearing_profile',hear_profile};
                    out = verhulst2018(insig,fs,fc_flags{:},fc_kv{:},afb_flags{:},'anfL','anfM','anfH',an_kv{:},an_flags{:},nomfb_flags{:});
                    fs_an = out.fs_abr;

                    num_tot = numH+numM+numL;
                    an_summed = out.an_summed(:,idx_cf)/num_tot;
                    anfH = out.anfH;
                    anfM = out.anfM;
                    anfL = out.anfL;

                    c = [];
                    c.an_summed = an_summed;
                    c.fs_an = fs_an;
                    c.anfH  = anfH;
                    c.anfM  = anfM;
                    c.anfL  = anfL;
                    amt_cache('set',fname,c);
                end
                YL = [-20 780];
            case {'bruce2018'}
                if bRun
                    kv = {'numH',numH,'numM',numM,'numL',numL,'psthbinwidth_mr',psth_binwidth, ...
                        'nrep',nrep,'specificSRautoTiming'};
                    out = bruce2018(insig,fs,fc,kv{:},afb_flags{:});
                    an_summed = out.meanrate;
                    anfH=out.meanrate_HSR; anfM = out.meanrate_MSR; anfL = out.meanrate_LSR;
                    out_psth.psth = out.psth;
                    out_psth.psth_HSR = out.psth_HSR;
                    out_psth.psth_MSR = out.psth_MSR;
                    out_psth.psth_LSR = out.psth_LSR;
                    out_psth.psth_binwidth = psth_binwidth;
                    fs_an = fs;

                    c = [];
                    c.an_summed = an_summed;
                    c.fs_an = fs_an;
                    c.anfH  = anfH;
                    c.anfM  = anfM;
                    c.anfL  = anfL;
                    c.out_psth = out_psth;
                    amt_cache('set',fname,c);
                end
                YL = [-20 780];
            case 'king2019'
                if bRun
                    fc_kv = {'flow',fc,'fhigh',fc,'basef',fc,'dboffset',dBFS};
                    out = king2019(insig,fs,fc_kv{:},afb_flags{:},afb_kv{:}, ...
                        ihc_flags{:},an_flags{:},nomfb_flags{:});
                    fs_an = fs;
                    an_summed = out;
                    anfH = out;
                    anfM = []; % 0*out;
                    anfL = []; % 0*out;

                    c = [];
                    c.fs_an = fs_an;
                    c.an_summed = an_summed;
                    c.anfH  = anfH;
                    c.anfM  = anfM;
                    c.anfL  = anfL;
                    amt_cache('set',fname,c);
                end

                factor = 1000;
                if factor == 1000
                    text4ylabel='(a.u. x 10^-3)';
                else
                    text4ylabel='(a.u.)';
                end
                YL = 1e-3*[-1.1 1.1];

            case 'relanoiborra2019'
                if bRun
                    fc_kv = {'flow',fc,'fhigh',fc,'basef',fc};
                    [~,~,out] = relanoiborra2019_featureextraction(insig, fs,fc_kv{:},afb_flags{:});
                    fs_an = fs;
                    an_summed = out;
                    anfH = out;
                    anfM = []; % 0*out;
                    anfL = []; % 0*out;

                    c = [];
                    c.an_summed = an_summed;
                    c.fs_an = fs_an;
                    c.anfH  = anfH;
                    c.anfM  = anfM;
                    c.anfL  = anfL;
                    c.out_psth = out;
                    amt_cache('set',fname,c);
                end
                text4ylabel='(MU)';
                YL = [-230 1450];

            case 'osses2021'
                if bRun
                    fc_kv = {'flow',fc,'fhigh',fc,'basef',fc,'dboffset',dBFS};
                    out = osses2021(insig,fs,fc_kv{:},afb_flags{:}, ...
                        ihc_flags{:},an_flags{:},nomfb_flags{:});
                    fs_an = fs;
                    an_summed = out;
                    anfH = out;
                    anfM = []; % 0*out;
                    anfL = []; % 0*out;

                    c = [];
                    c.an_summed = an_summed;
                    c.fs_an = fs_an;
                    c.anfH  = anfH;
                    c.anfM  = anfM;
                    c.anfL  = anfL;
                    amt_cache('set',fname,c);
                end
                text4ylabel='(MU)';
                YL = [-230 1450];
        end

        t_anf = (1:length(an_summed))'/fs_an;

        L_bin = round(psth_bin4plot*fs_an); % samples
        L_overlap = round((percent4plot/100)*L_bin);
        t_anf = buffer(t_anf, L_bin, L_overlap,'nodelay');
        t_anf = t_anf(1,:);

        anf = buffer(an_summed, L_bin, L_overlap,'nodelay');
        anf = mean(anf); % one value per bin

        if flags.do_plot
            offx = 2.5;
            if flags.do_fig12a
                XL = [0-offx 50+offx];
            elseif flags.do_fig12b
                XL = [350-offx 400+offx];
            end

            switch models{k}
                case {'dau1997','relanoiborra2019','osses2021'}
                    if k == 1
                        fig_functional = figure;
                        figure_handle_ax(end+1) = gca;
                        figure_handle(end+1) = gcf;
                        figure_name{end+1}   = ['fig12-Adaptation-' models{k}];
                        hold on, grid on;
                    else
                        figure(fig_functional);
                    end
                    if strcmp(models{k},'osses2021')
                        Style = '--';
                    else
                        Style = '-';
                    end

                case {'zilany2014','bruce2018'}
                    figure;

                    figure_handle_ax(end+1) = gca;
                    figure_handle(end+1) = gcf;
                    figure_name{end+1}   = ['fig12-Adaptation-' models{k}];
                    hold on, grid on;
                    Style = '-';

                case {'verhulst2015','verhulst2018'}
                    if ~exist('fig_verhulst','var')
                        fig_verhulst = figure;

                        figure_handle_ax(end+1) = gca;
                        figure_handle(end+1) = gcf;
                        figure_name{end+1}   = ['fig12-Adaptation-' models{k}];
                        hold on, grid on;
                    else
                        figure(fig_verhulst);
                    end
                    Style = '-';

                case 'king2019'
                    figure;
                    figure_handle_ax(end+1) = gca;
                    figure_handle(end+1) = gcf;
                    figure_name{end+1}   = ['fig12-Adaptation-' models{k}];
                    hold on, grid on;

                    Style = '-';
            end
            plot(t_anf*1000,factor*anf,Style,'Color',Colours{k},'LineWidth',2); grid on, hold on
        end

        switch models{k}
            case {'zilany2014','bruce2018'}
                psth_count = out_psth.psth;
                t_psth = (1:length(psth_count))*out_psth.psth_binwidth;
                if flags.do_plot
                    stairs(t_psth*1000,psth_count,Style,'Color',0.5*Colours{k},'LineWidth',2); grid on, hold on
                end

                outs_anf(k).psth_count = psth_count;
                outs_anf(k).t_psth = t_psth;

                if flags.do_plot
                    % Trick to add labels later:
                    plot(t_psth(1)*1000,psth_count(1),'Color',0.5*Colours{k},'LineWidth',2);
                end
        end
            
        if flags.do_plot
            xlim(XL);
            ylim(factor*YL);
            xlabel('Time (ms)');
            ylabel(text4ylabel);

            Pos = get(gcf,'Position');
            Pos(3:4) = [360 400];
            set(gcf,'Position',Pos);
        end


        outs_anf(k).t_anf = t_anf;
        outs_anf(k).fs_an = fs_an;
        outs_anf(k).anf  = anf;
        outs_anf(k).anfH = anfH;
        outs_anf(k).anfM = anfM;
        outs_anf(k).anfL = anfL;
    end
    % hl = legend(text4leg,'Location','SouthEast');
    % set(hl,'FontSize',8);

    data.models = models;
    data.insig = insig;
    data.fs    = fs;
    data.outs_anf = outs_anf;
    if flags.do_plot
        if flags.do_fig12a
            data.figure_flag   = 'do_fig12a';
        elseif flags.do_fig12b
            data.figure_flag   = 'do_fig12b';
        end
        data.figure_handle = figure_handle;
        data.figure_name   = figure_name;
        data.figure_handle_ax = figure_handle_ax;
    end
end

%% ------ FIG 13 Osses et al. 2021 ----------------------------------------
if flags.do_fig13

    dur = 300e-3;
    rampdur = 10e-3;
    lvl = 50;

    basef = 1000;
    baseaud = freqtoaud(basef);

    N_freqs     = 7;
    erb_step    = 3;
    freqs       =  audtofreq(baseaud-N_freqs+1:1:baseaud);
    freqs4tones = freqs(1:erb_step:end);
    starting_phases = [0.7145 4.2320 2.2943]; % Three fixed (frozen) phases, without using fixed seeds
    [insig,f2test,starting_phases] =  local_Profile_Analysis(dur, rampdur, lvl, fs,freqs4tones,starting_phases);
    insig = insig(:);

    ti = 220; % ms
    tf = 260;
    ti_e = ti+10;
    %%%

    for k = 1:length(models)
        out = [];
        suff_str = [];

        %%% Loading flags and keyvals
        [fg,kv] = local_get_flags(models{k});
        fc_flags = fg.fc_flags;
        afb_flags = fg.afb_flags;
        afb_keyvals = kv.afb_keyvals;
        ihc_flags = fg.ihc_flags;
        an_flags = fg.an_flags;
        nomfb_flags = fg.nomfb_flags;
        an_kv    = kv.an_keyvals;
        %%%

        fname = ['fig13_profile_' models{k}];
        c = amt_cache('get',fname,flags.cachemode);

        if ~isempty(c)
            bRun  = 0;

            out   = c.out;
            fc    = c.fc;
            fs_an = c.fs_an;
        else
            bRun = 1;
            out = [];
            amt_disp(['Calculating ' models{k} '...']);
        end

        fc_green = local_m2hz(401);

        for j = 1:N_freqs
            bin_nrs(j) = find(fc_green>freqs(j),1,'last');
        end
        % bin_nrs = bin_nrs(1:erb_step:end);
        CFs = fc_green(bin_nrs);

        if ~isempty(afb_keyvals)
            idx = find(strcmp(afb_keyvals,'compression_n'));
            if ~isempty(idx)
                compression_n = afb_keyvals{idx+1};
            else
                error('mfb_keyvals should contain BMF')
            end
        end
        if bRun && ~isempty(an_kv)
            % Loading common parameters
            switch models{k}
                case {'verhulst2015','verhulst2018','zilany2014','bruce2018'}
                    if strcmp(an_kv{1},'numH')
                        numH = an_kv{2};
                    else
                        error('an_kv should contain numH')
                    end
             
                    if strcmp(an_kv{3},'numM')
                        numM = an_kv{4};
                    else
                        error('an_kv should contain numM')
                    end
             
                    if strcmp(an_kv{5},'numL')
                        numL = an_kv{6};
                    else
                        error('an_kv should contain numM')
                    end
            end 
            switch models{k}
                case 'bruce2018'
                    idx = find(strcmp(an_kv,'psth_binwidth'));
                    if ~isempty(idx)
                        % Only true for bruce2018:
                        psth_binwidth = an_kv{idx+1};
                    end
                    nrep = 100;
            end

        end
        
        extra_after = '';
        switch models{k}
            case 'dau1997'
                if bRun
                    for j = 1:length(CFs)
                        fc_kv = {'basef',CFs(j),'flow',CFs(j),'fhigh',CFs(j),'dboffset',dBFS};
                        out_tmp = dau1997(insig,fs,afb_flags{:},ihc_flags{:}, ...
                            an_flags{:},nomfb_flags{:},fc_kv{:});
                        out(:,j) = out_tmp;
                        fc(j) = freqs(j);
                    end

                    fs_an = fs;

                    c = [];
                    c.out   = out;
                    c.fc    = fc;
                    c.fs_an = fs_an;
                    amt_cache('set',fname,c);
                end
                fc_ref = fc;
                ZL = [-240 1000];
                ymax = 1600;
                yscale = ymax/2;
                ymax_str = sprintf('%.0f MU',yscale);
                extra_str = 'a) ';

            case 'zilany2014'
                if bRun
                    % Only one repetition because we will use the mean-rate output.
                    kv={'fiberType',4,'numH',numH,'numM',numM,'numL',numL,'nrep',1};
                    out = zilany2014(insig,fs,fc_ref,kv{:},afb_flags{:});
                    out = out(1:length(insig),:);

                    fs_an = fs;

                    c = [];
                    c.out   = out;
                    c.out_description = 'AN mean-rate output';
                    c.fc    = fc_ref;
                    c.fs_an = fs_an;
                    amt_cache('set',fname,c);
                end

                ymax = 1000;
                yscale = ymax/2;
                ymax_str = sprintf('%.0f spikes/s',yscale);
                extra_str = 'c) ';
                extra_after = ', mean rate';

            case 'verhulst2015'
                if bRun
                    out_tmp = verhulst2015(insig,fs,fc_flags{:},afb_flags{:},an_flags{:},an_kv{:},nomfb_flags{:});

                    fc = fc_green(bin_nrs);
                    fs_an = out_tmp.fs_an;

                    numTot = numL+numM+numH;
                    out = out_tmp.an_summed(:,bin_nrs)/numTot;

                    c = [];
                    c.out   = out;
                    c.fc    = fc;
                    c.fs_an = fs_an;
                    amt_cache('set',fname,c);
                end

                extra_str = 'e) ';
                ymax = 1000;
                yscale = ymax/2;
                ymax_str = sprintf('%.0f spikes/s',yscale);

            case 'verhulst2018'
                if bRun
                    out_tmp = verhulst2018(insig,fs,fc_flags{:},afb_flags{:},an_flags{:},an_kv{:},nomfb_flags{:});

                    fc = fc_green(bin_nrs);
                    fs_an = out_tmp.fs_an;
                    
                    numTot = numL+numM+numH;
                    out = out_tmp.an_summed(:,bin_nrs)/numTot;

                    c = [];
                    c.out   = out;
                    c.fc    = fc;
                    c.fs_an = fs_an;
                    amt_cache('set',fname,c);
                end

                extra_str = 'f) ';
                ymax = 1000;
                yscale = ymax/2;
                ymax_str = sprintf('%.0f spikes/s',yscale);

            case 'bruce2018'

                bUse_PSTH = 1;
                if bRun
                    fc = fc_ref;
                    kv = {'numH',numH,'numM',numM,'numL',numL, ...
                        'psthbinwidth_mr',psth_binwidth,'nrep',nrep,'specificSRautoTiming'};
                    for j = 1:length(fc_ref)
                        output = bruce2018(insig,fs,fc_ref(j),kv{:},afb_flags{:});
                        out(:,j) = output.psth;
                    end
                    dt = psth_binwidth;
                    fs_an = 1/dt;
                    Nf = round(dur*fs_an);
                    out = out(1:Nf,:);

                    c = [];
                    c.out   = out;
                    c.out_description = 'PSTH output';
                    c.fc    = fc;
                    c.fs_an = fs_an;
                    amt_cache('set',fname,c);
                end

                if bUse_PSTH == 0
                    ymax = 300;
                else
                    ymax = 1000;
                end
                yscale = ymax/2;
                ymax_str = sprintf('%.0f spikes/s',yscale);
                extra_str = 'd) ';
                extra_after = ', PSTH';

            case 'king2019'

                switch compression_n
                    case 0.3
                        ymax = 2.5e-3; % all channels compressed
                end

                if bRun
                    fc = [];
                    for j = 1:length(fc_ref)
                        fc_kv = {'basef',fc_ref(j),'flow',fc_ref(j),'fhigh',fc_ref(j),'dboffset',dBFS};

                        [out_tmp,fc(j)] = king2019(insig,fs,afb_flags{:},ihc_flags{:}, ...
                            an_flags{:},nomfb_flags{:},fc_kv{:},'compression_n',compression_n);
                        out(:,j) = out_tmp;
                    end
                    fs_an = fs;

                    c = [];
                    c.out   = out;
                    c.fc    = fc;
                    c.fs_an = fs_an;
                    amt_cache('set',fname,c);
                end

                extra_str = 'b) ';
                yscale = ymax/2;
                ymax_str = sprintf('%.4f a.u',yscale);

            case 'relanoiborra2019'
                if bRun
                  fc = [];
                  for j = 1:length(fc_ref)
                    fc_kv = {'basef',fc_ref(j),'flow',fc_ref(j),'fhigh',fc_ref(j),'dboffset',dBFS};
                    [~,~,out_tmp,fc(j)] = relanoiborra2019_featureextraction(insig, fs, ...
                        fc_kv{:}, afb_flags{:});
                    out(:,j) = out_tmp;
                  end
                  % out = out(:,1:2*erb_step:end);
                  % fc  = fc(1:2*erb_step:end);
                  fs_an = fs;

                  c = [];
                  c.out   = out;
                  c.fc    = fc;
                  c.fs_an = fs_an;
                  amt_cache('set',fname,c);
                end

                ZL = [-240 1000];
                ymax = 1600;
                yscale = ymax/2;
                ymax_str = sprintf('%.0f MU',yscale);
                extra_str = 'g) ';

            case 'osses2021'
                if bRun
                    for j = 1:length(fc_ref)
                        fc_kv = {'basef',fc_ref(j),'flow',fc_ref(j),'fhigh',fc_ref(j),'dboffset',dBFS};
                        out_tmp = osses2021(insig,fs,afb_flags{:},ihc_flags{:}, ...
                            an_flags{:},nomfb_flags{:},fc_kv{:});
                        out(:,j) = out_tmp;
                        fc(j) = fc_ref(j);
                    end

                    fs_an = fs;

                    c = [];
                    c.out   = out;
                    c.fc    = fc;
                    c.fs_an = fs_an;
                    amt_cache('set',fname,c);
                end

                extra_str = 'h) ';
                ZL = [-240 1000];
                ymax = 1600;
                yscale = ymax/2;
                ymax_str = sprintf('%.0f MU',yscale);
        end

        idxi_=round(ti*1e-3*fs_an);
        idxf_=round(tf*1e-3*fs_an);

        offy = [];

        t_ms   = 1000*(1:size(insig,1))/fs;
        tan_ms = 1000*(1:size(out,1))/fs_an;

        bw_erb = audfiltbw(fc);

        if flags.do_plot
            figure;
        end
        for i = 1:length(fc)

            me = mean(out(idxi_:idxf_,i)); % looks for peaks above the avg in the steady section
            [env_here,idx_here] = local_Get_envelope2(out(:,i),me);

            tan_ms_here = tan_ms(idx_here);
            idxs = find(tan_ms >= ti-10 & tan_ms <= tf-7.5); % 10 ms before ti
            
            if flags.do_plot
                offy(i) = (i-1);
                plot(tan_ms(idxs),out(idxs,i)/ymax + offy(i),'Color',Colours{k}); hold on, grid on
            end

            idxs = find(tan_ms_here >= ti-10 & tan_ms_here <= tf-7.5);

            if flags.do_plot
                plot(tan_ms_here(idxs),env_here(idxs)/ymax + offy(i),'Color',[0.6 0.6 0.6],'LineWidth',2);
            end

            env_mean(k,i) = mean(env_here(idxs));
            env_std(k,i)  = std(env_here(idxs));
            env_scale(k) = yscale;            
        end

        factor_std = 30; % scaling (arbitrary) applied to the std of each model/channel
        if flags.do_plot
            %%% YAxis: scale
            tf_here=tf-7.5-.5;
            plot(tf_here*[1 1],[0 0.5],'k-','LineWidth',2);
            plot([tf_here-2 tf_here],0.5*[1 1],'k-','LineWidth',2);
            plot([tf_here-2 tf_here],    [0 0],'k-','LineWidth',2);
            text(tf_here-6,-0.35,ymax_str,'FontWeight','Bold');

            %%%
            ylabel('Simulated CF_n (Hz)')
            xlim([ti-2.5 tf]);
            ylim([-.5 7.5])
            YL = get(gca,'YLim');

            tf_i_here = tf+2.5-10; % '0' for envelope metric (vertical black line)

            env_std_here = env_std(k,:)/(env_scale(k)/factor_std);
            ampl_here = env_mean(k,:)/ymax + offy;
            plot(tf_i_here+env_std_here, ampl_here,'o--','Color',local_rgb('Maroon'),'LineWidth',2,'MarkerFaceColor','w');
            plot(tf_i_here*[1 1],YL,'k-');

            tf_f_here = tf_i_here + 2*3;
            y_here = 6.8;
            plot([tf_i_here tf_f_here],y_here*[1 1],'-','Color',local_rgb('Maroon'),'LineWidth',2);

            amp_val = 0;
            for ii=1:4
                tf_f_here = tf_i_here + amp_val;
                plot(tf_f_here*[1 1],[y_here-.2 y_here],'-','Color',local_rgb('Maroon'),'LineWidth',2);
                % text(tf_f_here,y_here+.1,num2str(amp_val),'FontSize',8);
                amp_val = amp_val+2;
            end
        
            if k == 1
                x_fc   = (1:length(fc));

                YTL = [];
                for j = 1:length(CFs)
                    YTL{j} = [num2str(round(CFs(j))) ' Hz'];
                end
            end
            set(gca,'YTick',x_fc-1);
            set(gca,'YTickLabel',YTL);

            set(gca,'XTick',ti:10:tf-10);
            % set(gca,'XTickLabel',[])

            ht = text(0.15,0.97,[extra_str models{k} extra_after],'Units','Normalized','FontSize',MediumFontSize);

            Pos = get(gcf,'Position');
            Pos(3:4) = [360 300]; % Previous size: [500 300];
            set(gcf,'Position',Pos);

            figure_handle(end+1) = gcf;
            figure_name{end+1} = ['fig13-profile-' models{k} suff_str];
        end
    end
    
    data.env_mean = env_mean;
    data.env_std  = env_std;
    data.env_scale = env_scale;
    data.env_scale_std = env_scale/factor_std;

    data.insig = insig;
    data.starting_phases = starting_phases;
    data.fs    = fs;
    data.models = models;
    
    if flags.do_plot
        data.figure_flag = 'do_fig13';
        data.figure_handle = figure_handle;
        data.figure_name   = figure_name;
    end
end

%% ------ FIG 14 ----------------------------------------------------------
if flags.do_fig14a || flags.do_fig14b   
    % Adapted from g20190618_investigating_IC_CN.m
    %%% Stimulus generation:
    if flags.do_fig14a
        L=30; % dB SPL
    elseif flags.do_fig14b
        L=70; % dB SPL
    end
    dur= 300e-3; % ms
    t  = 0:1/fs:dur-1/fs;
    t  = t(:); % column vector
    fc = 1000; % Hz
    dur_samples = length(t);

    fmod_step = 5;
    fmods = 10:fmod_step:130; % 250;
    N_signals = length(fmods);

    insig = zeros(dur_samples,N_signals); % Memory allocation

    % Up/down cosine ramp (fixed)
    dur_ramp_ms = 5;
    dur_ramp = round((dur_ramp_ms*1e-3)*fs); % duration ramp in samples

    rp    = ones(dur_samples,1);
    rp(1:dur_ramp) = rampup(dur_ramp);
    rp(end-dur_ramp+1:end) = rampdown(dur_ramp);

    %Stimulus
    for i = 1:length(fmods)
        fmod = fmods(i);
        mod_index=1;
        carrier=sin(2*pi*fc.*t);
        modulator=mod_index*cos(2*pi*fmod.*t+pi);
        insig_tmp=(1+modulator).*carrier;

        insig_tmp = scaletodbspl(insig_tmp,L,dBFS);

        insig(:,i) = rp.*insig_tmp;
    end

    ti = 190e-3;
    tf = 290e-3;
    for k = 1:N_models

        %%% Loading flags and keyvals
        [fg,kv] = local_get_flags(models{k});
        fc_flags = fg.fc_flags;
        afb_flags = fg.afb_flags;
        afb_keyvals = kv.afb_keyvals;
        ihc_flags = fg.ihc_flags;
        an_flags = fg.an_flags;
        an_keyvals = kv.an_keyvals;
        mfb_flags = fg.mfb_flags;
        mfb_keyvals = kv.mfb_keyvals;
        %%%

        switch L
            case 70
                suff_lvl='';
            case 30
                suff_lvl = '-30-dB';
            otherwise
                error('This function has been tested with levels of L=30 or 70 dB SPL only')
                % suff_lvl = sprintf('-%.0f-dB',L);
        end
        fname = ['fig14_MTF-' models{k} suff_lvl];
        c = amt_cache('get',fname,flags.cachemode);

        if ~isempty(c)
            bRun  = 0;

            outsig = c.outsig;
            subfs  = c.subfs;
            mfc    = c.mfc;
        else
            bRun = 1;
            outsig = [];
        end

        mfc_target = 100;
        if bRun && ~isempty(afb_keyvals)
            idx = find(strcmp(afb_keyvals,'compression_n'));
            if ~isempty(idx)
                compression_n = afb_keyvals{idx+1};
            else
                error('mfb_keyvals should contain BMF')
            end
        end
        if bRun && ~isempty(an_keyvals)
            % Loading common parameters
            switch models{k}
                case {'verhulst2015','verhulst2018','zilany2014','bruce2018'}
                    if strcmp(an_keyvals{1},'numH')
                        numH = an_keyvals{2};
                    else
                        error('an_keyvals should contain numH')
                    end
             
                    if strcmp(an_keyvals{3},'numM')
                        numM = an_keyvals{4};
                    else
                        error('an_keyvals should contain numM')
                    end
             
                    if strcmp(an_keyvals{5},'numL')
                        numL = an_keyvals{6};
                    else
                        error('an_keyvals should contain numM')
                    end
            end
             
            switch models{k}
                case {'zilany2014','bruce2018'}
                    if ~isempty(mfb_keyvals)
                        idx = find(strcmp(mfb_keyvals,'BMF'));
                        if ~isempty(idx)
                            BMF = mfb_keyvals{idx+1};
                            mfb_keyvals(idx:idx+1) = []; % removed from the keyvals
                        else
                            error('mfb_keyvals should contain BMF')
                        end
                    end
            
                    switch models{k}
                        case 'bruce2018'
                            idx = find(strcmp(an_keyvals,'psth_binwidth'));
                            if ~isempty(idx)
                                % Only true for bruce2018:
                                psth_binwidth = an_keyvals{idx+1};
                            end
                            nrep = 100;
                    end
            end
        end
        switch models{k}
            case 'dau1997'
                if bRun
                    afb_kv = {'basef',fc,'flow',fc,'fhigh',fc,'dboffset',dBFS};
                    for j = 1:N_signals
                        [out_tmp,fc,mfc_here] = dau1997(insig(:,j),fs,afb_kv{:},afb_flags{:}, ...
                            ihc_flags{:},an_flags{:},mfb_flags{:},mfb_keyvals{:});
                        mfc_idx = find(mfc_here<mfc_target,1,'last');
                        outsig(:,j) = out_tmp{1}(:,mfc_idx);

                        if j == 1
                            mfc = mfc_here(mfc_idx);

                            if size(out_tmp{1},1) == size(insig,1)
                                % Makers sure that no other default subfs is used
                                subfs = fs;
                            end
                        end
                    end

                    c = [];
                    c.outsig = outsig;
                    c.subfs  = subfs;
                    c.mfc    = mfc;
                    amt_cache('set',fname,c);
                end

            case 'zilany2014'
                if bRun
                    insigL = size(insig,1);
                    kv={'fiberType',4,'numH',numH,'numM',numM,'numL',numL,'nrep',1}; % 1 repetition is enough (no PSTH here)
                    for j = 1:N_signals
                        mean_rate = zilany2014(insig(:,j),fs,fc,kv{:},afb_flags{:});

%                         if j == 1
%                             L = size(insig,1);
%                         end
                        mean_rate = mean_rate(1:insigL);
                        [out_ic,~,~,par] = carney2015(mean_rate,BMF,fs,mfb_keyvals{:});
                        outsig(:,j) = out_ic(1:insigL); % column array

                        if j == 1
                            subfs = fs;
                            mfc = BMF;
                            amt_disp(sprintf('Zilany2014: CN: tau_exc=%f, tan_inh=%f, D=%f, S=%f',par.tau_ex_cn,par.tau_inh_cn,par.cn_delay,par.Sinh_cn));
                            amt_disp(sprintf('Zilany2014: IC: tau_exc=%f, tan_inh=%f, D=%f, S=%f',par.tau_ex_ic,par.tau_inh_ic,par.ic_delay_inh,par.Sinh_ic));
                        end
                    end

                    c = [];
                    c.outsig = outsig;
                    c.subfs  = subfs;
                    c.mfc    = mfc;
                    amt_cache('set',fname,c);
                end

            case 'verhulst2015'
                if bRun
                    fc_green = local_m2hz(401);
                    idx = find(fc_green>fc,1,'last');

                    out = verhulst2015(insig,fs,fc_flags{:},afb_flags{:},ihc_flags{:},an_flags{:}, ...
                        kv.an_keyvals{:},mfb_flags{:}); % no_ihc just reduce the data load here
                    subfs = out(1).fs_abr;

                    %%%
                    numTot = numH + numM + numL;
                    %%%

                    for j = 1:N_signals
                        outsig(:,j) = out(j).ic(:,idx)/numTot;
                    end
                    mfc = [];

                    c = [];
                    c.outsig = outsig;
                    c.subfs  = subfs;
                    c.mfc    = mfc;
                    amt_cache('set',fname,c);
                end

            case 'verhulst2018'
                if bRun
                    fc_green = local_m2hz(401);
                    idx = find(fc_green>fc,1,'last');

                    out = verhulst2018(insig,fs,fc_flags{:},afb_flags{:},ihc_flags{:},an_flags{:}, ...
                        kv.an_keyvals{:},mfb_flags{:},'no_ihc','no_an','no_cn');
                    subfs = out(1).fs_abr;

                    %%%
                    numTot = numH + numM + numL;
                    %%%
                    
                    for j = 1:N_signals
                        outsig(:,j) = out(j).ic(:,idx)/numTot;
                    end
                    mfc = [];

                    c = [];
                    c.outsig = outsig;
                    c.subfs  = subfs;
                    c.mfc    = mfc;
                    amt_cache('set',fname,c);
                end

            case 'bruce2018'
                if bRun
                    fc_kv = {};
                    kv = {'numH',numH,'numM',numM,'numL',numL,'nrep',nrep, ...
                        'specificSRautoTiming','psthbinwidth_mr',psth_binwidth};
                    for j = 1:N_signals

                        out_tmp = bruce2018(insig(:,j),fs,fc,fc_kv{:},kv{:},afb_flags{:});
                        if j == 1
                            dt = psth_binwidth;
                            fs_an = 1/dt;
                            insigL = round(dur*fs_an);
                        end
                        out_psth = out_tmp.psth(1:insigL);

                        outsig(:,j) = carney2015(out_psth,BMF,fs_an,mfb_keyvals{:});

                        if j == 1
                            subfs = fs_an;
                            mfc = BMF;
                        end
                    end

                    c = [];
                    c.outsig = outsig;
                    c.subfs  = subfs;
                    c.mfc    = mfc;
                    amt_cache('set',fname,c);
                end

            case 'relanoiborra2019'
                if bRun
                    afb_kv = {'basef',fc,'flow',fc,'fhigh',fc};
                    for j = 1:N_signals
                        [out_tmp,mfc_here,~,~] = relanoiborra2019_featureextraction(insig(:,j), fs, afb_kv{:}, afb_flags{:});

                        mfc_idx = find(mfc_here<mfc_target,1,'last');
                        outsig(:,j) = squeeze(out_tmp(:,1,mfc_idx));
                        subfs = fs;
                        if j == 1
                            mfc = mfc_here(mfc_idx);
                        end
                    end
                    c = [];
                    c.outsig = outsig;
                    c.subfs  = subfs;
                    c.mfc    = mfc;
                    amt_cache('set',fname,c);
                end

            case 'king2019'
                if bRun
                    for j = 1:N_signals
                        afb_kv = {'basef',fc,'flow',fc,'fhigh',fc,'compression_n',compression_n,'dboffset',dBFS,'subfs',fs};
                        [out_tmp,fc,mfc_here,extras] = king2019(insig(:,j),fs,afb_kv{:});
                        mfc_idx = find(mfc_here<mfc_target,1,'last');

                        out_tmp = squeeze(out_tmp);
                        mfc = mfc_here(mfc_idx);
                        outsig(:,j) = out_tmp(:,mfc_idx);
                        if j == 1
                            subfs = extras.subfs;
                        end
                    end
                    c = [];
                    c.outsig = outsig;
                    c.subfs  = subfs;
                    c.mfc    = mfc;
                    amt_cache('set',fname,c);
                end

            case 'osses2021'
                if bRun
                    afb_kv = {'basef',fc,'flow',fc,'fhigh',fc,'dboffset',dBFS};
                    for j = 1:N_signals
                        [out_tmp,fc,mfc_here] = osses2021(insig(:,j),fs,afb_kv{:},afb_flags{:},ihc_flags{:},an_flags{:},mfb_flags{:});
                        mfc_idx = find(mfc_here<mfc_target,1,'last');
                        outsig(:,j) = out_tmp{1}(:,mfc_idx);
                        subfs = fs;
                        if j == 1
                            mfc = mfc_here(mfc_idx);
                        end
                    end

                    c = [];
                    c.outsig = outsig;
                    c.subfs  = subfs;
                    c.mfc    = mfc;
                    amt_cache('set',fname,c);
                end

        end

        if ~isempty(outsig)

            switch models{k}
                case {'dau1997','relanoiborra2019','osses2021'}
                    YL = [-100 1000];
                case {'zilany2014','verhulst2015','verhulst2018','bruce2018'}
                    YL = [-200 500];
                case 'king2019'
                    YL = [-.1 1]*1e-4;
            end

            idxi = round(ti*subfs)+1;
            idxf = round(tf*subfs);

            % Maximum values:
            [vals_ma,idx_ma] = max(outsig(idxi:idxf,:));

            MFic(k,1:N_signals)      = vals_ma; % figure; plot(fmods,vals); hold on
            MFic_norm(k,1:N_signals) = vals_ma / max(abs(vals_ma));

            switch L
                case 30
                    MF_label = 'a) AM tones at 30 dB';
                case 70
                    MF_label = 'b) AM tones at 70 dB';
            end
        end
    end

    if flags.do_plot
        %%% Custom figure format here:
        for k = 1:N_models % little adjustments in format for the coming figures
            switch models{k}
                case 'dau1997'
                    Marker = [Markers{k} '-'];
                    LW = LineWidth(k);
                case 'zilany2014'
                    Marker = [Markers{k} LineStyle{k}];
                    LW = LineWidth(k);
                case 'bruce2018'
                    Marker = LineStyle{k};
                    LW = 4;
                case 'king2019'
                    Marker = LineStyle{k};
                    LW = 4;
                otherwise
                    Marker = LineStyle{k};
                    LW = LineWidth(k);
            end

            if k == 1
                figure;
            end
            opts = {Marker,'Color',Colours{k},'LineWidth',LW,'MarkerFaceColor',Colours{k}}; % ,'MarkerFaceColor','w','Markers',Markers{k},'MarkerSize',MarkersSize(k)}
            plot(fmods,MFic_norm(k,:),opts{:}); grid on, hold on
        end

        xlabel('Modulation frequency (Hz)')
        ylabel('On-CF response (Normalised)')
        xlim([min(fmods)-5 max(fmods)+5]);
        ylim([-0.05 1.15])

        Pos = get(gcf,'Position');
        Pos(3:4) = [400 360]; % [550 360];
        set(gcf,'Position',Pos);

        XTL = [];
        for i = 1:length(fmods)
            if mod(i,2) == 1
                XTL{i} = num2str(fmods(i));
            else
                XTL{i} = '';
            end
        end

        set(gca,'XTick'     ,fmods); 
        set(gca,'XTickLabel',XTL); 

        set(gca,'YTick',-1:.1:1);
        figure_handle(end+1) = gcf;
        figure_name{end+1} = ['fig14-modulation-strength' suff_lvl];

        title(MF_label);
    end
    
    %%% End: figure

    data.figure_flag = 'do_fig14';
    data.fmods = fmods;
    data.MFic = MFic;
    data.MFic_norm = MFic_norm;
    data.models = models;
    if flags.do_plot
        data.figure_handle = figure_handle;
        data.figure_name   = figure_name;
    end
end

%% ------ FIG 15 ----------------------------------------------------------
if flags.do_fig15
    % Local script: g20191107_rerun_Verhulst2018a.m

    % Generating click
    lvl = 70; % lvl_dBnHL + 30;
    p0  = 2e-5;
    dur_click = 1; % s
    N_samples = dur_click*fs;
    
    t=(0:1/fs:dur_click);
    click_duration = (100e-6*fs); % 100 us click (Osses and Verhulst 2019 used 80 us)

    N_length = fs/10;
    sil_skip = N_length-click_duration; % 90.1 ms

    samples_click      = [];
    samples_click_even = [];

    Nr_clicks = floor(fs*dur_click/(sil_skip+click_duration));

    idx_offset = round(10e-3*fs); % click starts 10 ms after time 0
    for i = 1:Nr_clicks
        start_sample = (i-1)*N_length+idx_offset;
        if mod(i,2) == 1
            samples_click      = [samples_click start_sample+click_duration:start_sample+click_duration+click_duration];
        else
            samples_click_even = [samples_click_even start_sample+click_duration:start_sample+click_duration+click_duration];
        end
    end

    A = 2*sqrt(2)*p0*10^(lvl/20);
    insig = zeros(N_samples,1);
    insig(samples_click)      =  A;
    insig(samples_click_even) = -A;

    for k = 1:N_models
        %%% Loading flags and keyvals
        [fg,kv] = local_get_flags(models{k});
        fc_flags = fg.fc_flags;
        afb_flags = fg.afb_flags;
        ihc_flags = fg.ihc_flags;
        an_flags = fg.an_flags;
        afb_keyvals = kv.afb_keyvals;
        an_keyvals = kv.an_keyvals;
        mfb_flags = fg.mfb_flags;
        % nomfb_flags = fg.nomfb_flags;
        mfb_keyvals = kv.mfb_keyvals;
        %%%

        fname = ['fig15_click-' models{k}];
        c = amt_cache('get',fname,flags.cachemode);

        if ~isempty(c)
            bRun  = 0;

            out   = c.out;
            subfs = c.subfs;
            mfc   = c.mfc;
        else
            bRun = 1;
            out = [];
        end

        mfc_target = 100;
        if bRun && ~isempty(afb_keyvals)
            idx = find(strcmp(afb_keyvals,'compression_n'));
            if ~isempty(idx)
                compression_n = afb_keyvals{idx+1};
            else
                error('mfb_keyvals should contain BMF')
            end
        end
        if bRun && ~isempty(an_keyvals)
            % Loading common parameters
            switch models{k}
                case {'verhulst2015','verhulst2018','zilany2014','bruce2018'}
                    if strcmp(an_keyvals{1},'numH')
                        numH = an_keyvals{2};
                    else
                        error('an_keyvals should contain numH')
                    end
                    
                    if strcmp(an_keyvals{3},'numM')
                        numM = an_keyvals{4};
                    else
                        error('an_keyvals should contain numM')
                    end
                    
                    if strcmp(an_keyvals{5},'numL')
                        numL = an_keyvals{6};
                    else
                        error('an_keyvals should contain numM')
                    end
            end
            
            switch models{k}
                case {'zilany2014','bruce2018'}
                    if ~isempty(mfb_keyvals)
                        idx = find(strcmp(mfb_keyvals,'BMF'));
                        if ~isempty(idx)
                            BMF = mfb_keyvals{idx+1};
                            mfb_keyvals(idx:idx+1) = []; % removed from the keyvals
                        else
                            error('mfb_keyvals should contain BMF')
                        end
                    end
                    
                    switch models{k}
                        case 'bruce2018'
                            idx = find(strcmp(an_keyvals,'psth_binwidth'));
                            if ~isempty(idx)
                                % Only true for bruce2018:
                                psth_binwidth = an_keyvals{idx+1};
                            end
                            nrep = 100;
                    end
                    idx_step = 8; % used in zilany2014 and bruce2018
            end
        end
        
        switch models{k}
            case 'dau1997'
                if bRun
                    subfs = 20000; % Hz, fs to be used in the modulation filter bank
                    afb_kv = {'subfs',subfs,'dboffset',dBFS,'basef',[]};

                    t_start = tic;
                    [out_tmp,fc_here,mfc_here] = dau1997(insig,fs,afb_kv{:},afb_flags{:},ihc_flags{:},an_flags{:},mfb_flags{:});
                    out.t_elapsed = toc(t_start);

                    mfc_idx = find(mfc_here<mfc_target,1,'last');

                    Nr_fc = length(fc_here);
                    outsig = [];

                    fcs = [];
                    for j = 1:Nr_fc
                        if size(out_tmp{j},2) >= mfc_idx
                            fcs(end+1) = fc_here(j); % only adding the fcs of filters containing mfc_idx
                            outsig(:,end+1) = out_tmp{j}(:,mfc_idx);
                        end
                    end

                    out.fcs    = fcs;
                    out.oustig = outsig;
                    out.outsig_description = 'Modulation filter bank output';

                    mfc = mfc_here(mfc_idx);

                    c = [];
                    c.out   = out;
                    c.subfs = subfs;
                    c.mfc   = mfc;
                    amt_cache('set',fname,c);
                end

                factor = 1;
                unit = 'MU';
                YL  = [-20 560];
                YT = 0:50:500;

            case 'zilany2014'
                if bRun
                    fcs = local_m2hz(401);
                    idxs = length(fcs):-idx_step:1;
                    fcs = fcs(idxs);
                    fcs = fcs(fcs>125);

                    kv={'fiberType',4,'numH',numH,'numM',numM,'numL',numL,'nrep',1};

                    for j = 1:length(fcs)
                        amt_disp(sprintf('zilany2014: band %.0f of %.0f\n',j,length(fcs)),'volatile');
                        t_start = tic;
                        mean_rate = zilany2014(insig,fs,fcs(j),afb_flags{:},kv{:});
                        if j == 1
                            L = size(insig,1);
                        end
                        mean_rate = mean_rate(1:L);
                        [ic_out,~,~,par] = carney2015(mean_rate, BMF,fs, mfb_keyvals{:});
                        t_elapsed = toc(t_start);

                        out.ic(:,j) = ic_out(1:L);
                        out.t_elapsed(j) = t_elapsed;
                        if j == 1
                            subfs = fs;
                            amt_disp(sprintf('Zilany2014: CN: tau_exc=%f, tan_inh=%f, D=%f, S=%f',par.tau_ex_cn,par.tau_inh_cn,par.cn_delay,par.Sinh_cn));
                            amt_disp(sprintf('Zilany2014: IC: tau_exc=%f, tan_inh=%f, D=%f, S=%f',par.tau_ex_ic,par.tau_inh_ic,par.ic_delay_inh,par.Sinh_ic));
                        end
                    end
                    out.fcs = fcs;
                    out.fcs_description = 'CFs that were tested';
                    mfc = [];

                    out.ic_all = sum(out.ic,2);
                    out = rmfield(out,'ic'); % to reduce size

                    c = [];
                    c.out   = out;
                    c.subfs = subfs;
                    c.mfc   = mfc;
                    amt_cache('set',fname,c);
                end

                factor = 1;
                unit = 'spikes/s';
                YL  = [-15 145];
                YT = 0:15:130;

            case 'verhulst2015'
                if bRun
                    t_start = tic;
                    out = verhulst2015(insig,fs,fc_flags{:},afb_flags{:},ihc_flags{:},an_flags{:}, ...
                        an_keyvals{:},mfb_flags{:});
                    t_elapsed = toc(t_start);

                    subfs = out(1).fs_abr;
                    out.t_elapsed = t_elapsed;
                    out = rmfield(out,{'fs_bm','ihc','an_summed','ic'});
                    mfc = [];

                    c = [];
                    c.out   = out;
                    c.subfs = subfs;
                    c.mfc   = mfc;
                    amt_cache('set',fname,c);
                end

                factor = 1e6; unit = '\muV';
                YL  = [-.14 .18];
                YT = -.3:.03:.3;

            case 'verhulst2018'
                if bRun
                    t_start = tic;
                    out = verhulst2018(insig,fs,fc_flags{:},afb_flags{:},ihc_flags{:},an_flags{:}, ...
                        kv.an_keyvals{:},mfb_flags{:},'no_ihc','no_an','no_cn');
                    t_elapsed = toc(t_start);

                    subfs = out(1).fs_abr;
                    out.t_elapsed = t_elapsed;
                    out = rmfield(out,{'fs_bm','ic'});
                    mfc = [];

                    c = [];
                    c.out   = out;
                    c.subfs = subfs;
                    c.mfc   = mfc;
                    amt_cache('set',fname,c);
                end

                factor = 1e6; unit = '\muV';
                YL  = [-.22 .32];
                YT = -.3:.05:.3;

            case 'bruce2018'
                if bRun
                    fcs = local_m2hz(401);
                    idxs = length(fcs):-idx_step:1;
                    fcs = fcs(idxs);
                    fcs = fcs(fcs>125);

                    kv = {'numH',numH,'numM',numM,'numL',numL,'nrep',nrep, ...
                        'specificSRautoTiming','psthbinwidth_mr',psth_binwidth};
                    kv_one_rep = {'numH',numH,'numM',numM,'numL',numL,'nrep',1,'specificSRautoTiming'};

                    dt = psth_binwidth;
                    subfs = 1/dt;
                    L = round(dur_click*subfs);

                    for j = 1:length(fcs)
                        amt_disp(sprintf('bruce2018: band %.0f of %.0f\n',j,length(fcs)));
                        t_start = tic;
                        %%% To measure computing time only:
                        out_tmp = bruce2018(insig,fs,fcs(j),afb_flags{:},kv_one_rep{:});
                        ic_out = carney2015(out_tmp.psth(1:L),BMF,subfs, mfb_keyvals{:});
                        t_elapsed = toc(t_start);

                        %%% To assess a reliable PSTH:
                        out_tmp = bruce2018(insig,fs,fcs(j),afb_flags{:},kv{:});
                        out_psth = out_tmp.psth(1:L);
                        ic_out = carney2015(out_psth,BMF,subfs,mfb_keyvals{:});
                        out.ic(1:L,j) = ic_out(1:L);
                        out.t_elapsed(j) = t_elapsed;
                    end

                    out.fcs = fcs;
                    out.fcs_description = 'CFs that were tested';
                    mfc = [];

                    out.ic_all = sum(out.ic,2);
                    out = rmfield(out,'ic'); % to reduce size

                    c = [];
                    c.out   = out;
                    c.subfs = subfs;
                    c.mfc   = mfc;
                    amt_cache('set',fname,c);
                end

                factor = 1;
                unit = 'spikes/s';
                YL  = [-15 145];
                YT = 0:15:130;

            case 'king2019'

                if bRun
                    subfs = 20000;
                    afb_kv = {'basef',[],'flow',80,'fhigh',8000, ... % arbitrary
                        'compression_n',compression_n,'dboffset',dBFS,'subfs',subfs, ...
                        'modbank_Nmod',10}; % 10 modulation filters instead of 5

                    t_start = tic;
                    [outsig,fcs,mfc_here] = king2019(insig,fs,afb_kv{:},afb_flags{:},ihc_flags{:},an_flags{:},mfb_flags{:});
                    out.t_elapsed = toc(t_start);

                    mfc_idx = find(mfc_here<mfc_target,1,'last');
                    outsig = squeeze(outsig(:,:,mfc_idx));

                    out.fcs    = fcs;
                    out.oustig = outsig;
                    out.outsig_description = 'Modulation filter bank output';

                    mfc = mfc_here(mfc_idx);

                    c = [];
                    c.out   = out;
                    c.subfs = subfs;
                    c.mfc   = mfc;
                    amt_cache('set',fname,c);
                end

                factor = 1e4;
                unit = 'a.u. x 10^{-4}';
                YL  = [-0.1 2.1];
                YT = [0:.2:2];

            case 'relanoiborra2019'

                if bRun
                    % afb_kv = {'erbspacebw','no_internalnoise'};
                    t_start = tic;
                    [out_tmp,mfc_here,~,fc_here] = relanoiborra2019_featureextraction(insig, fs, afb_flags{:});
                    out.t_elapsed = toc(t_start);

                    mfc_idx = find(mfc_here<mfc_target,1,'last');
                    mfc = mfc_here(mfc_idx);

                    Nr_fc = length(fc_here);
                    outsig = [];
                    fcs = [];
                    for j = 1:Nr_fc
                        if mfc < fc_here(j)/4 % limiting the bands (Verhey1999)
                            if size(out_tmp,3) >= mfc_idx
                                fcs(end+1) = fc_here(j); % only adding the fcs of filters containing mfc_idx
                                outsig(:,end+1) = out_tmp(:,j,mfc_idx);

                                % outsig_adapt(:,end+1) = out_tmp_an(:,j);
                            end
                        end
                    end

                    out.Nr_fc = Nr_fc;
                    out.fcs    = fcs;
                    out.oustig = outsig;
                    out.outsig_description = 'Modulation filter bank output';

                    subfs = fs;

                    c = [];
                    c.out   = out;
                    c.subfs = subfs;
                    c.mfc   = mfc;
                    amt_cache('set',fname,c);
                end

                factor = 1;
                unit = 'MU';
                YL  = [-20 500];
                YT = 0:50:450;

            case 'osses2021'
                if bRun
                    subfs = 20000;
                    afb_kv = {'subfs',subfs,'dboffset',dBFS};

                    t_start = tic;
                    [out_tmp,fc_here,mfc_here] = osses2021(insig,fs,afb_kv{:},afb_flags{:},ihc_flags{:},an_flags{:},mfb_flags{:});
                    out.t_elapsed = toc(t_start);

                    mfc_idx = find(mfc_here<mfc_target,1,'last');

                    Nr_fc = length(fc_here);
                    outsig = [];

                    fcs = [];
                    for j = 1:Nr_fc
                        if size(out_tmp{j},2) >= mfc_idx
                            fcs(end+1) = fc_here(j); % only adding the fcs of filters containing mfc_idx
                            outsig(:,end+1) = out_tmp{j}(:,mfc_idx);
                        end
                    end

                    out.fcs    = fcs;
                    out.oustig = outsig;
                    out.outsig_description = 'Modulation filter bank output';

                    mfc = mfc_here(mfc_idx);

                    c = [];
                    c.out   = out;
                    c.subfs = subfs;
                    c.mfc   = mfc;
                    amt_cache('set',fname,c);
                end

                factor = 1;
                unit = 'MU';
                YL  = [-20 230];
                YT = 0:20:220;
        end

        if ~isempty(out)

            w5 = [];
            switch models{k}
                case {'dau1997','king2019'}
                    if ~isfield(out,'Nr_fc')
                        Nr_fc = 31;
                    else
                        Nr_fc = out.Nr_fc;
                    end
                    num_cfs = length(out.fcs);
                    w5 = sum(out.oustig,2)/num_cfs;

                case {'relanoiborra2019','osses2021'}
                    num_cfs = length(out.fcs);
                    if num_cfs == 26
                        if ~isfield(out,'Nr_fc')
                            Nr_fc = 31;
                        else
                            Nr_fc = out.Nr_fc;
                        end
                    elseif num_cfs == 51
                        if ~isfield(out,'Nr_fc')
                            Nr_fc = 60;
                        else
                            Nr_fc = out.Nr_fc;
                        end
                    end
                    w5 = sum(out.oustig,2)/num_cfs;

                case {'zilany2014','bruce2018'}
                    num_cfs = length(out.fcs);
                    w5 = out.ic_all/num_cfs;

                    if ~isfield(out,'Nr_fc')
                        Nr_fc = num_cfs;
                    else
                        Nr_fc = out.Nr_fc;
                    end
                case {'verhulst2015','verhulst2018'}
                    if ~isfield(out,'Nr_fc')
                        Nr_fc = 401;
                    else
                        Nr_fc = out.Nr_fc;
                    end
                    w5 = out.w5;
            end

            % -------------------------------------------------------------
            % Time elapsed:
            if length(out.t_elapsed) ~= 1
                t_elapsed_here = sum(out.t_elapsed);
            else
                t_elapsed_here = out.t_elapsed;
            end
            t_elapsed(1,k) = t_elapsed_here;
            t_elapsed(2,k) = Nr_fc;
            t_elapsed(3,k) = t_elapsed(1,k)/Nr_fc;
            if k == 1
                t_elapsed_description = 'Row 1: Total time elapsed (s); Row 2: Total number of simulated channels; Row 3 = time per channel (Row 1/Row 2)';
            end
            if dispflag.do_verbose
                amt_disp(sprintf('Model: %s; number of channels=%.0f, nr of bands with a ''BMF=100 Hz''=%.0f',models{k},Nr_fc,num_cfs));
                amt_disp(sprintf('\t time elapsed=%.3f s (%.4f s per channel, %.0f channels)',t_elapsed(1,k),t_elapsed(3,k),t_elapsed(2,k)));
            end
            % END Time elapsed --------------------------------------------

            t_ms = 1000*(1:length(w5))/subfs;

            disp('Using Wave-V estimate only...')
            efr_here = w5;

            if flags.do_plot
                figure; % For each model there is one new figure
                figure_handle(end+1) = gcf;
                figure_name{end+1} = ['fig15-click-' models{k}];

                opts     = {'LineStyle','-' ,'Color',Colours{k},'LineWidth',LineWidth(k)};
                opts_neg = {'LineStyle','--','Color',local_rgb('Gray'),'LineWidth',3};
                offx = 10;
            end
            t_ms_orig = t_ms;
            %%% Getting the clicks to be plotted only
            idxs = [1 2 9 10]; % Click numbers to be plotted
            
            N = length(efr_here)/Nr_clicks;
            M = Nr_clicks;

            t_ms_orig = t_ms;
            t_ms     = reshape(t_ms    ,N,M);

            efr_here = reshape(efr_here,N,M);
            efr_here = efr_here(:,idxs);
            Nr_clicks_here = length(idxs);

            t_ms = t_ms(:,idxs);

            efr_here = efr_here(:); % truncated clicks back to a one-column array
            if Nr_clicks_here == 4
                efr_here = reshape(efr_here,N*2,Nr_clicks_here/2);
            end
            
            switch models{k}
                case {'zilany2014','bruce2018'}
                    t_ms_here = t_ms_orig(1:size(efr_here,1));
                    idx_pos = find(efr_here(:,2)>=0); % positive-valued amplitudes
                    idx_neg = find(efr_here(:,2)<0);  % negative-valued amplitudes

                    efr2plot = nan(size(efr_here(:,2)));
                    efr2plot(idx_pos) = efr_here(idx_pos,2);
                    if flags.do_plot
                        pl(k) = plot(t_ms_here,factor*efr2plot,opts{:}); hold on, grid on
                    end

                    efr2plot = nan(size(efr_here(:,2)));
                    efr2plot(idx_neg) = zeros(size(efr_here(idx_neg,2)));
                    if flags.do_plot
                        plot(t_ms_here,factor*efr2plot,opts_neg{:}); hold on, grid on
                    end

                otherwise
                    if flags.do_plot    
                        pl(k) = plot(t_ms_orig(1:size(efr_here,1)),factor*efr_here(:,2),opts{:}); hold on, grid on
                    end
            end

            efr_here = efr_here(:);
            
            if flags.do_plot
                xlabel('Time relative to Click A onset (ms)');
                title(models{k})
            
                step_time = 25; % 30
                XT = 10:step_time:210;
                XT_for_tick = XT;

                % XL = [-5 205];
                XTL = [];
                for ii = 1:length(XT)
                    if mod(ii,2) == 1
                        XTL{ii} = num2str(XT_for_tick(ii)-offx);
                    else
                        XTL{ii} = '';
                    end
                end

                set(gca,'XTick'     ,XT);
                set(gca,'XTickLabel',XTL);

                ylim(YL);
            
                YTL = [];
                for ii = 1:length(YT)
                    if mod(ii,2) == 1
                        YTL{ii} = num2str(YT(ii));
                    else
                        YTL{ii} = '';
                    end
                end

                set(gca,'YTick',YT);
                set(gca,'YTickLabel',YTL);

                ylabel(['(' unit ')']);

                Pos = get(gcf,'Position');
                Pos(3:4) = [380 300];
                set(gcf,'Position',Pos);
            end % flags.do_plot
        end % end ~ifempty(out)

        Nr_samples = length(efr_here);
        output_buf = reshape(efr_here,round(Nr_samples/Nr_clicks_here),Nr_clicks_here);

        Amp_max(k,:) = max(output_buf);
        Amp_min(k,:) = min(output_buf);
    end % end for k

    data.Amp_max = Amp_max;
    data.Amp_min = Amp_min;
    data.Amp_pp  = Amp_max - Amp_min;
    data.models  = models;

    if flags.do_plot
        data.figure_flag = 'do_fig15';
        data.figure_handle = figure_handle;
        data.figure_name   = figure_name;
    end
end

%%% End of file exp_osses2022.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of subfunctions:
function [flags,keyvals] = local_get_flags(model_str)

switch model_str
    case {'verhulst2015','verhulst2018'}
        kSR_H  = 68.5; % spikes/s, spontaneous rate: only Verhulst 2018
        kSR_M  = 10; % spikes/s, spontaneous rate: Verhulst 2015 and 2018
        kSR_L  = 1; % spikes/s, spontaneous rate: Verhulst 2015 and 2018
end

numH = 12;
numM =  4;
numL =  4;

bAMT = 1; % Is AMT v1.0

switch model_str
    case 'dau1997'
        fc_flags   = [];
        afb_flags = {};
        afb_keyvals = [];
        ihc_flags  = {'ihc','ihc_dau1996'};
        noihc_flags  = {'no_ihc'};
        an_flags   = {'adt','adt_dau1997'};
        noan_flags = {'no_adt','no_mfb'};
        an_keyvals = [];
        mfb_flags  = {'mfb','mfb_dau1997'};
        nomfb_flags  = {'no_mfb'};
        mfb_keyvals = {'subfs',[]}; % subfs will be the same as fs

    case 'osses2021'
        fc_flags   = [];
        afb_flags  = {'afb_osses2021'};
        afb_keyvals = [];
        ihc_flags  = {'ihc','ihc_breebaart2001'};
        noihc_flags  = {'no_ihc'};
        an_flags   = {'adt', 'adt_osses2021'};
        noan_flags = {'no_adt','no_mfb'};
        an_keyvals = [];
        mfb_flags  = {'mfb','mfb_jepsen2008'};
        nomfb_flags  = {'no_mfb'};
        mfb_keyvals = [];

    case 'relanoiborra2019'
        fc_flags   = [];
        afb_flags  = {'no_internalnoise','erbspacebw'};
        afb_keyvals = [];
        ihc_flags  = {'ihc','ihc_relanoiborra2019'};
        noihc_flags  = {'no_ihc'};
        an_flags   = {'adt', 'adt_relanoiborra2019'};
        noan_flags = {'no_an'};% ,'no_mfb'};
        an_keyvals = [];
        mfb_flags  = {'mfb','mfb_jepsen2008'};
        nomfb_flags  = {'no_mfb'};
        mfb_keyvals = [];

    case 'king2019'
        fc_flags   = [];
        afb_flags  = {'afb','compression_brokenstick'};
        afb_keyvals = {'compression_n',0.3};
        ihc_flags  = {'ihc'};
        noihc_flags  = {'no_ihc'};
        an_flags   = {'adt'};
        noan_flags = {'no_adt','no_mfb'};
        an_keyvals = [];
        mfb_flags  = {'mfb'};
        nomfb_flags  = {'no_mfb'};
        mfb_keyvals = [];

    case 'verhulst2015'
        fc_flags   = {'abr','debug'}; % abr = 401 frequency channels
        afb_flags  = {'no_outerear','middleear','no_v','no_oae','no_y'}; % these are the defaults
        afb_keyvals = [];
        ihc_flags  = {'ihc'}; % this is the default
        noihc_flags  = {'no_ihc'};
        an_flags   = {'an'};
        noan_flags = {'no_an','no_cn','no_ic'}; % if no_an, then no_cn and no_ic
        an_keyvals = {'numH',numH,'numM',numM,'numL',numL,'kSR_H',kSR_H,'kSR_M',kSR_M,'kSR_L',kSR_L};
        mfb_flags  = {'no_cn','ic'};
        nomfb_flags  = {'no_cn','no_ic'};
        mfb_keyvals = [];

    case 'verhulst2018'
        fc_flags   = {'abr','debug'}; % abr = 401 frequency channels
        afb_flags  = {'no_outerear','middleear','no_v','no_oae','no_y'}; % these are the defaults
        afb_keyvals = [];
        ihc_flags  = {'ihc'}; % this is the default
        noihc_flags  = {'no_ihc'};
        an_flags   = {'an'};
        noan_flags = {'no_an','no_cn','no_ic'};
        an_keyvals = {'numH',numH,'numM',numM,'numL',numL,'kSR_H',kSR_H,'kSR_M',kSR_M,'kSR_L',kSR_L};
        mfb_flags  = {'no_cn','ic'};
        nomfb_flags  = {'no_cn','no_ic'};
        mfb_keyvals = [];

    case 'zilany2014'
        species   = 'human';
        noiseType = 'fixedFGn';
        % noiseType = 'varFGn';
        fc_flags   = [];
        if bAMT
            afb_flags = {species};
        else
            afb_flags = {'middleear',species}; % this is the default, making sure that 'species' is always loaded
        end
        afb_keyvals = [];
        ihc_flags  = {'ihc'}; % this is the default
        noihc_flags  = {'no_ihc'};
        an_flags   = {'an',species,noiseType};
        if numH ~= 0
            an_flags(end+1) = {'anfH'};
        end
        if numM ~= 0
            an_flags(end+1) = {'anfM'};
        end
        if numL ~= 0
            an_flags(end+1) = {'anfL'};
        end
        noan_flags = {'no_an',species,noiseType};
        an_keyvals = {'numH',numH,'numM',numM,'numL',numL,'psth_binwidth',0.5e-3};
        mfb_flags  = {'cn','ic'};
        nomfb_flags  = {'no_cn','no_ic'};
        mfb_keyvals = {'BMF',90,'ic_delay_inh',0.0011};

    case 'bruce2018'
        species   = 'human';
        noiseType = 'fixedFGn';
        % noiseType = 'varFGn';

        fc_flags   = [];
        if bAMT
            afb_flags = {species};
        else
            afb_flags = {'middleear',species}; % this is the default, making sure that 'species' is always loaded
        end
        afb_keyvals = [];
        ihc_flags  = {'ihc'}; % this is the default
        noihc_flags  = {'no_ihc'};
        an_flags   = {'an',species,noiseType};
        if numH ~= 0
            an_flags(end+1) = {'anfH'};
        end
        if numM ~= 0
            an_flags(end+1) = {'anfM'};
        end
        if numL ~= 0
            an_flags(end+1) = {'anfL'};
        end
        noan_flags = {'no_an',species,noiseType};
        an_keyvals = {'numH',numH,'numM',numM,'numL',numL,'psth_binwidth',0.5e-3};
        mfb_flags  = {'cn','ic'};
        nomfb_flags  = {'no_cn','no_ic'};
        mfb_keyvals = {'BMF',90,'ic_delay_inh',0.0011};
end

flags = [];
flags.fc_flags   = fc_flags;
flags.afb_flags  = afb_flags;
flags.ihc_flags  = ihc_flags;
flags.noihc_flags  = noihc_flags;
flags.an_flags   = an_flags;
flags.noan_flags = noan_flags;
flags.mfb_flags  = mfb_flags;
flags.nomfb_flags= nomfb_flags;

keyvals = [];
keyvals.afb_keyvals = afb_keyvals;
keyvals.an_keyvals  = an_keyvals;
keyvals.mfb_keyvals = mfb_keyvals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rms_val = local_rmsdb(insig)

rms_val = 20*log10(rms(insig));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ac_target,dc_target,vrest] = local_get_ac_dc_osses2022(ihc_target,ac_reg,dc_reg)
% function [ac_target,dc_target,vrest] = local_get_ac_dc_osses2022(ihc_target,ac_reg,dc_reg)
%
% DC: is computed from the signal (mean value using the ac_reg samples), and
%     then this value is referenced to the resting potential
% AC: I compute it here as the peak-to-peak voltage, in line with Russel and
%     Palmer 1986

vrest     = median(ihc_target(dc_reg,:)); % median of the silent initial segment
dc_target = mean(ihc_target(ac_reg,:))-vrest;

v_min   = min(ihc_target(ac_reg,:));
v_max   = max(ihc_target(ac_reg,:));

dist_ac_pp = v_max-v_min; % peak to peak
ac_target  = dist_ac_pp; % approximation to the RMS of the AC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,x] = local_m2hz(Nr_sections,k)
% function [f,x] = local_m2hz(Nr_sections,k)
%
% 1. Description:
%       See Greenwood 1990, Equation 1 with A = 165.4, k = 0.85.
%       bApex_to_base is referenced to Helicotrema.
%
% Programmed by Alejandro Osses, WAVES, UGent, Belgium, 2018-2020

if nargin < 2
    k = 0.85; % As used in the models by Verhulst et al
end

if nargin == 0
    Nr_sections = 500;
end

cochleaLength = 0.035; % m
helicotremaLength = 0.001; % m
A = 165.4118; % For human specie
a = 61.765; %
bm_length = cochleaLength-helicotremaLength;

switch Nr_sections
    case {201,401}
        % Nr_sections = 500;
        f_range = [12000 112];
        step_bm = bm_length/500;
    otherwise
        f_range = [];
        step_bm = bm_length/Nr_sections;
end

A0 = 20682; % cochlear_model2018.py, L114
switch Nr_sections
    case 201
        x = step_bm  :2*step_bm:bm_length;

    case 401
        x = step_bm  :step_bm:bm_length;

    case 500
        x = step_bm/2:step_bm:bm_length;

    case 1000
        x = step_bm  :step_bm:bm_length;

    otherwise
        error('Check this function: %s',mfilename);
end

f = A0*(10.^(-a*x)) - A*k;

if ~isempty(f_range)
    idxi = find(f>= f_range(1),1,'last');
    idxf = find(f>= f_range(2),1,'last');

    f = f(idxi:idxf);
    x = x(idxi:idxf);
end
dispflag = amt_flags;
if nargout == 0
    if dispflag.do_verbose
        amt_disp(sprintf('The output CF derived from %.0f cochlear sections with frequencies between %.1f and %.1f [Hz]',length(f),f(1),f(end)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function col_rgb = local_rgb(colour)

switch colour
    case 'Gray'
        col_rgb = [0.5 0.5 0.5];
    case 'LightGray'
        col_rgb = [0.8242 0.8242 0.8242];
    case 'Green'
        col_rgb = [0 0.5 0];
    case 'LightSkyBlue'
        col_rgb = [0.5273 0.8047 0.9792];
    case 'Maroon'
        col_rgb = [0.5 0 0];
    case 'mediumorchid'
        col_rgb = [0.7266 0.3320 0.8242];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [env,idx_env] = local_Get_envelope2(insig,me)

method = 2;

switch method
    case 1
        error('Not used, tested on 7/4/2021')
    case 2

        env = [];
        idx_env = [];
        if nargin < 4
            me = mean(insig); % prctile(insig,75);
        end

        L1 = [1:length(insig); insig'];
        L2 = [1 length(insig); me me];

        P = local_InterX(L1,L2);
        N_P = size(P,2);

        for j = 1:2:N_P-2
            idx1 = round(P(1,j));
            idx2 = round(P(1,j+1));
            [env(end+1),idx_here] = max(insig(idx1:idx2));

            idx_env(end+1) = idx_here+idx1-1;
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BW_ERB, Q_ERB, Q03, Q10, extra]= local_Get_ERB_estimation_multi(impulse_ch, CF, fs, N_section)
dispflag = amt_flags;
N_samples = length(impulse_ch);
insig = reshape(impulse_ch,N_samples/N_section, N_section);

N = size(insig,1);
K = floor(N/2);         % K-non redundant points

freq_here = (1:K)*(fs/2)/K;
extra.freq_here = freq_here;
Amp_avg = zeros(K,1);

bUse_MA = 1;

for i = 1:N_section
    impulse_ch = insig(:,i);

    Amp = (2*abs(fft(impulse_ch,N))/N).^2; 
    Amp = Amp(:);
    Amp = Amp(1:K);

    Amp_raw(:,i) = Amp(1:K);

    Amp_avg = Amp_avg + Amp(1:K);
end

Amp_avg = Amp_avg / N_section;

Amp_raw = 10*log10(Amp_raw);
Amp_dB = 10*log10(Amp_avg);
max_dB = max(Amp_dB(2:end)); % starts at second sample to exclude the DC
Amp_dB = Amp_dB-max_dB;

extra.Amp_raw = Amp_raw;
extra.Amp_dB  = Amp_dB;
extra.max_dB  = max_dB;

if bUse_MA
    amt_disp(sprintf('\n\tApplying moving average\n'),'volatile');

    % N_moving = 15;
    Moving_step = 0.3; % ERB_N
    CF_for_step = audtofreq(freqtoaud(CF)-Moving_step);
    df = freq_here(2)-freq_here(1);
    N_moving = round((CF-CF_for_step)/df);
    N_moving = 2*N_moving+1;
    Amp_MA = local_moving_average(Amp_dB,N_moving);
end

if bUse_MA
    max_MA = max(Amp_MA(2:end));
    Amp_dB = Amp_MA-max_MA;
end
extra.Amp_dB_MA  = Amp_dB;
extra.max_MA = max_MA;

%%%
Amp_avg = 10.^((Amp_dB+max_dB)/10); % reverting the maximum
[M,Mn]  = max(Amp_avg);
E       = sum(Amp_avg);

BW_ERB  = (E./M)*fs/N;
Q_ERB   = CF./BW_ERB;
%%%

L1 = [freq_here; transpose(Amp_dB)];
L2 = [min(freq_here) max(freq_here); [-10 -10]];
P = local_InterX(L1, L2);
if ~isempty(P) && size(P,1)==2
    if size(P,2) > 2
        P = local_check_values(P);

        if size(P,2) > 2
            disp('')
        end
    end
    flow_10 = P(1,1);
    fhigh_10 = P(1,2);
    BW10 = P(1,2)-P(1,1);
    Q10  = CF/BW10;
else
    flow_10 = nan;
    fhigh_10 = nan;
    BW10 = nan;
    Q10  = nan;
end

L2 = [min(freq_here) max(freq_here); [-3 -3]];
P = local_InterX(L1, L2);
if ~isempty(P) && size(P,1)==2
    if size(P,2) > 2
        P = local_check_values(P);

        if size(P,2) > 2
            if dispflag.do_verbose
                amt_disp(sprintf('\tStill more than two cut-offs, keeping the ones with the largest BW (manual check: this can be the case for verhulst2018)...'));
            end
            [xx,idx] = max(diff(P(1,:)));

            P = P(:,idx:idx+1);
        end
    end
    try
        flow_03 = P(1,1);
        fhigh_03 = P(1,2);
    catch
        disp('')
    end
    BW03 = P(1,2)-P(1,1);
    Q03  = CF/BW03;

    if Q03 > 12
        disp('')
    end
else
    flow_03 = nan;
    fhigh_03 = nan;
    BW03 = nan;
    Q03  = nan;
end

extra.BW_ERB = BW_ERB;
extra.BW10 = BW10;
extra.BW03 = BW03;
extra.Q_ERB = Q_ERB;
extra.Q10 = Q10;
extra.Q03 = Q03;
extra.flow_03 = flow_03;
extra.flow_10 = flow_10;
extra.fhigh_03 = fhigh_03;
extra.fhigh_10 = fhigh_10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = local_check_values(P)

idxs2remove = find(P(1,:) < 20);
P(:,idxs2remove) = [];

tolerance_val = 40;
idxs2remove = find(diff(P(1,:)) < tolerance_val);

P(:,idxs2remove+1) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Amp_Ma = local_moving_average(Amp_dB,N_Ma)

N_i = (N_Ma-1)/2 + 1;
step_downup = (N_Ma-1)/2;

Amp_Ma = Amp_dB;
for i = N_i:length(Amp_dB)-N_i
    try
        Amp_Ma(i) = mean(Amp_dB(i-step_downup:i+step_downup));
    catch
        disp('')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = local_InterX(L1,L2)
%LOCAL_INTERX Intersection of curves
%   P = LOCAL_INTERX(L1,L2) returns the intersection points of two curves L1
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point.
%   Each factor of the 'C' arrays is essentially a matrix containing
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.
%
% This function is based on InterX.m (v3.0, 21 Sept. 2010) by author 'NS'
% that is available within the Mathworks FileExchange.

% Reordering the data before intersecting L1 and L2:
x1  = L1(1,:)';  x2 = L2(1,:);
y1  = L1(2,:)';  y2 = L2(2,:);
dx1 = diff(x1); dy1 = diff(y1);
dx2 = diff(x2); dy2 = diff(y2);

%...Determine 'signed distances'
S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);

C1 = le(local_D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
C2 = le(local_D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

%...Obtain the segments where an intersection is expected
[i,j] = find(C1 & C2);
if isempty(i),P = zeros(2,0);return; end;

%...Transpose and prepare for output
i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0

%...Solve system of eqs to get the common points
P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
            dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = local_D(x,y)
u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [waveform,freq_array,phases] = local_Profile_Analysis(dur, rampdur, stim_dB, fs, freq_array, phases)
% function [waveform,freq_array,phases] = local_Profile_Analysis(dur, rampdur, stim_dB, fs, freq_array, phases)
%
% This function produces a vector of log-spaced components based on Lentz (JASA 2005)
%Variables:
% dur = duration in sec  (e.g. 0.5 sec)
% ramp_time: onset/offset in seconds
% stim_dB = overall stimulus level in dB SPL (both inetervals are scaled to same dB SPL)
% fs = sampling rate (Hz)
%
% Code adapted from UR EAR function

L = round(dur * fs); % stimulus length (in samples)
t = (0:L-1)/fs; % time vector
waveform = 0;  % initialize variable

if nargin < 6
    phases = 2*pi * rand(1,length(freq_array));
end
for i = 1:length(freq_array)  % step through each frequency component in the complex
    f = freq_array(i);
    phase = phases(i); % Randomly vary the starting phase of each component
    waveform = waveform + cos(2 * pi * f * t + phase);
    % end
end

waveform = tukeywin(L, 2*rampdur/dur)' .* waveform; % gate
waveform = 20e-6 * 10.^(stim_dB/20) * waveform/rms(waveform);% convert signal into Pascals


