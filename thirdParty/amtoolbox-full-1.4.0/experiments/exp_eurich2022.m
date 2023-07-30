function exp_eurich2022(varargin)
%EXP_EURICH2022 Results from Eurich et al. (2022)
%   Usage: data = exp_eurich2022(flag)
%
%   EXP_EURICH2022(flag) reproduces figures of the study from
%   Eurich et al. (2022)
%
%   The following flags can be specified
%
%     'fig3'    Reproduces Figure 3;
%               Experimental data from van der Heijden and
%               Trahiotis (1999) (symbols). The continuous lines show the predictions of
%               the presented model including the across-channel incoherence interference.
%               The dashed lines show predictions for ODN excluding interference (single-
%               channel version), equivalent to Encke and Dietz (2022). Upper panel:
%               Detection thresholds with S0 target; lower panel: Spi target.
%
%     'fig4'    Reproduces Figure 4;
%               Experimental data from Marquardt and McAlpine
%               (2009) (symbols) and model predictions (lines). Detection thresholds are
%               given as function of the inner-band bandwidth. The inner band contains
%               delayed noise (triangles) or opposingly delayed noises (diamonds and bul-
%               lets) with a fixed ITD = 1 ms while the flanking bands have a constant IPD
%               of +pi/2 (upward triangle and diamond) and -pi/2, or vice versa (down-
%               ward triangle and bullet). Continuous and dashed lines again show predic-
%               tions with and without across-frequency incoherence interference,
%               respectively.
%
%
%     'fig5'    Reproduces Figure 5;
%               Symbols denote data from binaural detection experi-
%               ments with the configuration Np0p Sp as a function of the inner-band (N0)
%               bandwidth; continuous and dotted line: model prediction with and without
%               across-incoherence incoherence interference, respectively.
%
%
%   Further, cache flags (see amt_cache) and plot flags can be specified:
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'no_plot'  Don't plot, only return data.
%
%
%   Examples:
%   ---------
%
%   To display Figure 3 use :
%
%     exp_eurich2022('fig3');
%
%   To display Figure 4 use :
%
%     exp_eurich2022('fig4');
%
%   To display Figure 5 use :
%
%     exp_eurich2022('fig5');
%
%
%   See also: eurich2022, eurich2022_processing, eurich2022_decision,
%   sig_vanderheijden1999, sig_marquardt2009, sig_kolarik2010
%
%   References:
%     B. Eurich, J. Encke, S. D. Ewert, and M. Dietz. Lower interaural
%     coherence in off-signal bands impairs binaural detection. The Journal
%     of the Acoustical Society of America, 151(6):3927--3936, 06 2022.
%     [1]arXiv | [2]http ]
%     
%     References
%     
%     1. http://arxiv.org/abs/https://pubs.aip.org/asa/jasa/article-pdf/151/6/3927/16528275/3927\_1\_online.pdf
%     2. https://doi.org/10.1121/10.0011673
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_eurich2022.php


%   #Author: Bernhard Eurich (2023): Original code
%   #Author: Piotr Majdak (2023): Integration in the AMT 1.4

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% ------ Check input options --------------------------------------------


definput.import={'amt_cache'};

definput.flags.type = {'missingflag','fig3','fig4','fig5'};
definput.flags.plot = {'plot','no_plot'};

definput.keyvals.FontSize = 9;
definput.keyvals.MarkerSize = 3;
definput.import={'amt_cache'};

flags  = ltfatarghelper({'FontSize','MarkerSize'},definput,varargin);

%the error message:
if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
        sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end



%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%% ------ FIGURE 3 -----------------------------------------------------------
if flags.do_fig3
    %any calculations associated to Figure 3 in the publication
    if flags.do_plot
        fs             = 48e3;
        dur            = 300e-3; % s
        cos_rise_time  = 20e-3;
        bw             = 900;
        fc             = 500;
        spl            = 45.5; % noise spectrum level

        % variable stimulus parameters
        itd            = (0:0.125:4)*1e-3; % this is the tracking variable
        sphase         = [-1 1];
        tone_level     = 30:2:65; % target level / dB SPL
        noise_mode     = [1 2]; % 1: single delayed noise; 2: opposingly delayed noises
        repititions    = 1:50;

        tic
        % model parameters (default -- adjusted parameters for S0 condition are set in local_evaluation_fig3)
        mpar = arg_eurich2022;
        mpar.interference_sigma = [0 mpar.interference_sigma];
        mpar.idx_of_target_channel = 9; % holds if 1 filter per ERB; if empty, the model picks the channel offering largest sensitivity.


        filename = 'processed_average_fig3';%1
        processed_average_fig3 = amt_cache('get',filename,flags.cachemode); %2

        if isempty(processed_average_fig3)%3
            parfor i_itd = 1:length(itd)
                processed_average_fig3(i_itd,:,:,:,:,:,:) = local_evaluation_fig3(itd(i_itd),noise_mode,sphase,tone_level,repititions,fs,dur,cos_rise_time,bw,fc,spl,mpar);                
            end
            amt_cache('set',filename,processed_average_fig3);
        end



        % decision
        filename2 = 'dprime_fig3';%1
        dprime_fig3 = amt_cache('get',filename2,flags.cachemode);
        if isempty(dprime_fig3)
            for i_interference_sigma = 1:length(mpar.interference_sigma)
                for i_noise_mode = 1:length(noise_mode)
                    for i_itd = 1:length(itd)
                        for i_sphase = 1:length(sphase)
                            for i_tone_level = 1:length(tone_level)
                                for i_repititions = 1:length(repititions)
                                    dprime_fig3(i_interference_sigma,i_noise_mode,i_itd, i_sphase,i_tone_level) = eurich2022_decision(squeeze(processed_average_fig3(i_itd,i_interference_sigma,i_noise_mode, i_sphase,i_tone_level,:,:)), mpar);
                                end
                            end
                        end
                    end
                end
            end
            amt_cache('set',filename2,dprime_fig3);
        end

        % determine threshold

        % dprime scheint jetzt zu stimmen -- ab hier noch updaten bzgl. interference
        dprime_threshold = 0.78;

        % range of assumed linear psychometric function
        dprime0 = 0.5;
        dprime1 = 10^0.4;

        for i_interference_sigma = 1:length(mpar.interference_sigma)
            for i_noise_mode = 1:length(noise_mode)
                for i_itd = 1:length(itd)
                    for i_sphase = 1:length(sphase)

                        dprime_temp = squeeze(dprime_fig3(i_interference_sigma,i_noise_mode,i_itd,i_sphase,:));
                        idx0 = find(dprime_temp > dprime0 & dprime_temp < dprime1);

                        d_eval = squeeze(log10(dprime_temp(idx0)));

                        p0 = 20e-6;
                        level = tone_level(idx0);
                        %             level(:,c,b,a) = (p0 * 10.^(spar.dbspl_tone(idx0) / 20)) ./ val; %

                        %             fit regression line
                        d_eval(~isfinite(d_eval)) = NaN;
                        level1 = [ones(size(level')) level'];
                        slope = level1 \ d_eval;
                        d_fit = level1 * slope;
                        level_thresh(i_interference_sigma,i_noise_mode,i_itd, i_sphase) = (log10(dprime_threshold) - slope(1)) / slope(2);
                    end
                end
            end
        end

        % Plotting

        x.paperdata_pi=amt_load('eurich2022', 'vdHTr3b_1.mat');
        x.paperdata_0=amt_load('eurich2022', 'vdHTr3a_2.mat');

        data_raw_pi = x.paperdata_pi.paperdata_pi(:,2:3);
        data_raw_0  = x.paperdata_0.paperdata_0(:,2:3);

        data_pi(:,1) = data_raw_pi(~isnan(data_raw_pi(:,1)),1);
        data_pi(:,2) = data_raw_pi(~isnan(data_raw_pi(:,2)),2);

        data_0(:,1) = data_raw_0(~isnan(data_raw_0(:,1)),1);
        data_0(:,2) = data_raw_0(~isnan(data_raw_0(:,2)),2);

        itd_data = 0:0.125:4.125;


        col = [0 0 1; 0.5 0.5 0.5];%colororder;
        figure('Name','fig3', 'Position',[500 80 400 20/16*400]);

        subplot 211
        hold on;

        plot(itd_data, data_0(:,1), 's','Color',col(1,:));
        plot(itd_data, data_0(:,2), 'd','Color',col(2,:));


        % i_interference_sigma,i_noise_mode,i_itd, i_sphase
        plot(itd*1e3,squeeze(level_thresh(1,2,:,2)),'LineStyle',':','Color',col(1,:));
        plot(itd*1e3,squeeze(level_thresh(1,1,:,2)),'LineStyle',':','Color',col(2,:));
        plot(itd*1e3,squeeze(level_thresh(2,2,:,2)),'LineStyle','-','Color',col(1,:));
        plot(itd*1e3,squeeze(level_thresh(2,1,:,2)),'LineStyle','-','Color',col(2,:));


        ylim([45 61])
        ylabel('Threshold (dB)');
        xlim([-0.1 4.2])
        xticks([0 1 2 3 4])
        yticks([46:2:60])

        grid on;
        set(gca,'GridLineStyle',':','LineWidth',1)

        subplot 212

        hold on;


        plot(itd_data, data_pi(:,1), 's','Color',col(1,:));
        plot(itd_data, data_pi(:,2), 'd','Color',col(2,:));

        plot(itd*1e3,squeeze(level_thresh(1,2,:,1)),'LineStyle',':','Color',col(1,:));
        plot(itd*1e3,squeeze(level_thresh(1,1,:,1)),'LineStyle',':','Color',col(2,:));
        plot(itd*1e3,squeeze(level_thresh(2,2,:,1)),'LineStyle','-','Color',col(1,:));
        plot(itd*1e3,squeeze(level_thresh(2,1,:,1)),'LineStyle','-','Color',col(2,:));

        xlim([-0.1 4.2])
        xticks([0 1 2 3 4])
        yticks(46:2:60)

        ylim([45 61])

        grid on;
        set(gca,'GridLineStyle',':','LineWidth',1)

        ylabel('Threshold (dB)');

        xlabel('ITD (ms)')

        annotation('textbox',[.85 .57 .1 .1],'String','$S_0$','interpreter','Latex','EdgeColor',[1 1 1])
        annotation('textbox',[.85 .1,.1 .1],'String','$S_{\pi}$','interpreter','Latex','EdgeColor',[1 1 1])


        legend( 'SDN', 'ODN','Location','South','NumColumns',2);

        legend boxoff
    end
end



%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% ------ FIGURE 4 -----------------------------------------------------------
if flags.do_fig4
    %any calculations associated to Figure 4 in the publication
    if flags.do_plot
        fs             = 48e3;
        dur            = 300e-3; % s
        cos_rise_time  = 20e-3;
        bw             = 900;
        fc             = 500;
        spl            = 45.5; % noise spectrum level

        % variable stimulus parameters
        tone_level     = 40:2:70; % target level / dB SPL
        noise_mode     = [1 2]; % 1: single delayed noise; 2: opposingly delayed noises
        flanking_phase = [pi/2 -pi/2]; % the other flank gets the other sign
        innerbw        = [900 600 400 200 100 50 0];
        repititions    = 1:50;

        tic
        % model parameters
        mpar = arg_eurich2022;
        mpar.interference_sigma = [0 mpar.interference_sigma];
        mpar.idx_of_target_channel = 9;

        % adjust model parameters as optimized for predicting the data by
        % Marquardt & McAlpine
        mpar.binaural_internal_noise_sigma = 0.24;
        mpar.monaural_internal_noise_sigma = 0.4;
        mpar.rho_max = 0.89;
        mpar.interference_sigma = [0 0.65];


        %C.H.: example for using amt_cache---------------------------------
        filename = 'processed_average_fig4';%1
        processed_average_fig4 = amt_cache('get',filename,flags.cachemode); %2

        if isempty(processed_average_fig4)%3
            parfor i_innerbw = 1:length(innerbw)
                processed_average_fig4(i_innerbw,:,:,:,:,:,:) = ...
                    local_evaluation_fig4(innerbw(i_innerbw),noise_mode,flanking_phase,tone_level,repititions,fs,dur,cos_rise_time,bw,fc,spl,mpar);
            end
            amt_cache('set',filename,processed_average_fig4);%4
        end

        %% decision
        filename2 = 'dprime_fig4';%1
        dprime_fig4 = amt_cache('get',filename2,flags.cachemode);
        if isempty(dprime_fig4)
            for i_interference_sigma = 1:length(mpar.interference_sigma)
                for i_noise_mode = 1:length(noise_mode)
                    for i_innerbw = 1:length(innerbw)
                        for i_flanking_phase = 1:length(flanking_phase)
                            for i_tone_level = 1:length(tone_level)
                                dprime_fig4(i_interference_sigma,i_noise_mode,i_innerbw, i_flanking_phase,i_tone_level) = eurich2022_decision(squeeze(processed_average_fig4(i_innerbw,i_interference_sigma,i_noise_mode,i_flanking_phase,i_tone_level,:,:)), mpar);
                            end
                        end
                    end
                end
            end
            amt_cache('set',filename2,dprime_fig4);
        end

        %% determine threshold

        dprime_threshold = 1.14;

        % range of assumed linear psychometric function
        dprime0 = 0.5;
        dprime1 = 10^0.4;

        for i_interference_sigma = 1:length(mpar.interference_sigma)
            for i_noise_mode = 1:length(noise_mode)
                for i_innerbw = 1:length(innerbw)
                    for i_flanking_phase = 1:length(flanking_phase)

                        dprime_temp = squeeze(dprime_fig4(i_interference_sigma,i_noise_mode,i_innerbw,i_flanking_phase,:));
                        idx0 = find(dprime_temp > dprime0 & dprime_temp < dprime1);

                        d_eval = squeeze(log10(dprime_temp(idx0)));

                        p0 = 20e-6;
                        level = tone_level(idx0);
                        %             level(:,c,b,a) = (p0 * 10.^(spar.dbspl_tone(idx0) / 20)) ./ val; %

                        %             fit regression line
                        d_eval(~isfinite(d_eval)) = NaN;
                        level1 = [ones(size(level')) level'];
                        slope = level1 \ d_eval;
                        d_fit = level1 * slope;
                        level_thresh(i_interference_sigma,i_noise_mode,i_innerbw, i_flanking_phase) = (log10(dprime_threshold) - slope(1)) / slope(2);
                    end
                end
            end
        end


        % Plotting
        x.paperdata = amt_load('eurich2022', 'MqMc_2c.mat');


        data_raw = x.paperdata.paperdata(:,2:5);

        legstr = {'$+\pi/2$, SDN, $-\pi/2$'; '$+\pi/2$, ODN, $-\pi/2$';'$-\pi/2$, SDN, $+\pi/2$'; '$-\pi/2$, ODN, $+\pi/2$';'diotic'};

        col = [0 0 1; 0.5 0.5 0.5; 0.5 0.8 0.3; 0.1 .75 0.85];%col = colororder;

        datadiotic = 61.95;
        preddiotic = 62.1;


        f = figure('Position',[60   100  309.9213  193.7008]);
        hold on;

        % Predictions

        % with interference
        p1 = [fliplr(squeeze(level_thresh(2,2,:,1)))];
        p2 = [fliplr(squeeze(level_thresh(2,1,:,1)))];
        p3 = [fliplr(squeeze(level_thresh(2,2,:,2)))];
        p4 = [fliplr(squeeze(level_thresh(2,1,:,2)))];

        %without interference
        p1_no = [fliplr(squeeze(level_thresh(1,2,:,1)))];
        p2_no = [fliplr(squeeze(level_thresh(1,1,:,1)))];
        p3_no = [fliplr(squeeze(level_thresh(1,2,:,2)))];
        p4_no = [fliplr(squeeze(level_thresh(1,1,:,2)))];




        h1(1) = plot(p1,'-','Color',col(1,:),'DisplayName','none','LineWidth',2);
        hold on;
        h1(2) = plot(p2,'-','Color',col(2,:),'DisplayName','none','LineWidth',2);
        h1(3) = plot(p3,'-','Color',col(3,:),'DisplayName','none','LineWidth',2);
        h1(4) = plot(p4,'-','Color',col(4,:),'DisplayName','none','LineWidth',2);

        h1(1) = plot(p1_no,'-','Color',col(1,:),'DisplayName','none','LineStyle',':','LineWidth',2);
        h1(2) = plot(p2_no,'-','Color',col(2,:),'DisplayName','none','LineStyle',':','LineWidth',2);
        h1(3) = plot(p3_no,'-','Color',col(3,:),'DisplayName','none','LineStyle',':','LineWidth',2);
        h1(4) = plot(p4_no,'-','Color',col(4,:),'DisplayName','none','LineStyle',':','LineWidth',2);

        xticks(1:length(innerbw))

        hold on;
        plot(7,preddiotic,':','Color',col(1,:))
        plot(7,preddiotic,':','Color',col(2,:))
        plot(7,preddiotic,':','Color',col(3,:))
        plot(7,preddiotic,':','Color',col(4,:))


        xticklabels({'full','600','400','200','100','50','0'})

        % h1(5) = plot(predictions.spar.innerbw(end),62,'Color',col(1,:))

        xlabel('Inner-band Bandwidth (Hz)')
        % ylabel('Threshold / dB SPL');

        % Data
        markerstyles = {'^','v','d','s'};


        paperdata2 = data_raw;%(:,2:end);
        for n = 1:size(paperdata2,2)

            data(:,n) = paperdata2(~isnan(paperdata2(:,n)),n);
            k(n) = plot(flipud(data(:,n)),'LineStyle','none','Marker',markerstyles{n},'Color',col(n,:),'MarkerFaceColor',col(n,:),'LineWidth',2);

        end

        diotic = plot(7,datadiotic,'o','Color',[0.3 0.3 0.3],'LineWidth',2);

        ylim([50 63])
        xlim([0.5 7.3])
        ylabel('Threshold (dB)');


        grid on;
        set(gca,'GridLineStyle',':','LineWidth',0.7)


        lg = legend([k(1),k(3),diotic,k(2),k(4)], legstr([1 2 5 3 4]),'Location','Northwest', 'NumColumns',2,'FontSize',8,'interpreter','latex');
        legend boxoff
        lg.Position(1) = 0.06;
        lg.Position(2) = 0.7;
        lg.ItemTokenSize(1) = 10;
    end
end


%% ------ FIGURE 5 -----------------------------------------------------------
if flags.do_fig5
    %any calculations associated to Figure 5 in the publication
    if flags.do_plot
        fs             = 48e3;
        dur            = 300e-3; % s
        cos_rise_time  = 20e-3;
        bw             = [1000];
        fc             = 500;
        spl            = 45.5; % noise spectrum level
        noise_mode     = [1]; % 1: single delayed noise; 2: opposingly delayed noises (not used here)
        flanking_phase = [pi]; % IPD of flanking bands

        % variable stimulus parameters
        tone_level     = [40:2:70]; % target level / dB SPL

        innerbw = [0 50 100 200 300 400 600 800 1000];
        repititions    = 1:100;

        % model parameters
        mpar = arg_eurich2022;
        mpar.interference_sigma = [0 mpar.interference_sigma];
        mpar.idx_of_target_channel = 9;
        mpar.rho_max = 0.91;


        %C.H.: example for using amt_cache---------------------------------
        filename = 'processed_average_fig5';%1
        processed_average_fig5 = amt_cache('get',filename,flags.cachemode); %2

        if isempty(processed_average_fig5)%3
            parfor i_innerbw = 1:length(innerbw)
                processed_average_fig5(i_innerbw,:,:,:,:,:,:) = ...
                    local_evaluation_fig5(innerbw(i_innerbw),noise_mode,flanking_phase,tone_level,repititions,fs,dur,cos_rise_time,bw,fc,spl,mpar);
            end
            amt_cache('set',filename,processed_average_fig5);%4
        end



        % decision
        filename2 = 'dprime_fig5';%1
        dprime_fig5 = amt_cache('get',filename2,flags.cachemode);
        if isempty(dprime_fig5)
            for i_interference_sigma = 1:length(mpar.interference_sigma)
                for i_noise_mode = 1:length(noise_mode)
                    for i_innerbw = 1:length(innerbw)
                        for i_flanking_phase = 1:length(flanking_phase)
                            for i_tone_level = 1:length(tone_level)
                                dprime_fig5(i_interference_sigma,i_noise_mode,i_innerbw, i_flanking_phase,i_tone_level) = eurich2022_decision(squeeze(processed_average_fig5(i_innerbw,i_interference_sigma,i_noise_mode,i_flanking_phase,i_tone_level,:,:)), mpar);
                            end
                        end
                    end
                end
            end
            amt_cache('set',filename2,dprime_fig5);
        end

        % determine threshold
        dprime_threshold = 0.78;

        % range of assumed linear psychometric function
        dprime0 = 0.5;
        dprime1 = 10^0.4;

        for i_interference_sigma = 1:length(mpar.interference_sigma)
            for i_noise_mode = 1:length(noise_mode)
                for i_innerbw = 1:length(innerbw)
                    for i_flanking_phase = 1:length(flanking_phase)

                        dprime_temp = squeeze(dprime_fig5(i_interference_sigma,i_noise_mode,i_innerbw,i_flanking_phase,:));
                        idx0 = find(dprime_temp > dprime0 & dprime_temp < dprime1);

                        d_eval = squeeze(log10(dprime_temp(idx0)));

                        p0 = 20e-6;
                        level = tone_level(idx0);
                        %             level(:,c,b,a) = (p0 * 10.^(spar.dbspl_tone(idx0) / 20)) ./ val; %

                        %             fit regression line
                        d_eval(~isfinite(d_eval)) = NaN;
                        level1 = [ones(size(level')) level'];
                        slope = level1 \ d_eval;
                        d_fit = level1 * slope;
                        level_thresh(i_interference_sigma,i_noise_mode,i_innerbw, i_flanking_phase) = (log10(dprime_threshold) - slope(1)) / slope(2);
                    end
                end
            end
        end

        % Plotting

        % represent simulations as relative to NoSpi
        thresh_level_med_rel = level_thresh - level_thresh(:,:,1);


        % read literature data
        x.paperdata_H98_2 = amt_load('eurich2022', 'Holube98_2_ll.mat');
        x.paperdata_H98_3 = amt_load('eurich2022', 'Holube98_3_ll.mat');
        x.paperdata_SG66  = amt_load('eurich2022', 'Sondhi66_6.mat');
        x.paperdata_KC10  = amt_load('eurich2022', 'KolCul10_4_500.mat');

        % common x axis: inner band bandwidth
        paperdata_y = [0 0 0 100 100 100 141 141 141 200 200 200 283 283 283 400 400 400]';

        % prepare data Kolarik & Culling 2010
        data_KC10 = x.paperdata_KC10.paperdata_KC10;
        KolCul10_1 = data_KC10(~isnan(data_KC10(:,3)),3);
        data_KC10(:,1) = paperdata_y;


        % prepare data by Sondhi & Guttman 1966 and Holube et al. 1998
        innerbw_KC10 = [0 100 141 200 283 400]';
        innerbw_SoGu = [0 10 40 100 200 300 bw-20];
        innerbw_Hol = [0 20 50 100 150 200 400 600 800 bw-20];

        data_SG66 = NaN(size(x.paperdata_SG66.paperdata_SG66));
        data_Hol2 = NaN(size(x.paperdata_H98_2.paperdata_H98_2));
        data_Hol3 = NaN(size(x.paperdata_H98_3.paperdata_H98_3));


        for s = 2:size(data_SG66,2)
            rawSo = x.paperdata_SG66.paperdata_SG66(~isnan(x.paperdata_SG66.paperdata_SG66(:,s)),s);
            data_SG66(1:length(rawSo),s) = rawSo;
        end
        data_SG66 = data_SG66(1:7,2:end);

        for s = 2:size(data_Hol3,2)
            rawH98_2 = x.paperdata_H98_2.paperdata_H98_2(~isnan(x.paperdata_H98_2.paperdata_H98_2(:,s)),s);
            rawH98_3 = x.paperdata_H98_3.paperdata_H98_3(~isnan(x.paperdata_H98_3.paperdata_H98_3(:,s)),s);

            data_Hol2(1:length(rawH98_2),s) = rawH98_2;
            data_Hol3(1:length(rawH98_3),s) = rawH98_3;
        end

        data_Hol2 = data_Hol2(1:10,2:end);
        data_Hol3 = data_Hol3(1:10,2:end);


        % Plot
        f = figure('Position',[60   100  309.9213  193.7008]);
        col = [0 0 1; 0.5 0.5 0.5; 0.5 0.8 0.3; 0.1 .75 0.85];%col = colororder;

        hold on;


        h0(1) = plot(innerbw_KC10, KolCul10_1(2:end), 'o','Color',col(1,:),'LineWidth',2);
        ho(1).Marker = 'o';

        innerbw_SoGu(end) = 1000;
        h4 = plot(innerbw_SoGu,-data_SG66(:,2),'d','LineWidth',2);

        innerbw_Hol(end) = 1000;
        h1 = plot(innerbw_Hol,data_Hol2(:,1),'^','LineWidth',2);
        h1.Color = col(3,:);
        h1.Marker = '^';

        h2 = plot(innerbw_Hol,data_Hol3(:,1),'^','LineWidth',2);
        h2.Color = col(4,:);
        h2.Marker = 'v';

        %         broken_innerbw = [innerbw(1:end-1) 900 innerbw(end)];
        broken_koldata(1,:) =  [squeeze(thresh_level_med_rel(1,:,1:end-1)); NaN;  squeeze(thresh_level_med_rel(1,:,end))];


        %         h1(1) =
        plot(innerbw, squeeze(thresh_level_med_rel(1,:,:)),'color',[0 0 0],'Linestyle',':','LineWidth',2);
        %         h1(2) =
        plot(innerbw, squeeze(thresh_level_med_rel(2,:,:)),'color',[0 0 0],'LineWidth',2);


        % plot break
        %         x_breakline = [870:0.1:885.5];
        %         y_breakline = [-2.18,-11.80];
        %         h1(3) = plot(x_breakline, linspace(y_breakline(1),y_breakline(2),length(x_breakline)),'DisplayName','None','Color','black','LineWidth',0.5);
        %         h1(4) = plot(x_breakline+24, linspace(y_breakline(1),y_breakline(2),length(x_breakline)),'DisplayName','None','Color','black','LineWidth',0.5);

        xticks([innerbw([1 4 6 7 8]) 1000]);
        xticklabels({'0','200','400','600','800','full'})
        ylim([-16 3]);


        grid on;
        set(gca,'GridLineStyle',':','LineWidth',1)
        xlabel('Inner-band Bandwidth (Hz)');
        ylabel('Threshold Shift (dB)');


        lg =legend('Kolarik \& Culling 2010','Sondhi \& Guttman 1966','Holube et al. 1998, subj. RN','Holube et al. 1998, subj. IH','Location','NorthEast','NumColumns',1,...
            'FontSize',8,'Interpreter','Latex');
        legend boxoff
        lg.Position(1) = 0.4;
        lg.Position(2) = 0.65;
        lg.ItemTokenSize(1) = 10;

    end
end
end

function processed_average = local_evaluation_fig3(itd,noise_mode,sphase,tone_level,repititions,fs,dur,cos_rise_time,bw,fc,spl,mpar)

all_interference_sigma = mpar.interference_sigma;
for i_interference_sigma = 1:length(mpar.interference_sigma)
    for i_noise_mode = 1:length(noise_mode)
        for i_sphase = 1:length(sphase)
            for i_tone_level = 1:length(tone_level)
                for i_repititions = 1:length(repititions)


                    mpar.interference_sigma = all_interference_sigma(i_interference_sigma);
                    c_itd            = itd;
                    c_sphase         = sphase(i_sphase);
                    c_noise_mode     = noise_mode(i_noise_mode);
                    c_repitition     = repititions(i_repititions);


                    % stimulus noise only
                    c_tone_level = -inf;
                    mRef = sig_vanderheijden1999(c_itd,c_sphase,c_noise_mode,c_tone_level,fs,dur,cos_rise_time,bw,fc,spl);

                    % stimulus tone in noise
                    c_tone_level     = tone_level(i_tone_level);
                    mTest = sig_vanderheijden1999(c_itd,c_sphase,c_noise_mode,c_tone_level,fs,dur,cos_rise_time,bw,fc,spl);

                    % as this model cannot distinguish between an N0Spi and
                    % an NpiS0 threshold, model parameters are adjusted in the S0 condition for
                    % precise predictions.
                    if c_sphase == 1
                        mpar.binaural_internal_noise_sigma = 0.17;
                        mpar.rho_max = 0.86;
                        if mpar.interference_sigma > 0
                            mpar.interference_sigma = 0.65;
                        end
                        %else
                        %   mpar_temp = arg_eurich2022;

                    end


                    % run model for this condition
                    [processed(:,:,:,i_repititions)] = eurich2022(mRef,mTest,mpar);


                end


                processed_average(i_interference_sigma,i_noise_mode, i_sphase,i_tone_level,:,:) = mean(processed,4);

            end
        end
    end
end
end


function processed_average = local_evaluation_fig4(innerbw,noise_mode,flanking_phase,tone_level,repititions,fs,dur,cos_rise_time,bw,fc,spl,mpar)

all_interference_sigma = mpar.interference_sigma;
for i_interference_sigma = 1:length(mpar.interference_sigma)
    for i_noise_mode = 1:length(noise_mode)
        for i_flanking_phase = 1:length(flanking_phase)
            for i_tone_level = 1:length(tone_level)
                for i_repititions = 1:length(repititions)


                    mpar.interference_sigma = all_interference_sigma(i_interference_sigma);
                    c_innerbw        = innerbw;
                    c_itd            = 1e-3;
                    c_sphase         = 1;
                    c_noise_mode     = noise_mode(i_noise_mode);
                    c_flanking_phase = flanking_phase(i_flanking_phase);
                    c_repitition     = repititions(i_repititions);


                    % stimulus noise only
                    c_tone_level = -inf;
                    mRef = sig_marquardt2009(c_innerbw,c_flanking_phase,c_noise_mode,c_tone_level,c_itd,c_sphase,fs,dur,cos_rise_time,bw,fc,spl);
                    %                                                      c_itd,c_sphase,c_noise_mode,c_tone_level,fs,dur,cos_rise_time,bw,fc,spl
                    % stimulus tone in noise
                    c_tone_level     = tone_level(i_tone_level);
                    mTest = sig_marquardt2009(c_innerbw,c_flanking_phase,c_noise_mode,c_tone_level,c_itd,c_sphase,fs,dur,cos_rise_time,bw,fc,spl);

                    % run model for this condition

                    [processed(:,:,:,i_repititions)] = eurich2022(mRef,mTest,mpar);


                end


                processed_average(i_interference_sigma,i_noise_mode, i_flanking_phase,i_tone_level,:,:) = mean(processed,4);

            end
        end
    end
end
end

function processed_average = local_evaluation_fig5(innerbw,noise_mode,flanking_phase,tone_level,repititions,fs,dur,cos_rise_time,bw,fc,spl,mpar)

all_interference_sigma = mpar.interference_sigma;
for i_interference_sigma = 1:length(mpar.interference_sigma)
    for i_noise_mode = 1:length(noise_mode)
        for i_flanking_phase = 1:length(flanking_phase)
            for i_tone_level = 1:length(tone_level)
                for i_repititions = 1:length(repititions)


                    mpar.interference_sigma = all_interference_sigma(i_interference_sigma);
                    c_innerbw        = innerbw;
                    c_itd            = 0;
                    c_sphase         = -1;
                    c_noise_mode     = noise_mode(i_noise_mode);
                    c_flanking_phase = flanking_phase(i_flanking_phase);
                    c_repitition     = repititions(i_repititions);


                    % stimulus noise only
                    c_tone_level = -inf;
                    mRef = sig_kolarik2010(c_innerbw,c_flanking_phase,c_noise_mode,c_tone_level,c_itd,c_sphase,fs,dur,cos_rise_time,bw,fc,spl);
                    %                                                      c_itd,c_sphase,c_noise_mode,c_tone_level,fs,dur,cos_rise_time,bw,fc,spl
                    % stimulus tone in noise
                    c_tone_level     = tone_level(i_tone_level);
                    mTest = sig_kolarik2010(c_innerbw,c_flanking_phase,c_noise_mode,c_tone_level,c_itd,c_sphase,fs,dur,cos_rise_time,bw,fc,spl);

                    % run model for this condition

                    [processed(:,:,:,i_repititions)] = eurich2022(mRef,mTest,mpar);


                end


                processed_average(i_interference_sigma,i_noise_mode, i_flanking_phase,i_tone_level,:,:) = mean(processed,4);

            end
        end
    end
end
end
