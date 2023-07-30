function varargout = exp_barumerli2023(varargin)
%EXP_BARUMERLI2023 Experiments and results of Barumerli et al. (2023)
%   Usage: [] = exp_barumerli2023(flag) 
%
%   EXP_BARUMERLI2023(flag) reproduces figures and results of the study  
%   from Barumerli et al. (2023). Note: in the paper, we refer to magnitue
%   profiles (MP) and gradient profiles (GP). This new naming has been
%   introduced in a later stage. The code instead uses DTF instead of MP
%   and PGE instead of PG.
%
%   The following flags can be specified
%
%     'tab2'    Report data in Tab.2:
%               The metrics LE, PE and QE computed for each fitted model 
%               are compared to the actual performances of five listeners in
%               Majdak et al. (2010). Three different versions are tested 
%               each one with a different features space. 
%               The first considers with only binaural cues, 
%               the second combines binaural cues with spectral amplitudes
%               and the third relies on binaural cues and spectral gradients. 
% 
%     'gain_test'    Report results on prior contribution:
%                    This experiment prints the results of the analysis which
%                    employs the polar gain metric as in Ege et al. (2018) to
%                    quantify the contribution of the prior distribution. 
% 
%     'fig2'    Reproduce Fig.2:
%               Binarual feature examples. Visual comparison 
%               between the interaural cues for one subject in 
%               the frontal and horizontal plane.
%
%     'fig3'    Reproduce Fig.3:
%               Monaural feature examples. Visual comparison 
%               between the spectral cues for one subject in 
%               the median plane.
% 
%     'fig4'    Reproduce Fig.4:
%               Prediction examples. Example of the model inferring the sound
%               direction. Particularly, the plot shows the prediction of
%               the bayesian observer based on the posterior distribution.
%               Finally, such estimate is corrupted by motor noise required
%               by the listener to provide a response.
% 
%     'fig5'    Reproduce Fig.5:
%               Example of the fitted model based on spectral gradients.
%               The actual data for both lateral and polar dimensions
%               are from subject NH16 from Majdak et al. (2010). 
%               Moreover, model predictions, black dots, are: 
%               (a) based only on the likelihood function (i.e. inference driven 
%               only by sensory evidence) as in Reijniers et al. (2014) 
%               (b) Bayesian inference with both prior belief and sensory evidence. 
%               (c) Full model (i.e. Bayesian inference and sensorimotor stage. 
% 
%     'fig6'    Reproduce Fig.6:
%               Comparison between fitted models and actual data for 
%               five listeners in Majdak et al (2010).  
%               Actual (gray) and predicted (black) are the sound-localization 
%               performance metrics obtained by models based on 
%               spectral amplitues or gradients. 
%               Each row reports a different metrics: 
%               the first is about the Lateral Error (LE) as function of the
%               lateral angle, the second and the third show the 
%               Polar Error (PE) and Quadrant Error (QE), respectively as 
%               a function of the polar angle, calculated for all targets 
%               within the lateral interval [-30, 30]deg.
%
%     'fig7'    Reproduce Fig.7: see flag exp_middlebrooks1999
%
%     'fig8'    Reproduce Fig.8: see flag exp_macpherson2003
%               
%     'exp_middlebrooks1999'    Reproduce Fig.7:
%                               Predicted localization performance obtained for 
%                               the individual (Own) and non-individual (Other) HRTFs 
%                               with models based on two feature spaces.
%                               Additionally, predictions from Reijniers et al. (2014) and 
%                               Baumgartner et al. (2014) as well as actual data from 
%                               the original experiment Middlebrooks (1999) are shown. 
% 
%     'exp_macpherson2003'      Reproduce Fig.8:
%                               Effect of the spectral ripples on sound localization 
%                               performance by means of the polar error metric. 
%                               Top and bottom left panels show differences to 
%                               the reference condition in the right-most bottom panel 
%                               which reports the polar errors obtained with 
%                               broadband noise without spectral ripples. 
%                               All panels show, in addition to predictions from our models, 
%                               predictions from Reijniers et al (2014) and 
%                               Baumgartner et al. (2014) as well as actual data 
%                               from the original experiment Macpherson and Middlebrooks (2003).
% 
%   Further, cache flags (see amt_cache) and plot flags can be specified:
% 
%     'plot'    Plot the output of the experiment. This is the default.
% 
%     'no_plot'  Don't plot, only return data.
% 
%     'test'    Run one iteration for the experiment for testing code.
% 
%     'redo'    Recompute all results (it can take a while)
% 
%     'redo_fast' Recumpute all results but with less iterations. Cached files are not changed.
% 
%   Requirements: 
% 
%   1) SOFA API v1 or higher from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in auxdata/barumerli2023
% 
%   3) Statistics Toolbox and Computer Vision Toolbox for Matlab
% 
%   Examples:
% 
%   To display Fig.5 use :
% 
%     exp_barumerli2023('fig5');
% 
%   To display exp_middlebrooks1999 use :
% 
%     exp_barumerli2023('exp_middlebrooks1999');
% 
%   To display exp_macpherson2003 use :
% 
%     exp_barumerli2023('exp_macpherson2003');
%
% 
%   References:
%     R. Barumerli, P. Majdak, M. Geronazzo, D. Meijer, F. Avanzini, and
%     R. Baumgartner. A Bayesian model for human directional localization of
%     broadband static sound sources. Acta Acust., 7:12, 2023. [1]http ]
%     
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%     P. Majdak, M. J. Goupell, and B. Laback. 3-D localization of virtual
%     sound sources: Effects of visual environment, pointing method and
%     training. Atten Percept Psycho, 72:454--469, 2010.
%     
%     J. C. Middlebrooks. Virtual localization improved by scaling
%     nonindividualized external-ear transfer functions in frequency. J.
%     Acoust. Soc. Am., 106:1493--1510, 1999.
%     
%     E. A. Macpherson and J. C. Middlebrooks. Vertical-plane sound
%     localization probed with ripple-spectrum noise. J. Acoust. Soc. Am.,
%     114:430--445, 2003.
%     
%     J. Reijniers, D. Vanderleist, C. Jin, C. S., and H. Peremans. An
%     ideal-observer model of human sound localization. Biological
%     Cybernetics, 108:169--181, 2014.
%     
%     References
%     
%     1. https://doi.org/10.1051/aacus/2023006
%     
%
%   See also: demo_barumerli2023 barumerli2023 barumerli2023_featureextraction barumerli2023_metrics
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_barumerli2023.php



%   AUTHOR: Roberto Barumerli
%   Information Engineering Dept., University of Padova, Italy, 2021
%   Acoustics Research Institute, OeAW, Wien, Austria, 2023

%% ------ Check input options ---------------------------------------------

definput.import = {'amt_cache'};
definput.keyvals.MarkerSize = 6;
definput.keyvals.FontSize = 9;

definput.flags.type = {'missingflag', 'tab1', 'tab2', 'gain_test', 'fig2', 'fig3', 'fig4', 'fig5', 'fig6', 'fig7', 'fig8', 'fig9', 'exp_middlebrooks1999', 'exp_macpherson2003'};
definput.flags.plot = {'plot', 'no_plot'};

definput.flags.redo = {'no_redo_fast','redo_fast', 'redo'};
definput.flags.test = {'no_test','test'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},...
             definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.', ...
      upper(mfilename),flagnames);
end

%% ------ tab1 - fitted parameters
if flags.do_tab1
    data_majdak = data_majdak2010('Learn_M');
    data_majdak([1:5]) = [];

    calibrations = amt_load('barumerli2023', 'barumerli2023_calibration.mat');
    calibrations = calibrations.cache.value;
    
    if size(calibrations.sigma, 1) ~= length(data_majdak)
        warning('sigma values not enough for the provided subejcts')
    end

    num_calib = size(calibrations.combination,1);
    num_sub = size(calibrations.sigma, 1);
    
    for c = 1:num_calib
        fprintf("CALIBRATION %s\n", calibrations.combination{c,:})
        fprintf("ID   & PRIOR & ILD  & MON  & MOTOR\n")
        for s = 1:num_sub % subjects
            sigma_ild = calibrations.sigma(s,c).values(2);
            sigma_mon = calibrations.sigma(s,c).values(3);
            sigma_motor = calibrations.sigma(s,c).values(4);
            sigma_prior = calibrations.sigma(s,c).values(5);
            
            fprintf("%s & %.2f & %.2f & %.2f & %.2f\n", data_majdak(s).id, sigma_prior, sigma_ild, sigma_mon, sigma_motor)
        end
    end
end


%% ------ tab2 - fitted models and predicted perfomances
if flags.do_tab2
    data_majdak = data_majdak2010('Learn_M');
    data_majdak([1:5]) = [];

    calibrations = amt_load('barumerli2023', 'barumerli2023_calibration.mat');
    calibrations = calibrations.cache.value;
    
    if size(calibrations.sigma, 1) ~= length(data_majdak)
        warning('sigma values not enough for the provided subejcts')
    end

    tab2 = amt_cache('get', 'barumerli2023_tab2',flags.cachemode);
    num_exp = 300;
    num_calib = size(calibrations.combination,1);
    num_sub = size(calibrations.sigma, 1);
    
    if flags.do_redo_fast
        num_exp = 20;
        tab2 = [];
    end
    
    % Preallocation
    if isempty(tab2)
        tab2 = repmat(struct('err', ...
            struct([])),length(data_majdak), size(calibrations.combination, 1));

        for s = 1:num_sub % subjects
            fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n SUBJECT %s\n', data_majdak(s).id)

            sofa = amt_load('barumerli2023',['ARI_' upper(calibrations.name{s}) '_hrtf_M_dtf 256.sofa']);
            comb = calibrations.combination;
            for c = 1:num_calib% feature space
                fprintf('\nCOMBINATION %i ', c)
                fprintf('%s ', comb{c,:})

                [template_par, target_par] = barumerli2023_featureextraction(sofa, ...
                                        calibrations.combination{c,1});


                sigma_l = calibrations.sigma(s,c).values(1);
                sigma_l2 = calibrations.sigma(s,c).values(2);
                sigma_mon = calibrations.sigma(s,c).values(3);
                sigma_m = calibrations.sigma(s,c).values(4);
                sigma_prior = calibrations.sigma(s,c).values(5);

                if calibrations.sigma(s,c).values(3) == 0
                    sigma_mon= [];
                end

                m = barumerli2023('template', template_par, ...
                                    'target', target_par, ...
                                    'num_exp', num_exp, ...
                                    'sigma_itd', sigma_l, ...
                                    'sigma_ild', sigma_l2, ...
                                    'sigma_spectral', sigma_mon, ...
                                    'sigma_motor', sigma_m,...
                                    'MAP',...
                                    'sigma_prior', sigma_prior);

                tab2(s,c).err = barumerli2023_metrics(m,'middle_metrics');
            end
        end
        
        if ~flags.do_redo_fast
            amt_cache('set', 'barumerli2023_tab2',tab2);
        end
    end
    % metrics
    amt_disp(sprintf('CHECK CONDITIONs'))
    metric_calib = zeros(num_calib,3);
    for s = 1:num_sub
        amt_disp(sprintf('\nSUBJECT %i', s))
        mtx = data_majdak(s).mtx;
        real = barumerli2023_metrics(mtx,'middle_metrics');
        real_metric(s,1) = real.rmsL;
        real_metric(s,2) = real.rmsP;
        real_metric(s,3) = real.querr;
        for c = 1:num_calib
            % plot
            amt_disp(sprintf('\tFEATURE SPACE %s ', calibrations.combination{c,:}))
            le_tau = abs(tab2(s,c).err.rmsL - real.rmsL)/real.rmsL;
            pe_tau = abs(tab2(s,c).err.rmsP - real.rmsP)/real.rmsP;
            qe_tau_rau = abs(rau(tab2(s,c).err.querr, 1, 'PC') - rau(real.querr, 1, 'PC'))/rau(real.querr, 1, 'PC');
            fprintf('le_tau %.2f pe_tau %.2f querr_tau %.2f\n', le_tau, pe_tau, qe_tau_rau);

            metric_calib(c,1) = metric_calib(c,1) + tab2(s,c).err.rmsL;
            metric_calib(c,2) = metric_calib(c,2) + tab2(s,c).err.rmsP;
            metric_calib(c,3) = metric_calib(c,3) + tab2(s,c).err.querr;
        end
    end
    
    metric_calib = metric_calib./num_sub;
    
    % std stuff
    for s = 1:num_sub
        for c = 1:num_calib
            metric_calib_std(c,1,s) = tab2(s,c).err.rmsL;
            metric_calib_std(c,2,s) = tab2(s,c).err.rmsP;
            metric_calib_std(c,3,s) = tab2(s,c).err.querr;
        end
    end
    
    amt_disp(sprintf('\n############ TAB2 - AVERAGE OVER SUBJECTS'))
    for c = 1:size(calibrations.combination, 1)
        amt_disp(sprintf('\tFEATURE SPACE %s ', calibrations.combination{c,:}))
        fprintf('le %.2f pe %.2f querr %.2f\n', metric_calib(c,:));
        amt_disp(sprintf('le_std %.2f, pe_std %.2f, qe_std %.2f ',...
            std(metric_calib_std(c,:,:), 0, 3)));
    end
    
    amt_disp(sprintf('\n############ TAB2 - REAL SUBJECTS'))
        amt_disp(sprintf('le %.2f, pe %.2f, qe %.2f ',...
            mean(squeeze(real_metric(:,1))), mean(squeeze(real_metric(:,2))), mean(squeeze(real_metric(:,3)))));
        amt_disp(sprintf('le_std %.2f, pe_std %.2f, qe_std %.2f ',...
            std(squeeze(real_metric(:,1))), std(squeeze(real_metric(:,2))), std(squeeze(real_metric(:,3)))));
    
    
end

%% ------ gain metric to evaluate prior contribution
if flags.do_gain_test
    data_majdak = data_majdak2010('Learn_M');
    data_majdak([1:5]) = [];

    calibrations = amt_load('barumerli2023', 'barumerli2023_calibration.mat');
    calibrations = calibrations.cache.value;
    
    % Preallocation
    gain_test = amt_cache('get', 'barumerli2023_gain_test',flags.cachemode);

    % Preallocation
    if isempty(gain_test)    
        for s = 1:length(data_majdak)
            sofa = amt_load('barumerli2023',['ARI_' upper(data_majdak(s).id) '_hrtf_M_dtf 256.sofa']);
            mtx = data_majdak(s).mtx;

            [template, target] = barumerli2023_featureextraction(sofa, ...
                                    'pge', ...
                                    'targ_az', mtx(:, 1), ...
                                    'targ_el', mtx(:, 2));

            m_noprior = barumerli2023('template', template, ...
                                     'target', target, ...
                                     'num_exp', 5, ...
                                     'sigma_itd',calibrations.sigma(s,2).values(1), ...
                                     'sigma_ild', calibrations.sigma(s,2).values(2), ...
                                     'sigma_spectral', calibrations.sigma(s,2).values(3), ...
                                     'sigma_prior', [],...
                                     'sigma_motor', []);

            m_prior = barumerli2023('template', template, ...
                                     'target', target, ...
                                     'num_exp', 5, ...
                                     'sigma_itd',calibrations.sigma(s,2).values(1), ...
                                     'sigma_ild', calibrations.sigma(s,2).values(2), ...
                                     'sigma_spectral', calibrations.sigma(s,2).values(3), ...
                                     'sigma_prior', calibrations.sigma(s,2).values(5),...
                                     'sigma_motor', []);

            m_prior_motor = barumerli2023('template', template, ...
                                     'target', target, ...
                                     'num_exp', 5, ...
                                     'sigma_itd',calibrations.sigma(s,2).values(1), ...
                                     'sigma_ild', calibrations.sigma(s,2).values(2), ...
                                     'sigma_spectral', calibrations.sigma(s,2).values(3), ...
                                     'sigma_prior', calibrations.sigma(s,2).values(5),...
                                     'sigma_motor', calibrations.sigma(s,2).values(4));

            results  = {mtx, m_noprior, m_prior, m_prior_motor};

            metrics = {'gainPfront'};

            for m = 1:length(results)
                for sp = 1:length(metrics)
                    gain_test(s, m, sp) = localizationerror(results{m}, metrics{sp});
                end
            end
        end
        
        amt_cache('set', 'barumerli2023_gain_test', gain_test);
    end
    
    gain_test = squeeze(gain_test(:,:,1));
    res = mean(gain_test, 1);
    fprintf("Gains frontal plane averaged over 5 subjects\n")
    fprintf("Real:\t\t%.2f\n", res(1))
    fprintf("Only Likelihood:%.2f\n", res(2))
    fprintf("Full Model:\t%.2f\n", res(4))
end

%% ------ figure with binaural features
if flags.do_fig2
    if flags.do_redo
        fig2 = [];
    else
        fig2 = amt_cache('get','barumerli2023_fig2',flags.cachemode);
    end
    
    if isempty(fig2)
        sofa = amt_load('barumerli2023', 'ARI_NH12_hrtf_M_dtf 256.sofa');
        template = barumerli2023_featureextraction(sofa);
        c = template.coords.return_positions('horizontal-polar');

        % select directions and values with polar angle less than 2 (eye-level)
        % in the paper the line is smoother since I have been using more points
        % for interpolation
        pos = find(abs(c(:, 2)) < 2);
        itd = template.itd(pos);
        ild = template.ild(pos);

        % sort them for plotting
        [pos, idx] = sort(c(pos,1));
        itd = -itd(idx);
        ild = ild(idx);
        fig2.pos = pos;
        fig2.itd = itd;
        fig2.ild = ild;

        amt_cache('set','barumerli2023_fig2',fig2);
    end
    
    if flags.do_plot
        % figure 
        FontSize = kv.FontSize;
        figure('Units', 'points', 'Position', [200 200 245 150]);
%         colororder([[55,126,184]./255-0.2; [77,175,74]./255-0.2])

        %yyaxis left
        plot(fig2.pos, fig2.itd)
        ylabel('ITD [JND]')
        yticks(sort([-15:7.5:15]))
        hold on 
        xticks([-90:45:90])
        grid on
        yyaxis right
        plot(fig2.pos, fig2.ild)
        ylabel('ILD [dB]')
        yticks([-30:15:30])

        xlabel('Lateral angle [deg]')
        set(gca, 'FontSize', FontSize)
    end
end

%% ------ figure with monaural features
if flags.do_fig3
    if flags.do_redo
        fig3 = [];
    else
        fig3 = amt_cache('get','barumerli2023_fig3',flags.cachemode);
    end
               
    if isempty(fig3)
        data_majdak = data_majdak2010('Learn_M');
        data_majdak([1:5]) = [];

        calibrations = amt_load('barumerli2023', 'barumerli2023_calibration.mat');
        calibrations = calibrations.cache.value;

        sofa = amt_load('barumerli2023', ['ARI_' upper(data_majdak(1).id) '_hrtf_M_dtf 256.sofa']);

        % Preprocessing source information for both directions
        [dtf] = barumerli2023_featureextraction(sofa, 'template','dtf');
        [pge] = barumerli2023_featureextraction(sofa, 'template','pge');
        
        fig3.dtf = dtf;
        fig3.pge = pge;

        amt_cache('set','barumerli2023_fig3',fig3);
    end

    if flags.do_plot
        fig = figure('Units', 'points', 'Position', [200 200 245 300]);
%         (Nh, Nw, gap, marg_h, marg_w)
        [ha, ~] = tight_subplot(2,1, [.06 0.08],[.1 .05],[.15 .15]);

        hor = fig3.dtf.coords.return_positions('horizontal-polar');
        med_plane_idx = find(abs(hor(:,1))<1);% & hor(:,2)<90 & hor(:,2)>=0;
        
        t = hor(med_plane_idx,2);
        
        [t_sort, t_idx] = sort(t);
         t_sort(1) = -30;
         t_sort(end) = 210;

        f_dtf = fig3.dtf.fc./1e3;
        
        sp = 1;
        axes(ha(sp));
        X = fig3.dtf.monaural(med_plane_idx(t_idx), 1:length(fig3.dtf.fc));
        image(f_dtf,t_sort,X,'CDataMapping','scaled')
        shading flat
        title('Magnitude [dB]', 'FontWeight', 'normal')
        view(0,90);
        cdtf = colorbar('location', 'eastoutside', 'colormap', colormap(hot));
%         cpge.Label.String = 'Magnitude [dB]';
%         cdtf.Label.Rotation = 0;
%         cdtf.Label.Position = [0.689523637294769 -23.4999957765852 0];

        sp = 2;
        axes(ha(sp));
        X = fig3.pge.monaural(med_plane_idx(t_idx), 1:length(fig3.pge.fc));
        image(f_dtf(1:end-1),t_sort,X,'CDataMapping','scaled')
        shading flat
        view(0,90);
        title('Gradient [dB]', 'FontWeight', 'normal')

        cpge = colorbar('location', 'eastoutside', 'colormap', colormap(hot));
%         cpge.Label.String = 'Gradient [dB]';
%         cpge.Label.Rotation = 0;
%         cpge.Label.Position = [0.689523662839617 6.77143143245152 0];
         

       
        %% LABELs and AXIS
        set(ha, 'YDir','normal')
        set(ha, 'FontSize',kv.FontSize)
        set(ha, 'YLim', [-30, 210])

        set(ha, 'YTick', [0:90:180])
        set(ha, 'XScale', 'log')
        set(ha, 'XTick', [0, 1, 5, 10])
        set(ha, 'XTickLabel', {'0.1', '1', '5', '10'})
        set(ha, 'XLim', [0.7, f_dtf(end-1)])
        
        set(ha(1), 'XTick', [])
%         set(ha([1, 2]), 'XTick', [])
        
        set(ha([1]), 'CLim', [-35 -25])
        set(ha([2]), 'CLim', [0 6])

        ha(2).XLabel.String  = 'Frequency [kHz]';
        ha(1).YLabel.String  = 'Polar angle [deg]';
        ha(2).YLabel.String  = 'Polar angle [deg]';
        
        set(cdtf, 'Position',   [0.879510204081633 0.5795 0.022 0.35]);
        set(cpge, 'Position',   [0.879510204081633 0.123 0.022 0.35]);
        
%          saveas(fig, 'new_plots/monaural_features.eps', 'epsc')
%         print(fig, 'fig4', '-dpng', '-r600')

    end
end

%% ------ fig.4 from paper - model estimation
if flags.do_fig4
    data_majdak = data_majdak2010('Learn_M');
    data_majdak([1:5]) = [];

    calibrations = amt_load('barumerli2023', 'barumerli2023_calibration.mat');
    calibrations = calibrations.cache.value;

    fig4 = amt_cache('get', 'barumerli2023_fig4',flags.cachemode);
    
    % select one point from real data
    az_target = data_majdak(5).mtx(177, 1);
    el_target = data_majdak(5).mtx(177, 2);
    
    if isempty(fig4)   
        s = 5;
        sofa = amt_load('barumerli2023',['ARI_' upper(data_majdak(s).id) '_hrtf_M_dtf 256.sofa']);

        [template, target] = barumerli2023_featureextraction(sofa, ...
                                'pge', ...
                                'targ_az', az_target, ...
                                'targ_el', el_target); % defined in spherical coordinates

        fig4.template = template;
        fig4.target = target;

        [fig4.m, fig4.doa, fig4.target_coords] = ...
                            barumerli2023('template', template, ...
                                 'target', target, ...
                                 'num_exp', 1, ...
                                 'sigma_itd', calibrations.sigma(s,2).values(1), ...
                                 'sigma_ild', calibrations.sigma(s,2).values(2), ...
                                 'sigma_spectral', calibrations.sigma(s,2).values(3), ...
                                 'sigma_prior', calibrations.sigma(s,2).values(5),...
                                 'sigma_motor', calibrations.sigma(s,2).values(4));
          amt_cache('set', 'barumerli2023_fig4', fig4);
    end

    fig = figure('Units', 'points', 'Position', [100 100 245 300]);

    temp_c = fig4.template.coords.return_positions('cartesian');
    
    % load full sphere
    dirs = amt_load('barumerli2023','dirs.mat');
    dirs = dirs.cache.value;
    % pad posterior otherwise artifacts in the plot
    posterior = [fig4.doa.posterior, zeros(1,500)]; 
    
    [~, cbar] = plot_reijniers2014(dirs, max(log10(posterior), -10), 'FontSize', kv.FontSize);
    set(cbar, 'colormap', colormap(flipud(gray)))
    set(cbar, 'Location', 'northoutside')
    ctitle = get(cbar, 'Title');
    set(ctitle, 'String', "$log(p(\varphi|t))$")
    set(ctitle, 'Interpreter', "latex")
    set(cbar, 'Position', [0.526530612244898 0.73 0.438819197403734 0.0266666666666666])

    %% real direction
    [x,y,x_sign] = lambert_area_projection(deg2rad(az_target),deg2rad(el_target));
        
    q0 = plot(x+(1-x_sign),y,'x');
%     plot(x+(1-x_sign),y,'x');
    q0.MarkerSize = 10;
    q0.LineWidth = 1.5;
    q0.MarkerEdgeColor = [228,0,28]./255;
    
    %% subject estimate
    % select one point from real data
    az_sbj = deg2rad(data_majdak(5).mtx(177, 3));
    el_sbj = deg2rad(data_majdak(5).mtx(177, 4));
    [x,y,x_sign] = lambert_area_projection(az_sbj,el_sbj);
        
    q1 = plot(x+(1-x_sign),y,'+');
    q1.MarkerSize = 10;
    q1.LineWidth = 1.5;
    q1.Color = [55,126,184]./255;
    
    %% estimated direction by bayesian observer
    [~, i] = max(posterior);
    [az,el]=cart2sph(temp_c(i,1),temp_c(i,2),temp_c(i,3));
    [x,y,x_sign] = lambert_area_projection(az,el);

    q2 = plot(x+(1-x_sign),y,'+');
    q2.MarkerSize = 10;
    q2.LineWidth = 1.5;
    q2.Color = [77,175,74]./255;
    
    %% final estimate corrupted by sensorymotor scatter
    est = fig4.doa.estimations;
    [az,el]=cart2sph(est(:,:,1),est(:,:,2),est(:,:,3));
    
    % lambert equal area projection
    [x,y,x_sign] = lambert_area_projection(az,el);

    q3 = plot(x+(1-x_sign),y,'+');
    q3.MarkerSize = 10;
    q3.LineWidth = 1.5;
    q3.Color = [255,127,0]./255; % [0 0.4470 0.7410];

    l = legend(["", repmat("", 1, 29), ...
                "Source direction", "Human response", "Model estimate", "Model response"], ...
                'Interpreter', 'latex', ...
                'Position',[0.0303890612964734 0.691896928826252 0.433163763552296 0.188103071173748]);

%     saveas(fig, 'new_plots/single_trial.eps', 'epsc')
%     print(fig, 'fig5', '-dpng', '-r600')
end

%% -------------- fig 5 paper - individual perforances
if flags.do_fig5
    data_majdak = data_majdak2010('Learn_M');
    data_majdak([1:5]) = [];

    calibrations = amt_load('barumerli2023', 'barumerli2023_calibration.mat');
    calibrations = calibrations.cache.value;
   
    if size(calibrations.sigma, 1) ~= length(data_majdak)
        warning('sigma values not enough for the provided subejcts')
    end

    fig5 = amt_cache('get', 'barumerli2023_fig5',flags.cachemode);

    % Preallocation
    if isempty(fig5)
        for s = 1:length(data_majdak)
            sofa = amt_load('barumerli2023', ['ARI_' upper(data_majdak(s).id) '_hrtf_M_dtf 256.sofa']);
            mtx = data_majdak(s).mtx;

            fprintf("SUBJECT %s\n", data_majdak(s).id);

            %% DTF
            calibs_dtf = calibrations.sigma(s,1); % select dtf

            calibs.sigma = calibs_dtf.values(1,1:3);
            calibs.motor_sigma = calibs_dtf.values(1,4);
            calibs.prior = calibs_dtf.values(1,5);
            [template, target] = barumerli2023_featureextraction(sofa, ...
                                                'dtf', ...
                                                'targ_az', mtx(:, 1), ...
                                                'targ_el', mtx(:, 2));

            m_motor_dtf = barumerli2023('template', template, ...
                                         'target', target, ...
                                         'num_exp', 300, ...
                                         'sigma_itd', calibs.sigma(1), ...
                                         'sigma_ild', calibs.sigma(2), ...
                                         'sigma_spectral', calibs.sigma(3), ...
                                         'sigma_prior', calibs.prior,...
                                         'sigma_motor', calibs.motor_sigma);

            %% PGE  
            calibs_pge = calibrations.sigma(s,2); % select dtf

            calibs.sigma = calibs_pge.values(1,1:3);
            calibs.motor_sigma = calibs_pge.values(1,4);
            calibs.prior = calibs_pge.values(1,5);
            [template, target] = barumerli2023_featureextraction(sofa, ...
                            'pge', ...
                            'targ_az', mtx(:, 1), ...
                            'targ_el', mtx(:, 2));

            m_motor_pge = barumerli2023('template', template, ...
                             'target', target, ...
                             'num_exp', 300, ...
                             'sigma_itd', calibs.sigma(1), ...
                             'sigma_ild', calibs.sigma(2), ...
                             'sigma_spectral', calibs.sigma(3), ...
                             'MAP',...
                             'sigma_prior', calibs.prior,...
                             'sigma_motor', calibs.motor_sigma);

            
            results(s,:) = {m_motor_dtf, m_motor_pge};
        end
        
        amt_cache('set', 'barumerli2023_fig5', results);
    end

    %% COMPUTE RMS 
    lat = [-90, -40, 0, 40, 90];
    lat_label = [-65, -20, 20, 65];

    pol =     [-30 30 150 210];
    pol_label = [0 90 180];

    % check 
    rmsL_dtf = zeros(length(data_majdak), length(lat)-1);
    rmsL_pge = zeros(length(data_majdak), length(lat)-1);
    rmsL_real = zeros(length(data_majdak), length(lat)-1);

    rmsP_dtf = zeros(1, length(pol)-1);
    rmsP_pge = zeros(1, length(pol)-1);
    rmsP_real = zeros(1, length(pol)-1);

    querr_dtf= zeros(1, length(pol)-1);
    querr_pge= zeros(1, length(pol)-1);
    querr_real = zeros(1, length(pol)-1);

    for s = 1:length(data_majdak)
        results = fig5(s,:);
        
        m_dtf = results{1};
        m_pge = results{2};
        mtx = data_majdak(s).mtx;
        
        for i=1:length(lat)-1
            % real
            mtx_temp = mtx((mtx(:, 5) > lat(i) & (mtx(:, 5) < lat(i+1))),:);
            rmsL_real(s,i) = localizationerror(mtx_temp, 'rmsL');

            % simulation
            m_temp = m_pge((m_pge(:, 5) > lat(i) & (m_pge(:, 5) < lat(i+1))),:);
            rmsL_pge(s,i) = localizationerror(m_temp, 'rmsL');
            m_temp = m_dtf((m_dtf(:, 5) > lat(i) & (m_dtf(:, 5) < lat(i+1))),:);
            rmsL_dtf(s,i) = localizationerror(m_temp, 'rmsL');
        end

        mtx(abs(mtx(7,:)) <= 30,:)=[];
        m_pge(abs(m_pge(7,:)) <= 30,:)=[];
        m_dtf(abs(m_dtf(7,:)) <= 30,:)=[];
        for i=1:length(pol)-1
            m_temp = mtx((mtx(:, 6) > pol(i) & (mtx(:, 6) < pol(i+1))),:);
            rmsP_real(s,i) = localizationerror(m_temp, 'rmsPmedianlocal');
            querr_real(s,i) = localizationerror(m_temp, 'querrMiddlebrooks');

            m_temp = m_pge((m_pge(:, 6) > pol(i) & (m_pge(:, 6) < pol(i+1))),:);
            rmsP_pge(s,i) = localizationerror(m_temp, 'rmsPmedianlocal');    
            querr_pge(s,i) = localizationerror(m_temp, 'querrMiddlebrooks');    

            m_temp = m_dtf((m_dtf(:, 6) > pol(i) & (m_dtf(:, 6) < pol(i+1))),:);
            rmsP_dtf(s,i) = localizationerror(m_temp, 'rmsPmedianlocal');    
            querr_dtf(s,i) = localizationerror(m_temp, 'querrMiddlebrooks');    
        end
    end
    
    
    %% FIGURE
    FontSize = 9;
%     fig = figure('Units', 'Points', 'Position', [1e3 1e3 510 200]);
    fig = figure('Units', 'Points', 'Position', [1 1 510 300]);
    [ha, ~] = tight_subplot(3, length(data_majdak), [.1 0.01],[.15 0.08],[.09 .01]);

    xangle = 0;
    sp = 1;
    
    Size = 36;
    
    for s = 1:length(data_majdak)
        marker_dtf = 's';
        marker_pge = 'v';

        % lateral
        axes(ha(sp))
        ha(sp).Position(2) = 0.6675;
        ax = gca;
        scatter(lat_label, rmsL_real(s,:), Size, 0.8*[1 1 1], 'filled')
        hold on
        scatter(lat_label, rmsL_dtf(s,:), Size, 0.2*[1 1 1], marker_dtf)
        scatter(lat_label, rmsL_pge(s,:), Size, 0.2*[1 1 1], marker_pge)
        set(gca, 'YLim', [0 20], 'Ytick', [0,10,20], 'Xtick', lat_label, 'Xlim', [-90, 90],'FontSize',FontSize)
        title({upper(data_majdak(s).id)})
        grid on

        % polar
        axes(ha(sp+length(data_majdak)))
        ha(sp+length(data_majdak)).Position(4) = ha(sp+length(data_majdak)).Position(4);
        ha(sp+length(data_majdak)).Position(2) = 0.37;
        scatter(pol_label, rmsP_real(s,:), Size, 0.8*[1 1 1], 'filled')
        hold on
        scatter(pol_label, rmsP_dtf(s,:), Size, 0.2*[1 1 1], marker_dtf)
        scatter(pol_label, rmsP_pge(s,:), Size, 0.2*[1 1 1], marker_pge)
        set(gca, 'Xtick', pol_label, 'Xlim', [-90, 270], 'Ylim', [0 60], 'Ytick', [0 30 60],'FontSize',FontSize)
        xtickangle(xangle)
        grid on

        % querr
        axes(ha(sp+2*length(data_majdak)))
        ha(sp+2*length(data_majdak)).Position(2) = 0.12;

        scatter(pol_label, querr_real(s,:), Size, 0.8*[1 1 1], 'filled')
        hold on
        scatter(pol_label, querr_dtf(s,:), Size, 0.2*[1 1 1], marker_dtf)
        scatter(pol_label, querr_pge(s,:), Size, 0.2*[1 1 1], marker_pge)
        set(gca, 'Xtick', pol_label, 'Xlim', [-90, 270], 'Ylim', [-5 40], 'Ytick', [0, 20,40], 'FontSize',FontSize)
        xtickangle(xangle)
        grid on

        sp = sp + 1;
    end

    set(ha([2:5, 7:10, 12:15]),'YTickLabel','')

    set(ha([6:10]),'XTickLabel','')
    for i=6:10
        ha(i).Title.String  = [];
    end

    ha(3).XLabel.String  = '      Lateral angle [deg]';
    ha(3).XLabel.FontSize = FontSize;

    
    ha(13).XLabel.String  = '      Polar angle [deg]';
    ha(13).XLabel.FontSize = FontSize;

    for i=1:5
        ha(i).TitleFontWeight = 'normal';
    end

    x_pos = -9;
    ha(1).YLabel.String  = {'Lateral';'error [deg]'};
%     ha(1).YLabel.Position = ha(1).YLabel.Position + [x_pos, 0, 0];
    ax.YLabel.FontSize = FontSize;
    ha(6).YLabel.String  = {'Polar';'error [deg]'};
%     ha(6).YLabel.Position = ha(6).YLabel.Position + [x_pos, 0, 0];
    ax.YLabel.FontSize = FontSize;
    ha(11).YLabel.String  = {'Quadrant';'error [%]'};
%     ha(11).YLabel.Position = ha(11).YLabel.Position + [x_pos, 0, 0];

    ax.YLabel.FontSize = FontSize;

    legend({'Actual data', 'MP variant', 'GP variant'}, ...
        'Orientation','horizontal', 'Interpreter', 'latex', ...
        'Position', [0.330346204909391 0.92194871928753 0.412033794094542 0.0652307678919574])

%       saveas(fig, 'new_plots/model_estimations_sectors.eps', 'epsc')
%         print(fig, 'fig6', '-dpng', '-r600')

end

%% ------ fig6 from paper - model stages
if flags.do_fig6
    data_majdak = data_majdak2010('Learn_M');
    data_majdak([1:5]) = [];

    calibrations = amt_load('barumerli2023', 'barumerli2023_calibration.mat');
    calibrations = calibrations.cache.value;
    
    if size(calibrations.sigma, 1) ~= length(data_majdak)
        warning('sigma values not enough for the provided subejcts')
    end

    fig6 = amt_cache('get', 'barumerli2023_fig6',flags.cachemode);

    % Preallocation
    if isempty(fig6)
        s = 3; % select subject 3

        calibrations = calibrations.sigma(s,2); % select pge
        
        sofa = amt_load('barumerli2023',['ARI_' upper(data_majdak(s).id) '_hrtf_M_dtf 256.sofa']);
        mtx = data_majdak(s).mtx;

        [template, target] = barumerli2023_featureextraction(sofa, ...
                                'pge', ...
                                'targ_az', mtx(:, 1), ...
                                'targ_el', mtx(:, 2));


        m_noprior = barumerli2023('template', template, ...
                                 'target', target, ...
                                 'num_exp', 1, ...
                                 'sigma_itd',calibrations.values(1), ...
                                 'sigma_ild', calibrations.values(2), ...
                                 'sigma_spectral', calibrations.values(3), ...
                                 'sigma_prior', [],...
                                 'sigma_motor', []);

        m_prior = barumerli2023('template', template, ...
                                 'target', target, ...
                                 'num_exp', 1, ...
                                 'sigma_itd',calibrations.values(1), ...
                                 'sigma_ild', calibrations.values(2), ...
                                 'sigma_spectral', calibrations.values(3), ...
                                 'sigma_prior', calibrations.values(5),...
                                 'sigma_motor', []);

        m_prior_motor = barumerli2023('template', template, ...
                                 'target', target, ...
                                 'num_exp', 1, ...
                                 'sigma_itd',calibrations.values(1), ...
                                 'sigma_ild', calibrations.values(2), ...
                                 'sigma_spectral', calibrations.values(3), ...
                                 'sigma_prior', calibrations.values(5),...
                                 'sigma_motor', calibrations.values(4));
                             
        fig6  = {m_noprior, m_prior, m_prior_motor, mtx};
        
        amt_cache('set', 'barumerli2023_fig6', fig6);
    end

    m_noprior = fig6{1,1};
    m_prior = fig6{1,2};
    m_prior_motor = fig6{1,3};
    mtx = fig6{1,4};
    
    %% FIGURE
    Size = 10;
    fig = figure('Units', 'points', 'Position', [10 10 245 400]);
    [ha, ~] = tight_subplot(3, 2, [.07 0.1], [.07 .07],[.13 .04]);

    %% SUBPLOTs
    sp = 1;
    
    axes(ha(sp));
    plot([-100 100], [-100 100], '-', 'Color', [1 1 1 1]*0.5)
    hold on
    p1=scatter(mtx(:, 5), mtx(:, 7), Size, 0.6*[1 1 1], 'filled');
    p2=scatter(m_noprior(:, 5), m_noprior(:, 7),  Size, 0.2*[1 1 1]); 

    sp = sp + 1;
    axes(ha(sp));
    plot([-90 270], [-90 270], '-', 'Color', [1 1 1 1]*0.5)
    hold on
    plot([-90 270], [-90 270]+90, '--', 'Color', [0 0 1 1]*0.5)
    plot([-90 270], [-90 270]-90, '--', 'Color', [0 0 1 1]*0.5)
    scatter(mtx(:, 6), mtx(:, 8), Size, 0.6*[1 1 1], 'filled');
    m_noprior(~((m_noprior(:, 6) > -30) & (m_noprior(:, 6) < 210)),:)=[];
    scatter(m_noprior(:, 6), m_noprior(:, 8),  Size, 0.2*[1 1 1]); 
%     grid on

    
    %% PRIOR
    sp = sp + 1;
    axes(ha(sp));
    plot([-100 100], [-100 100], '-', 'Color', [1 1 1 1]*0.5)
    hold on
    scatter(mtx(:, 5), mtx(:, 7), Size, 0.6*[1 1 1], 'filled');
    scatter(m_prior(:, 5), m_prior(:, 7), Size, 0.2*[1 1 1]); 
%     grid on  
    
    sp = sp + 1;
    axes(ha(sp));
    plot([-90 270], [-90 270], '-', 'Color', [1 1 1 1]*0.5)
    hold on
    plot([-90 270], [-90 270]+90, '--', 'Color', [0 0 1 1]*0.5)
    plot([-90 270], [-90 270]-90, '--', 'Color', [0 0 1 1]*0.5)
    scatter(mtx(:, 6), mtx(:, 8), Size, 0.6*[1 1 1], 'filled');
    hold on
    m_prior(~((m_prior(:, 6) > -30) & (m_prior(:, 6) < 210)),:)=[];
    scatter(m_prior(:, 6), m_prior(:, 8), Size, 0.2*[1 1 1]); 
%     grid on
    
    %% FULL MODEL
    sp = sp + 1;
    axes(ha(sp));
    plot([-100 100], [-100 100], '-', 'Color', [1 1 1 1]*0.5)
    hold on
    scatter(mtx(:, 5), mtx(:, 7), Size, 0.6*[1 1 1], 'filled');
    scatter(m_prior_motor(:, 5), m_prior_motor(:, 7), Size, 0.2*[1 1 1]); 
%     grid on
    
    sp = sp + 1;
    axes(ha(sp));
    plot([-90 270], [-90 270], '-', 'Color', [1 1 1 1]*0.5)
    hold on
    plot([-90 270], [-90 270]+90, '--', 'Color', [0 0 1 1]*0.5)
    plot([-90 270], [-90 270]-90, '--', 'Color', [0 0 1 1]*0.5)
    scatter(mtx(:, 6), mtx(:, 8), Size, 0.6*[1 1 1], 'filled');
    hold on
    m_prior_motor(~((m_prior_motor(:, 6) > -30) & (m_prior_motor(:, 6) < 210)),:)=[];
    scatter(m_prior_motor(:, 6), m_prior_motor(:, 8), Size, 0.2*[1 1 1]);  

%     grid on
    
    %% FRONTAL POSITION 
    for i=1:6
        axes(ha(i));
        plot(0,0,'r+')
    end
    
    %% LABELs and AXIS
    for i=1:2:5
        axis(ha(i), 'equal')
        set(ha(i), 'YTick',[-90, 0, 90], 'XTick',[-90, 0, 90],...
            'XLim', [-100, 100], 'YLim', [-100, 100])
    end 

    for i=2:2:6
        axis(ha(i), 'equal')
        set(ha(i), 'YTick',[-90, 0, 90, 180, 270],...
            'XTick',[-90, 0, 90, 180, 270], 'XLim', [-90, 270], 'YLim', [-90, 270])
        xtickangle(ha(i),0) 
    end 

    set(ha([1:4]),'XTickLabel','')
    set(ha,'FontSize',kv.FontSize-1)

    ha(1).Title.String = "Lateral";
    ha(2).Title.String = "Polar";

    ha(1).TitleFontWeight = 'normal';
    ha(2).TitleFontWeight = 'normal';

    ha(1).TitleFontSizeMultiplier = 1;
    ha(2).TitleFontSizeMultiplier = 1;
    
    ha(5).XLabel.String  = '      Target angle [deg]';
    ha(6).XLabel.String  = '      Target angle [deg]';

    ha(1).YLabel.String  = 'Response angle [deg]';
    ha(3).YLabel.String  = 'Response angle [deg]';
    ha(5).YLabel.String  = 'Response angle [deg]';

    titles = {'a) Sensory evidence only', ...
                'b) including prior beliefs', ...
                'c) including response noise'};

    yt = [0.95, 0.63, 0.316];

    for i=1:3
        annotation(fig,'textbox',...
        [0.1 yt(i) 0.8 0.05],...
        'String',titles(i),...
        'LineStyle','none', ...
        'FontWeight', 'normal', ...
        'FontSize', kv.FontSize + 1, ...
        'HorizontalAlignment', 'left');
    end
    
    legend([p1 p2], {'Actual', 'Simulated'}, ...
        'Orientation','vertical', ...
        'Interpreter', 'latex', ...
        'Position',[0.652990351141813 0.62463162317031 0.305810393087726 0.0533707851774237]);
    
%     saveas(fig, 'new_plots/model_stages.eps', 'epsc')
%     print(fig, 'fig7', '-dpng', '-r300')
end

%% ------ middlebrooks -------------------------------------
if flags.do_fig7 || flags.do_exp_middlebrooks1999
    data_majdak = data_majdak2010('Learn_M');
    data_majdak([1:5]) = [];

    exp_middlebrooks = [];

    if ~flags.do_redo
        exp_middlebrooks = amt_cache('get', ...
            'exp_middlebrooks1999',flags.cachemode);
    end
    
    calibrations = amt_load('barumerli2023', 'barumerli2023_calibration.mat');
    calibrations = calibrations.cache.value;
    if size(calibrations.sigma, 1) ~= length(calibrations.name)
        error('sigma values not enough for the provided subejcts')
    end
    
    % remove feature spaces without monaural cues
    monaural_none_idx = find(strcmp(calibrations.combination, 'monaural_none'));
    if monaural_none_idx > 0
        calibrations.combination(monaural_none_idx) = [];
        calibrations.sigma(:,monaural_none_idx) = [];
    end 
    
    % setting 
    sbj_num = length(calibrations.name);
    cal_num = size(calibrations.combination, 1);
    num_exp = 50;
    
    if flags.do_redo_fast
        num_exp = 5;
        exp_middlebrooks = [];
    end
    
    if flags.do_test
        num_exp = 1;
        cal_num = 1;
        sbj_num = 1;
        exp_middlebrooks = [];
    end
    
    if isempty(exp_middlebrooks)
        % preprocess templates for each user
        amt_disp('Processing subjects'' templates');

        for s = 1:sbj_num
            amt_disp(['Pre-processing subject #' num2str(s)]);

            sofa = amt_load('barumerli2023',['ARI_' upper(calibrations.name{s}) '_hrtf_M_dtf 256.sofa']);

            for c = 1:cal_num
                [template(c,s), target(c,s)] = ...
                    barumerli2023_featureextraction(sofa, ...
                                            calibrations.combination{c,1});
            end
        end
        
        % preallocation for results
        amt_disp('Allocating memory for results');
        estimations = struct('m', []);
        estimations = repmat(estimations, cal_num, ...
            sbj_num, sbj_num); % all vs all
        
        for c = 1:cal_num
            amt_disp(sprintf('Combination #%i', c));
            for s = 1:sbj_num
                amt_disp(sprintf('\tSubject #%i', s));
                
                assert(length(calibrations.sigma(s,c).values) == 5, ...
                    'something is wrong with the calibration file')
                sigma_l = calibrations.sigma(s,c).values(1);
                sigma_l2 = calibrations.sigma(s,c).values(2);
                sigma_mon = calibrations.sigma(s,c).values(3);
                sigma_m = calibrations.sigma(s,c).values(4);
                sigma_prior = calibrations.sigma(s,c).values(5);
                
                for j = 1:sbj_num
                    amt_disp(num2str(j));
                    estimations(c, s, j).m = ... 
                        barumerli2023('template', template(c, s),...
                                            'target', target(c, j), ...
                                            'num_exp', num_exp, ...
                                            'sigma_itd', sigma_l, ...
                                            'sigma_ild', sigma_l2, ...
                                            'sigma_spectral', sigma_mon,...
                                            'sigma_motor', sigma_m, ...
                                            'sigma_prior', sigma_prior);
                end
            end
        end

        % compute metrics
        for c = 1:size(estimations, 1)
            for i = 1:size(estimations, 2)
                for j = 1:size(estimations, 3)
                    metrics(c, i, j) = barumerli2023_metrics(estimations(c, i, j).m, 'middle_metrics'); 
                end
            end
        end
        
        exp_middlebrooks = metrics;
        
        if ~flags.do_redo_fast && ~flags.do_test
            amt_cache('set','exp_middlebrooks1999',exp_middlebrooks);
        end
    end
    
    metrics_all = exp_middlebrooks;
    metrics_all(strcmp(calibrations.combination(:,1), 'monaural_none'), :, :) = [];
    num_calib = size(metrics_all, 1);
        quants = [0,0.05,0.25,0.5,0.75,0.95,1];
    
    % iterate over calibrations
    for c=1:num_calib
        metrics = squeeze(metrics_all(c,:,:));
        
        % aggregate metrics
        ns = size(metrics,1);
        own = logical(eye(ns));
        other = not(own);

        % code similar to baumgartner2014 - fig9
        le_own(c,1).quantiles = quantile([metrics(own).rmsL], quants);
        lb_own(c,1).quantiles = quantile([metrics(own).accL], quants);
        qe_own(c,1).quantiles = quantile([metrics(own).querr], quants);
        pe_own(c,1).quantiles = quantile([metrics(own).rmsP], quants);
        pb_own(c,1).quantiles = quantile([metrics(own).accP], quants);
        le_own(c,1).mean = mean([metrics(own).rmsL]);
        lb_own(c,1).mean = mean([metrics(own).accL]);
        qe_own(c,1).mean = mean([metrics(own).querr]);
        pe_own(c,1).mean = mean([metrics(own).rmsP]);
        pb_own(c,1).mean = mean([metrics(own).accP]);

        le_other(c,1).quantiles = quantile([metrics(other).rmsL], quants);
        lb_other(c,1).quantiles = quantile([metrics(other).accL], quants);
        qe_other(c,1).quantiles = quantile([metrics(other).querr], quants);
        pe_other(c,1).quantiles = quantile([metrics(other).rmsP], quants);
        pb_other(c,1).quantiles = quantile([metrics(other).accP], quants);
        le_other(c,1).mean = mean([metrics(other).rmsL]);
        lb_other(c,1).mean = mean([metrics(other).accL]);
        qe_other(c,1).mean = mean([metrics(other).querr]);
        pe_other(c,1).mean = mean([metrics(other).rmsP]);
        pb_other(c,1).mean = mean([metrics(other).accP]);
    end
    
    % load reference data
    data_middle = data_middlebrooks1999;

    % baumgartner data
    data_baum_temp = exp_baumgartner2014('fig9', 'no_plot');
    data_baum.qe_pool = data_baum_temp(1).qe;
    data_baum.pe_pool = data_baum_temp(1).pe;
    data_baum.pb_pool = data_baum_temp(1).pb;

    ns = size(data_baum.pe_pool,1);
    own = eye(ns) == 1;
    other = not(own);
    data_baum.pb_pool = abs(data_baum.pb_pool);
    data_baum.qe_own.quantiles = quantile(data_baum.qe_pool(own),quants);
    data_baum.pe_own.quantiles = quantile(data_baum.pe_pool(own),quants);
    data_baum.pb_own.quantiles = quantile(data_baum.pb_pool(own),quants);
    data_baum.qe_own.mean = mean(data_baum.qe_pool(own));
    data_baum.pe_own.mean = mean(data_baum.pe_pool(own));
    data_baum.pb_own.mean = mean(data_baum.pb_pool(own));

    data_baum.qe_other.quantiles = quantile(data_baum.qe_pool(other),quants);
    data_baum.pe_other.quantiles = quantile(data_baum.pe_pool(other),quants);
    data_baum.pb_other.quantiles = quantile(data_baum.pb_pool(other),quants);
    data_baum.qe_other.mean = mean(data_baum.qe_pool(other));
    data_baum.pe_other.mean = mean(data_baum.pe_pool(other));
    data_baum.pb_other.mean = mean(data_baum.pb_pool(other));
    
    % reijniers
    data_reij_temp = exp_reijniers2014('fig2_barumerli2020forum', 'no_plot');

    ns = size(data_reij_temp,1);
    own = logical(eye(ns));
    other = not(own);
    data_reij.le_own.quantiles = quantile([data_reij_temp(own).rmsL], quants);
    data_reij.lb_own.quantiles = quantile([data_reij_temp(own).accL], quants);
    data_reij.qe_own.quantiles = quantile([data_reij_temp(own).querr], quants);
    data_reij.pe_own.quantiles = quantile([data_reij_temp(own).rmsP], quants);
    data_reij.pb_own.quantiles = quantile([data_reij_temp(own).accP], quants);
    data_reij.le_own.mean = mean([data_reij_temp(own).rmsL]);
    data_reij.lb_own.mean = mean([data_reij_temp(own).accL]);
    data_reij.qe_own.mean = mean([data_reij_temp(own).querr]);
    data_reij.pe_own.mean = mean([data_reij_temp(own).rmsP]);
    data_reij.pb_own.mean = mean([data_reij_temp(own).accP]);

    data_reij.le_other.quantiles = quantile([data_reij_temp(other).rmsL], quants);
    data_reij.lb_other.quantiles = quantile([data_reij_temp(other).accL], quants);
    data_reij.qe_other.quantiles = quantile([data_reij_temp(other).querr], quants);
    data_reij.pe_other.quantiles = quantile([data_reij_temp(other).rmsP], quants);
    data_reij.pb_other.quantiles = quantile([data_reij_temp(other).accP], quants);
    data_reij.le_other.mean = mean([data_reij_temp(other).rmsL]);
    data_reij.lb_other.mean = mean([data_reij_temp(other).accL]);
    data_reij.qe_other.mean = mean([data_reij_temp(other).querr]);
    data_reij.pe_other.mean = mean([data_reij_temp(other).rmsP]);
    data_reij.pb_other.mean = mean([data_reij_temp(other).accP]);
    
    % plot
    if flags.do_plot
%         calib_plot_order = [3,1,2]; %[lat, dtf, pge]
        calib_plot_order = [1,2]; %[dtf, pge]
        % spacing
        dx = 0.11;
        % multiplier for horizontal shift
        middle_off = 2;
        cdist_init = -1;
        reij_off = -3+1;
        baum_off = -2+1;
        
        Marker = 's-';
        LineColor = [[0.9290, 0.6940, 0.1250]; ...
            [0.4940, 0.1840, 0.5560]; ...
            [0.4660, 0.6740, 0.1880]; ...
            [0.3010, 0.7450, 0.9330]; ...
            [0.6350, 0.0780, 0.1840]];
        data_middle.Marker = 'ko-';
        data_middle.LineColor = 'k';%[1 1 1]*0.3;
        
%         data_majdak.Marker = 'b^-';
%         data_majdak.LineColor = 'b';%[0 0 1]*0.3;
        
        data_baum.Marker = 'd-';
        data_baum.LineColor = [0.8500 0.3250 0.0980];
        
        data_reij.Marker = 'v-';
        data_reij.LineColor = [0 0.4470 0.7410];
        
        mFig = figure;
        mFig.Units = 'points';
        mFig.Position = [0, 0, 510, 300];
        tile_left = 0.02;%[left bottom width height]
        tile_width = 0.32;

        %% SUBPLOT 1
        sp = 1; 
        subplot(1, 3, sp);
        ax=gca;
        ax.OuterPosition(1) = tile_left;
        ax.OuterPosition(3) = tile_width;
        ax.Position(3) = 0.26;
        
        % reference
        local_middlebroxplot(ax, 1-middle_off*dx,data_middle.le_own, data_middle.Marker, kv.MarkerSize, data_middle.LineColor, data_middle.LineColor);
        local_middlebroxplot(ax, 2-middle_off*dx,data_middle.le_other, data_middle.Marker, kv.MarkerSize, data_middle.LineColor, data_middle.LineColor);
        
        % baseline
%         local_middlebroxplot(ax, 1-majdak_off*dx,data_majdak.le, data_majdak.Marker, kv.MarkerSize, data_majdak.LineColor, data_majdak.LineColor);
        
        % simulation
        cdist = cdist_init;
        for c=calib_plot_order
            local_middlebroxplot(ax, 1+cdist*dx,le_own(c,1), Marker, kv.MarkerSize, LineColor(c,:), 'w');
            local_middlebroxplot(ax, 2+cdist*dx,le_other(c,1), Marker, kv.MarkerSize, LineColor(c,:), 'w');
            cdist = cdist+1;
        end
        
        % reijniers2014
        local_middlebroxplot(ax, 1-reij_off*dx,data_reij.le_own, data_reij.Marker, kv.MarkerSize, data_reij.LineColor, 'w');
        local_middlebroxplot(ax, 2-reij_off*dx,data_reij.le_other, data_reij.Marker, kv.MarkerSize, data_reij.LineColor, 'w');
%         
        ylabel('Lateral error [deg]','FontSize',kv.FontSize)
        set(gca,'YLim',[0 45],'YTick', 0:10:40,'XLim',[0.5 2.5],...
          'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
            'TickLength',2*get(gca,'TickLength'))
        
        %% SUBPLOT 2
        sp = sp +1;
        subplot(1, 3, sp);
        ax=gca;
        ax.OuterPosition(1) = tile_left*2 + tile_width;
        ax.OuterPosition(3) = tile_width;
        ax.Position(3) = 0.26;

        % reference
        local_middlebroxplot(ax, 1-middle_off*dx, data_middle.pe_own, data_middle.Marker, kv.MarkerSize, data_middle.LineColor, data_middle.LineColor);
        local_middlebroxplot(ax, 2-middle_off*dx, data_middle.pe_other, data_middle.Marker, kv.MarkerSize, data_middle.LineColor, data_middle.LineColor);
        
        % baseline
%         local_middlebroxplot(ax, 1-majdak_off*dx, data_majdak.pe, data_majdak.Marker, kv.MarkerSize, data_majdak.LineColor, data_majdak.LineColor);

        % simulations
        cdist = cdist_init;
        for c=calib_plot_order
            local_middlebroxplot(ax, 1+cdist*dx, pe_own(c,1), Marker, kv.MarkerSize, LineColor(c,:),'w');
            local_middlebroxplot(ax, 2+cdist*dx, pe_other(c,1), Marker, kv.MarkerSize, LineColor(c,:),'w');
            cdist = cdist + 1;
        end
        
        % reijniers2014
        local_middlebroxplot(ax, 1-reij_off*dx, data_reij.pe_own,data_reij.Marker, kv.MarkerSize, data_reij.LineColor,'w');
        local_middlebroxplot(ax, 2-reij_off*dx, data_reij.pe_other,data_reij.Marker, kv.MarkerSize, data_reij.LineColor,'w');

        % baumgartner2014
        local_middlebroxplot(ax, 1-baum_off*dx, data_baum.pe_own,data_baum.Marker, kv.MarkerSize, data_baum.LineColor,'w');
        local_middlebroxplot(ax, 2-baum_off*dx, data_baum.pe_other,data_baum.Marker, kv.MarkerSize, data_baum.LineColor,'w');
        
        ylabel('Polar error [deg]','FontSize',kv.FontSize)
        set(gca,'YLim',[0 65],'YTick', 0:10:60,'XLim',[0.5 2.5],...
          'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
            'TickLength',2*get(gca,'TickLength'))

         %% SUBPLOT 3
        sp = sp +1;
        sp_ref = subplot(1, 3, sp);
        ax=gca;
        
        ax.OuterPosition(1) = tile_left*3 + tile_width*2;
        ax.OuterPosition(3) = tile_width;
                ax.Position(3) = 0.26;

        % reference
        middle = local_middlebroxplot(ax, 1-middle_off*dx,data_middle.qe_own, data_middle.Marker, kv.MarkerSize, data_middle.LineColor, data_middle.LineColor);
        local_middlebroxplot(ax, 2-middle_off*dx,data_middle.qe_other, data_middle.Marker, kv.MarkerSize, data_middle.LineColor, data_middle.LineColor);

        % baseline
%         baseline = local_middlebroxplot(ax, 1-majdak_off*dx,data_majdak.qe, data_majdak.Marker, kv.MarkerSize, data_majdak.LineColor, data_majdak.LineColor);
        
        % simulations
        cdist = cdist_init;
        for c=calib_plot_order
            baru(1,c) = local_middlebroxplot(ax, 1+cdist*dx,qe_own(c,1), Marker, kv.MarkerSize, LineColor(c,:),'w');
            local_middlebroxplot(ax, 2+cdist*dx,qe_other(c,1), Marker, kv.MarkerSize, LineColor(c,:),'w');
            cdist = cdist + 1;
        end
        
        % reijniers2014
        reij = local_middlebroxplot(ax, 1-reij_off*dx,data_reij.qe_own, data_reij.Marker,kv.MarkerSize, data_reij.LineColor,'w');
        local_middlebroxplot(ax, 2-reij_off*dx,data_reij.qe_other, data_reij.Marker,kv.MarkerSize, data_reij.LineColor,'w');

        % baumgartner2014
        baum = local_middlebroxplot(ax, 1-baum_off*dx,data_baum.qe_own, data_baum.Marker,kv.MarkerSize, data_baum.LineColor,'w');
        local_middlebroxplot(ax, 2-baum_off*dx,data_baum.qe_other, data_baum.Marker,kv.MarkerSize, data_baum.LineColor,'w');

        ylabel('Quadrant error [%]','FontSize',kv.FontSize)
        set(gca,'YLim',[-5 55],'YTick', 0:10:50,'XLim',[0.5 2.5],...
          'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
            'TickLength',2*get(gca,'TickLength'))
        
        for c=calib_plot_order
            switch lower(calibrations.combination{c,1})
                case 'dtf'
                    labels{1,c} = 'MP variant';
                case 'pge'
                    labels{1,c} = 'GP variant';
            end
        end

        leg = legend(sp_ref, [middle, baru(calib_plot_order), baum, reij], horzcat({'Actual data'},labels, {'baumgartner2014'}, {'reijniers2014'}));
        leg.FontSize = kv.FontSize - 2;
        leg.Units = 'centimeters';
        leg.Interpreter = 'latex';
        leg.Units = 'normalized';
        leg.Position = [0.0854369511609941 0.711535203689433 0.211155508622859 0.187531802247802];
        
%         saveas(mFig, 'new_plots/middlebrooks.eps', 'epsc')
%         print(mFig, 'fig8', '-dpng', '-r600')

    end
end

%% ------ exp_macpherson -------------------------------------
if flags.do_fig8 || flags.do_exp_macpherson2003
    exp_macpherson = [];

    if ~flags.do_redo
        exp_macpherson = amt_cache('get', ...
            'exp_macpherson2003',flags.cachemode);
    end
    
    calibrations = amt_load('barumerli2023', 'barumerli2023_calibration.mat');
    calibrations = calibrations.cache.value;
    
    if size(calibrations.sigma, 1) ~= length(calibrations.name)
        error('sigma values not enough for the provided subejcts')
    end
    
    % remove features spaces without monaural features
    monaural_none_idx = find(strcmp(calibrations.combination, 'monaural_none'));
    if monaural_none_idx > 0
        calibrations.combination(monaural_none_idx) = [];
        calibrations.sigma(:,monaural_none_idx) = [];
    end 
    
    % Settings
    num_exp = 50;
    num_sbj = length(calibrations.name);
    num_calib = size(calibrations.combination, 1);
    
    if flags.do_redo_fast
        exp_macpherson = [];
        num_exp = 2;
    end
    
    if flags.do_test
        exp_macpherson = [];
        num_exp = 1;
        num_calib = 1;
        num_sbj = 1; 
    end
    
    if isempty(exp_macpherson)
        sofa = amt_load('barumerli2023',['ARI_' upper(calibrations.name{1}) '_hrtf_M_dtf 256.sofa']);
        % generate stimulus
        % copyed from exp_baumgartner2014/do_fig10
        density = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8]; % ripples/oct
        depth =   10:10:40;        % ripple depth (peak-to-trough) in dB
        
        % 250-ms bursts, 20-ms raised-cosine fade in/out, flat from 0.6-16kHz
        fs = sofa.Data.SamplingRate;
        flow = 1e3;   % lower corner frequency of ripple modification in Hz
        fhigh = 16e3; % upper corner frequency of ripple modification in Hz
        Nf = 2^10;    % # Frequency bins

        f = 0:fs/2/Nf:fs/2;	% frequency bins
        id600 = find(f<=600,1,'last'); % index of 600 Hz (lower corner frequency of stimulus energy)
        idlow = find(f<=flow,1,'last'); % index of flow (ripples)
        idhigh = find(f>=fhigh,1,'first');  % index of fhigh (ripples)
        N600low = idlow - id600 +1;   % # bins without ripple modification
        Nlowhigh = idhigh - idlow +1; % # bins with ripple modification     % 
        O = log2(f(idlow:idhigh)/1e3);   % freq. trafo. to achieve equal ripple density in log. freq. scale

        % Raised-cosine '(i.e., cos^2)' ramp 1/8 octave wide
        fup = f(idlow)*2^(1/8);       % upper corner frequency of ramp upwards 
        idup = find(f<=fup,1,'last');
        Nup = idup-idlow+1;
        rampup = cos(-pi/2:pi/2/(Nup-1):0).^2;
        fdown = f(idhigh)*2^(-1/8);  % lower corner frequency of ramp downwards
        iddown = find(f>=fdown,1,'first');
        Ndown = idhigh-iddown+1;
        rampdown = cos(0:pi/2/(Ndown-1):pi/2).^2;
        ramp = [rampup ones(1,Nlowhigh-Nup-Ndown) rampdown];
        ramp = [-inf*ones(1,id600-1) zeros(1,N600low) ramp -inf*ones(1,Nf - idhigh)];

        % Ripples of Experiment I
        Sexp1 = zeros(Nf+1,length(density),2);  % 3rd dim: 1:0-phase 2:pi-phase
        Sexp1(idlow:idhigh,:,1) = (40/2* sin(2*pi*density'*O+ 0))';  % depth: 40dB, 0-phase
        Sexp1(idlow:idhigh,:,2) = (40/2* sin(2*pi*density'*O+pi))';  % depth: 40dB, pi-phase
        Sexp1 = repmat(ramp',[1,length(density),2]) .* Sexp1;
        Sexp1 = [Sexp1;Sexp1(Nf:-1:2,:,:)];
        Sexp1(isnan(Sexp1)) = -100;
        sexp1 = ifftreal(10.^(Sexp1/20),2*Nf);
        sexp1 = circshift(sexp1,Nf);  % IR corresponding to ripple modification
        sexp1 = squeeze(sexp1(:,:,1));
        % Ripples of Experiment II
        Sexp2 = zeros(Nf+1,length(depth),2);  % 3rd dim: 1:0-phase 2:pi-phase
        Sexp2(idlow:idhigh,:,1) = (depth(:)/2*sin(2*pi*1*O+ 0))';  % density: 1 ripple/oct, 0-phase
        Sexp2(idlow:idhigh,:,2) = (depth(:)/2*sin(2*pi*1*O+pi))';  % density: 1 ripple/oct, pi-phase
        Sexp2 = repmat(ramp',[1,length(depth),2]) .* Sexp2;
        Sexp2 = [Sexp2;Sexp2(Nf-1:-1:2,:,:)];
        Sexp2(isnan(Sexp2)) = -100;
        sexp2 = ifftreal(10.^(Sexp2/20),2*Nf);
        sexp2 = circshift(sexp2,Nf);  % IR corresponding to ripple modification
        sexp2 = squeeze(sexp2(:,:,1));

        if flags.do_test
             density(2:end) = []; % ripples/oct
            depth(2:end) = []; 
        end
        % preprocess templates for each user
        for i = 1:num_sbj
            amt_disp(['Processing subject ', num2str(i)]);

            for c = 1:num_calib
                amt_disp(['Pre-processing calibration #' num2str(c)]);
    
                sofa = amt_load('barumerli2023',['ARI_' upper(calibrations.name{i}) '_hrtf_M_dtf 256.sofa']);
                % extract directions
                % filter targets' coordinates
                % convert from spherical to horizontal-polar coordinates
                horpolar_coords = barumerli2023_coordinates(sofa).return_positions('horizontal-polar');
  
                % polar in [60, 120]
                % lateral = 0
                idx = find(((horpolar_coords(:, 2) >= -60 ...
                                & horpolar_coords(:, 2) <= 60) ...
                                | (horpolar_coords(:, 2) >= 120 & ...
                                        horpolar_coords(:, 2) <= 240)) ...
                                & (horpolar_coords(:, 1) <= 30 & horpolar_coords(:, 1) >= -30));

                amt_disp(['Pre-processing subject #' num2str(i)]);
                [template(c, i), target_flat(c, i)] = ...
                    barumerli2023_featureextraction(sofa, ...
                        'targ_az', sofa.SourcePosition(idx, 1), ...
                        'targ_el', sofa.SourcePosition(idx, 2), ...
                        calibrations.combination{c,1});

                amt_disp('Densities conditions');
                for j = 1:length(density)
                    target_exp1(c, i, j) = ...
                        barumerli2023_featureextraction(sofa, 'target', 'source', ...
                        'source_ir', squeeze(sexp1(:, j)), 'source_fs', fs, ...
                        'targ_az', sofa.SourcePosition(idx, 1), ...
                        'targ_el', sofa.SourcePosition(idx, 2), ...
                                            calibrations.combination{c,1});
                end

                amt_disp('Depth conditions');
                for j = 1:length(depth)
                    target_exp2(c, i, j) = ...
                        barumerli2023_featureextraction(sofa, 'target', 'source', ...
                        'source_ir', squeeze(sexp2(:, j)), 'source_fs', fs, ...
                        'targ_az', sofa.SourcePosition(idx, 1), ...
                        'targ_el', sofa.SourcePosition(idx, 2), ...
                        calibrations.combination{c,1});
                end
            end
        end
        
        % preallocation for results
        amt_disp('Allocating memory for results');
        estimations = struct('m',[]);
        est_expflat = repmat(estimations, num_calib, num_sbj); 
        est_exp1 = repmat(estimations, num_calib, ...
            num_sbj,length(density)); 
        est_exp2 = repmat(estimations, num_calib, ...
            num_sbj,length(depth)); 

        % data for prior computation
        data_majdak = data_majdak2010('Learn_M');
        data_majdak([1:5]) = [];
        
        % simulations
        for i = 1:num_sbj
            for c = 1:num_calib
                amt_disp(sprintf('\tCalibration #%i', c));

                assert(length(calibrations.sigma(i,c).values) == 5, 'something is wrong with the calibration file')
                sigma_l = calibrations.sigma(i,c).values(1);
                sigma_l2 = calibrations.sigma(i,c).values(2);
                sigma_mon = calibrations.sigma(i,c).values(3);
                sigma_m = calibrations.sigma(i,c).values(4);
                sigma_prior = calibrations.sigma(i,c).values(5);
                
                amt_disp(sprintf('\tSubject #%i', i));
                % flat spectrum estimations
                est_expflat(c, i, 1).m = barumerli2023('template', template(c, i), 'target', target_flat(c, i), ...
                                    'num_exp', num_exp, ...
                                    'sigma_itd', sigma_l, ...
                                    'sigma_ild', sigma_l2, ...
                                    'sigma_spectral', sigma_mon,...
                                    'sigma_motor', sigma_m, ...
                                    'sigma_prior', sigma_prior);

                % rippled estimations
                for j = 1:length(density)
                    est_exp1(c, i, j).m = barumerli2023('template', template(c, i), ...
                                    'target', target_exp1(c, i, j), ...
                                    'num_exp', num_exp, ...
                                    'sigma_itd', sigma_l, ...
                                    'sigma_ild', sigma_l2, ...
                                    'sigma_spectral', sigma_mon,...
                                    'sigma_motor', sigma_m, ...
                                    'sigma_prior', sigma_prior);
                end

                for j =1:length(depth)
                    est_exp2(c, i, j).m = barumerli2023('template', template(c, i), ...
                                    'target', target_exp2(c, i, j), ...
                                    'num_exp', num_exp, ...
                                    'sigma_itd', sigma_l, ...
                                    'sigma_ild', sigma_l2, ...
                                    'sigma_spectral', sigma_mon,...
                                    'sigma_motor', sigma_m, ...
                                    'sigma_prior', sigma_prior);
                end
            end        
        end

        % metrics
        % allocate memory for results
        % aggregate over different lateral angles
        pe_exp1 = zeros(num_calib, num_sbj, length(density));
        pe_exp2 = zeros(num_calib, num_sbj, length(depth));
        pe_flat = zeros(num_calib, num_sbj, 1);

        for c = 1:num_calib
            for i = 1:num_sbj
                % compute iterative regression (see Macpherson paper and localizationerror.m)   
                [f,r] = localizationerror(est_expflat(c,i).m, 'sirpMacpherson2000');
                
                pe_flat(c,i) = localizationerror(est_expflat(c,i).m, f, r, 'perMacpherson2003');

                for j = 1:length(density)
                    pe_exp1(c, i, j) = localizationerror(est_exp1(c, i, j).m, f, r, 'perMacpherson2003');
                end

                for j = 1:length(depth)
                    pe_exp2(c, i, j) = localizationerror(est_exp2(c ,i, j).m, f, r, 'perMacpherson2003');
                end
            end
        end
        
        % save cache
        exp_macpherson.pe_flat = pe_flat;
        exp_macpherson.pe_exp1 = pe_exp1;
        exp_macpherson.pe_exp2 = pe_exp2;

        if ~flags.do_redo_fast && ~flags.do_test
            amt_cache('set','exp_macpherson2003', exp_macpherson);
        end
    end
    
    % Original data:
    data = data_macpherson2003;

    % Reijniers2014's data
    data_reij = exp_reijniers2014('fig4_barumerli2020forum', 'no_plot');

    % Baumgartner2014's data
    % varargout{1} = {pe_exp1,pe_exp2,pe_flat,noDCN};
    data_baum_temp = exp_baumgartner2014('fig10', 'no_plot');
    data_baum.pe_exp1 = data_baum_temp{1,1};
    data_baum.pe_exp2 = data_baum_temp{1,2};
    data_baum.pe_flat = data_baum_temp{1,3};

    % Phase condition handling
    % average across the phase condition
    % real data
    data.pe_exp1 = mean(data.pe_exp1,3);
    data.pe_exp2 = mean(data.pe_exp2,3);
    % baumgartner data 
    data_baum.pe_exp1 = mean(data_baum.pe_exp1,3);
    data_baum.pe_exp2 = mean(data_baum.pe_exp2,3);
    idphase = 1;

    % Increase
    % reijniers2014
    data_reij.pe_exp1 = data_reij.pe_exp1 - repmat(data_reij.pe_flat(:), 1, size(data_reij.pe_exp1, 2));
    data_reij.pe_exp2 = data_reij.pe_exp2 - repmat(data_reij.pe_flat(:), 1, size(data_reij.pe_exp2, 2));
    % baumgartner data
    data_baum.pe_exp1 = data_baum.pe_exp1 - repmat(data_baum.pe_flat(:),1,size(data_baum.pe_exp1,2));
    data_baum.pe_exp2 = data_baum.pe_exp2 - repmat(data_baum.pe_flat(:),1,size(data_baum.pe_exp2,2));

    % Statistics
    % real data
    data.quart_pe_flat = quantile(data.pe_flat,[.25 .50 .75]);
    data.quart_pe_exp1 = quantile(data.pe_exp1,[.25 .50 .75]);
    data.quart_pe_exp2 = quantile(data.pe_exp2,[.25 .50 .75]);

    % reijniers data
    data_reij.quart_pe_flat = quantile(data_reij.pe_flat,[.25 .50 .75]);
    data_reij.quart_pe_exp1 = quantile(data_reij.pe_exp1,[.25 .50 .75]);
    data_reij.quart_pe_exp2 = quantile(data_reij.pe_exp2,[.25 .50 .75]);
    
    % baumgartner data
    data_baum.quart_pe_flat = quantile(data_baum.pe_flat,[.25 .50 .75]);
    data_baum.quart_pe_exp1 = quantile(data_baum.pe_exp1,[.25 .50 .75]);
    data_baum.quart_pe_exp2 = quantile(data_baum.pe_exp2,[.25 .50 .75]);
    
    for c = 1:num_calib
        % simulations data
        pe_flat = exp_macpherson.pe_flat(c,:);
        pe_exp1 = squeeze(exp_macpherson.pe_exp1(c,:,:));
        pe_exp2 = squeeze(exp_macpherson.pe_exp2(c,:,:));

        % simulations
        pe_exp1 = pe_exp1 - repmat(pe_flat(:), 1, size(pe_exp1, 2) );
        pe_exp2 = pe_exp2 - repmat(pe_flat(:), 1, size(pe_exp2, 2) );
        
        % simulations
        quart_pe_flat(c,:) = quantile(pe_flat,[.25 .50 .75]);
        quart_pe_exp1(c,:,:) = quantile(pe_exp1,[.25 .50 .75]);
        quart_pe_exp2(c,:,:) = quantile(pe_exp2,[.25 .50 .75]);
    end


    % plot
    if flags.do_plot
        calib_plot_order = [1,2]; % lat, dtf, pge
        
        dx = 1.05;
        FontSize = kv.FontSize;
        MarkerSize = kv.MarkerSize;

        LineColor = [[0.9290, 0.6940, 0.1250]; ...
                    [0.4940, 0.1840, 0.5560]; ...
                    [0.4660, 0.6740, 0.1880]; ...
                    [0.3010, 0.7450, 0.9330]; ...
                    [0.6350, 0.0780, 0.1840]];
        data.Marker = 'ko-';
        data.LineColor = [1 1 1]*0;

        data_reij.Marker = 'v-';
        data_reij.LineColor = [0 0.4470 0.7410];
        
        data_baum.Marker = 'd-';
        data_baum.LineColor = [0.8500 0.3250 0.0980];

        % Exp1
        mFig = figure;
        mFig.Units = 'points';
        mFig.Position = [0 0 500 265];
%         sgtitle(sprintf('%s ', calibrations.combination{c,:}), 'Interpreter', 'none')

        subplot(2,8,1:8)
        mach = errorbar(data.density,data.quart_pe_exp1(2,:,idphase),...
        data.quart_pe_exp1(2,:,idphase) - data.quart_pe_exp1(1,:,idphase),...
        data.quart_pe_exp1(3,:,idphase) - data.quart_pe_exp1(2,:,idphase),...
        'o-','MarkerSize',MarkerSize, 'Color', data.LineColor, ...
        'MarkerFaceColor', data.LineColor);
        hold on
        for c = calib_plot_order
            baru(1,c) = errorbar(data.density/dx,squeeze(quart_pe_exp1(c, 2,:)),...
            squeeze(quart_pe_exp1(c,2,:) - quart_pe_exp1(c,1,:)),...
            squeeze(quart_pe_exp1(c,3,:) - quart_pe_exp1(c,2,:)),...
            's-','MarkerSize',MarkerSize, 'Color', LineColor(c,:),'MarkerFaceColor','w');
        end
        hold on
        reij = errorbar(data.density*dx,data_reij.quart_pe_exp1(2,:),...
        data_reij.quart_pe_exp1(2,:) - data_reij.quart_pe_exp1(1,:),...
        data_reij.quart_pe_exp1(3,:) - data_reij.quart_pe_exp1(2,:),...
        'v--','MarkerSize',MarkerSize, 'Color', data_reij.LineColor,'MarkerFaceColor','w');

        baum = errorbar(data.density*dx,data_baum.quart_pe_exp1(2,:,idphase),...
        data_baum.quart_pe_exp1(2,:,idphase) - data_baum.quart_pe_exp1(1,:,idphase),...
        data_baum.quart_pe_exp1(3,:,idphase) - data_baum.quart_pe_exp1(2,:,idphase),...
        'd--','MarkerSize',MarkerSize, 'Color', data_baum.LineColor,'MarkerFaceColor','w');

        set(gca,'XScale','log','YMinorTick','on')
        set(gca,'XLim',[0.25/1.2 8*1.2],'XTick',data.density,'YLim',[-16 70],'FontSize',FontSize)
        xlabel('Ripple Density [ripples/octave]','FontSize',FontSize)
        ylabel({'Increase in';'Polar Error Rate [%]'},'FontSize',FontSize)

        % Exp2
        sp_ref = subplot(2,8,9:13);
        errorbar(data.depth,data.quart_pe_exp2(2,:,idphase),...
        data.quart_pe_exp2(2,:,idphase) - data.quart_pe_exp2(1,:,idphase),...
        data.quart_pe_exp2(3,:,idphase) - data.quart_pe_exp2(2,:,idphase),...
        'o-','MarkerSize',MarkerSize, 'Color', data.LineColor, ...
        'MarkerFaceColor', data.LineColor);
        hold on
        for c = calib_plot_order
            errorbar(data.depth-0.5,squeeze(quart_pe_exp2(c, 2,:)),...
            squeeze(quart_pe_exp2(c,2,:) - quart_pe_exp2(c,1,:)),...
            squeeze(quart_pe_exp2(c,3,:) - quart_pe_exp2(c,2,:)),...
            's-','MarkerSize',MarkerSize, 'Color', LineColor(c,:),'MarkerFaceColor','w');
        end
        hold on
        errorbar(data.depth+1,data_reij.quart_pe_exp2(2,:),...
        data_reij.quart_pe_exp2(2,:) - data_reij.quart_pe_exp2(1,:),...
        data_reij.quart_pe_exp2(3,:) - data_reij.quart_pe_exp2(2,:),...
        'v--','MarkerSize',MarkerSize, 'Color', data_reij.LineColor,'MarkerFaceColor','w');

        errorbar(data.depth+1,data_baum.quart_pe_exp2(2,:,idphase),...
        data_baum.quart_pe_exp2(2,:,idphase) - data_baum.quart_pe_exp2(1,:,idphase),...
        data_baum.quart_pe_exp2(3,:,idphase) - data_baum.quart_pe_exp2(2,:,idphase),...
        'd--','MarkerSize',MarkerSize, 'Color', data_baum.LineColor,'MarkerFaceColor','w');

        set(gca,'XLim',[data.depth(1)-5 data.depth(end)+5],'XTick',data.depth,...
        'YLim',[-16 70],'YMinorTick','on','FontSize',FontSize)
        xlabel('Ripple Depth [dB]','FontSize',FontSize)
        ylabel({'Increase in';'Polar Error Rate [%]'},'FontSize',FontSize)
        ytick = get(gca,'YTick');
        ticklength = get(gca,'TickLength');

        % Baseline
        subplot(2,8,14:16)
        errorbar(-1,data.quart_pe_flat(2),...
        data.quart_pe_flat(2) - data.quart_pe_flat(1),...
        data.quart_pe_flat(3) - data.quart_pe_flat(2),...
        'o-','MarkerSize',MarkerSize, 'Color', data.LineColor, ...
        'MarkerFaceColor', data.LineColor);
        hold on
        for c = calib_plot_order
            errorbar(-0.5 + (c-1)*0.5,quart_pe_flat(c,2),...
                quart_pe_flat(c,2) - quart_pe_flat(c,1),...
                quart_pe_flat(c,3) - quart_pe_flat(c,2),...
                's-','MarkerSize',MarkerSize, 'Color', LineColor(c,:),'MarkerFaceColor','w');
        end
        hold on
        errorbar(0.5,data_baum.quart_pe_flat(2),...
        data_baum.quart_pe_flat(2) - data_baum.quart_pe_flat(1),...
        data_baum.quart_pe_flat(3) - data_baum.quart_pe_flat(2),...
        'd--','MarkerSize',MarkerSize, 'Color', data_baum.LineColor,'MarkerFaceColor','w');

        errorbar(1,data_reij.quart_pe_flat(2),...
        data_reij.quart_pe_flat(2) - data_reij.quart_pe_flat(1),...
        data_reij.quart_pe_flat(3) - data_reij.quart_pe_flat(2),...
        'd--','MarkerSize',MarkerSize, 'Color', data_reij.LineColor,'MarkerFaceColor','w');

        set(gca,'XLim',[-3 3],'XTick',0,'XTickLabel',{'Baseline'},...
        'YLim',[-15 59],'YTick',ytick,'TickLength',3*ticklength,...
        'FontSize',FontSize,'YAxisLocation','right')
        xlabel(' ','FontSize',FontSize)
        ylabel({'Polar Error Rate [%]'},'FontSize',FontSize)

        %legend
        for c = calib_plot_order
            switch lower(calibrations.combination{c,1})
                case 'dtf'
                    labels{1,c} = 'MP variant';
                case 'pge'
                    labels{1,c} = 'GP variant';
            end
        end
        
        leg = legend(sp_ref, [mach, baru(calib_plot_order), baum, reij], horzcat({'Actual data'}, labels(calib_plot_order), {'baumgartner2014'}, {'reijniers2014'}));
        leg.FontSize = FontSize - 2;
        leg.Units = 'normalized';
        leg.Interpreter = 'latex';
        leg.Orientation = 'horizontal';
        leg.Position = [0.147112074200069 0.944044062817491 0.739274312186318 0.0430594892744974];
        % Overall correlation between actual and predicted median values
        for c=calib_plot_order
            m_pe_pred = [squeeze(quart_pe_exp1(c,2,:))' squeeze(quart_pe_exp2(c,2,:))'];
            m_pe_actual = [data.quart_pe_exp1(2,:) data.quart_pe_exp2(2,:)];
            r = corrcoef(m_pe_pred,m_pe_actual);
            r_sqr = r(2);

            amt_disp('Correlation between actual and predicted median values (15 conditions):')
            amt_disp(sprintf('%s: r = %0.2f', labels{1,c}, r_sqr))
        end
        
%          saveas(mFig, 'new_plots/macpherson.eps', 'epsc')
%          print(mFig, 'fig9', '-dpng', '-r600')

    end
end

function hg = local_middlebroxplot(ax, x, data, Marker, MarkerSize, LineColor, FaceColor)
    lilen = 0.05; % length of horizontal lines
    
    hb=[];
    % Symbols
    i=1; hb(i) = plot(ax, x, data.quantiles(1),'x','MarkerSize',MarkerSize, 'MarkerEdgeColor', LineColor, 'MarkerFaceColor', LineColor); % min
    hold on
    i=i+1; hb(i) = plot(ax, x,data.quantiles(7),'x','MarkerSize',MarkerSize, 'MarkerEdgeColor', LineColor, 'MarkerFaceColor', LineColor); % max

    % Horizontal lines
    i=i+1; hb(i:(i+1)) = line(ax, x+0.5*[-lilen,lilen],repmat(data.quantiles(2),2),'Color',LineColor); % lower whisker
    i=i+2; hb(i:(i+1)) = line(ax, x+[-lilen,lilen],repmat(data.quantiles(3),2),'Color',LineColor); % 25% Quartile
    i=i+2; hb(i:(i+1)) = line(ax, x+[-lilen,lilen],repmat(data.quantiles(4),2),'Color',LineColor); % Median
    i=i+2; hb(i:(i+1)) = line(ax, x+[-lilen,lilen],repmat(data.quantiles(5),2),'Color',LineColor); % 75% Quartile
    i=i+2; hb(i:(i+1)) = line(ax, x+0.5*[-lilen,lilen],repmat(data.quantiles(6),2),'Color',LineColor); % upper whisker

    % Vertical lines
    i=i+2; hb(i:(i+1)) = line(ax, [x,x],data.quantiles(2:3),'Color',LineColor); % connector lower whisker
    i=i+2; hb(i:(i+1)) = line(ax, [x,x],data.quantiles(5:6),'Color',LineColor); % connector upper whisker
    i=i+2; hb(i:(i+1)) = line(ax, [x,x]-lilen,data.quantiles([3,5]),'Color',LineColor); % left box edge
    i=i+2; hb(i:(i+1)) = line(ax, [x,x]+lilen,data.quantiles([3,5]),'Color',LineColor); % left box edge
    
    % middle value
    i=i+1; hb(i) = plot(x,data.mean, Marker,'MarkerSize', MarkerSize, 'MarkerFaceColor', FaceColor, 'MarkerEdgeColor',LineColor);
    
    % create a group to avoid issues with the legend
    % https://stackoverflow.com/questions/12894652/matlab-how-to-make-a-custom-legend
    hg = hggroup;
    set(hb,'Parent',hg);
    set(get(get(hg,'Annotation'),'LegendInformation'),...
      'IconDisplayStyle','off');
  
  
function [x,y,x_sign] = lambert_area_projection(az, el)
    % az and el has to be in radiants!!!!
    % lambert equal area projection 
    x_sign = 1; 
    if abs(az) > pi/2; az = az - pi; x_sign = -1; end
    k = sqrt(2 ./ (eps + 1  + (cos(el) .* cos(az))));
    x = x_sign * k * 1 .* cos(el) .* sin(az) ./ sqrt(2); % ./sqrt(2) normalizing
    y = k * 1 .* sin(el) ./ sqrt(2);
    
function rau=rau(X,N,opt)

% RAU   rationalized arcsine transform
% RAU(X,N) transforms the number of correct responses X to the
% rationalized arcsine (rau). N gives the number of repetitions.
% 
% This function allows to use ANOVA statistics with percent correct scores
% because: 1) RAUs are normally distributed; 2) mean and variance of RAUs
% are not correlated with eachother; and 3) likelihood that a score will
% increase/decrease will remain constant over the range.
% 
% RAU=RAU(X,N,opt) defines one of the following options:
%  'Pc'  ... X is given in percent correct scores (0..100%)
%  'X'   ... X is given in the number of correct responses (default)
%
% The formula are based on Sherbecoe and Studebaker,
% Int. J. of Audiology 2004; 43; 442-448
%
% See also IRAU.

% 30.8.2007, Piotr Majdak
%

if exist('opt','var')
  if strcmp(upper(opt),'PC')
    X=X/100*N;
  elseif strcmp(upper(opt),'X')
  else
    error('OPT must be Pc (=percent correct) or X (=number of correct responses)');
  end
end
th=asin(sqrt(X/(N+1)))+asin(sqrt((X+1)/(N+1)));
rau=146/pi*(th)-23;

function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering

% Pekka Kumpulainen (2021). tight_subplot(Nh, Nw, gap, marg_h, marg_w) 
% (https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w),
% MATLAB Central File Exchange. Retrieved June 18, 2021. 

% This function is covered by Copyright (c) 2016, Pekka Kumpulainen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1; 
    gap = [gap gap];
end
if numel(marg_w)==1; 
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1; 
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh; 

% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);

