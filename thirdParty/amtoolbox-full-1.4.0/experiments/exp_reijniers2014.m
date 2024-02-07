function varargout = exp_reijniers2014(varargin)
%EXP_REIJNIERS2014 Experiments of Reijniers et al. (2014)
%   Usage: [] = exp_reijniers2014(flag) 
%
%   EXP_REIJNIERS2014(flag) reproduces figures of the study from 
%   Reijniers et al. (2014). Note that HRTFs from the ARI database
%   are used for Fig. 4-6, which is most probably different to the HRTFs 
%   used for the calculations in the paper. For Fig. 2 and 3 HRTFs from the
%   Symare database are used, which are also most probably non original.
%
%
%   The following flags can be specified
%   
%     'fig2a'          Reproduce Fig. 2 (a):
%                      The results for simulated trial of a locali-
%                      zation experiment for subject 1 [see also Fig. 3(a)] 
%                      in front Theta = (-34,45)deg (a)  is shown. The
%                      figures show the simulated posterior angular
%                      probabilities and the templates corresponding to the
%                      original direction Theta (solid black line) and
%                      estimated direction Theta_est (dotted black line).
%                      The simulated auditory input is shown as the solid 
%                      red line.
%
%     'fig2b'          Reproduce Fig. 2 (b):
%                      The results for simulated trial of a locali-
%                      zation experiment for subject 1 [see also Fig. 3(a)] 
%                      in back Theta = (14,171)deg (b) is shown. The figures 
%                      show the simulated posterior angular probabilities
%                      and the templates corresponding to the original 
%                      direction Theta (solid black line) and estimated
%                      direction Theta_est (dotted black line). The
%                      simulated auditory input is shown as the solid red
%                      line.
%
%     'fig3'           Reproduce Fig. 3:
%                      The simulated mean spherical error as function
%                      of the source position for 3 different subjects. The
%                      average was taken over 500 localization trials for  
%                      each source position.
%
%     'fig4'           Reproduce Fig. 4(a):                    
%                      The mean localization performance when the
%                      simulation results are averaged over 100 subjects.
%                      Superimposed arrows indicate the size and direction 
%                      of local population response biases for different
%                      source positions.  
%
%     'fig5'           Reproduce Fig. 5:
%                      Sensitivity analysis of the Bayesian localization
%                      model. The simulated mean spherical error is shown 
%                      as function of the source position, when each of the
%                      model parameters is varied separately. The average 
%                      was taken over 500 localization trials for each 
%                      source position and 100 subjects were pooled.
%                      (a) The model with input parameters as described in
%                      the main text. The standard deviation is doubled,
%                      respectively for (b) the noise on the ITD, (c) the 
%                      internal noise and the variation on (d) the source
%                      spectrum. In (e), the bandwidth of the source is 
%                      reduced to [300Hz-8kHz].  
%
%     'fig6'           Reproduce Fig. 6:
%                      The mean spherical error for different values of the
%                      SNR. SNR=75dB corresponds to the control situation,
%                      see Fig.5(a), as the magnitude in all frequency 
%                      channels is above the system noise level.
%
%     'tab1_barumerli2020aes'           Reproduce Tab. 1:
%                      Comparison between actual (majdak2010) and simulated,
%                      performances by relying on different perceptual
%                      metrics (middlebrooks1999b). The variable 'multiplier'
%                      allows to tune the internal noise.
%
%     'fig2_barumerli2020forum'    Reproduce Fig.2 of Barumerli et al. (2020):
%                                  comparison between virtual estimations
%                                  and real data (see
%                                  data_middlebrooks1999()) for individual
%                                  and non-individual DTFs
%
%     'fig3_barumerli2020forum'    Reproduce Fig.3 of Barumerli et al. (2020):
%                                  estimations' evaluation with band limited spectra.
%                                  See Best et al. 2005 and Fig 11 baumgartner2014
%
%     'fig4_barumerli2020forum'    Reproduce Fig.4 of Barumerli et al. (2020):
%                                  estimations' evaluation with rippled
%                                  spectra. See Macpherson and Middlebrooks 2003
%                                  and Fig 10 baumgartner2014
%
%   Further, cache flags (see amt_cache) and plot flags can be specified
%   (Warning: Enforcing recalculation of cached data might take several
%   hours).
%     'no_plot'        Do not compute plots. Flag for cluster execution.
%
%     'interp'         Plot scattered interpolated data (default).
%   
%     'scatter'        Plot only discrete scattered data instead of inter-
%                      polated scattered data.
%
%     'redo_fast'      Quickly recalculate data for figures without using
%                      cache. To speed up computation the amount of locali-
%                      zation trials for Fig. 3 is reduced to 20 and the 
%                      amount of subjects for Fig. 4-6 is reduced to 12.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API v1.1 or higher from 
%      http://sourceforge.net/projects/sofacoustics for Matlab (e.g. in 
%      thirdparty/SOFA)
%
%   Examples:
%   ---------
%
%   To display right part of Fig. 2 (a) use :
%
%     exp_reijniers2014('fig2a');
%
%   To display Fig. 3 and do a quick recalcultaion use :
%
%     exp_reijniers2014('fig3','redo_fast');
%
%
%
%   See also: reijniers2014 plot_reijniers2014 reijniers2014_featureextraction
%   
%   References:
%     R. Barumerli, P. Majdak, R. Baumgartner, J. Reijniers, M. Geronazzo,
%     and F. Avanzini. Predicting directional sound-localization of human
%     listeners in both horizontal and vertical dimensions. In Audio
%     Engineering Society Convention 148. Audio Engineering Society, 2020.
%     
%     R. Barumerli, P. Majdak, R. Baumgartner, M. Geronazzo, and F. Avanzini.
%     Evaluation of a human sound localization model based on bayesian
%     inference. In Forum Acusticum, 2020.
%     
%     J. Reijniers, D. Vanderleist, C. Jin, C. S., and H. Peremans. An
%     ideal-observer model of human sound localization. Biological
%     Cybernetics, 108:169--181, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_reijniers2014.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   #Author: Jonas Reijniers 
%   #Author: Roberto Barumerli (2019): Modified and adapted for amtoolbox
%   #Author: Michael Sattler (2019): Modified and adapted for amtoolbox
%   #Author: Clara Hollomey (2021): Modified and adapted for AMT 1.0


%% ------ Check input options ---------------------------------------------
definput.import = {'amt_cache'};
definput.keyvals.MarkerSize = 6;
definput.keyvals.FontSize = 12;

definput.flags.type = {'missingflag', 'fig2a','fig2b', ...
                        'fig3','fig4','fig5','fig6', ...
                        'tab1_barumerli2020aes', ...
                        'fig2_barumerli2020forum', ...
                        'fig3_barumerli2020forum',...
                        'fig4_barumerli2020forum'};
definput.flags.plot = {'plot', 'no_plot'};
definput.flags.plot_type = {'interp','scatter'};
definput.flags.redo = {'no_redo_fast','redo_fast'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},...
             definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.', ...
      upper(mfilename),flagnames);
end


%% ------ fig2a -----------------------------------------------------------
if flags.do_fig2a  
    if flags.do_redo_fast || flags.do_redo
        num_exp = 1;
        fig2a = [];
    else
        num_exp = 1; % definedin case of flag 'redo' 
        fig2a = amt_cache('get','fig2a',flags.cachemode);
    end
               
    if isempty(fig2a)
        % Load SOFA file 
        SOFA_obj=amt_load('reijniers2014','HA01.sofa');

        % Get index of closest directions to Theta shown in Fig. 2(a) 
        az = -34;
        el = 45;
        % Preprocessing source information for both directions
        [template, target] = reijniers2014_featureextraction(SOFA_obj, ... 
                                        'targ_az', az, 'targ_el', el);
        
        [fig2a.doa, fig2a.params] = ...
            reijniers2014(template, target, 'num_exp', num_exp);
        amt_cache('set','fig2a',fig2a);
    end

    if flags.do_plot
        % plot left side of Fig. 2 (a)   
        plot_reijniers2014(fig2a.params.template_coords, ...
                           squeeze(fig2a.params.post_prob), ...
                           'target', fig2a.doa.real, ...
                           flags.plot_type,flags.type);
        title('Fig. 2(a), posterior probability density');  

        % plot right side of fig 2(a)
        figure('NumberTitle', 'off', 'Name', 'Fig. 2 (a) templates');
        set(gcf, 'Position',  [100,150, 1100, 550]);
        s1 = subplot(1,3,1);
        set(s1, 'position', [0.07, 0.11, 0.05, 0.8] );
        plot(0, fig2a.params.T_target(fig2a.params.Tidx.itd),'k.', ...
             0, fig2a.params.T_template(fig2a.params.Tidx.itd),'b.', ...
             0, fig2a.params.X(fig2a.params.Tidx.itd),'r.','MarkerSize',20);    
        set(gca,'Color',[0.99,0.95,0.93]);
        ylim([-10,10]);
        xlim([-0.2,0.2]);
        ax1 = gca;
        title('ITD','FontSize', 14);
        ylabel('itd (jnd)','FontSize', 13);
        set(ax1,'xticklabel',[]);

        s2 = subplot(1,3,2);
        set(s2, 'position', [0.19, 0.11, 0.35, 0.8] );
        plot(fig2a.params.freq_channels,...
             fig2a.params.T_template(fig2a.params.est_idx, fig2a.params.Tidx.Hm),'k', ...
             fig2a.params.freq_channels,fig2a.params.T_target(fig2a.params.Tidx.Hm),'k:', ...
             fig2a.params.freq_channels,squeeze(fig2a.params.X(1,1,fig2a.params.Tidx.Hm)),'r', ...
             fig2a.params.freq_channels,squeeze(fig2a.params.X(1,1,fig2a.params.Tidx.Hm)),'r.', ...
             'MarkerSize',20,'LineWidth',1.5);
        legend({'$T_{-}(\Theta)$','$T_{-}(\widehat{\Theta})$', ... 
                '$X_{-}(\Theta)$'},'Interpreter','latex','FontSize',13);
        title('spectral difference','FontSize', 14);
        ylabel('logmagnitude (dB)','FontSize',13);
        xlabel('log_{10}(f) (Hz)','Interpreter','tex','FontSize',13);
        ylim([-20,30]);
        xlim([min(fig2a.params.freq_channels) max(fig2a.params.freq_channels)]);
        set(gca,'Xscale','log')

        s3 = subplot(1,3,3); 
        set(s3, 'position', [0.61, 0.11, 0.35, 0.8] );
        plot(fig2a.params.freq_channels,...
             fig2a.params.T_template(fig2a.params.est_idx, fig2a.params.Tidx.Hp),'k', ...
             fig2a.params.freq_channels,fig2a.params.T_target(fig2a.params.Tidx.Hp),'k:', ...
             fig2a.params.freq_channels,squeeze(fig2a.params.X(1,1,fig2a.params.Tidx.Hp)),'r', ...
             fig2a.params.freq_channels,squeeze(fig2a.params.X(1,1,fig2a.params.Tidx.Hp)),'r.', ...
             'MarkerSize',20,'LineWidth',1.5);
        legend({'$T_{+}(\Theta)$','$T_{+}(\widehat{\Theta})$',  ...
                '$X_{+}(\Theta)$'},'Interpreter','latex','FontSize',13);
        title('spectral sum','FontSize', 14);
        ylabel('logmagnitude (dB)','FontSize',13);
        xlabel('log_{10}(f) (Hz)','Interpreter','tex','FontSize',13);
        ylim([-20,30]);
        xlim([min(fig2a.params.freq_channels) max(fig2a.params.freq_channels)]);
        set(gca,'Xscale','log');
    end
end

%% ------ fig2b -----------------------------------------------------------
if flags.do_fig2b  
    if flags.do_redo_fast || flags.do_redo
        num_exp = 1;
        fig2b = [];
    else
        num_exp = 1; 
        fig2b = amt_cache('get','fig2b',flags.cachemode);
    end
               
    if isempty(fig2b)
        % Load SOFA file 
        SOFA_obj=amt_load('reijniers2014','HA01.sofa');

        % Get index of closest directions to Theta shown in Fig. 2(a) 
        az = 14;
        el = 171;
       
        % Preprocessing source information for both directions
        [template, target] = reijniers2014_featureextraction(SOFA_obj, ... 
                            'targ_az', az, 'targ_el', el);
        
        % Model
        [fig2b.doa, fig2b.params] = ...
            reijniers2014(template, target, 'num_exp', num_exp);
        amt_cache('set','fig2b',fig2b);
    end
    
    if ~flags.do_no_plot
        % plot left side of Fig. 2b   
        plot_reijniers2014(fig2b.params.template_coords, ...
                           squeeze(fig2b.params.post_prob), ...
                           'target', fig2b.doa.real, ...
                           flags.plot_type,flags.type);
        title('Fig. 2(b), posterior probability density');    
        % plot right side of fig 2(b)
        figure('NumberTitle', 'off', 'Name', 'Fig. 2 (b) templates');
        set(gcf, 'Position',  [100,150, 1100, 550]);
        s1 = subplot(1,3,1);
        set(s1, 'position', [0.07, 0.11, 0.05, 0.8] );
        plot(0, fig2b.params.T_target(fig2b.params.Tidx.itd),'k.', ...
             0, fig2b.params.T_template(fig2b.params.Tidx.itd),'b.', ...
             0, fig2b.params.X(fig2b.params.Tidx.itd),'r.','MarkerSize',20); 
        set(gca,'Color',[0.99,0.95,0.93]);
        ylim([-10,10]);
        xlim([-0.2,0.2]);
        ax1 = gca;
        title('ITD','FontSize', 14);
        ylabel('itd (jnd)','FontSize', 13);
        set(ax1,'xticklabel',[]);

        s2 = subplot(1,3,2);
        set(s2, 'position', [0.19, 0.11, 0.35, 0.8] );
        plot(fig2b.params.freq_channels,...
             fig2b.params.T_template(fig2b.params.est_idx, fig2b.params.Tidx.Hm),'k', ...
             fig2b.params.freq_channels,fig2b.params.T_target(fig2b.params.Tidx.Hm),'k:', ...
             fig2b.params.freq_channels,squeeze(fig2b.params.X(1,1,fig2b.params.Tidx.Hm)),'r', ...
             fig2b.params.freq_channels,squeeze(fig2b.params.X(1,1,fig2b.params.Tidx.Hm)),'r.', ...
             'MarkerSize',20,'LineWidth',1.5);
        legend({'$T_{-}(\Theta)$','$T_{-}(\widehat{\Theta})$',...
                '$X_{-}(\Theta)$'},'Interpreter','latex','FontSize',13);
        title('spectral difference','FontSize', 14);
        ylabel('logmagnitude (dB)','FontSize',13);
        xlabel('log_{10}(f) (Hz)','Interpreter','tex','FontSize',13);
        ylim([-20,30]);
        xlim([min(fig2b.params.freq_channels) max(fig2b.params.freq_channels)]);
        set(gca,'Xscale','log');

        s3 = subplot(1,3,3); 
        set(s3, 'position', [0.61, 0.11, 0.35, 0.8] );
        plot(fig2b.params.freq_channels,...
             fig2b.params.T_template(fig2b.params.est_idx, fig2b.params.Tidx.Hp),'k', ...
             fig2b.params.freq_channels,fig2b.params.T_target(fig2b.params.Tidx.Hp),'k:', ...
             fig2b.params.freq_channels,squeeze(fig2b.params.X(1,1,fig2b.params.Tidx.Hp)),'r', ...
             fig2b.params.freq_channels,squeeze(fig2b.params.X(1,1,fig2b.params.Tidx.Hp)),'r.', ...
             'MarkerSize',20,'LineWidth',1.5);
        legend({'$T_{+}(\Theta)$','$T_{+}(\widehat{\Theta})$', ... 
                '$X_{+}(\Theta)$'},'Interpreter','latex','FontSize',13);
        title('spectral sum','FontSize', 14);
        ylabel('logmagnitude (dB)','FontSize',13);
        xlabel('log_{10}(f) (Hz)','Interpreter','tex','FontSize',13);
        ylim([-20,30]);
        xlim([min(fig2b.params.freq_channels) max(fig2b.params.freq_channels)]);
        set(gca,'Xscale','log');
    end
end

%% ------ fig3 ------------------------------------------------------------
%if flags.do_fig3
%    if flags.do_redo_fast
%        num_exp = 20;
%        fig3 = [];
%    else
%        num_exp = 500;
%        fig3 = amt_cache('get','fig3',flags.cachemode);
%    end
%        
%    if isempty(fig3)
%        filenames={'HA03.sofa'; 'HA01.sofa';'HA06.sofa'};
%        for i = 1:length(filenames)
%            path = ['db://reijniers2014/',filenames{i}];
%            SOFA_obj(i)=amt_load('reijniers2014',filenames{i});
%        end
%        % Preallocation
%        fig3 = struct('metrics', struct([]), 'err',struct([]), 'doa_real',struct([]));
%        fig3 = repmat(fig3,length(filenames),1);
%        
%        amt_disp('Processing three subjects in parallel...');
%        
%        for i = 1:length(filenames)
%            amt_disp(['Processing subject #' num2str(i)]);
%            % Get directions from SOFA files
%            [template, target] = reijniers2014_featureextraction(SOFA_obj(i));
%            doa = reijniers2014(template, target,'num_exp',num_exp);
%            
%            fig3(i).err = reijniers2014_metrics(doa);
%            fig3(i).doa_real = doa.real;
%            fig3(i).metrics = reijniers2014_metrics(doa, 'middle_metrics'); 
%        end
%        
%        if ~flags.do_redo_fast
%             amt_cache('set','fig3',fig3);
%        end
%    end
    
%    if flags.do_plot
%        % plot for each subject
%        for i = 1:length(fig3)
%            [~, cbar] = plot_reijniers2014(fig3(i).doa_real, ...
%                                    fig3(i).err, ...
%                                    flags.plot_type, flags.type);
%            title(sprintf('Fig. 3(a), Subject %i', i));
%            caxis([0,20]);
%            cbar.Label.String = 'Mean spherical error [^\circ]';
%            cbar.Label.FontSize = 12;    
%            ct1=get(cbar,'TickLabels');
%            for k=1:numel(ct1)
%                ct1{k}=sprintf('%s',ct1{k});
%            end
%            set(cbar,'xticklabel',ct1);
%        end
%    end
%    
%    for i = 1:length(fig3)
%        amt_disp(sprintf('\nSUBJECT %i', i),'documentation');
%        amt_disp(sprintf('\t\t\tSIM'),'documentation');
%        amt_disp(sprintf('lateral_bias [deg]:\t%0.2f', ...
%            mean([fig3(i).metrics.accL])),'documentation');
%        amt_disp(sprintf('lateral_rms_error [deg]:\t%0.2f', ...
%            mean([fig3(i).metrics.rmsL])),'documentation');
%        amt_disp(sprintf('elevation_bias [deg]:\t%0.2f', ...
%            mean([fig3(i).metrics.accE])),'documentation');
%        amt_disp(sprintf('local_rms_polar [deg]:\t%0.2f', ...
%            mean([fig3(i).metrics.rmsP])),'documentation');
%        amt_disp(sprintf('quadrant_err [%%]:\t%0.2f', ...
%            mean([fig3(i).metrics.querr])),'documentation');
%    end
%end

%% ------ fig4 ------------------------------------------------------------
if flags.do_fig4
    if flags.do_redo_fast
        num_exp = 20;
        num_sub = 12;
        fig4 = [];
    else
        num_exp = 500;
        num_sub = 100;
        fig4 = amt_cache('get','fig4',flags.cachemode);
    end
    
    if isempty(fig4)
        offset = 13; % start at 14th file in folder 
        % Get names of all used hrtfs
        x=amt_load('reijniers2014','hrtfnames.mat');  
        % Load SOFA files  
        for i=1:num_sub             
            SOFA_obj(i)=amt_load('reijniers2014',x.hrtfnames{i+offset});
        end
        % Preallocation
        fig4 = struct('exp',struct(), ...
                      'err',struct([]), ...
                      'bias',struct([]), ...
                      'doa_real',struct([]));
        fig4 = repmat(fig4,num_sub,1);
        
        % Compute mean spherical error for all subjects 
        amt_disp(['Processing ' num2str(num_sub) ...
                  ' subjects in parallel...']);
        for i = 1:num_sub
            amt_disp(['Processing subject #' num2str(i)]);
            % Get directions from SOFA files
            [template, target] = reijniers2014_featureextraction(SOFA_obj(i));
            doa = reijniers2014(template, target,'num_exp',num_exp);
            
            [fig4(i).err, fig4(i).bias] = reijniers2014_metrics(doa);
            fig4(i).doa_real = doa.real;
            fig4(i).exp = reijniers2014_metrics(doa, 'middle_metrics');
        end
        
        if ~flags.do_redo_fast
            amt_cache('set','fig4',fig4);
        end
    end
 
    % remove ill formed hrtf
    if(length(fig4(51).err) ~= length(fig4(1).err))
        fig4(51) = [];
    end
    if flags.do_plot
        % Calculate averages
        mean_error = zeros(size(fig4(1).doa_real, 1), 1);
        mean_bias = zeros(size(fig4(1).doa_real, 1), 3);
        
        for k = 1:length(fig4)
            mean_error = mean_error + fig4(k).err;
            mean_bias = mean_bias + fig4(k).bias;
        end
        
        mean_error = abs(mean_error/length(fig4)); 
        mean_bias = mean_bias/length(fig4); 
        
        [~, cbar] = plot_reijniers2014(fig4(1).doa_real, ... % assuming same dirs
                                mean_error,'bias', ...
                                mean_bias,flags.plot_type,flags.type);                        
        title('Fig. 4(a), simulation');
        caxis([0,35]);
        cbar.Label.String = 'Mean spherical error [^\circ]';
        cbar.Label.FontSize = 12;    
        ct1=get(cbar,'TickLabels');
        for k=1:numel(ct1)
            ct1{k}=sprintf('%s',ct1{k});
        end
        set(cbar,'xticklabel',ct1);
    end
    
    metrics =  struct2cell(fig4);
    metrics = cell2mat(metrics(1,:));
    amt_disp(sprintf('\t\t\tSIM'))
    amt_disp(sprintf('lateral_bias [dg]:\t%0.2f', ...
        mean([metrics.accL], 'all')),'documentation');
    amt_disp(sprintf('lateral_rms_error [deg]:\t%0.2f', ...
        mean([metrics.rmsL], 'all')),'documentation');
    amt_disp(sprintf('elevation_bias [deg]:\t%0.2f', ...
        mean([metrics.accE], 'all')),'documentation');
    amt_disp(sprintf('local_rms_polar [deg]:\t%0.2f', ...
        mean([metrics.rmsP], 'all')),'documentation');
    amt_disp(sprintf('quadrant_err [%%]:\t%0.2f', ...
        mean([metrics.querr], 'all')),'documentation');
end

%% ------ fig5 ------------------------------------------------------------
if flags.do_fig5
    if flags.do_redo_fast
        num_exp = 20;
        num_sub = 12;
        fig5 = [];
    else
        num_exp = 500;
        num_sub = 100;
        fig5 = amt_cache('get','fig5',flags.cachemode);
    end
    
    if isempty(fig5) 
        offset = 13; % start at 14th file in folder 
        % Get names of all used hrtfs
        x=amt_load('reijniers2014','hrtfnames.mat');  
        % Load SOFA files  
        for i=1:num_sub             
            SOFA_obj(i)=amt_load('reijniers2014',x.hrtfnames{i+offset});
        end
        
        % preallocation
        fig5 = struct('a',struct('err',[], 'exp', struct([])), ...
                      'b',struct('err',[], 'exp', struct([])), ...
                      'c',struct('err',[], 'exp', struct([])), ...
                      'd',struct('err',[], 'exp', struct([])), ...
                      'e',struct('err',[], 'exp', struct([])), ...
                      'doa_real', []);
        fig5 = repmat(fig5,num_sub,1);

        amt_disp(['Processing ' num2str(num_sub) ...
                  ' subjects in parallel...']);
        % Compute mean sphercial error for all subjects
        for i = 1:num_sub       
            %% Preprocessing source and template information 
            amt_disp(['Processing subject #' num2str(i)]);
            [template, target] = reijniers2014_featureextraction(SOFA_obj(i));

            % Model
            % Fig. 5 (a)
            amt_disp('Doing control');
            doa = reijniers2014(template, target, ...
                                    'num_exp',num_exp);
            fig5(i).a.err = reijniers2014_metrics(doa);
            fig5(i).a.exp = reijniers2014_metrics(doa, 'middle_metrics');
            
            % Fig. 5 (b)
            amt_disp('Doing sigma_itd doubled');
            doa = reijniers2014(template, target, ...
                                    'sig_itd',0.569*2, 'num_exp',num_exp);                         
            fig5(i).b.err = reijniers2014_metrics(doa);
            fig5(i).b.exp = reijniers2014_metrics(doa, 'middle_metrics');

            % Fig. 5 (c)
            amt_disp('Doing sigma_I doubled');
            doa = reijniers2014(template, target, ...
                                    'sig_I',3.5*2,'num_exp',num_exp); 
            fig5(i).c.err = reijniers2014_metrics(doa);
            fig5(i).c.exp = reijniers2014_metrics(doa, 'middle_metrics');

            % Fig. 5 (d)
            amt_disp('Doing sigma_S doubled');
            doa = reijniers2014(template, target, ...
                                    'sig_S',3.5*2,'num_exp',num_exp);                         
            fig5(i).d.err = reijniers2014_metrics(doa);
            fig5(i).d.exp = reijniers2014_metrics(doa, 'middle_metrics');

            % Fig. 5 (e)
            amt_disp('Doing LPF');
            [templpf, targlpf] = reijniers2014_featureextraction(SOFA_obj(i), ...
                                     'fb_low', 3e2, 'fb_high', 8e3);
            doa = reijniers2014(templpf, targlpf, ...
                                     'num_exp',num_exp);
            fig5(i).e.err = reijniers2014_metrics(doa);
            fig5(i).e.exp = reijniers2014_metrics(doa, 'middle_metrics');

            fig5(i).doa_real = doa.real;
        end
                    
        if ~flags.do_redo_fast
            amt_cache('set','fig5',fig5);
        end
    end
    
    graphs_lb = {'a', 'b', 'c', 'd', 'e'};

    if flags.do_plot
        % remove anomalous subject
        fig5(51) = [];
        num_sub = length(fig5);
        % Calculate averages and plot
        mean_error = zeros(length(fig5(1).a.err), length(graphs_lb));

        titles = {'Fig. 5 (a), control', ...
                  'Fig. 5 (b), $2\sigma_{itd}$', ...
                  'Fig. 5 (c), $2\sigma_{I}$',...
                  'Fig. 5 (d), $2\sigma_{S}$',...
                  'Fig. 5 (e), LPF'};

        for j =1:length(graphs_lb)
            for k = 1:num_sub
                mean_error(:,j) = mean_error(:,j) + fig5(k).(graphs_lb{j}).err;
            end

            mean_error(:,j) = abs((mean_error(:,j)/num_sub)); 

            [~, cbar] = plot_reijniers2014(fig5(1).doa_real, ...
                                mean_error(:,j),flags.plot_type,flags.type);                        
            title(titles{j},'Interpreter','Latex');
            caxis([0,20]);
            cbar.Label.String = 'Mean spherical error [^\circ]';
            cbar.Label.FontSize = 12;    
            ct1=get(cbar,'TickLabels');
            for k=1:numel(ct1)
                ct1{k}=sprintf('%s',ct1{k});
            end
            set(cbar,'xticklabel',ct1);
        end
    end
    
    accL = zeros(length(fig5(1).a.err), length(graphs_lb));
    rmsL = zeros(length(fig5(1).a.err), length(graphs_lb));
    accE = zeros(length(fig5(1).a.err), length(graphs_lb));
    rmsP = zeros(length(fig5(1).a.err), length(graphs_lb));
    querr = zeros(length(fig5(1).a.err), length(graphs_lb));
    
    for j =1:length(graphs_lb)
        for k = 1:num_sub
            accL(k,j) = accL(k,j) + fig5(k).(graphs_lb{j}).exp.accL;
            rmsL(k,j) = rmsL(k,j) + fig5(k).(graphs_lb{j}).exp.rmsL;
            accE(k,j) = accE(k,j) + fig5(k).(graphs_lb{j}).exp.accE;
            rmsP(k,j) = rmsP(k,j) + fig5(k).(graphs_lb{j}).exp.rmsP;
            querr(k,j) = querr(k,j) + fig5(k).(graphs_lb{j}).exp.querr;
        end

        amt_disp(sprintf('\n EXP %s', graphs_lb{j}),'documentation');
        amt_disp(sprintf('\t\t\tSIM'),'documentation');
        amt_disp(sprintf('lateral_bias [deg]:\t%0.3f', ...
            mean(accL(:,j), 'all')),'documentation');
        amt_disp(sprintf('lateral_rms_error [deg]:\t%0.3f', ...
            mean(rmsL(:,j), 'all')),'documentation');
        amt_disp(sprintf('elevation_bias [deg]:\t%0.3f', ...
            mean(accE(:,j), 'all')),'documentation');
        amt_disp(sprintf('local_rms_polar [deg]:\t%0.3f', ...
            mean(rmsP(:,j), 'all')),'documentation');
        amt_disp(sprintf('quadrant_err [%%]:\t%0.3f', ...
            mean(querr(:,j), 'all')),'documentation');
    end
end

%% ------ fig6 ------------------------------------------------------------
if flags.do_fig6 
    if flags.do_redo_fast
        num_exp = 20;
        num_sub = 12;
        fig6 = [];
    else
        num_exp = 500;
        num_sub = 100;
        fig6 = amt_cache('get','fig6',flags.cachemode);
    end
    
    SNRs = [75, 50, 40, 30];
    graphs_lb = {'a', 'b', 'c', 'd'};

    if isempty(fig6) 
        offset = 13; % start at 14th file in folder 
        % Get names of all used hrtfs
        x=amt_load('reijniers2014','hrtfnames.mat');
        % Load SOFA files
        for i=1:num_sub   
            SOFA_obj(i)=amt_load('reijniers2014',x.hrtfnames{i+offset});
        end

        % preallocation
        fig6 = struct('a',struct('err',[], 'exp', struct([])), ...
                      'b',struct('err',[], 'exp', struct([])), ...
                      'c',struct('err',[], 'exp', struct([])), ...
                      'd',struct('err',[], 'exp', struct([])), ...
                      'doa_real', []);
        fig6 = repmat(fig6,num_sub,1);
        
        % Compute mean sphercial error for all subjects 
        amt_disp(['Processing ' num2str(num_sub) ...
                  ' subjects in parallel...']);
        for i = 1:num_sub
            amt_disp(['Processing subject #' num2str(i)]);
            [template, target] = reijniers2014_featureextraction(SOFA_obj(i));

            for j = 1:length(SNRs)
                amt_disp(sprintf('Doing SNR=%idB', SNRs(j)));
                doa = reijniers2014(template, target, ...
                                    'SNR',SNRs(j),'num_exp',num_exp);
                fig6(i).(graphs_lb{j}).err = reijniers2014_metrics(doa);
                fig6(i).(graphs_lb{j}).exp = reijniers2014_metrics(doa, 'middle_metrics');
                fig6(i).doa_real = doa.real;
            end
        end

        if ~flags.do_redo_fast
            amt_cache('set','fig6',fig6);
        end
    end
    
    if flags.do_plot
        num_sub = length(fig6);
        % Calculate averages and plot
        mean_error = zeros(length(fig6(1).doa_real), length(graphs_lb));

        for j =1:length(graphs_lb)
            for k = 1:num_sub
                mean_error(:,j) = mean_error(:,j) + fig6(k).(graphs_lb{j}).err;
            end
            mean_error(:,j) = abs((mean_error(:,j)/num_sub)); 

            [~, cbar] = plot_reijniers2014(fig6(1).doa_real, ...
                                mean_error(:,j),flags.plot_type,flags.type);                        
            title(sprintf('Fig.6 (%c), SNR=%idB',graphs_lb{j}, SNRs(j)),...
                        'Interpreter','Latex');
            caxis([0,20]);
            cbar.Label.String = 'Mean spherical error [^circ]';
            cbar.Label.FontSize = 12;    
            ct1=get(cbar,'TickLabels');
            for k=1:numel(ct1)
                ct1{k}=sprintf('%s',ct1{k});
            end
            set(cbar,'xticklabel',ct1);
        end     
    end
    
    accL = zeros(length(fig6(1).a.err), length(graphs_lb));
    rmsL = zeros(length(fig6(1).a.err), length(graphs_lb));
    accE = zeros(length(fig6(1).a.err), length(graphs_lb));
    rmsP = zeros(length(fig6(1).a.err), length(graphs_lb));
    querr = zeros(length(fig6(1).a.err), length(graphs_lb));
    
    for j =1:length(graphs_lb)
        for k = 1:num_sub
            accL(:,j) = accL(:,j) + fig6(k).(graphs_lb{j}).exp.accL;
            rmsL(:,j) = rmsL(:,j) + fig6(k).(graphs_lb{j}).exp.rmsL;
            accE(:,j) = accE(:,j) + fig6(k).(graphs_lb{j}).exp.accE;
            rmsP(:,j) = rmsP(:,j) + fig6(k).(graphs_lb{j}).exp.rmsP;
            querr(:,j) = querr(:,j) + fig6(k).(graphs_lb{j}).exp.querr;
        end

        amt_disp(sprintf('\n EXP %s', graphs_lb{j}),'documentation');
        amt_disp(sprintf('\t\t\tSIM'),'documentation');
        amt_disp(sprintf('lateral_bias [deg]:\t%0.2f', ...
            mean([accL], 'all')),'documentation');
        amt_disp(sprintf('lateral_rms_error [deg]:\t%0.2f', ...
            mean([rmsL], 'all')),'documentation');
        amt_disp(sprintf('elevation_bias [deg]:\t%0.2f', ...
            mean([accE], 'all')),'documentation');
        amt_disp(sprintf('local_rms_polar [deg]:\t%0.2f', ...
            mean([rmsP], 'all')),'documentation');
        amt_disp(sprintf('quadrant_err [%%]:\t%0.2f', ...
            mean([querr], 'all')),'documentation');
    end
end

%% ------ tab1_barumerli2020aes ----------------------------------------------
if flags.do_tab1_barumerli2020aes
    data_baseline = data_majdak2010('Learn_M');
    % remove subjects with no data
    data_baseline(1:5) = [];

    % uncertainties tuner 
    multiplier = 3;
    amt_disp(sprintf('Uncertanties multiplied by a factor: %i\n', multiplier));
    
    tab1_barumerli2020aes = [];
    
    if ~flags.do_redo
        tab1_barumerli2020aes = amt_cache('get','tab1_barumerli2020aes',flags.cachemode);
    end
    
    if isempty(tab1_barumerli2020aes)
        amt_disp('Loading SOFA files');

        for i = 1:length(data_baseline)
            data_baseline(i).sofa = amt_load('reijniers2014',['dtf_' lower(data_baseline(i).id) '.sofa']);
        end

        % Preallocation
        tab1_barumerli2020aes = struct('doa', struct([]));
        tab1_barumerli2020aes = repmat(tab1_barumerli2020aes,length(data_baseline),1);

        for i = 1:length(data_baseline)
            amt_disp(['Processing subject #' num2str(i)]);
            % Get directions from SOFA files
            targ_az = data_baseline(i).mtx(:,1);
            targ_el = data_baseline(i).mtx(:,2);
            % preprocessing
            [template, target] = ...
                reijniers2014_featureextraction(data_baseline(i).sofa, ...
                'targ_az', targ_az, 'targ_el', targ_el);

            % model esecution
            [tab1_barumerli2020aes(i).doa] = ... 
                reijniers2014(template, target, 'num_exp', 1, ...
                            'sig_itd', 0.569*multiplier, ...
                            'sig_I', 3.5*multiplier, ...
                            'sig_S', 3.5*multiplier, ...
                            'sig', 5*multiplier);
        end

        amt_cache('set','tab1_barumerli2020aes',tab1_barumerli2020aes);
    end
        
    for i = 1:length(tab1_barumerli2020aes)
        % lateral_bias 
        exp(i).accL = reijniers2014_metrics(tab1_barumerli2020aes(i).doa, 'accL'); 
        % lateral_rms_error
        exp(i).rmsL = reijniers2014_metrics(tab1_barumerli2020aes(i).doa, 'rmsL'); 
        % elevation_bias
        exp(i).accE = reijniers2014_metrics(tab1_barumerli2020aes(i).doa, 'accE'); 
        % local_rms_polar
        exp(i).rmsP = ... 
            reijniers2014_metrics(tab1_barumerli2020aes(i).doa, 'rmsPmedianlocal'); 
        % quadrant_err
        exp(i).querr = ...
            reijniers2014_metrics(tab1_barumerli2020aes(i).doa, 'querrMiddlebrooks'); 

        real(i).accL = localizationerror(data_baseline(i).mtx, 'accL');
        real(i).rmsL = localizationerror(data_baseline(i).mtx, 'rmsL');
        real(i).accE = localizationerror(data_baseline(i).mtx, 'accE');
        real(i).rmsP = ...
            localizationerror(data_baseline(i).mtx, 'rmsPmedianlocal');
        real(i).querr = ...
            localizationerror(data_baseline(i).mtx, 'querrMiddlebrooks');
    end
    
    amt_disp(sprintf('\t\t\tSIM\t\tREAL'),'documentation');
    amt_disp(sprintf('lateral_bias [deg]:\t%0.2f\t\t%0.2f', ...
        mean([exp.accL]), mean([real.accL])),'documentation');
    amt_disp(sprintf('lateral_rms_error [deg]:\t%0.2f\t\t%0.2f', ...
        mean([exp.rmsL]), mean([real.rmsL])),'documentation');
    amt_disp(sprintf('elevation_bias [deg]:\t%0.2f\t\t%0.2f', ...
        mean([exp.accE]), mean([real.accE])),'documentation');
    amt_disp(sprintf('local_rms_polar [deg]:\t%0.2f\t\t%0.2f', ...
        mean([exp.rmsP]), mean([real.rmsP])),'documentation');
    amt_disp(sprintf('quadrant_err [%%]:\t%0.2f\t\t%0.2f', ...
        mean([exp.querr]), mean([real.querr])),'documentation');
end


%% ------ fig2_barumerli2020forum -------------------------------------
if flags.do_fig2_barumerli2020forum
    fig2_barumerli2020forum = [];

    if ~flags.do_redo
        fig2_barumerli2020forum = amt_cache('get', ...
            'fig2_barumerli2020forum',flags.cachemode);
    end
    
    if isempty(fig2_barumerli2020forum)
        % load 23 DTFs from baumgartner2014         
        amt_disp('Loading SOFA files');
        sbj_dtf = data_baumgartner2014('pool','cached');

        % preprocess templates for each user
        amt_disp('Processing subjects'' templates');
        
        for i = 1:length(sbj_dtf)
            amt_disp(['Pre-processing subject #' num2str(i)]);
            [sbj_template(i), sbj_target(i)] = reijniers2014_featureextraction(sbj_dtf(i).Obj);
        end
        
        % preallocation for results
        amt_disp('Allocating memory for results');
        estimations = struct('doa', struct([]));
        estimations = repmat(estimations, ...
            length(sbj_dtf),length(sbj_dtf)); % all vs all

        for i = 1:length(sbj_dtf)
            amt_disp(['Processing subject #' num2str(i)]);
            for j = 1:length(sbj_dtf)
                amt_disp(num2str(j));
                %TODO: select points in the sphere grid
                estimations(i, j).doa = ... 
                    reijniers2014(sbj_template(i), sbj_target(j), 'num_exp', 1);
            end
        end

        % compute metrics
        for i = 1:size(estimations, 1)
            for j = 1:size(estimations, 2)
                % lateral bias
                metrics(i, j).accL = reijniers2014_metrics(estimations(i, j).doa, 'accL'); 
                % lateral_rms_error
                metrics(i, j).rmsL = reijniers2014_metrics(estimations(i, j).doa, 'rmsL'); 
                % elevation_bias
                metrics(i, j).accE = reijniers2014_metrics(estimations(i, j).doa, 'accE'); 
                % polar bias
                metrics(i, j).accP = ... 
                    reijniers2014_metrics(estimations(i, j).doa, 'accP'); 
                % local rms polar
                metrics(i, j).rmsP = ... 
                    reijniers2014_metrics(estimations(i, j).doa, 'rmsPmedianlocal'); 
                % quadrant_err
                metrics(i, j).querr = ...
                    reijniers2014_metrics(estimations(i, j).doa, 'querrMiddlebrooks'); 
            end
        end
        
        fig2_barumerli2020forum = metrics;
        
        amt_cache('set','fig2_barumerli2020forum',fig2_barumerli2020forum);
    end
    
    metrics = fig2_barumerli2020forum;
    varargout{1} = metrics;

    % aggregate metrics
    ns = size(metrics,1);
    own = logical(eye(ns));
    other = not(own);
    quants = [0,0.05,0.25,0.5,0.75,0.95,1];
    % code similar to baumgartner2014 - fig9
    le_own.quantiles = quantile([metrics(own).rmsL], quants);
    lb_own.quantiles = quantile([metrics(own).accL], quants);
    qe_own.quantiles = quantile([metrics(own).querr], quants);
    pe_own.quantiles = quantile([metrics(own).rmsP], quants);
    pb_own.quantiles = quantile([metrics(own).accP], quants);
    le_own.mean = mean([metrics(own).rmsL]);
    lb_own.mean = mean([metrics(own).accL]);
    qe_own.mean = mean([metrics(own).querr]);
    pe_own.mean = mean([metrics(own).rmsP]);
    pb_own.mean = mean([metrics(own).accP]);

    le_other.quantiles = quantile([metrics(other).rmsL], quants);
    lb_other.quantiles = quantile([metrics(other).accL], quants);
    qe_other.quantiles = quantile([metrics(other).querr], quants);
    pe_other.quantiles = quantile([metrics(other).rmsP], quants);
    pb_other.quantiles = quantile([metrics(other).accP], quants);
    le_other.mean = mean([metrics(other).rmsL]);
    lb_other.mean = mean([metrics(other).accL]);
    qe_other.mean = mean([metrics(other).querr]);
    pe_other.mean = mean([metrics(other).rmsP]);
    pb_other.mean = mean([metrics(other).accP]);
    
    data = data_middlebrooks1999;
    
    % baumgartner data
    data_baum_temp = exp_baumgartner2014('fig9', 'no_plot');
    data_baum.qe_pool = data_baum_temp(1).qe;
    data_baum.pe_pool = data_baum_temp(1).pe;
    data_baum.pb_pool = data_baum_temp(1).pb;
    
    ns = size(data_baum.pe_pool,1);
    own = eye(ns) == 1;
    other = not(own);
    data_baum.pb_pool = abs(data_baum.pb_pool);
    data_baum.qe_own.quantiles = quantile(data_baum.qe_pool(own),[0,0.05,0.25,0.5,0.75,0.95,1]);
    data_baum.pe_own.quantiles = quantile(data_baum.pe_pool(own),[0,0.05,0.25,0.5,0.75,0.95,1]);
    data_baum.pb_own.quantiles = quantile(data_baum.pb_pool(own),[0,0.05,0.25,0.5,0.75,0.95,1]);
    data_baum.qe_own.mean = mean(data_baum.qe_pool(own));
    data_baum.pe_own.mean = mean(data_baum.pe_pool(own));
    data_baum.pb_own.mean = mean(data_baum.pb_pool(own));

    data_baum.qe_other.quantiles = quantile(data_baum.qe_pool(other),[0,0.05,0.25,0.5,0.75,0.95,1]);
    data_baum.pe_other.quantiles = quantile(data_baum.pe_pool(other),[0,0.05,0.25,0.5,0.75,0.95,1]);
    data_baum.pb_other.quantiles = quantile(data_baum.pb_pool(other),[0,0.05,0.25,0.5,0.75,0.95,1]);
    data_baum.qe_other.mean = mean(data_baum.qe_pool(other));
    data_baum.pe_other.mean = mean(data_baum.pe_pool(other));
    data_baum.pb_other.mean = mean(data_baum.pb_pool(other));
  
    % plot
    if flags.do_plot
        dx = 0.22;
        dx_lat = 0.15;
        Marker = 's-';
        LineColor = [0 0.4470 0.7410];
        data.Marker = 'ko-';
        data.LineColor = [1 1 1]*0.3;
        data_baum.Marker = 'd-';
        data_baum.LineColor = [0.8500 0.3250 0.0980];
        
        mFig = figure;
        mFig.Units = 'centimeters';
        mFig.Position = [5,5,35,10];
        subplot(1, 5, 1)
        middlebroxplot(1-dx_lat,data.le_own.quantiles,kv.MarkerSize, data.LineColor)
        plot(1-dx_lat,data.le_own.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.LineColor)
        middlebroxplot(1+dx_lat,le_own.quantiles,kv.MarkerSize, LineColor)
        plot(1+dx_lat,le_own.mean,Marker,'MarkerSize', kv.MarkerSize, 'MarkerFaceColor', LineColor, 'MarkerEdgeColor', LineColor)
        
        middlebroxplot(2-dx_lat,data.le_other.quantiles,kv.MarkerSize, data.LineColor)
        plot(2-dx_lat,data.le_other.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.LineColor)
        middlebroxplot(2+dx_lat,le_other.quantiles,kv.MarkerSize, LineColor)
        plot(2+dx_lat,le_other.mean,Marker,'MarkerSize',kv.MarkerSize, 'MarkerFaceColor', LineColor, 'MarkerEdgeColor', LineColor)
        
        ylabel('RMS Lateral Error [deg]','FontSize',kv.FontSize)
        set(gca,'YLim',[-10 60],'XLim',[0.5 2.5],...
          'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
            'TickLength',2*get(gca,'TickLength'))

        subplot(1, 5, 2)
        middlebroxplot(1-dx_lat,data.lb_own.quantiles,kv.MarkerSize, data.LineColor)
        plot(1-dx_lat,data.lb_own.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.LineColor)
        middlebroxplot(1+dx_lat,lb_own.quantiles,kv.MarkerSize, LineColor)
        plot(1+dx_lat,lb_own.mean,Marker,'MarkerSize',kv.MarkerSize, 'MarkerFaceColor', LineColor, 'MarkerEdgeColor', LineColor)
        
        middlebroxplot(2-dx_lat,data.lb_other.quantiles,kv.MarkerSize, data.LineColor)
        plot(2-dx_lat,data.lb_other.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.LineColor)
        middlebroxplot(2+dx_lat,lb_other.quantiles,kv.MarkerSize, LineColor)
        plot(2+dx_lat,lb_other.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor', LineColor, 'MarkerEdgeColor', LineColor)
        ylabel('Magnitude Lateral Bias [deg]','FontSize',kv.FontSize)
        set(gca,'YLim',[-10 60],'XLim',[0.5 2.5],...
          'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
            'TickLength',2*get(gca,'TickLength'))

        subplot(1, 5, 3)
        middlebroxplot(1-dx,data.qe_own.quantiles,kv.MarkerSize, data.LineColor)
        midd = plot(1-dx,data.qe_own.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.LineColor);
        middlebroxplot(1,qe_own.quantiles,kv.MarkerSize, LineColor)
        reij = plot(1,qe_own.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor', LineColor, 'MarkerEdgeColor', LineColor);
        middlebroxplot(1+dx,data_baum.qe_own.quantiles,kv.MarkerSize, data_baum.LineColor)
        baum = plot(1+dx,data_baum.qe_own.mean,data_baum.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor', data_baum.LineColor, 'MarkerEdgeColor', data_baum.LineColor);
        
        middlebroxplot(2-dx,data.qe_other.quantiles,kv.MarkerSize, data.LineColor)
        plot(2-dx,data.qe_other.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.LineColor)
        middlebroxplot(2,qe_other.quantiles,kv.MarkerSize, LineColor)
        plot(2,qe_other.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor', LineColor, 'MarkerEdgeColor', LineColor)
        middlebroxplot(2+dx,data_baum.qe_other.quantiles,kv.MarkerSize, data_baum.LineColor)
        plot(2+dx,data_baum.qe_other.mean,data_baum.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor', data_baum.LineColor, 'MarkerEdgeColor', data_baum.LineColor)
        
        ylabel('Quadrant Errors (%)','FontSize',kv.FontSize)
        set(gca,'YLim',[0 50],'XLim',[0.5 2.5],...
          'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
            'TickLength',2*get(gca,'TickLength'))
        leg = legend([midd, reij, baum], {'Actual', 'SA', 'SP'}, 'Location', 'none');
        leg.Units = 'normalized';
        leg.Position = [0.4745,0.831536390309064,0.097524379329044,0.159029645257883];

        subplot(1, 5, 4)
        middlebroxplot(1-dx,data.pe_own.quantiles,kv.MarkerSize, data.LineColor)
        plot(1-dx,data.pe_own.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.LineColor)
        middlebroxplot(1,pe_own.quantiles,kv.MarkerSize, LineColor)
        plot(1,pe_own.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',LineColor, 'MarkerEdgeColor', LineColor)
        middlebroxplot(1+dx,data_baum.pe_own.quantiles,kv.MarkerSize, data_baum.LineColor)
        plot(1+dx,data_baum.pe_own.mean,data_baum.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data_baum.LineColor, 'MarkerEdgeColor', data_baum.LineColor)
        
        middlebroxplot(2-dx,data.pe_other.quantiles,kv.MarkerSize, data.LineColor)
        plot(2-dx,data.pe_other.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.LineColor)
        middlebroxplot(2,pe_other.quantiles,kv.MarkerSize, LineColor)
        plot(2,pe_other.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',LineColor, 'MarkerEdgeColor', LineColor)
        middlebroxplot(2+dx,data_baum.pe_other.quantiles,kv.MarkerSize, data_baum.LineColor)
        plot(2+dx,data_baum.pe_other.mean,data_baum.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data_baum.LineColor, 'MarkerEdgeColor', data_baum.LineColor)
        
        ylabel('Local Polar RMS Error (deg)','FontSize',kv.FontSize)
        set(gca,'YLim',[-10 60],'XLim',[0.5 2.5],...
          'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
            'TickLength',2*get(gca,'TickLength'))

        subplot(1, 5, 5)
        middlebroxplot(1-dx,data.pb_own.quantiles,kv.MarkerSize, data.LineColor)
        plot(1-dx,data.pb_own.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.LineColor)
        middlebroxplot(1,pb_own.quantiles,kv.MarkerSize, LineColor)
        plot(1,pb_own.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',LineColor, 'MarkerEdgeColor', LineColor)
        middlebroxplot(1+dx,data_baum.pb_own.quantiles,kv.MarkerSize, data_baum.LineColor)
        plot(1+dx,data_baum.pb_own.mean,data_baum.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data_baum.LineColor, 'MarkerEdgeColor', data_baum.LineColor)
        
        middlebroxplot(2-dx,data.pb_other.quantiles,kv.MarkerSize, data.LineColor)
        plot(2-dx,data.pb_other.mean,data.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',data.LineColor)
        middlebroxplot(2,pb_other.quantiles,kv.MarkerSize, LineColor)
        plot(2,pb_other.mean,Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor',LineColor, 'MarkerEdgeColor', LineColor)
        middlebroxplot(2+dx,data_baum.pb_other.quantiles,kv.MarkerSize, data_baum.LineColor)
        plot(2+dx,data_baum.pb_other.mean,data_baum.Marker,'MarkerSize',kv.MarkerSize,'MarkerFaceColor', data_baum.LineColor, 'MarkerEdgeColor', data_baum.LineColor)
        
        ylabel('Magnitude of Elevation Bias (deg)','FontSize',kv.FontSize)
        set(gca,'YLim',[-10 60],'XLim',[0.5 2.5],...
          'XTick',1:2,'XTickLabel',{'Own' 'Other'},'FontSize',kv.FontSize,...
            'TickLength',2*get(gca,'TickLength'))
    end
end

%% ------ fig3_barumerli2020forum -------------------------------------
if flags.do_fig3_barumerli2020forum
    fig3_barumerli2020forum = [];

    if ~flags.do_redo
        fig3_barumerli2020forum = amt_cache('get', ...
            'fig3_barumerli2020forum',flags.cachemode);
    end
    
    if isempty(fig3_barumerli2020forum)
        % Settings
        num_exp = 1;

        % Load Data
        % Speech Samples from Harvard Word list
        speechsample = amt_cache('get','../experiments%2Fexp_baumgartner2014.m/best2005speechSamples');
        samples_num = length(speechsample);
        % FIR Low-pass filters at 8kHz
        % Brick-wall (aka sinc-filter): fir1(200,1/3) -> -60 dB
        x=amt_load('baumgartner2014','highfreqatten_filters.mat');
        LP{1} = [1 zeros(1,100)];
        LP{2} = x.fir20db;
        LP{3} = x.fir40db;
        LP{4} = x.fir60db;

        % Model Data
        sbj_dtf = data_baumgartner2014('pool', 'cached');

        estimations = struct('doa', struct([]));
        est_expflat = repmat(estimations, length(sbj_dtf), 1); 

        est_expLP = repmat(estimations, ...
            length(sbj_dtf),length(LP),length(speechsample)); 
        
        % Simulations
        amt_disp('Processing subjects HRIR');
        for i = 1:length(sbj_dtf)
            amt_disp(['Pre-processing subject #' num2str(i)]);

            % noise stimulus
            [sbj_template(i), sbj_target_flat(i)] = ...
                reijniers2014_featureextraction(sbj_dtf(i).Obj);

            % flat spectrum estimations
            est_expflat(i, 1).doa = ... 
                reijniers2014(sbj_template(i), sbj_target_flat(i), 'num_exp', num_exp);

            amt_disp('Computing localization');
            for f = 1:length(LP)
                amt_disp(sprintf('Filter %i\n', f))
                for s = 1:samples_num
                    stim = filter(LP{f},1,speechsample{s});

                    [~, trgt] = ...
                            reijniers2014_featureextraction(sbj_dtf(i).Obj, ...
                            'source_ir', stim);
                    est_expLP(i, f, s).doa = ... 
                            reijniers2014(sbj_template(i), trgt, 'num_exp', num_exp);
                end
            end
        end

        % metrics
        % allocate memory for results
        % aggregate over different lateral angles
        ale_expflat = zeros(length(sbj_dtf),1);
        ale_expLP = zeros(length(sbj_dtf),length(LP), samples_num);
        ape_expflat = zeros(length(sbj_dtf),1);
        ape_expLP = zeros(length(sbj_dtf),length(LP), samples_num);
        qe_expflat = zeros(length(sbj_dtf),1);
        qe_expLP = zeros(length(sbj_dtf),length(LP), samples_num);

        for i = 1:length(sbj_dtf)
            amt_disp(['Computing metrics #' num2str(i)]);

            % flat spectrum estimations
            ale_expflat(i,1) = reijniers2014_metrics(est_expflat(i, 1).doa, 'accabsL');
            ape_expflat(i,1) = reijniers2014_metrics(est_expflat(i, 1).doa, 'accabsP');
            qe_expflat(i,1) = reijniers2014_metrics(est_expflat(i, 1).doa, 'querr');

            for f = 1:length(LP)
                for s = 1:samples_num
                    ale_expLP(i, f, s) = reijniers2014_metrics(est_expLP(i, f, s).doa, 'accabsL');
                    ape_expLP(i, f, s) = reijniers2014_metrics(est_expLP(i, f, s).doa, 'accabsP');
                    qe_expLP(i, f, s) = reijniers2014_metrics(est_expLP(i, f, s).doa, 'querr');
                end
            end
        end

        % save cache
        fig3_barumerli2020forum.ale_expflat = ale_expflat;
        fig3_barumerli2020forum.ale_expLP   = ale_expLP;
        fig3_barumerli2020forum.ape_expflat = ape_expflat;
        fig3_barumerli2020forum.ape_expLP   = ape_expLP;        
        fig3_barumerli2020forum.qe_expflat  = qe_expflat;
        fig3_barumerli2020forum.qe_expLP    = qe_expLP;
        
        amt_cache('set','fig3_barumerli2020forum', fig3_barumerli2020forum);
    end
    varargout{1} = [fig3_barumerli2020forum];
    
    ale_expflat = fig3_barumerli2020forum.ale_expflat;
    ale_expLP = fig3_barumerli2020forum.ale_expLP;
    ape_expflat = fig3_barumerli2020forum.ape_expflat;
    ape_expLP = fig3_barumerli2020forum.ape_expLP;        
    qe_expflat = fig3_barumerli2020forum.qe_expflat;
    qe_expLP = fig3_barumerli2020forum.qe_expLP;
    
    % load real and baum2014 data
    data = data_best2005;

    % DCN enabled
    [baum_temp, ~] = exp_baumgartner2014('fig11', 'no_plot'); 
    data_baum.ape = zeros(size(ape_expLP));
    data_baum.qe = zeros(size(qe_expLP));
    
    for i = 1:size(ape_expLP,3)
        data_baum.ape(:,:,i) = transpose(baum_temp{1}(:,:,i));
        data_baum.qe(:,:,i) = transpose(baum_temp{2}(:,:,i));
    end

    data_baum.ape_expflat = baum_temp{3};
    data_baum.qe_expflat = baum_temp{4};
     
    % Pool Samples
    ale_pooled = mean(ale_expLP,3);
    ape_pooled = mean(ape_expLP,3);
    qe_pooled = mean(qe_expLP,3);
    data_baum.ape_pooled = mean(data_baum.ape,3);
    data_baum.qe_pooled = mean(data_baum.qe,3);
    
    % Confidence Intervals or standard errors
    % reijniers model 
    df_speech = size(ale_expLP,1)-1;
    tquant_speech = 1;
    seale_speech = std(ale_pooled,0,1)*tquant_speech/(df_speech+1);
    df_noise = size(ale_expflat,1)-1;
    tquant_noise = 1;
    seale_noise = std(ale_expflat)*tquant_noise/(df_noise+1);
    seale = [seale_noise, seale_speech];
    
    df_speech = size(ape_expLP,1)-1;
    tquant_speech = 1;
    seape_speech = std(ape_pooled,0,1)*tquant_speech/(df_speech+1);
    df_noise = size(ape_expflat,1)-1;
    tquant_noise = 1;
    seape_noise = std(ape_expflat)*tquant_noise/(df_noise+1);
    seape = [seape_noise, seape_speech];
    
    % baumgartner2014 model
    df_speech = size(data_baum.ape,1)-1;
    tquant_speech = 1;
    seape_speech = std(data_baum.ape_pooled,0,1)*tquant_speech/(df_speech+1);
    df_noise = size(data_baum.ape_expflat,1)-1;
    tquant_noise = 1;
    seape_noise = std(data_baum.ape_expflat)*tquant_noise/(df_noise+1);
    data_baum.seape = [seape_noise, seape_speech];
    
    % averages
    % reij2014
    ale = mean([ale_expflat, ale_pooled],1);
    ape = mean([ape_expflat, ape_pooled],1);
    qe = mean([qe_expflat, qe_pooled],1);

    % baum2014
    data_baum.ape = mean([data_baum.ape_expflat', data_baum.ape_pooled],1);
    data_baum.qe = mean([data_baum.qe_expflat', data_baum.qe_pooled],1);

    if flags.do_plot
        MarkerSize = kv.MarkerSize;
        FontSize = kv.FontSize;
        LineColor = [0 0.4470 0.7410];
        Marker = 's-';
        data.Marker = 'o-';
        data.LineColor = [1 1 1]*0.3;
        data_baum.Marker = 'd--';
        data_baum.LineColor = [0.8500 0.3250 0.0980];
        
        mFig = figure;
        mFig.Units = 'centimeters';
        mFig.Position = [5,5,13.5,15];
        
        dx = 0;
        xticks = 0:size(ale_pooled,2);
        subplot(3,1,1)
        plot(xticks-dx,data.ale, data.Marker,'Color', data.LineColor, 'MarkerFaceColor',data.LineColor,'MarkerSize',MarkerSize)
        hold on
        plot(xticks-dx,ale, Marker,'Color', LineColor, 'MarkerFaceColor',LineColor,'MarkerSize',MarkerSize)
        ylabel('Lateral Error (deg)','FontSize',FontSize)
        set(gca,'XTick',xticks,'XTickLabel',[],'FontSize',FontSize)
        set(gca,'XLim',[-0.5 4.5],'YLim',[0 75],'YMinorTick','on')

        subplot(3,1,2)
        plot(xticks,data.ape,data.Marker,'Color', data.LineColor, 'MarkerFaceColor',data.LineColor,'MarkerSize',MarkerSize)
        hold on
        plot(xticks-dx,ape, Marker,'Color', LineColor, 'MarkerFaceColor',LineColor,'MarkerSize',MarkerSize)
        plot(xticks,data_baum.ape, data_baum.Marker,'Color', data_baum.LineColor, 'MarkerFaceColor',data_baum.LineColor,'MarkerSize',MarkerSize)
        ylabel('Polar Error (deg)','FontSize',FontSize)
        set(gca,'XTick',xticks,'XTickLabel',[],'FontSize',FontSize)
        set(gca,'XLim',[-0.5 4.5],'YLim',[0 75],'YMinorTick','on')

        subplot(3,1,3)
        best = plot(xticks([1 2 5]),data.qe([1 2 5]), data.Marker,'Color', data.LineColor, 'MarkerFaceColor',data.LineColor,'MarkerSize',MarkerSize);
        hold on
        reij = plot(xticks-dx,qe, Marker,'Color', LineColor, 'MarkerFaceColor',LineColor,'MarkerSize',MarkerSize);
        baum = plot(xticks,data_baum.qe, data_baum.Marker,'Color', data_baum.LineColor, 'MarkerFaceColor',data_baum.LineColor,'MarkerSize',MarkerSize);
        ylabel('Quadrant Err. (%)','FontSize',FontSize)
        set(gca,'XTick',xticks,'XTickLabel',data.meta,'FontSize',FontSize,...
        'XLim',[-0.5 4.5],'YLim',[-3 54],'YMinorTick','on')
        
        leg = legend([best, reij, baum], {'Actual', 'SA', 'SP'});
        leg.FontSize = FontSize - 1;
        leg.Units = 'centimeters';
        leg.Position = [9.281,12.762,3.466,1.561];
    end
end

if flags.do_fig4_barumerli2020forum
    fig4_barumerli2020forum = [];

    if ~flags.do_redo
        fig4_barumerli2020forum = amt_cache('get', ...
            'fig4_barumerli2020forum',flags.cachemode);
    end
    
    if isempty(fig4_barumerli2020forum)
        % load 23 DTFs from baumgartner2014         
        amt_disp('Loading SOFA files');
        sbj_dtf = data_baumgartner2014('pool','global');        
        num_exp = 5;
        
        % generate stimulus
        % copyed from exp_baumgartner2014/do_fig10
        density = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8]; % ripples/oct
        depth =   10:10:40;        % ripple depth (peak-to-trough) in dB
        
        % 250-ms bursts, 20-ms raised-cosine fade in/out, flat from 0.6-16kHz
        fs = sbj_dtf(1).Obj.Data.SamplingRate;
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

        % Raised-cosine "(i.e., cos^2)" ramp 1/8 octave wide
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
        % Ripples of Experiment II
        Sexp2 = zeros(Nf+1,length(depth),2);  % 3rd dim: 1:0-phase 2:pi-phase
        Sexp2(idlow:idhigh,:,1) = (depth(:)/2*sin(2*pi*1*O+ 0))';  % density: 1 ripple/oct, 0-phase
        Sexp2(idlow:idhigh,:,2) = (depth(:)/2*sin(2*pi*1*O+pi))';  % density: 1 ripple/oct, pi-phase
        Sexp2 = repmat(ramp',[1,length(depth),2]) .* Sexp2;
        Sexp2 = [Sexp2;Sexp2(Nf-1:-1:2,:,:)];
        Sexp2(isnan(Sexp2)) = -100;
        sexp2 = ifftreal(10.^(Sexp2/20),2*Nf);
        sexp2 = circshift(sexp2,Nf);  % IR corresponding to ripple modification
                        
        % preprocess templates for each user
        amt_disp('Processing subjects HRIR');
        for i = 1:length(sbj_dtf)
            % extract directions
            % filter targets' coordinates
            % convert from spherical to horizontal-polar coordinates
            %horpolar_coords = SOFAconvertCoordinates(...
            %                    sbj_dtf(i).Obj.SourcePosition, ...
            %                    'spherical', 'horizontal-polar');
            horpolar_coords = zeros(size(sbj_dtf(i).Obj.SourcePosition));

            [horpolar_coords(:,2),horpolar_coords(:,1),horpolar_coords(:,3)]=...
                cart2sph(sbj_dtf(i).Obj.SourcePosition(:,1),sbj_dtf(i).Obj.SourcePosition(:,3),-sbj_dtf(i).Obj.SourcePosition(:,2));
            
            horpolar_coords(:,1:2)=rad2deg(horpolar_coords(:,1:2));
            horpolar_coords(:,1)=-horpolar_coords(:,1);
            horpolar_coords(:,2)=mod(horpolar_coords(:,2),360);
            % polar in [60, 120]
            % lateral = 0
            idx = find(((horpolar_coords(:, 2) >= 300 ...
                            | horpolar_coords(:, 2) <= 60) ...
                            | (horpolar_coords(:, 2) >= 120 & ...
                                    horpolar_coords(:, 2) <= 240)) ...
                            & (horpolar_coords(:, 1) <= 30 & horpolar_coords(:, 1) >= -30));

            amt_disp(['Pre-processing subject #' num2str(i)]);
            [sbj_template(i), sbj_target_flat(i)] = ...
                reijniers2014_featureextraction(sbj_dtf(i).Obj, ...
                    'targ_az', sbj_dtf(i).Obj.SourcePosition(idx, 1), ...
                    'targ_el', sbj_dtf(i).Obj.SourcePosition(idx, 2));
            
            amt_disp('Densities conditions');
            for j = 1:length(density)
                [~, sbj_target_exp1(i, j)] = ...
                    reijniers2014_featureextraction(sbj_dtf(i).Obj, ...
                    'source_ir', squeeze(sexp1(:, j, :)), ...
                    'targ_az', sbj_dtf(i).Obj.SourcePosition(idx, 1), ...
                    'targ_el', sbj_dtf(i).Obj.SourcePosition(idx, 2));
            end
            
            amt_disp('Depth conditions');
            for j = 1:length(depth)
                [~, sbj_target_exp2(i, j)] = ...
                    reijniers2014_featureextraction(sbj_dtf(i).Obj, ...
                    'source_ir', squeeze(sexp2(:, j, :)), ...
                    'targ_az', sbj_dtf(i).Obj.SourcePosition(idx, 1), ...
                    'targ_el', sbj_dtf(i).Obj.SourcePosition(idx, 2));
            end
        end
        
        % preallocation for results
        amt_disp('Allocating memory for results');
        estimations = struct('doa', struct([]));
        est_expflat = repmat(estimations, length(sbj_dtf)); 
        est_exp1 = repmat(estimations, ...
            length(sbj_dtf),length(density)); 
        est_exp2 = repmat(estimations, ...
            length(sbj_dtf),length(depth)); 

        % simulations
        for i = 1:length(sbj_dtf)
            amt_disp(['Processing subject #' num2str(i)]);
            % flat spectrum estimations
            est_expflat(i, 1).doa = ... 
                reijniers2014(sbj_template(i), sbj_target_flat(i), 'num_exp', num_exp);
            
            % rippled estimations
            for j = 1:length(density)
                est_exp1(i, j).doa = ... 
                    reijniers2014(sbj_template(i), sbj_target_exp1(i, j), 'num_exp', num_exp);
            end
            
            for j =1:length(depth)
                est_exp2(i, j).doa = ... 
                reijniers2014(sbj_template(i), sbj_target_exp2(i, j), 'num_exp', num_exp);
            end
        end        

        % metrics
        % allocate memory for results
        % aggregate over different lateral angles
        pe_exp1 = zeros(length(sbj_dtf),length(density));
        pe_exp2 = zeros(length(sbj_dtf),length(depth));
        pe_flat = zeros(length(sbj_dtf), 1);

        for i = 1:length(sbj_dtf)
            % compute regression (see paper)   
            [f,r] = reijniers2014_metrics(est_expflat(i).doa,'sirpMacpherson2000');
            pe_flat(i) = reijniers2014_metrics(est_expflat(i).doa,f,r,'perMacpherson2003');

            for j = 1:size(est_exp1, 2)
                pe_exp1(i, j) = reijniers2014_metrics(est_exp1(i, j).doa, f, r, 'perMacpherson2003');
            end
            
            for j = 1:size(est_exp2, 2)
                pe_exp2(i, j) = reijniers2014_metrics(est_exp2(i, j).doa, f, r, 'perMacpherson2003');
            end
        end

        
        % save cache
        fig4_barumerli2020forum.pe_flat = pe_flat;
        fig4_barumerli2020forum.pe_exp1 = pe_exp1;
        fig4_barumerli2020forum.pe_exp2 = pe_exp2;
        amt_cache('set','fig4_barumerli2020forum', fig4_barumerli2020forum);
    end
    
    varargout{1} = fig4_barumerli2020forum;

    % simulations data
    pe_flat = fig4_barumerli2020forum.pe_flat;
    pe_exp1 = fig4_barumerli2020forum.pe_exp1;
    pe_exp2 = fig4_barumerli2020forum.pe_exp2;
    
    % Original data:
    data = data_macpherson2003;
    
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
    % simulations
    %pe_exp1 = pe_exp1 - pe_flat(:);
    pe_exp1 = pe_exp1 -    repmat(pe_flat(:), 1, size(pe_exp1, 2));
    %pe_exp2 = pe_exp2 - pe_flat(:);
    pe_exp2 = pe_exp2 -    repmat(pe_flat(:), 1, size(pe_exp2, 2));
    % baumgartner data
    data_baum.pe_exp1 = data_baum.pe_exp1 - repmat(data_baum.pe_flat(:),1,size(data_baum.pe_exp1,2));
    data_baum.pe_exp2 = data_baum.pe_exp2 - repmat(data_baum.pe_flat(:),1,size(data_baum.pe_exp2,2));

    % Statistics
    % simulations
    quart_pe_flat = quantile(pe_flat,[.25 .50 .75]);
    quart_pe_exp1 = quantile(pe_exp1,[.25 .50 .75]);
    quart_pe_exp2 = quantile(pe_exp2,[.25 .50 .75]);
    
    % real data
    data.quart_pe_flat = quantile(data.pe_flat,[.25 .50 .75]);
    data.quart_pe_exp1 = quantile(data.pe_exp1,[.25 .50 .75]);
    data.quart_pe_exp2 = quantile(data.pe_exp2,[.25 .50 .75]);

    % baumgartner data
    data_baum.quart_pe_flat = quantile(data_baum.pe_flat,[.25 .50 .75]);
    data_baum.quart_pe_exp1 = quantile(data_baum.pe_exp1,[.25 .50 .75]);
    data_baum.quart_pe_exp2 = quantile(data_baum.pe_exp2,[.25 .50 .75]);

    % plot
    if flags.do_plot
        dx = 1.05;
        FontSize = kv.FontSize;
        MarkerSize = kv.MarkerSize;

        LineColor = [0 0.4470 0.7410];
        data.Marker = 'ko-';
        data.LineColor = [1 1 1]*0.3;
        data_baum.Marker = 'd-';
        data_baum.LineColor = [0.8500 0.3250 0.0980];
        
        % Exp1
        mFig = figure;
        mFig.Units = 'centimeters';
        mFig.Position = [5,5,13.5,12];
        subplot(2,8,1:8)
        mach = errorbar(data.density,data.quart_pe_exp1(2,:,idphase),...
        data.quart_pe_exp1(2,:,idphase) - data.quart_pe_exp1(1,:,idphase),...
        data.quart_pe_exp1(3,:,idphase) - data.quart_pe_exp1(2,:,idphase),...
        'o-','MarkerSize',MarkerSize, 'Color', data.LineColor, ...
        'MarkerFaceColor', data.LineColor);
        hold on
        reij = errorbar(data.density/dx,quart_pe_exp1(2,:,idphase),...
        quart_pe_exp1(2,:,idphase) - quart_pe_exp1(1,:,idphase),...
        quart_pe_exp1(3,:,idphase) - quart_pe_exp1(2,:,idphase),...
        's-','MarkerSize',MarkerSize, 'Color', LineColor,'MarkerFaceColor',LineColor);
        hold on
        baum = errorbar(data.density*dx,data_baum.quart_pe_exp1(2,:,idphase),...
        data_baum.quart_pe_exp1(2,:,idphase) - data_baum.quart_pe_exp1(1,:,idphase),...
        data_baum.quart_pe_exp1(3,:,idphase) - data_baum.quart_pe_exp1(2,:,idphase),...
        'd--','MarkerSize',MarkerSize, 'Color', data_baum.LineColor,'MarkerFaceColor',data_baum.LineColor);
        
        set(gca,'XScale','log','YMinorTick','on')
        set(gca,'XLim',[0.25/1.2 8*1.2],'XTick',data.density,'YLim',[-16 59],'FontSize',FontSize)
        xlabel('Ripple Density (ripples/octave)','FontSize',FontSize)
        ylabel({'Increase in';'Polar Error Rate (%)'},'FontSize',FontSize)
        
        % Exp2
        subplot(2,8,9:13)
        errorbar(data.depth,data.quart_pe_exp2(2,:,idphase),...
        data.quart_pe_exp2(2,:,idphase) - data.quart_pe_exp2(1,:,idphase),...
        data.quart_pe_exp2(3,:,idphase) - data.quart_pe_exp2(2,:,idphase),...
        'o-','MarkerSize',MarkerSize, 'Color', data.LineColor, ...
        'MarkerFaceColor', data.LineColor);
        hold on
        errorbar(data.depth-1,quart_pe_exp2(2,:,idphase),...
        quart_pe_exp2(2,:,idphase) - quart_pe_exp2(1,:,idphase),...
        quart_pe_exp2(3,:,idphase) - quart_pe_exp2(2,:,idphase),...
        's-','MarkerSize',MarkerSize, 'Color', LineColor,'MarkerFaceColor',LineColor);
        hold on
        errorbar(data.depth+1,data_baum.quart_pe_exp2(2,:,idphase),...
        data_baum.quart_pe_exp2(2,:,idphase) - data_baum.quart_pe_exp2(1,:,idphase),...
        data_baum.quart_pe_exp2(3,:,idphase) - data_baum.quart_pe_exp2(2,:,idphase),...
        'd--','MarkerSize',MarkerSize, 'Color', data_baum.LineColor,'MarkerFaceColor',data_baum.LineColor);

        set(gca,'XLim',[data.depth(1)-5 data.depth(end)+5],'XTick',data.depth,...
        'YLim',[-16 59],'YMinorTick','on','FontSize',FontSize)
        xlabel('Ripple Depth (dB)','FontSize',FontSize)
        ylabel({'Increase in';'Polar Error Rate (%)'},'FontSize',FontSize)
        ytick = get(gca,'YTick');
        ticklength = get(gca,'TickLength');

        % Baseline
        subplot(2,8,14:15)
        errorbar(-0.5,data.quart_pe_flat(2),...
        data.quart_pe_flat(2) - data.quart_pe_flat(1),...
        data.quart_pe_flat(3) - data.quart_pe_flat(2),...
        'o-','MarkerSize',MarkerSize, 'Color', data.LineColor, ...
        'MarkerFaceColor', data.LineColor);
        hold on
        errorbar(0,quart_pe_flat(2),...
        quart_pe_flat(2) - quart_pe_flat(1),...
        quart_pe_flat(3) - quart_pe_flat(2),...
        's-','MarkerSize',MarkerSize, 'Color', LineColor,'MarkerFaceColor',LineColor);
        hold on
        errorbar(0.5,data_baum.quart_pe_flat(2),...
        data_baum.quart_pe_flat(2) - data_baum.quart_pe_flat(1),...
        data_baum.quart_pe_flat(3) - data_baum.quart_pe_flat(2),...
        'd--','MarkerSize',MarkerSize, 'Color', data_baum.LineColor,'MarkerFaceColor',data_baum.LineColor);

        set(gca,'XLim',[-3 3],'XTick',0,'XTickLabel',{'Baseline'},...
        'YLim',[-15 59],'YTick',ytick,'TickLength',3*ticklength,...
        'FontSize',FontSize,'YAxisLocation','right')
        xlabel(' ','FontSize',FontSize)
        ylabel({'Polar Error Rate (%)'},'FontSize',FontSize)

        %legend
        leg = legend([mach, reij, baum], {'Actual', 'SA', 'SP'});
        leg.FontSize = FontSize - 1;
        leg.Units = 'centimeters';
        leg.Position = [9.281,10,3.466,1.561];
        
        % Overall correlation between actual and predicted median values
        m_pe_pred = [quart_pe_exp1(2,:,idphase) quart_pe_exp2(2,:,idphase)];
        m_pe_actual = [data.quart_pe_exp1(2,:,idphase) data.quart_pe_exp2(2,:,idphase)];
        r = corrcoef(m_pe_pred,m_pe_actual);
        r_sqr = r(2);

        amt_disp('Correlation between actual and predicted median values (15 conditions):','documentation');
        amt_disp(['w/  PSGE: r = ' num2str(r_sqr,'%0.2f')],'documentation');
    end
    
end

function middlebroxplot(x,quantiles,MarkerSize,LineColor)
    lilen = 0.1; % length of horizontal lines

    % Symbols
    plot(x,quantiles(1),'x','MarkerSize',MarkerSize, 'MarkerEdgeColor', LineColor) % min
    hold on
    plot(x,quantiles(7),'x','MarkerSize',MarkerSize, 'MarkerEdgeColor', LineColor) % max

    % Horizontal lines
    line(x+0.5*[-lilen,lilen],repmat(quantiles(2),2),'Color',LineColor) % lower whisker
    line(x+[-lilen,lilen],repmat(quantiles(3),2),'Color',LineColor) % 25% Quartile
    line(x+[-lilen,lilen],repmat(quantiles(4),2),'Color',LineColor) % Median
    line(x+[-lilen,lilen],repmat(quantiles(5),2),'Color',LineColor) % 75% Quartile
    line(x+0.5*[-lilen,lilen],repmat(quantiles(6),2),'Color',LineColor) % upper whisker

    % Vertical lines
    line([x,x],quantiles(2:3),'Color',LineColor) % connector lower whisker
    line([x,x],quantiles(5:6),'Color',LineColor) % connector upper whisker
    line([x,x]-lilen,quantiles([3,5]),'Color',LineColor) % left box edge
    line([x,x]+lilen,quantiles([3,5]),'Color',LineColor) % left box edge


