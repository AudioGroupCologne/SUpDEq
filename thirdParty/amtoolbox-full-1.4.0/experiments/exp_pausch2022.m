function fig_handle = exp_pausch2022(varargin)
%EXP_PAUSCH2022 - Reproduce the figures from Pausch et al. (2022)
%
%   Usage: fig_handle = exp_pausch2022(fig)
%
%   The fig flag may be one of:
%
%     'fig5'          : Comparison of measurement-based mean ITDs between 
%                       datasets in the horizontal plane, evaluated in a 
%                       frequency range of 0.5-1.5 kHz and averaged across 
%                       participants with standard deviations sigma (grey 
%                       areas). Dashed and dotted lines show the mean ITD 
%                       differences in front (F) and rear (R) HARTFs, 
%                       respectively, compared to the mean ITDs in HRTFs.
%                       scenarios, represent differences to scenario A.
%     'fig9'         :  a) Measurement-based (light grey) and b), c) model-based 
%                       broadband ITD estimations, colour-coded in grey, black 
%                       and blue, respectively, evaluated for directions in 
%                       the horizontal plane for HRTF and front (F) and 
%                       rear (R) HARTF datasets. Deviations from measurement-
%                       based ITDs are shown as dashed-dotted (Models~1 and~3) 
%                       or dotted lines (Model~2+), with standard deviations 
%                       as shaded areas. d) Scatter plots comparing measurement-
%                       based and model-based ITD maxima, fitted by linear 
%                       regression lines. Box plots show medians and IQRs of 
%                       differences in e) ITD maxima and f) arguments of the 
%                       ITD maxima, with whiskers covering 1.5 times the 
%                       IQR, and outliers displayed as crosses. Horizontal 
%                       black lines indicate non-significant (n.s.) mean 
%                       differences at the 95% confidence level.
%
%
%   Output parameters:
%
%     fig_handle      : Figure handle [matlab.ui.Figure].
%
%
%   Requirements:
%   -------------
%
%   1) SOFA API v0.4.3 or higher from http://sourceforge.net/projects/sofacoustics 
%      for Matlab (in e.g. thirdparty/SOFA)
%
%   2) Data in auxdata/pausch2022 (downloaded on the fly)
%
%
%   Examples:
%   ---------
%
%   To display results of Fig. 5:
%
%     exp_pausch2022('fig5','plot');  
%
%   To display results of Fig. 9:
%     exp_pausch2022('fig9');
%
%
%   References:
%     F. Pausch, S. Doma, and J. Fels. Hybrid multi-harmonic model for the
%     prediction of interaural time differences in individual behind-the-ear
%     hearing-aid-related transfer functions. Acta Acust., 6:34, 2022.
%     [1]http ]
%     
%     References
%     
%     1. https://doi.org/10.1051/aacus/2022020
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_pausch2022.php


%   #Author: Florian Pausch (2022): integrated in the AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Parse flags and keyvals

definput.import={'pausch2022'};
definput.flags.fig = {'missingflag','fig5','fig9'};
%definput.flags.plot = {'plot', 'no_plot'}
[flags,kv] = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.fig{2:end-2}),...
             sprintf('%s or %s',definput.flags.fig{end-1},definput.flags.fig{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% Figure 3
% if flags.do_fig3
%     amt_disp([mfilename,': Not yet implemented.'])
%     fig_handle = [];
% end

%% Figure 5
if flags.do_fig5

    % load directional datasets and estimate ITDs
    [~,itd_hrtf] = data_pausch2022('hrtf');
    [~,itd_hartf_front] = data_pausch2022(char('hartf_front'));
    [~,itd_hartf_rear] = data_pausch2022('hartf_rear');

    % plot data
    font_size = 18;
    line_width = 3.5;
    plot_size = [400,400,650,360];

    phi_vec = 0:2.5:357.5;
    phi_vec = phi_vec(:);
    
    fig_handle = figure('Position',plot_size);
    
    patch([phi_vec;flipud(phi_vec)],...
        [(itd_hrtf(1).itd_hor_mean-itd_hrtf(1).itd_hor_std)*1e6;...
        flipud(itd_hrtf(1).itd_hor_mean+itd_hrtf(1).itd_hor_std)*1e6],...
        kv.col_shad_mod2plus,'EdgeColor','none');
    hold on
    alpha(kv.alpha)
    ph = patch([phi_vec;flipud(phi_vec)],...
        [(itd_hartf_front(1).itd_hor_mean-itd_hartf_front(1).itd_hor_std)*1e6;...
        flipud(itd_hartf_front(1).itd_hor_mean+itd_hartf_front(1).itd_hor_std)*1e6],...
        kv.col_shad_mod1,'EdgeColor','none');
    alpha(kv.alpha)
    patch([phi_vec;flipud(phi_vec)],...
        [(itd_hartf_rear(1).itd_hor_mean-itd_hartf_rear(1).itd_hor_std)*1e6;...
        flipud(itd_hartf_rear(1).itd_hor_mean+itd_hartf_rear(1).itd_hor_std)*1e6],...
        kv.col_shad_meas,'EdgeColor','none');
    alpha(kv.alpha)

    l1 = plot(phi_vec, itd_hrtf(1).itd_hor_mean*1e6,'Color',kv.col_mod2plus,'LineWidth',line_width);
    hold on
    l2 = plot(phi_vec, itd_hartf_front(1).itd_hor_mean*1e6,'Color',kv.col_shad_mod2plus,'LineWidth',line_width);
    l3 = plot(phi_vec, (itd_hrtf(1).itd_hor_mean-itd_hartf_front(1).itd_hor_mean)*1e6,'Color',kv.col_mod1,'LineWidth',line_width,'LineStyle','-.');
    l4 = plot(phi_vec, itd_hartf_rear(1).itd_hor_mean*1e6,'Color',kv.col_meas,'LineWidth',line_width);
    l5 = plot(phi_vec, (itd_hrtf(1).itd_hor_mean-itd_hartf_rear(1).itd_hor_mean)*1e6,'Color',kv.col_shad_mod1,'LineWidth',line_width,'LineStyle',':');

    leg = legend([l1,l2,l4,l3,l5,ph],{'HRTF','HARTF (F)','HARTF (R)','Deviation (F)','Deviation (R)','$\sigma$'},'fontsize',font_size-2);
    %set(leg,'NumColumns', 2,'ItemTokenSize',[14,18],'Location','southwest','interpreter','latex')
    set(leg,'ItemTokenSize',[14,18],'Location','southwest','interpreter','latex')

    xlim([0,357.5])
    grid on
    box on
    set(gca,...
        'xlim',[phi_vec(1),phi_vec(end)],...
        'ylim',[-800,800],...
        'xTick', 0:30:360,...
        'yTick', -800:200:800,...
        'ticklabelinterpreter','latex',...
        'fontsize',font_size)
    xlabel('Azimuth (deg)','fontsize',font_size+1,'interpreter','latex')
    ylabel('ITD ($\mu$s)','fontsize',font_size,'interpreter','latex')
    xtickangle(0)

end

%% Figure 8a
% if flags.do_fig8a
%     amt_disp([mfilename,': Not yet implemented.'])
%     fig_handle = [];
% end

%% Figure 8b
% if flags.do_fig8b
%     amt_disp([mfilename,': Not yet implemented.'])
%     fig_handle = [];
% end

%% Figure 9
if flags.do_fig9

    models = {'kuhn','woodworth_ext','pausch'};
    stf = {'hrtf','hartf_front','hartf_rear'};
    title_stf = {'HRTF', 'HARTF (F)', 'HARTF (R)'};

    num_stf = numel(stf);
    num_mod = numel(models);

    % load directional datasets and estimate ITDs
    data_itd = cell(numel(stf),1);
    for idx_stf = 1:numel(stf)
        [~,data_itd{idx_stf}] = data_pausch2022(stf{idx_stf});
    end

    % load individual features
    features = data_pausch2022('features');

    % extract feature subsets as required per type_stf and type_mod
    feature_stf = cell2mat(features(2:end,4:8));
    feature_stf(:,:,2) = cell2mat(features(2:end,[4:6,33,34]));
    feature_stf(:,:,3) = feature_stf(:,:,2);

    feature_mod2_aux = cell2mat( features(2:end,[8,29,30,37]) ); % x5, d9, d10, Theta3

    % evaluated azimuth directions
    phi_vec = 0:2.5:180;
    phi_vec = phi_vec(:);
    num_dir = numel(phi_vec);

    num_part = size(data_itd{1},2);

    % extract modelled broadband ITDs (and max{ITD} and arg max_phi{ITD})
    itd_mod = zeros(num_dir,num_part,num_stf,num_mod);
    itd_mod_max = zeros(num_part,num_stf,num_mod);
    itd_mod_arg_max_phi = itd_mod_max;

    for idx_stf=1:num_stf
        for idx_part=1:num_part

            for idx_mod=1:num_mod

                switch idx_mod
                    case 1 % kuhn
                        [itd_mod(:,idx_part,idx_stf,idx_mod),...
                         itd_mod_max(idx_part,idx_stf,idx_mod),...
                         itd_mod_arg_max_phi(idx_part,idx_stf,idx_mod)] = ...
                                pausch2022(feature_stf(idx_part,:,idx_stf), ...
                                           models{idx_mod}, ...
                                           stf{idx_stf});
                    case 2 % woodworth_ext
                        [itd_mod(:,idx_part,idx_stf,idx_mod),...
                         itd_mod_max(idx_part,idx_stf,idx_mod),...
                         itd_mod_arg_max_phi(idx_part,idx_stf,idx_mod)] = ...
                                pausch2022(feature_stf(idx_part,:,idx_stf), ...
                                           models{idx_mod},...
                                           stf{idx_stf},...
                                           'x5',feature_mod2_aux(idx_part,1),...
                                           'd9',feature_mod2_aux(idx_part,2),...
                                           'd10',feature_mod2_aux(idx_part,3),...
                                           'Theta3',feature_mod2_aux(idx_part,4));
                    otherwise % pausch
                        [itd_mod(:,idx_part,idx_stf,idx_mod),...
                         itd_mod_max(idx_part,idx_stf,idx_mod),...
                         itd_mod_arg_max_phi(idx_part,idx_stf,idx_mod)] = ...
                                pausch2022(feature_stf(idx_part,:,idx_stf), ...
                                           models{idx_mod},...
                                           stf{idx_stf});
                end

                itd_mod(:,idx_part,idx_stf,idx_mod) = itd_mod(:,idx_part,idx_stf,idx_mod)*1e6;
                itd_mod_max(idx_part,idx_stf,idx_mod) = itd_mod_max(idx_part,idx_stf,idx_mod)*1e6;
                itd_mod_arg_max_phi(idx_part,idx_stf,idx_mod) = phi_vec(itd_mod_arg_max_phi(idx_part,idx_stf,idx_mod));

            end
        end
    end

    % extract measured broadband ITDs (and max{ITD} and arg max_phi{ITD})
    itd_meas = zeros(num_dir,num_part,num_stf);
    itd_meas_max = zeros(num_part,num_stf);
    itd_meas_arg_max_phi = itd_meas_max;

    for idx_stf=1:num_stf
        for idx_part=1:num_part

            itd_meas(:,idx_part,idx_stf) = data_itd{idx_stf}(idx_part).itd_hor(1:numel(phi_vec))*1e6;
            [itd_meas_max(idx_part,idx_stf), itd_meas_arg_max_phi(idx_part,idx_stf)] = ...
                max(itd_meas(:,idx_part,idx_stf));
            itd_meas_arg_max_phi(idx_part,idx_stf) = phi_vec(itd_meas_arg_max_phi(idx_part,idx_stf));

        end
    end

    % calculate direction-dependent mean +- std of ITD deviations across participants
    delta_itd_mean = squeeze( mean(itd_mod-itd_meas,2) );
    delta_itd_std = squeeze( std(itd_mod-itd_meas,0,2) );

    % calculate coefficients of the linear regression line fitting modelled to measured ITD maxima
    reg_coef = zeros(num_stf,2,num_mod);
    for idx_mod = 1:num_mod
        for idx_stf = 1:num_stf
            temp = fitlm(itd_meas_max(:,idx_stf),itd_mod_max(:,idx_stf,idx_mod),'linear');
            reg_coef(idx_stf,:,idx_mod) = temp.Coefficients.Estimate;
        end
    end

    % plot settings
    ylabel_meas = 'ITD$_{\textnormal{meas}}$ ($\mu s$)';

    xlabel_mod = 'Azimuth (deg)';
    ylabel_mod = {'ITD$_{\textnormal{Model\,1}}$ ($\mu s$)',...
        'ITD$_{\textnormal{Model\,2+}}$ ($\mu s$)',...
        'ITD$_{\textnormal{Model\,3}}$ ($\mu s)$'};

    xlabel_itd_max = 'max\{ITD$_{\textnormal{meas}}$\} ($\mu s$)';
    ylabel_itd_max = 'max\{ITD$_{\textnormal{mod}}$\} ($\mu s$)';

    ylabel_delta_itd_max = '$\Delta$max\{ITD\} ($\mu s$)';

    xlabel_delta_arg_max_phi = 'Model';
    ylabel_delta_arg_max_phi = '$\Delta$arg\,max\{ITD\} (deg)';

    shade_mod = {kv.col_shad_mod1, kv.col_shad_mod2plus, kv.col_shad_mod3};
    marker_mod = {'o','x','.'};
    marker_size = 10;

    fsize = 7.5;
    lwidth = 2;

    ymin_data = -200;
    ymax_data = 850;
    ystep_data = 200;

    xminax = 500;
    xmaxax = 850;
    yminax = xminax;
    ymaxax = xmaxax;
    xstep = 100;
    ystep = xstep;

    yminbox_ITDmax = -230;
    ymaxbox_ITDmax = 136;

    yminbox_ITDmaxphi = -20;
    ymaxbox_ITDmaxphi = 30;

    % horizontal lines denoting significant differences
    lineoff = 4.75;

    line_HRTF_ITDmaxphi_Delta12.xrange = [1 2];
    line_HRTF_ITDmaxphi_Delta12.y = 21.5;

    line_HRTF_ITDmaxphi_Delta23.xrange = [2 3];
    line_HRTF_ITDmaxphi_Delta23.y = line_HRTF_ITDmaxphi_Delta12.y+0.7*lineoff;

    line_rHARTF_ITDmaxphi_Delta23.xrange = [2 3];
    line_rHARTF_ITDmaxphi_Delta23.y = line_HRTF_ITDmaxphi_Delta12.y;

    % labels a)...d)
    xpos_anno1 = -0.45;
    ypos_anno1 = 1.05;

    % create figure
    fig_handle = figure;
    set(gcf,'Position',[281,208,427,788])

    for idx_tile = 1:18
         subplot(6,3,idx_tile)

         % itd_meas
        if ismember(idx_tile,1:3)
            plot(phi_vec,itd_meas(:,:,idx_tile),'color',kv.col_meas)
            if idx_tile==1
                text(xpos_anno1,ypos_anno1,'$\textbf{a)}$',...
                    'interpreter','latex','fontsize',fsize,'units','normalized')
                ylabel(ylabel_meas,'fontsize',fsize,'interpreter','latex')
            end
            grid on
            hold on
            title(title_stf{idx_tile},'interpreter','latex')
            set(gca,...
                'xlim',[phi_vec(1) phi_vec(end)],...
                'ylim',[ymin_data-50,ymax_data],...
                'xtick',0:30:180,...
                'ytick',ymin_data:ystep_data:ymax_data,...
                'ticklabelinterpreter','latex',...
                'fontsize',fsize)
            axis square
            box on

        % itd_mod1, itd_mod2plus, ITD deviations
        elseif ismember(idx_tile,4:6)

            fillarea1 = [delta_itd_mean(:,idx_tile-3,1)+delta_itd_std(:,idx_tile-3,1); ...
                flipud(delta_itd_mean(:,idx_tile-3,1)-delta_itd_std(:,idx_tile-3,1))];
            fill([phi_vec; flipud(phi_vec)], fillarea1, kv.col_mod1+0.4, 'edgecolor','none');
            hold on
            fillarea2 = [delta_itd_mean(:,idx_tile-3,2)+delta_itd_std(:,idx_tile-3,2); ...
                flipud(delta_itd_mean(:,idx_tile-3,2)-delta_itd_std(:,idx_tile-3,2))];
            fill([phi_vec; flipud(phi_vec)], fillarea2, kv.col_mod2plus+0.5, 'edgecolor','none');

            plot(phi_vec,delta_itd_mean(:,idx_tile-3,1),'color',kv.col_mod1,'linestyle','-.','linewidth',lwidth/2)
            plot(phi_vec,delta_itd_mean(:,idx_tile-3,2),'color',kv.col_mod2plus,'linestyle',':','linewidth',lwidth/2)
            plot(phi_vec,itd_mod(:,:,idx_tile-3,1),'color',kv.col_mod1)
            plot(phi_vec,itd_mod(:,:,idx_tile-3,2),'color',kv.col_mod2plus)

            if idx_tile==4
                text(xpos_anno1,ypos_anno1,'$\textbf{b)}$',...
                    'interpreter','latex','fontsize',fsize,'units','normalized')
                ylabel(ylabel_mod{idx_tile-2},'fontsize',fsize,'interpreter','latex')
                text(-.49,0.5,ylabel_mod{1},'fontsize',fsize,'interpreter','latex',...
                    'Rotation',90,'horizontalalignment','center',...
                    'color',kv.col_mod1,'units','normalized')
            end

            grid on
            set(gca,...
                'xlim',[phi_vec(1) phi_vec(end)],...
                'ylim',[ymin_data-50,ymax_data],...
                'xtick',0:30:180,...
                'ytick',ymin_data:ystep_data:ymax_data,...
                'ticklabelinterpreter','latex',...
                'fontsize',fsize)
            axis square
            box on

        % itd_mod3, ITD deviations
        elseif ismember(idx_tile,7:9)

            fillarea3 = [delta_itd_mean(:,idx_tile-6,3)+delta_itd_std(:,idx_tile-6,3); ...
                flipud(delta_itd_mean(:,idx_tile-6,3)-delta_itd_std(:,idx_tile-6,3))];
            fill([phi_vec; flipud(phi_vec)], fillarea3, kv.col_shad_mod3, 'edgecolor','none');
            hold on
            plot(phi_vec,delta_itd_mean(:,idx_tile-6,3),'color',kv.col_mod3,'linestyle','-.','linewidth',lwidth/2)
            plot(phi_vec,itd_mod(:,:,idx_tile-6,3),'color',kv.col_mod3)

            if idx_tile==7
                text(xpos_anno1,ypos_anno1,'$\textbf{c)}$',...
                    'interpreter','latex','fontsize',fsize,'units','normalized')
                ylabel(ylabel_mod{idx_tile-4},'fontsize',fsize,'interpreter','latex')
            end

            grid on
            set(gca,...
                'xlim',[phi_vec(1) phi_vec(end)],...
                'ylim',[ymin_data-50,ymax_data],...
                'xtick',0:30:180,...
                'ytick',ymin_data:ystep_data:ymax_data,...
                'ticklabelinterpreter','latex',...
                'fontsize',fsize)
            xlabel(xlabel_mod,'fontsize',fsize,'interpreter','latex')
            axis square
            box on

        % max{ITD}, scatter plots / regression lines
        elseif ismember(idx_tile,10:12)

            warning('off','MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout')

            plot(linspace(500,xmaxax,num_part),linspace(500,xmaxax,num_part),'linewidth',0.5,'linestyle',':','color','k');
            hold on
            scatter(itd_meas_max(:,idx_tile-9),itd_mod_max(:,idx_tile-9,1),marker_size,shade_mod{1},marker_mod{1})
            scatter(itd_meas_max(:,idx_tile-9),itd_mod_max(:,idx_tile-9,2),marker_size,shade_mod{2},marker_mod{2})
            scatter(itd_meas_max(:,idx_tile-9),itd_mod_max(:,idx_tile-9,3),marker_size,shade_mod{3},marker_mod{3})

            grid minor

            plot(itd_meas_max(:,idx_tile-9),reg_coef(idx_tile-9,2,1)*itd_meas_max(:,idx_tile-9) + ...
                reg_coef(idx_tile-9,1,1),'linewidth',lwidth,'color',kv.col_mod1);
            plot(itd_meas_max(:,idx_tile-9),reg_coef(idx_tile-9,2,2)*itd_meas_max(:,idx_tile-9) + ...
                reg_coef(idx_tile-9,1,2),'linewidth',lwidth,'color',kv.col_mod2plus);
            plot(itd_meas_max(:,idx_tile-9),reg_coef(idx_tile-9,2,3)*itd_meas_max(:,idx_tile-9) + ...
                reg_coef(idx_tile-9,1,3),'linewidth',lwidth,'color',kv.col_mod3);

            if idx_tile==10
                text(xpos_anno1,ypos_anno1,'$\textbf{d)}$',...
                    'interpreter','latex','fontsize',fsize,'units','normalized')
                ylabel(ylabel_itd_max,'fontsize',fsize,'interpreter','latex')
            end

            set(gca,'ticklabelinterpreter','latex','fontsize',fsize,...
                'xlim',[xminax xmaxax],'ylim',[yminax ymaxax],...
                'xtick',xminax:xstep:xmaxax,'ytick',yminax:ystep:ymaxax)
            xtickangle(0)
            box on
            axis square
            xlabel(xlabel_itd_max,'fontsize',fsize,'interpreter','latex')
            hA=get(gca);
            hA.XAxis.MinorTickValues = xminax:(xstep/2):xmaxax;
            hA.XAxis.MinorTick = 'on';
            hA.YAxis.MinorTickValues = yminax:(ystep/2):ymaxax;
            hA.YAxis.MinorTick = 'on';
            hA.MinorGrid = 'on';

        % Delta max{ITD}
        elseif ismember(idx_tile,13:15)

            boxplot([itd_mod_max(:,idx_tile-12,1)-itd_meas_max(:,idx_tile-12), ...
                itd_mod_max(:,idx_tile-12,2)-itd_meas_max(:,idx_tile-12), ...
                itd_mod_max(:,idx_tile-12,3)-itd_meas_max(:,idx_tile-12)],...
                'Notch','on','Labels',{'1','2+','3'},...
                'colors',[kv.col_mod1; kv.col_mod2plus; kv.col_mod3])

            set(findobj(gca,'type','line'),'linew',1.25)
            h=findobj(gca,'tag','Outliers');
%            set(h(3),'MarkerEdgeColor',kv.col_mod1);
%            set(h(2),'MarkerEdgeColor',kv.col_mod2plus);
%            set(h(1),'MarkerEdgeColor',kv.col_mod3);

            if idx_tile==13
                text(xpos_anno1,ypos_anno1,'$\textbf{e)}$',...
                    'interpreter','latex','fontsize',fsize,'units','normalized')
                ylabel(ylabel_delta_itd_max,'fontsize',fsize,'interpreter','latex')
            end

            grid on
            box on
            set(gca,'ticklabelinterpreter','latex',...
                'fontsize',fsize,...
                'yminorgrid','on',...
                'ylim',[yminbox_ITDmax ymaxbox_ITDmax],...
                'ytick',-200:50:250)
            axis square

        % Delta arg_max_phi{ITD}
        else

            boxplot([itd_mod_arg_max_phi(:,idx_tile-15,1)-itd_meas_arg_max_phi(:,idx_tile-15), ...
                itd_mod_arg_max_phi(:,idx_tile-15,2)-itd_meas_arg_max_phi(:,idx_tile-15), ...
                itd_mod_arg_max_phi(:,idx_tile-15,3)-itd_meas_arg_max_phi(:,idx_tile-15)],...
                'Notch','on','Labels',{'1','2+','3'},...
                'colors',[kv.col_mod1; kv.col_mod2plus; kv.col_mod3])

            set(findobj(gca,'type','line'),'linew',1.25)
            h=findobj(gca,'tag','Outliers');
            set(h(3),'MarkerEdgeColor',kv.col_mod1);
            set(h(2),'MarkerEdgeColor',kv.col_mod2plus);
            set(h(1),'MarkerEdgeColor',kv.col_mod3);

            if idx_tile==16
                text(xpos_anno1,ypos_anno1,'$\textbf{f)}$',...
                    'interpreter','latex','fontsize',fsize,'units','normalized')
                ylabel(ylabel_delta_arg_max_phi,'fontsize',fsize,'interpreter','latex')

                line(line_HRTF_ITDmaxphi_Delta12.xrange,...
                    [line_HRTF_ITDmaxphi_Delta12.y line_HRTF_ITDmaxphi_Delta12.y],...
                    'color','k','linewidth',lwidth/2)
                text(mean(line_HRTF_ITDmaxphi_Delta12.xrange),...
                    line_HRTF_ITDmaxphi_Delta12.y+2.5,'n.s.',...
                    'interpreter','latex',...
                    'fontsize',fsize,...
                    'horizontalalignment','center')

                line(line_HRTF_ITDmaxphi_Delta23.xrange,...
                    [line_HRTF_ITDmaxphi_Delta23.y line_HRTF_ITDmaxphi_Delta23.y],...
                    'color','k','linewidth',lwidth/2)
                text(mean(line_HRTF_ITDmaxphi_Delta23.xrange),...
                    line_HRTF_ITDmaxphi_Delta23.y+2.5,'n.s.',...
                    'interpreter','latex',...
                    'fontsize',fsize,...
                    'horizontalalignment','center')

            elseif idx_tile==18

                line(line_rHARTF_ITDmaxphi_Delta23.xrange,...
                    [line_rHARTF_ITDmaxphi_Delta23.y line_rHARTF_ITDmaxphi_Delta23.y],...
                    'color','k','linewidth',lwidth/2)
                text(mean(line_rHARTF_ITDmaxphi_Delta23.xrange),...
                    line_rHARTF_ITDmaxphi_Delta23.y+2.5,'n.s.',...
                    'interpreter','latex',...
                    'fontsize',fsize,...
                    'horizontalalignment','center')

            end

            xlabel(xlabel_delta_arg_max_phi,'fontsize',fsize,'interpreter','latex')

            grid on
            box on
            set(gca,'ticklabelinterpreter','latex',...
                'fontsize',fsize,...
                'ylim',[yminbox_ITDmaxphi ymaxbox_ITDmaxphi],...
                'ytick',-40:10:40,...
                'yminorgrid','on')
            axis square

        end

    end

end

%% Figure 10
% if flags.do_fig10
%     amt_disp([mfilename,': Not yet implemented.'])
%     fig_handle = [];
% end

%% Figure 11
% if flags.do_fig11
%     amt_disp([mfilename,': Not yet implemented.'])
%     fig_handle = [];
% end

%% no_fig
% if flags.do_no_fig
%     amt_disp([mfilename,': No figure selected for plotting.'])
%     fig_handle = [];
% end




