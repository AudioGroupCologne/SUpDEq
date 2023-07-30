function [FIG,RESULTS_struct] = plot_bischof2023(RESULTS_struct, flags)
%PLOT_BISCHOF2023 plot results from Bischof et al. 2023
%   Usage: [FIG,RESULTS_struct] = plot_bischof2023(RESULTS_struct, flags)
%
%
%   Input parameters:
%        RESULTS_struct    : structure with all predicted SNRs, BMLDs and
%                            better ear SNRs for the experiment reported in
%                            Bischof et al. (2023)
%
%   Output parameters:
%       FIG             :  structure with all figure handles of the
%                          reproduced Figures
%       RESULTS_struct  :  structure with all predicted SNRs, BMLDs and
%                          better ear SNRs for the experiment reported in
%                          Bischof et al. (2023)
%
%
%   PLOT_BISCHOF2023
%   reproduces all figures of the paper related to the model and the
%   underlying experiment (Figure 3, Figure 4, Figure 7 & Figure 8) based on
%   the model predictions from bischof2023 and returns a structure containing
%   all figure handles and the results structure containing the model
%   predictions as well as the experimental data for further analysis
%
%
%   See also: bischof2023_filterbank data_bischof2023 exp_bischof2023
%
%
%   References:
%     N. Bischof, P. Aublin, and B. Seeber. Fast processing models effects of
%     reflections on binaural unmasking. Acta Acustica, 2023.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_bischof2023.php


%   #Author: Norbert F. Bischof (2023)
%   #Author: Pierre G. Aublin
%   #Author: Bernhard Seeber (2023)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.

%% run exp_bischof2023 if OUT_struct is not yet in workspace
if nargin < 2
    %do all
    flags.do_fig3 = 1;
    flags.do_fig4 = 1;
    flags.do_fig7 = 1;
    flags.do_fig8 = 1;
end

if nargin < 1
    RESULTS_struct = exp_bischof2023(0);
end
expdata = RESULTS_struct.expdata;
RESULTS_struct = rmfield(RESULTS_struct,'expdata');

reflection_win = expdata.bischof2023_reflection_win;
L_interf = expdata.bischof2023_interf_level;
%% %%%%%%%%%%%%% reproduce Figure 3 for Bischof et al. 2023 %%%%%%%%%%%%%%%%
if flags.do_fig3
    FontName = 'Arial'; FontSize = 10;
    LineWidth = 2; MarkerSize = 6;

    % Sizes in cm for full figure, including axes etc
    fig_x = 0; fig_y = 0;
    fig_w = 8;
    fig_h = 8;
    e = 2e-3;   % Epsilon
    % Define sizes for plots and spaces around (more or less arbitrary, since
    % the sum of all spaces  needed will be normalized -> Only the ratios
    % between values are important here.
    x_left = 0.1; x_plot = 0.5; x_right = 0.01;
    y_top = 0.05; y_plot = 1; y_bottom = 0.17;

    % Normalize to 1 - adapt this to the plot configurations you want
    x_tot = x_left + x_plot + x_right;
    x_left = x_left/x_tot; x_plot = x_plot/x_tot; x_right = x_right/x_tot;
    y_tot = y_top + y_plot + y_bottom;
    y_top = y_top/y_tot; y_plot = y_plot/y_tot; y_bottom = y_bottom/y_tot;

    FIG.fig3 = figure('NumberTitle','off','Name','Bischof et al. (2023) Figure 3','MenuBar','none','Toolbar','none',...
                      'Units','centimeters','InnerPosition',[fig_x, fig_y, fig_w, fig_h],'PaperSize',[fig_w,fig_h]);
    ax_big = axes('Units','normalized','Position',[0-e 0-e 1+e 1+e],'Visible','on','XTick',[],'YTick',[],'Color','none'); % Big axis for overall x and y labels
    ax{1} = axes('Units','normalized','Position',[x_left y_bottom x_plot y_plot]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_alpha05_S0 = [14.4,19.2,43.22,72,144,240,480];
    x_alpha01_S0 = [14.8,19.73,44.4,74,148,246,493];
    x_alpha05_S60 = [15.2,20.3,45.6,76,152,253,507];
    x_alpha01_S60 = [15.6,20.8,46.8,78,156,260,520];

    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[L_interf L_interf],...
                            '--','Color',[0.8,0.8,0.8],'LineWidth',LineWidth);
    hold on;
    plot_early_alpha01_S0 = errorbar(x_alpha01_S0,...
                            expdata.bischof2023_exp_MEDIAN(:,2),...
                            expdata.bischof2023_exp_QUART25(:,2),...
                            expdata.bischof2023_exp_QUART75(:,2),...
                            '-^','Color',[0.8,0,0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
    plot_early_alpha01_S60 = errorbar(x_alpha01_S60,...
                            expdata.bischof2023_exp_MEDIAN(:,4),...
                            expdata.bischof2023_exp_QUART25(:,4),...
                            expdata.bischof2023_exp_QUART75(:,4),...
                            ':^','Color',[0.8,0,0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);

    plot_early_alpha05_S0 = errorbar(x_alpha05_S0,...
                            expdata.bischof2023_exp_MEDIAN(:,1),...
                            expdata.bischof2023_exp_QUART25(:,1),...
                            expdata.bischof2023_exp_QUART75(:,1),...
                            '-o','Color',[0,0,0.8],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
    plot_early_alpha05_S60 = errorbar(x_alpha05_S60,...
                            expdata.bischof2023_exp_MEDIAN(:,3),...
                            expdata.bischof2023_exp_QUART25(:,3),...
                            expdata.bischof2023_exp_QUART75(:,3),...
                            ':o','Color',[0,0,0.8],'LineWidth',LineWidth,'MarkerSize',MarkerSize);

    leg = legend([plot_early_alpha05_S0,plot_early_alpha01_S0,plot_early_alpha05_S60,plot_early_alpha01_S60],...
        {'\alpha_{0.5}S_{0}','\alpha_{0.1}S_{0}','\alpha_{0.5}S_{60}','\alpha_{0.1}S_{60}'},...
        'Position',[0.672525119941162,0.677433997675249,0.304820420514755,0.27528356263163]);
    grid on;
    xlabel('RIR truncation time (ms)','FontName',FontName,'FontSize',FontSize);
    ylabel('Threshold of HCT in noise (dB SPL)','FontName',FontName,'FontSize',FontSize);


    set(gca,'FontSize',FontSize,'FontName',FontName);
    set(gca,'XLim',[12 600],'XTick',reflection_win,'XTickLabel',reflection_win);
    set(gca,'YLim',[35,61],'YTick',35:5:60,'YTickLabel',35:5:60);
    set(gca,'XScale','log');

end
%% %%%%%%%%%%%% reproduce Figure 4 for Bischof et al. 2023 %%%%%%%%%%%%%%%%
if flags.do_fig4
    FontName = 'Arial'; FontSize = 10;
    LineWidth = 2; MarkerSize = 6;

    % Sizes in cm for full figure, including axes etc
    fig_x = 0; fig_y = 0;
    fig_w = 8;
    fig_h = 8;
    e = 2e-3;   % Epsilon
    % Define sizes for plots and spaces around (more or less arbitrary, since
    % the sum of all spaces  needed will be normalized -> Only the ratios
    % between values are important here.
    x_left = 0.1; x_plot = 0.5; x_right = 0.01;
    y_top = 0.05; y_plot = 1; y_bottom = 0.17;

    % Normalize to 1 - adapt this to the plot configurations you want
    x_tot = x_left + x_plot + x_right;
    x_left = x_left/x_tot; x_plot = x_plot/x_tot; x_right = x_right/x_tot;
    y_tot = y_top + y_plot + y_bottom;
    y_top = y_top/y_tot; y_plot = y_plot/y_tot; y_bottom = y_bottom/y_tot;

    FIG.fig4 = figure('NumberTitle','off','Name','Bischof et al. (2023) Figure 4','MenuBar','none','Toolbar','none',...
                      'Units','centimeters','InnerPosition',[fig_x, fig_y, fig_w, fig_h],'PaperSize',[fig_w,fig_h]);
    ax_big = axes('Units','normalized','Position',[0-e 0-e 1+e 1+e],'Visible','on','XTick',[],'YTick',[],'Color','none'); % Big axis for overall x and y labels
    ax{1} = axes('Units','normalized','Position',[x_left y_bottom x_plot y_plot]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_alpha05_S0 = [14.8,19.7,44.4,74,148,246,493];
    x_alpha01_S0 = [15.2,20.3,45.6,76,152,253,506];

    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[L_interf L_interf],...
                            '--','Color',[0.8,0.8,0.8],'LineWidth',LineWidth);
    hold on;
    plot_late_alpha01_S0 = errorbar(x_alpha01_S0,...
                            expdata.bischof2023_exp_MEDIAN(:,6),...
                            expdata.bischof2023_exp_QUART25(:,6),...
                            expdata.bischof2023_exp_QUART75(:,6),...
                            '-^','Color',[0.8,0,0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
    plot_late_alpha05_S0 = errorbar(x_alpha05_S0,...
                            expdata.bischof2023_exp_MEDIAN(:,5),...
                            expdata.bischof2023_exp_QUART25(:,5),...
                            expdata.bischof2023_exp_QUART75(:,5),...
                            '-o','Color',[0,0,0.8],'LineWidth',LineWidth,'MarkerSize',MarkerSize);

    leg = legend([plot_late_alpha05_S0,plot_late_alpha01_S0],...
        {'\alpha_{0.5}S_{0}','\alpha_{0.1}S_{0}'},'Position',[0.17057937774746,0.809048274050716,0.28780718742111,0.143160381728427]);
    grid on;
    xlabel('RIR cutting time (ms)','FontName',FontName,'FontSize',FontSize);
    ylabel('Threshold of HCT in noise (dB SPL)','FontName',FontName,'FontSize',FontSize);


    set(gca,'FontSize',FontSize,'FontName',FontName);
    set(gca,'XLim',[12 600],'XTick',reflection_win,'XTickLabel',reflection_win);
    set(gca,'YLim',[35,61],'YTick',35:5:60,'YTickLabel',35:5:60);
    set(gca,'XScale','log');    
end
%% %%%%%%%%%%%% reproduce Figure 7 for Bischof et al. 2023 %%%%%%%%%%%%%%%%
if flags.do_fig7
    FontName = 'Arial'; FontSize = 10;
    LineWidth = 2;

    % Sizes in cm for full figure, including axes etc
    fig_x = 0; fig_y = 0;
    fig_w = 16;
    fig_h = 22;
    e = 2e-3;   % Epsilon
    % Define sizes for plots and spaces around (more or less arbitrary, since
    % the sum of all spaces  needed will be normalized -> Only the ratios
    % between values are important here.
    x_left = 0.1; x_space = 0.004; x_plot = 0.5; x_right = 0.01;
    y_top = 0.1; y_space = 0.01; y_space2 = 0.25; y_plot = 1; y_bottom = 0.25;

    % Normalize to 1 - adapt this to the plot configurations you want
    x_tot = x_left + x_plot + x_space + x_plot + x_right;
    x_left = x_left/x_tot; x_space = x_space/x_tot; x_plot = x_plot/x_tot; x_right = x_right/x_tot;
    y_tot = y_top + y_plot + y_space + y_plot + y_space2 + y_plot + y_bottom;
    y_top = y_top/y_tot; y_space = y_space/y_tot; y_space2 = y_space2/y_tot; y_plot = y_plot/y_tot; y_bottom = y_bottom/y_tot;

    FIG.fig8 = figure('NumberTitle','off','Name','Bischof et al. (2023) Figure 7','MenuBar','none','Toolbar','none',...
                      'Units','centimeters','InnerPosition',[fig_x, fig_y, fig_w, fig_h],'PaperSize',[fig_w,fig_h]);
    ax_big = axes('Units','normalized','Position',[0-e 0-e 1+e 1+e],'Visible','on','XTick',[],'YTick',[],'Color','none'); % Big axis for overall x and y labels
    ax{1} = axes('Units','normalized','Position',[x_left                y_bottom+(y_plot+y_space2+y_plot+y_space) x_plot y_plot]);
    ax{2} = axes('Units','normalized','Position',[x_left+x_plot+x_space y_bottom+(y_plot+y_space2+y_plot+y_space) x_plot y_plot]);
    ax{3} = axes('Units','normalized','Position',[x_left                y_bottom+(y_plot+y_space2) x_plot y_plot]);
    ax{4} = axes('Units','normalized','Position',[x_left+x_plot+x_space y_bottom+y_plot+y_space2 x_plot y_plot]);
    ax{5} = axes('Units','normalized','Position',[x_left                y_bottom                    x_plot y_plot]);
    ax{6} = axes('Units','normalized','Position',[x_left+x_plot+x_space y_bottom                    x_plot y_plot]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_data_fast = RESULTS_struct.pred_SNR_fast_bischof2023;
    set(FIG.fig8,'CurrentAxes',ax{1})
    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[L_interf L_interf],...
        '--','Color',[0.8,0.8,0.8],'LineWidth',LineWidth,'Parent',ax{1});
    hold on;
    plot_bmld_early_alpha01_S0 = errorbar(reflection_win,...
        expdata.bischof2023_exp_MEDIAN(:,2),...
        expdata.bischof2023_exp_QUART25(:,2),...
        expdata.bischof2023_exp_QUART75(:,2),...
        '-^','Color',[0.8,0,0],'LineWidth',LineWidth);
    FAST = plot(reflection_win,...
        L_interf-plot_data_fast(:,2),...
        '-.s','Color',[0,142,44]./255,'LineWidth',LineWidth);
    set(ax{1},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',[],...
        'YLim',[30 64],'YTick',30:5:60,'YTickLabel',30:5:60,'FontSize',FontSize,'FontName',FontName);
    grid on;
    legend([plot_bmld_early_alpha01_S0,FAST],{'Exp. data','DynBU_{fast}'},'Location','southwest');
    text(410,32,'a','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')
    title('Room absorption coefficient \alpha_{0.1}');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{2})
    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[L_interf L_interf],...
        '--','Color',[0.8,0.8,0.8],'LineWidth',LineWidth);
    hold on;
    plot_bmld_early_alpha01_S0 = errorbar(reflection_win,...
        expdata.bischof2023_exp_MEDIAN(:,1),...
        expdata.bischof2023_exp_QUART25(:,1),...
        expdata.bischof2023_exp_QUART75(:,1),...
        '-o','Color',[0,0,0.8],'LineWidth',LineWidth);
    plot_est_early_alpha05_S0_fast = plot(reflection_win,...
        L_interf-plot_data_fast(:,1),...
        '-.s','Color',[0,142,44]./255,'LineWidth',LineWidth);
    set(ax{2},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',[],...
        'YLim',[30 64],'YTick',30:5:60,'YTickLabel',[],'FontSize',FontSize,'FontName',FontName);
    grid on;
    % title('BMLD \alpha_{0.5} S_{0} early');
    legend([plot_bmld_early_alpha01_S0,plot_est_early_alpha05_S0_fast],{'Exp. data','DynBU_{fast}'},'Location','southwest');
    text(410,32,'b','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')
    title('Room absorption coefficient \alpha_{0.5}');
    %colororder({'k','k'})
    text(15,63,'Target at 0° early reflections','HorizontalAlignment','center','VerticalAlignment','top','FontSize',10,'FontName','Arial','FontWeight','bold','BackgroundColor','w');
    set(ax{2},'YLim',[30 64],'YTick',30:5:60,'YTickLabel',[],'FontSize',FontSize,'FontName',FontName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{3})
    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[L_interf L_interf],...
        '--','Color',[0.8,0.8,0.8],'LineWidth',LineWidth);
    hold on;
    plot_bmld_early_alpha01_S0 = errorbar(reflection_win,...
        expdata.bischof2023_exp_MEDIAN(:,4),...
        expdata.bischof2023_exp_QUART25(:,4),...
        expdata.bischof2023_exp_QUART75(:,4),...
        '-^','Color',[0.8,0,0],'LineWidth',LineWidth);
    plot_est_early_alpha01_S0_fast = plot(reflection_win,...
        L_interf-plot_data_fast(:,4),...
        '-.s','Color',[0,142,44]./255,'LineWidth',LineWidth);
    set(ax{3},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',reflection_win,...
        'YLim',[30 64],'YTick',30:5:60,'YTickLabel',30:5:60,'FontSize',FontSize,'FontName',FontName);
    grid on;
    text(410,32,'c','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{4})
    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[L_interf L_interf],...
        '--','Color',[0.8,0.8,0.8],'LineWidth',LineWidth);
    hold on;
    plot_bmld_early_alpha01_S0 = errorbar(reflection_win,...
        expdata.bischof2023_exp_MEDIAN(:,3),...
        expdata.bischof2023_exp_QUART25(:,3),...
        expdata.bischof2023_exp_QUART75(:,3),...
        '-o','Color',[0,0,0.8],'LineWidth',LineWidth);
    plot_est_early_alpha05_S0_fast = plot(reflection_win,...
        L_interf-plot_data_fast(:,3),...
        '-.s','Color',[0,142,44]./255,'LineWidth',LineWidth);
    set(ax{4},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',reflection_win,...
        'YLim',[30 64],'YTick',30:5:60,'YTickLabel',[],'FontSize',FontSize,'FontName',FontName);
    grid on;
    text(410,32,'d','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')

    text(15,63,'Target at 60° early reflections','HorizontalAlignment','center','VerticalAlignment','top','FontSize',10,'FontName','Arial','FontWeight','bold','BackgroundColor','w');
    set(ax{4},'YLim',[30 64],'YTick',30:5:60,'YTickLabel',[],'FontSize',FontSize,'FontName',FontName);

    text(12,26,'RIR truncation time (ms)','HorizontalAlignment','center','VerticalAlignment','top','FontSize',11,'FontName','Arial','FontWeight','normal','BackgroundColor','w');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{5})
    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[L_interf L_interf],...
        '--','Color',[0.8,0.8,0.8],'LineWidth',LineWidth);
    hold on;
    plot_bmld_early_alpha01_S0 = errorbar(reflection_win,...
        expdata.bischof2023_exp_MEDIAN(:,6),...
        expdata.bischof2023_exp_QUART25(:,6),...
        expdata.bischof2023_exp_QUART75(:,6),...
        '-^','Color',[0.8,0,0],'LineWidth',LineWidth);
    plot_est_early_alpha05_S0_fast = plot(reflection_win,...
        L_interf-plot_data_fast(:,6),...
        '-.s','Color',[0,142,44]./255,'LineWidth',LineWidth);
    set(ax{5},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',reflection_win,...
        'YLim',[30 64],'YTick',30:5:60,'YTickLabel',30:5:60,'FontSize',FontSize,'FontName',FontName);
    grid on;
    text(410,32,'e','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{6})
    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[60 60],...
        '--','Color',[0.8,0.8,0.8],'LineWidth',LineWidth);
    hold on;
    plot_bmld_early_alpha01_S0 = errorbar(reflection_win,...
        expdata.bischof2023_exp_MEDIAN(:,5),...
        expdata.bischof2023_exp_QUART25(:,5),...
        expdata.bischof2023_exp_QUART75(:,5),...
        '-o','Color',[0,0,0.8],'LineWidth',LineWidth);
    plot_est_early_alpha05_S0_fast = plot(reflection_win,...
        L_interf-plot_data_fast(:,5),...
        '-.s','Color',[0,142,44]./255,'LineWidth',LineWidth);
    set(ax{6},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',{'',num2str(reflection_win(2)),num2str(reflection_win(3)),num2str(reflection_win(4)),num2str(reflection_win(5)),num2str(reflection_win(6)),num2str(reflection_win(7))},...
        'YLim',[30 64],'YTick',30:5:60,'YTickLabel',[],'FontSize',FontSize,'FontName',FontName);
    grid on;
    text(410,32,'f','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')
    text(15,63,'Target at 0° late reflections','HorizontalAlignment','center','VerticalAlignment','top','FontSize',10,'FontName','Arial','FontWeight','bold','BackgroundColor','w');
    set(ax{6},'YLim',[30 64],'YTick',30:5:60,'YTickLabel',[],'FontSize',FontSize,'FontName',FontName);

    text(12,26,'RIR cutting time (ms)','HorizontalAlignment','center','VerticalAlignment','top','FontSize',11,'FontName','Arial','FontWeight','normal','BackgroundColor','w');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax_big)
    % Do the labels
    ax_big.XLabel.String = [];
    ax_big.XLabel.HorizontalAlignment = 'center'; ax_big.XLabel.VerticalAlignment = 'bottom';
    ax_big.XLabel.Position(2) = 0.01;  % Flush with the lower bound of the figure
    ax_big.YLabel.String = 'Detection thresholds (dB SPL)';
    ax_big.YLabel.HorizontalAlignment = 'center'; ax_big.YLabel.VerticalAlignment = 'top';
    ax_big.YLabel.Position(1) = -0.0001; % Flush with the left bound of the figure
    ax_big.FontSize = 10;
    ax_big.FontName = FontName;
end
    %% %%%%%%%%%%%%% reproduce Figure 8 for Bischof et al. 2023 %%%%%%%%%%%%%%%%
if flags.do_fig8
    pred_BMLD_fast_0 = max(RESULTS_struct.pred_BMLD_fast_bischof2023,0);
    %% plot Figure 7 Bischof et al. (2022)
    FontName = 'Arial'; FontSize = 10;
    LineWidth = 2;

    % Sizes in cm for full figure, including axes etc
    fig_x = 0; fig_y = 0;
    fig_w = 16;
    fig_h = 12;
    e = 2e-3;   % Epsilon
    % Define sizes for plots and spaces around (more or less arbitrary, since
    % the sum of all spaces  needed will be normalized -> Only the ratios
    % between values are important here.
    x_left = 0.1; x_space = 0.004; x_plot = 0.5; x_right = 0.01;
    y_top = 0.2; y_space = 0.01; y_space2 = 0.35; y_plot = 1; y_bottom = 0.35;

    % Normalize to 1 - adapt this to the plot configurations you want
    x_tot = x_left + x_plot + x_space + x_plot + x_right;
    x_left = x_left/x_tot; x_space = x_space/x_tot; x_plot = x_plot/x_tot; x_right = x_right/x_tot;
    y_tot = y_top + y_plot + y_space + y_plot + y_space2 + y_plot + y_bottom;
    y_top = y_top/y_tot; y_space = y_space/y_tot; y_space2 = y_space2/y_tot; y_plot = y_plot/y_tot; y_bottom = y_bottom/y_tot;

    FIG.fig8 = figure('NumberTitle','off','Name','Bischof et al. (2022) Figure 8','MenuBar','none','Toolbar','none',...
                      'Units','centimeters','InnerPosition',[fig_x, fig_y, fig_w, fig_h],'PaperSize',[fig_w,fig_h]);
    ax_big = axes('Units','normalized','Position',[0-e 0-e 1+e 1+e],'Visible','on','XTick',[],'YTick',[],'Color','none'); % Big axis for overall x and y labels
    ax{1} = axes('Units','normalized','Position',[x_left                y_bottom+(y_plot+y_space2+y_plot+y_space) x_plot y_plot]);
    ax{2} = axes('Units','normalized','Position',[x_left+x_plot+x_space y_bottom+(y_plot+y_space2+y_plot+y_space) x_plot y_plot]);
    ax{3} = axes('Units','normalized','Position',[x_left                y_bottom+(y_plot+y_space2) x_plot y_plot]);
    ax{4} = axes('Units','normalized','Position',[x_left+x_plot+x_space y_bottom+y_plot+y_space2 x_plot y_plot]);
    ax{5} = axes('Units','normalized','Position',[x_left                y_bottom                    x_plot y_plot]);
    ax{6} = axes('Units','normalized','Position',[x_left+x_plot+x_space y_bottom                    x_plot y_plot]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{1})
    plot([reflection_win(1) reflection_win(end)],[0 0],...
        '--','Color','none','LineWidth',LineWidth,'Parent',ax{1});
    hold on;

    area_BetterEar_fast = fill([reflection_win(1) reflection_win(end) fliplr(reflection_win)],...
        [0 0 fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,2)')],...
        [0,142,44]./255);
        area_BetterEar_fast.EdgeColor = 'none'; area_BetterEar_fast.FaceAlpha = 0.8;

    area_BMLD_fast = fill([reflection_win fliplr(reflection_win)],...
        [RESULTS_struct.pred_BE_fast_bischof2023(:,2)' fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,2)'+ pred_BMLD_fast_0(:,2)')],...
        [0,142,44]./255);
        area_BMLD_fast.EdgeColor = 'none'; area_BMLD_fast.FaceAlpha = 0.3;

    FAST = plot(reflection_win,...
        RESULTS_struct.pred_SNR_fast_bischof2023(:,2),...
        '-s','Color',[0,0,0]./255,'LineWidth',LineWidth);

    set(ax{1},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',[],...
        'YLim',[0 40],'YTick',0:5:35,'YTickLabel',0:5:35,'FontSize',FontSize,'FontName',FontName);
    grid on;
    legend([FAST,area_BMLD_fast,area_BetterEar_fast],{'DynBU_{fast}','BMLD contribution','better-ear contribution'},'Position',[0.128344418245923,0.831434920760292,0.275283560378147,0.116813605618417]);
    text(13.5,32.5,'a','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')
    title('Room absorption coefficient \alpha_{0.1}','Position',[84.85335524290923,39.7,0]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{2})
    plot([reflection_win(1) reflection_win(end)],[0 0],...
        '--','Color','none','LineWidth',LineWidth);
    hold on;

    area_BetterEar_fast = fill([reflection_win(1) reflection_win(end) fliplr(reflection_win)],...
        [0 0 fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,1)')],...
        [0,142,44]./255);
        area_BetterEar_fast.EdgeColor = 'none'; area_BetterEar_fast.FaceAlpha = 0.8;

    area_BMLD_fast = fill([reflection_win fliplr(reflection_win)],...
        [RESULTS_struct.pred_BE_fast_bischof2023(:,1)' fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,1)'+pred_BMLD_fast_0(:,1)')],...
        [0,142,44]./255);
        area_BMLD_fast.EdgeColor = 'none'; area_BMLD_fast.FaceAlpha = 0.3;

    FAST = plot(reflection_win,...
        RESULTS_struct.pred_SNR_fast_bischof2023(:,1),...
        '-s','Color',[0,0,0]./255,'LineWidth',LineWidth);

    set(ax{2},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',[],...
        'YLim',[0 40],'YTick',0:5:35,'YTickLabel',[],'FontSize',FontSize,'FontName',FontName);
    grid on;
    % title('BMLD \alpha_{0.5} S_{0} early');
    text(410,32.5,'b','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')
    title('Room absorption coefficient \alpha_{0.5}','Position',[84.85335524290923,39.7,0]);
%    colororder({'k','k'})
    text(15,37,'Target at 0° early reflections','HorizontalAlignment','center','VerticalAlignment','top','FontSize',10,'FontName','Arial','FontWeight','bold','BackgroundColor','w');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{3})
    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[0 0],...
        '--','Color','none','LineWidth',LineWidth);
    hold on;

    area_BetterEar_fast = fill([reflection_win(1) reflection_win(end) fliplr(reflection_win)],...
        [0 0 fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,4)')],...
        [0,142,44]./255);
        area_BetterEar_fast.EdgeColor = 'none'; area_BetterEar_fast.FaceAlpha = 0.8;

    area_BMLD_fast = fill([reflection_win fliplr(reflection_win)],...
        [RESULTS_struct.pred_BE_fast_bischof2023(:,4)' fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,4)'+pred_BMLD_fast_0(:,4)')],...
        [0,142,44]./255);
        area_BMLD_fast.EdgeColor = 'none'; area_BMLD_fast.FaceAlpha = 0.3;

    FAST = plot(reflection_win,...
        RESULTS_struct.pred_SNR_fast_bischof2023(:,4),...
        '-s','Color',[0,0,0]./255,'LineWidth',LineWidth);

    set(ax{3},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',reflection_win,...
        'YLim',[0 40],'YTick',0:5:35,'YTickLabel',0:5:35,'FontSize',FontSize,'FontName',FontName);
    grid on;
    text(13.5,32.5,'c','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{4})
    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[0 0],...
        '--','Color','none','LineWidth',LineWidth);
    hold on;

    area_BetterEar_fast = fill([reflection_win(1) reflection_win(end) fliplr(reflection_win)],...
        [0 0 fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,3)')],...
        [0,142,44]./255);
        area_BetterEar_fast.EdgeColor = 'none'; area_BetterEar_fast.FaceAlpha = 0.8;

    area_BMLD_fast = fill([reflection_win fliplr(reflection_win)],...
        [RESULTS_struct.pred_BE_fast_bischof2023(:,3)' fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,3)'+pred_BMLD_fast_0(:,3)')],...
        [0,142,44]./255);
        area_BMLD_fast.EdgeColor = 'none'; area_BMLD_fast.FaceAlpha = 0.3;

    FAST = plot(reflection_win,...
        RESULTS_struct.pred_SNR_fast_bischof2023(:,3),...
        '-s','Color',[0,0,0]./255,'LineWidth',LineWidth);

    set(ax{4},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',reflection_win,...
        'YLim',[0 40],'YTick',0:5:35,'YTickLabel',[],'FontSize',FontSize,'FontName',FontName);
    grid on;
    text(410,32.5,'d','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')

    text(15,37,'Target at 60° early reflections','HorizontalAlignment','center','VerticalAlignment','top','FontSize',10,'FontName','Arial','FontWeight','bold','BackgroundColor','w');

    text(12,-8,'RIR truncation time (ms)','HorizontalAlignment','center','VerticalAlignment','top','FontSize',11,'FontName','Arial','FontWeight','normal','BackgroundColor','w');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{5})
    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[0 0],...
        '--','Color','none','LineWidth',LineWidth);
    hold on;

    area_BetterEar_fast = fill([reflection_win(1) reflection_win(end) fliplr(reflection_win)],...
        [0 0 fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,6)')],...
        [0,142,44]./255);
        area_BetterEar_fast.EdgeColor = 'none'; area_BetterEar_fast.FaceAlpha = 0.8;

    area_BMLD_fast = fill([reflection_win fliplr(reflection_win)],...
        [RESULTS_struct.pred_BE_fast_bischof2023(:,6)' fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,6)'+pred_BMLD_fast_0(:,6)')],...
        [0,142,44]./255);
        area_BMLD_fast.EdgeColor = 'none'; area_BMLD_fast.FaceAlpha = 0.3;

    FAST = plot(reflection_win,...
        RESULTS_struct.pred_SNR_fast_bischof2023(:,6),...
        '-s','Color',[0,0,0]./255,'LineWidth',LineWidth);

    set(ax{5},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',reflection_win,...
        'YLim',[0 40],'YTick',0:5:35,'YTickLabel',0:5:35,'FontSize',FontSize,'FontName',FontName);
    grid on;
    text(13.5,32.5,'e','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax{6})
    plot_noise_level = plot([reflection_win(1) reflection_win(end)],[0 0],...
        '--','Color','none','LineWidth',LineWidth);
    hold on;

    area_BetterEar_fast = fill([reflection_win(1) reflection_win(end) fliplr(reflection_win)],...
        [0 0 fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,5)')],...
        [0,142,44]./255);
        area_BetterEar_fast.EdgeColor = 'none'; area_BetterEar_fast.FaceAlpha = 0.8;

    area_BMLD_fast = fill([reflection_win fliplr(reflection_win)],...
        [RESULTS_struct.pred_BE_fast_bischof2023(:,5)' fliplr(RESULTS_struct.pred_BE_fast_bischof2023(:,5)'+pred_BMLD_fast_0(:,5)')],...
        [0,142,44]./255);
        area_BMLD_fast.EdgeColor = 'none'; area_BMLD_fast.FaceAlpha = 0.3;

    FAST = plot(reflection_win,...
        RESULTS_struct.pred_SNR_fast_bischof2023(:,5),...
        '-s','Color',[0,0,0]./255,'LineWidth',LineWidth);

    set(ax{6},'XScale','log',...
        'XLim',[12 600],'XTick',reflection_win,'XTickLabel',reflection_win,...
        'YLim',[0 40],'YTick',0:5:35,'YTickLabel',[],'FontSize',FontSize,'FontName',FontName);
    grid on;
    text(410,32.5,'f','HorizontalAlignment','left','Fontsize',12,'FontName','Arial','FontWeight','bold')
    text(15,37,'Target at 0° late reflections','HorizontalAlignment','center','VerticalAlignment','top','FontSize',10,'FontName','Arial','FontWeight','bold','BackgroundColor','w');

    text(12,-8,'RIR cutting time (ms)','HorizontalAlignment','center','VerticalAlignment','top','FontSize',11,'FontName','Arial','FontWeight','normal','BackgroundColor','w');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(FIG.fig8,'CurrentAxes',ax_big)
    % Do the labels
    ax_big.XLabel.String = [];
    ax_big.XLabel.HorizontalAlignment = 'center'; ax_big.XLabel.VerticalAlignment = 'bottom';
    ax_big.XLabel.Position(2) = 0.01;  % Flush with the lower bound of the figure
    ax_big.YLabel.String = 'Detection benefit (dB)';
    ax_big.YLabel.HorizontalAlignment = 'center'; ax_big.YLabel.VerticalAlignment = 'top';
    ax_big.YLabel.Position(1) = -0.0001; % Flush with the left bound of the figure
    ax_big.FontSize = 10;
    ax_big.FontName = FontName;
end
%%
RESULTS_struct.expdata = expdata;

