function plot_joergensen2011(dSRT, varargin)
%plot_joergensen2011  Plot Fig. 5 or 6 of Joergensen and Dau (2011)
%   Usage: plot_joergensen2011(dSRT,flag)
%
%   PLOT_JOERGENSEN2011(dSRT) plots the output of JOERGENSEN2011 in the
%   style of Fig. 5 or 6 of Joergensen and Dau (2011).
%
%   The flag may be one of:
%
%     'fig5' : The measured change in SRT (open squares), averaged across 6 
%              normal-hearing listeners, as a function of the reverberation 
%              time, T30. The mean SRT in the reference condition was -3 dB. 
%              Model predictions are indicated by the filled squares. The 
%              linear correlation coefficient (q) and RMSE is indicated in 
%              the upper left corner.
%     'fig6' : DSRT (left ordinate) as a function of the over-subtraction 
%              factor a for 4 normal-hearing listeners (open squares) and 
%              sEPSM predictions (filled squares). The right ordinate (with 
%              a reversed scale) shows the corresponding sSTI values as filled 
%              gray circles. These values are, however, not converted to DSRT 
%              values since these would be outside the left ordinate scale.
%
%   See also: joergensen2011, plot_joergensen2011
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/plot/plot_joergensen2011.php

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

% Define input flags
definput.flags.type = {'fig5','fig6'};

[flags,kv]  = ltfatarghelper({}, definput,varargin);

%% Plot figure 5 from Joergensen and Dau (2011)
if flags.do_fig5
    % Get measured dSRT
    [dSRTdata,SRTdata_std ] = data_joergensen2011('fig5');
    conditions = [0 0.4 0.7 1.3 2.3];
    
xmin = 0;
xmax = 6;
ymin = -1.5;
ymax = 12.5;
ytickmax = 12;
ytickmin = -1;
fnts = 14;
fig=figure;

sh = 400;
set(fig,'Position',[2*sh, 0.15*sh, 1.2*sh, 1*sh]);
 
    mark_sty = {'s','d','v','d','>','o'};
    mark_size = 10;%[6 6 6 7];
    mark_col(1,:) = [1 1 1]*0.01;
    mark_col(2,:) = [1 1 1]*0.4;
    mark_col(3,:) = [1 1 1];
    mark_col(4,:) = [1 1 1]*0.7;
    mark_col(4,:) = [1 1 1]*0.7;
  
    x = 1:5;
     
     RMSE = sqrt(mean(abs(dSRT- dSRTdata).^2));
    
     r_corr = corr(dSRT',dSRTdata', 'type', 'Pearson');
    
    offset = 0.15;
    h = plot(0,0,'color',[1 1 1]);
errorbar(x,dSRTdata,SRTdata_std,...
            'linestyle',    'none',...
            'linewidth',    1,...
            'color',           mark_col(1,:),...
            'marker',           char(mark_sty(1)),...
            'markerfacecolor',  [1 1 1],...
            'markersize',       mark_size); hold on

    idx_stim = 1;
    plot(x-offset,dSRT,...
            'linestyle',    'none',...
            'linewidth',    1,...
            'color',           mark_col(idx_stim,:),...
            'marker',           char(mark_sty(idx_stim)),...
            'markerfacecolor',  mark_col(idx_stim,:),...
            'markersize',       mark_size);hold on

    xlabel('Reverberation time T_{30} (s)','FontSize',fnts, 'FontName', 'Times');
    ylabel('\DeltaSRT (dB)','FontSize',fnts, 'FontName', 'Times');
    ylim([ymin ymax])
    xlim([xmin xmax])
    le = legend('Data','sEPSM');
    set(le,'box','off','fontsize',fnts,'location','southeast');
     text(.6,9,{[' \rho = ',num2str(r_corr,2)],['RMSE = ',num2str(RMSE,2), ' dB']},'fontsize',fnts,'FontName', 'Times');
    set(gca,'xTick',1:5 ,'xTickLabel',conditions,'ytick',ytickmin:ytickmax,'FontSize',fnts,'yticklabel',ytickmin:ytickmax);

set(gca, 'FontName', 'Times');
set(gcf, 'Color', 'w');
    
end

%% Plot figure 6 from Joergensen and Dau (2011)
if flags.do_fig6
    % Get measured dSRT
    [dSRTdata,SRTdata_std ] = data_joergensen2011('fig6');
    conditions = [0 0.5 1 2 4 8];
    xmin = 0;
    xmax = 7;
    ymin = -4.5;
    ymax = 4.5;
    ytickmax = 4;
    ytickmin = -4;
    
    fig=figure;
    
    sh = 400;
    set(fig,'Position',[2*sh, 0.15*sh, 1.2*sh, 1*sh]);
    
    
    line_sty = {'-','--','-.','-'};
    
    mark_sty = {'s','d','v','x','>','o'};
    mark_size = 10;%[6 6 6 7];
    mark_col(1,:) = [1 1 1]*0.01;
    mark_col(2,:) = [1 1 1]*0.4;
    mark_col(3,:) = [1 1 1]*0.3;
    mark_col(4,:) = [1 1 1]*0.5;
    mark_col(5,:) = [1 1 1]*0.2;
    fnts = 14;
    
    x = 1:6;
    
    offset = [0.21 -0.21 0.25 -0.25];
    RMSE = sqrt(mean(abs(dSRT- dSRTdata).^2));
    r_corr = corr(dSRT',dSRTdata', 'type', 'Pearson');
    
    h = plot(0,0,'color',[1 1 1]);
    
    
    errorbar(x,dSRTdata,SRTdata_std,...
        'linestyle',    'none',...
        'linewidth',    1,...
        'color',           mark_col(1,:),...
        'marker',           char(mark_sty(1)),...
        'markerfacecolor',  [1 1 1],...
        'markersize',       mark_size); hold on
    
    plot(x+offset(1),dSRT,...
        'linestyle',    'none',...
        'linewidth',    1,...
        'color',           mark_col(1,:),...
        'marker',           char(mark_sty(1)),...
        'markerfacecolor',  mark_col(1,:),...
        'markersize',       mark_size);hold on
    
    hold off
    
    axis([xmin xmax ymin ymax]);
    
    set(gca,'xTick',x,'fontsize',fnts);
    set(gca,'yTick',ytickmin:ytickmax);
    set(gca,'xTickLabel',conditions);
    set(gca,'yTickLabel',ytickmin:ytickmax);
    
    xlabel('Over-subtraction factor \alpha','Fontsize',16)
    ylabel('\DeltaSRT (dB)','Fontsize',fnts);
    
    text(.6,3.5,{[' \rho = ',num2str(floor(r_corr*100)/100,2)],['RMSE = ',num2str(RMSE,2), ' dB']},'fontsize',fnts,'FontName', 'Times-Roman');
    
    pos= get(gca,'Position');
    set(gca,'box','on','position', [pos(1)*.8 pos(2) pos(3)*1 pos(4)*1 ])
    
    le = legend('Data','sEPSM');
    set(le,'box','off','fontsize',14,'Location','southwest');
    
    
    
    set(gca, 'FontName', 'Times');
    set(gcf, 'Color', 'w');
    
end





