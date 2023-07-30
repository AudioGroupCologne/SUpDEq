function plot_joergensen2013(simSRTs, varargin)
%PLOT_JOERGENSEN2013  Plot SRTs in the style of Fig. 2 from Joergensen et al., (2013)
%
%   Usage: plot_joergensen2013(SRT,flag)
%
%   PLOT_JOERGENSEN2013(SRT) plots the output of
%   JOERGENSEN2013 in the style of Fig. 2 of Joergensen et al., (2013).
%
%   The flag is:
%
%     'fig2'  plots the data in the style of figure 2 of  Jørgensen, Ewert and Dau (2013)
%
%
%   Please cite Joergensen et al. (2013) if you use this model
%
%   See also: joergensen2013
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_joergensen2013.php


%   #StatusDoc: Submitted
%   #StatusCode: Submitted
%   #Verification: Untrusted
%   #Requirements: MATLAB M-Signal M-Stats
%   #Author: Peter L. Sondergaard (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% Define input flags
definput.flags.type = {'fig2'};

[flags,kv]  = ltfatarghelper({}, definput,varargin);
%%
    if flags.do_fig2
        % [SRT_data,STDs_data ] = data_joergensen2013('fig2');
        
        %SRTs{1} = [SRT_data.SRTs_Jetal2013(1) SRT_data.SRTs_FP1990(1) SRT_data.SRTs_Kjems2009(1) SRT_data.SRTs_Kjems2009(3:4)];
        %SRTs{2} = [SRT_data.SRTs_Kjems2009(2)  SRT_data.SRTs_Jetal2013(2:3) SRT_data.SRTs_FP1990(2:3)];
        %SRTs{3} =  SRT_data.SRTs_reverbJandD2011;
        %SRTs{4} = SRT_data.SRTs_specsubJandD2011;
        
        %STDs{1} = [ STDs_data.SRTs_Jetal2013_std(1) STDs_data.SRTs_FP1990_std(1) 0 0 0];
        %STDs{2} = [0  STDs_data.SRTs_Jetal2013_std(2:3) STDs_data.SRTs_FP1990_std(2:3)];
        %STDs{3} = STDs_data.SRTs_reverbJandD2011_std;
        %STDs{4} = STDs_data.SRTs_specsubJandD2011_std; 
                  
        sim{1} = [simSRTs.simSRTs_Jetal2013(1) simSRTs.simSRTs_FP1990(1) simSRTs.simSRTs_Kjems2009(1) simSRTs.simSRTs_Kjems2009(3:4)];
        sim{2} = [simSRTs.simSRTs_Kjems2009(2) simSRTs.simSRTs_Jetal2013(2:3) simSRTs.simSRTs_FP1990(2:3)];
        sim{3} = simSRTs.simSRTs_reverb;
        sim{4} = simSRTs.simSRTs_specsub;
        
              
        
        fig=figure;
        sh = 400;
        set(fig,'Position',[.2*sh, 0.15*sh, 1.5*sh, 1.5*sh]);
        xlabels = {'CLUE','PM','DAN','Car',' Bottle','Cafe','SAM','ISTS','SMN','RT'};
        
        fnts =14;
        fnts_txt  =12;
        fnts_le = 8;
        
        line_sty = {'-','--','-.','-'};
        
        mark_sty = {'s','>','d','d','d','d','s','s','>','>','v','v','d','v','v','v'};
        mark_size = [8 7];
        colors = {[1 1 1]*.001,[1 1 1]*.3,[1 1 1]*0.5, [1 1 1],[1 1 1]*.5, [1 1 1]*.4,[1 1 1]*.7};
        
        
        subplots1 = subplot(2,1,1);
        subplotPos1 = get(subplots1,'position');
        
        xmin = 0.5;
        xmax = 10.5;
        ymin = -30.5;
        ymax = 14;
        ytickmax = 14;
        ytickmin = -26;
        
        x = 1:10;
        %y = [SRTs{1} SRTs{2}];
        %ystd = [STDs{1} STDs{2}];
%         for k = 1:length(x)
%             errorbar(x(k),y(k),-ystd(k),ystd(k),...
%                 'linestyle',    'none',...
%                 'linewidth',    1,...
%                 'color',           colors{1},...
%                 'marker',           mark_sty{k},...
%                 'markerfacecolor',  [1 1 1],...
%                 'markersize',       mark_size(1));hold on
%         end
        
        y2 = [sim{1} sim{2}];
       offset =  -0.2;
        for k = 1:length(x)
            
            plot(x(k)+offset,y2(k),...
                'linestyle',    'none',...
                'linewidth',    1,...
                'color',           colors{1},...
                'marker',           'o',...
                'markerfacecolor',  colors{1},...
                'markersize',       mark_size(1));hold on
        end
        
        plot([5.5 5.5],[-40 20],'k-')
        plot([3.5 3.5],[-40 20],'k--')
        
        
%         RMSE(3) = sqrt(mean((y(10:end,1)-y(10:end,2)).^2));
        
        %creating legend
        markers = {'s','>','d','v','o'};
        colors = {[1 1 1]*.001,[1 1 1]*.3,[1 1 1]*0.5, [1 1 1]*.001,[1 1 1]*.5, [1 1 1]*.7,[1 1 1]*.7,[1 1 1]};
        for k = 1:4          
            hh1(k) = plot(100,100,'marker',markers{k},'color',colors{1}, 'markerfacecolor', colors{8},'markersize',  8,'linewidth', 1,'linestyle','none');
        end
         hh1(5) = plot(100,100,'marker',markers{5},'color',colors{1}, 'markerfacecolor', colors{1},'markersize',  8,'linewidth', 1,'linestyle','none');
        
        le = legend(hh1,'New data', 'Kjems et al. (2009)','Festen and Plomp (1990)', 'Jørgensen and Dau (2011)','Simulation');
        set(le,'box','on','fontsize',fnts_le);
        lepos = get(le,'position');
        set(le,'position',lepos + [0 0.07 0 0])
        hold off
        
        axis([xmin xmax ymin ymax]);
        
        set(gca,'fontsize',fnts);
        set(gca,'yTick',ytickmin:6:ytickmax,'yTickLabel',ytickmin:6:ytickmax,'fontsize',14,'fontname','times');
        set(gca,'xtick',x,'xTickLabel',xlabels,'fontsize',14,'fontname','times');
       
        ylabel1 = ylabel('SRT (dB)','Fontsize',fnts);
        
        text(6.4,-37.5,'Fluctuating interferer','Fontsize',fnts,'FontName', 'times')
        
        text(1.4,-37.5,'Stationary interferer','Fontsize',fnts,'FontName', 'times')
        
%         RMSE(1) = sqrt(mean((y(3:4)-y2(3:4)).^2));
%         RMSE(2) = sqrt(mean((y(5:9)-y2(5:9)).^2));
%         text(4,5.5,{'RMSE = ',[num2str(RMSE(1),2), ' dB']},'fontsize',fnts_txt,'FontName', 'times');
%         text(7.8,-25,{['RMSE = ',num2str(RMSE(2),2), ' dB']},'fontsize',fnts_txt,'FontName', 'times');
      
        text(1.8,7.6,'SSN','fontsize',14,'FontName', 'times','fontweight','bold');
        
        set(gca, 'FontName', 'times');
        
        ah=axes('position',get(gca,'position'),...
            'visible','off');
        
        subplots2 = subplot(2,1,2);
        subplotPos2 = get(subplots1,'position');
        xlabels = {'0','0.4','0.7','1.3','2.3','UNP','0.5','1','2','4','8'};
         markers = {'v','d'};
        xmin = 0.5;
        xmax = 11.5;
        ymin = -7.5;
        ymax = 9;
        ytickmax = 8;
        ytickmin = -6;
        
        x = 1:11;
        clear y;
        
        %y = [SRTs{3} SRTs{4}];
        %ystd = [STDs{3} STDs{4}];
%         for k = 1:length(x)
%             errorbar(x(k),y(k),-ystd(k),ystd(k),...
%                 'linestyle',    'none',...
%                 'linewidth',    1,...
%                 'color',           colors{1},...
%                 'marker',           markers{1},...
%                 'markerfacecolor',  [1 1 1],...
%                 'markersize',       mark_size(1));hold on
%         end
        
        y2 = [sim{3} sim{4}];
       offset =  -0.2;
        for k = 1:length(x)
            
            plot(x(k)+offset,y2(k),...
                'linestyle',    'none',...
                'linewidth',    1,...
                'color',           colors{1},...
                'marker',           'o',...
                'markerfacecolor',  colors{1},...
                'markersize',       mark_size(1));hold on
        end
        
        plot([5.5 5.5],[-30 15],'k-')
        
       
        colors = {[1 1 1],[1 1 1]*.3};
        for k = 1:2
            if k==1
                hh(k) = plot(100,100,'k','marker',markers{k}, 'markerfacecolor', colors{k},'markersize',  8,'linewidth', 1,'linestyle','none');
            else
                hh(k) = plot(100,100,'marker',markers{k}, 'color', colors{k},'markerfacecolor', colors{k},'markersize',  8,'linewidth', 1,'linestyle','none');
            end
        end
        
        hold off
        
        axis([xmin xmax ymin ymax]);
        
  set(gca,'fontsize',fnts);
        set(gca,'yTick',ytickmin:2:ytickmax,'yTickLabel',ytickmin:2:ytickmax,'fontsize',14,'fontname','times');
        set(gca,'xtick',x,'xTickLabel',xlabels,'fontsize',14,'fontname','times');
        ylabel1 = ylabel('SRT (dB)','Fontsize',fnts);
      
        text(1.2,-10,'Reverberation time (s)','Fontsize',fnts,'FontName', 'times')
        text(6.1,-10,'Spectral subtraction factor \alpha','Fontsize',fnts,'FontName', 'times')
        
%         RMSE(1) = sqrt(mean((y(2:5)-y2(2:5)).^2));
%         RMSE(2) = sqrt(mean((y(7:11)-y2(7:11)).^2));
%         text(2.5,-6,{['RMSE = ',num2str(RMSE(1),2), ' dB']},'fontsize',fnts_txt,'FontName', 'times');
%         text(8.6,7.2,{['RMSE = ',num2str(RMSE(2),2), ' dB']},'fontsize',fnts_txt,'FontName', 'times');

        set(gca, 'FontName', 'times');
        
      
    end;





