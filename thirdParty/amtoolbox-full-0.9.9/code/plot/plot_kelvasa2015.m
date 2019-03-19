function plot_kelvasa2015(data,varargin)
%PLOT_KELVASA2015 Plot figures from Kelvasa and Dietz (2015) 
%   Usage: plot_kelvasa2015('fig8a')
%
%   PLOT_KELVASA2015('fig8a','identifier','BTE','HRTFchannels',[3,4]);
%   plots paper figure 8a using the behind-the-ear
%   microphone channel('BTE')
%          
%   PLOT_KELVASA2015('fig8a',varargin);
%   plots figure data for fig8a with updated model
%   parameters specified as key/value pairs in 'varargin'.
%   NOTE: When re-producing figures with new paramters,
%   'identifier' must be set so previously computed data is not
%   
%        
%   The following flags can be specified;
%
%     'redo'    Recomputes data for specified figure
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'noplot'  Don't plot, only return data.
%
%     'fig5'    Reproduce Fig. 5.
%
%     'fig6'    Reproduce Fig. 6.
%
%     'fig8a'    Reproduce Fig. 8a.
%
%     'fig8b'    Reproduce Fig. 8b.
%
%     'fig9a'    Reproduce Fig. 9a.
%
%     'fig10'    Reproduce Fig. 10.
%
%     'fig12'    Reproduce Fig. 12.
%
%   References:
%     D. Kelvasa and M. Dietz. Auditory model-based sound direction
%     estimation with bilateral cochlear implants. Trends in Hearing,
%     19:2331216515616378, 2015.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/plot/plot_kelvasa2015.php

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

%   Authors: 
%            Daryl Kelvasa (daryl.kelvasa@uni-oldenburg.de) 2016
%            Mathias Dietz (mdietz@uwo.ca) 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve and compute model paramters
% Set flags
    
      definput.flags.type = {'missingflag','fig5','fig6','fig7','fig8a',...
                                           'fig8b','fig9a','fig10','fig12'};
      definput.flags.debugging = {'no_debug','debug'};
      definput.flags.plot = {'plot','no_plot'};
      definput.flags.plot_stage_fig = {'no_plot_stage_fig','plot_stage_fig'};
   
    % import default arguments from other functions
    definput.import={'kelvasa2015','amt_cache'};
                        
    [flags,kv]  = ltfatarghelper({},definput,varargin);
    
    if flags.do_missingflag
           flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
           sprintf('%s or %s',definput.flags.type{end-1},...
           definput.flags.type{end})];
           error('%s: You must specify one of the following flags: %s.',...
                 upper(mfilename),flagnames);end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Plot Figure 5 from Kelvasa and Dietz 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig5

binPos  = [ 10, 13, 16, 19, 22];
levels = 35:70;

%Set plot layout for consecutive plots
rows =  3;
axisHeight = 5;
axisWidth = 7;
xIndent = 4;
yIndent = 1.9;

%Set plot labels
colLabel{1} = {'SSN','\fontsize{18}noHRTF'};
colLabel{2} = {'SSN','\fontsize{18}HRTF'};
colLabel{3} = {'SSN','\fontsize{18}HRTF'};
colLabel{4} = {'Sp','\fontsize{18}noHRTF'};
colLabel{5} = {'Sp','\fontsize{18}HRTF'};
colLabel{6} = {'Sp','\fontsize{18}HRTF'};
colLabel{7} = {'PN','\fontsize{18}noHRTF'};
colLabel{8} = {'PN','\fontsize{18}HRTF'};
colLabel{9} = {'PN','\fontsize{18}HRTF'};
colLabel{10} = {'WN','\fontsize{18}noHRTF'};
colLabel{11} = {'WN','\fontsize{18}HRTF'};
colLabel{12} = {'WN','\fontsize{18}HRTF'};
                       
%Set colors for each bin plot
col = hsv(numel(binPos));    
col(2,:) = col(2,:) + [0.2 -0.5 0];

%Generate Main Figure
mainFig = figure('Units','centimeters','Position',[0 0 37 18.5], ...
    'PaperUnits','centimeters','PaperPosition',[0 0 37 18.5]);

for axisInd = 1 : 12
                                
%Generate axis
           cornerX = xIndent + axisWidth*(ceil(axisInd/rows)-1);
           if ismember(axisInd,[1,4,7,10]);  cornerY = yIndent;
           elseif ismember(axisInd,[2,5,8,11]);  cornerY = yIndent + axisHeight;
           else ismember(axisInd,[3,6,9,12]);  cornerY = yIndent + 2*axisHeight; 
           end
           
ax(axisInd)= axes          ('Parent',mainFig,...
                            'Units','centimeters',...
                            'ColorOrder',col,... 
                            'Fontsize',22,...
                            'XTick',45:10:65,...
                            'XTickLabel',[],...
                            'XLim',[35 75],...
                            'YTick',150:150:450,...
                            'YTickLabel',[],.....
                            'YLim',[0 500],...
                            'Position',...
                                [cornerX cornerY axisWidth axisHeight]);

                            hold(ax(axisInd),'on')           
                            axis(ax(axisInd),'manual')      
                     
dataInd = [1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8];                    
% Plot Data
if ismember(axisInd,[3,6,9,12])
    plotData = data(dataInd(axisInd)).currentPerBinPerLevel;
    plotData = plotData.*0.04;
else
    plotData = data(dataInd(axisInd)).spikeRatePerBinPerLevel;
end
plot(ax(axisInd),levels,plotData,'Linewidth',3)
    
%Labels, legend, grid  
               if ismember(axisInd,[3 6 9 12]); labPos = 40;
               else labPos = 410;end
               text(37,labPos,colLabel{axisInd},...
                        'Fontsize',22,'EdgeColor','black');
                    
               if ismember(axisInd,[1 4 7 10])
                    set(ax(axisInd),'XTickLabel',45:10:65);end
                              
               if ismember(axisInd,[3 6 9 12])
                   ylim(ax(axisInd),[0 50])
                    set(ax(axisInd),'YTick',[15,30,45]); 
               end
                
               if ismember(axisInd,[1 2])
                    set(ax(axisInd),'YTickLabel',[150,300,450]); 
               end
   
               if axisInd == 7
                    h = xlabel(ax(axisInd),'Signal Level (dB SPL)');
                    p = get(h,'Position');p(1) = p(1) - 20;
                    set(h,'position',p);
               end
                    
               if axisInd == 3
                    set(ax(axisInd),'YTickLabel',[15,30,45])
                    curTitle = sprintf('Total charge\n (nC)');
                    h =ylabel(ax(axisInd),curTitle,'Fontsize',22);
                    p = get(h,'position'); p = p + [-2.5 0 0];
                    set(h,'position',p);
               end
                    
          
               if axisInd == 1
                    ANTitle = sprintf('AN response rate\n(spikes/sec)');
                    h = ylabel(ax(axisInd),ANTitle,'Fontsize',22);
                    p = get(h,'position'); p = p + [- 1 300 0];
                    set(h,'position',p);
               end
                
end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Plot Figure 6 from Kelvasa and Dietz 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig6

%Set plot layout for consecutive plots
rows =  2;
columns = 3; 
axisHeight = 4.8;
axisWidth = 9;
xIndent = 4.5;
yIndent = 1.9;

%Needed parameters
azis = data(1).azis;
binPos  = [ 10, 13, 16, 19, 22];

%Set plot labels
xAxisLab = 'Azimuthal angle';
yAxisLab1 = 'Right - Left';
yAxisLab3 = 're 0^\circ';
yAxisLab2 = sprintf('AN response rate differences\n(spikes/sec)');  
colLabel{3} = sprintf('SSN\n65dB');
colLabel{2} = sprintf('SSN\n55dB');
colLabel{1} = sprintf('SSN\n45dB');
azisLabels = [num2str((0:30:90)'),repmat('^\circ',numel(0:30:90),1)];
%label = [num2str(binPos'),repmat(('mm'),5,1)];                         

%Set colors for each bin plot
col = hsv(numel(binPos));    
col(2,:) = col(2,:) + [0.2 -0.5 0];

%Generate Main Figure
mainFig = figure('units','centimeters','position',[0 0 38 13],...
    'PaperUnits','centimeters','PaperPosition',[0 0 38 13]);

for axisInd = 1 : rows*columns
%Generate axis
           cornerX = xIndent + axisWidth*(ceil(axisInd/rows)-1) + ...
                        ceil(axisInd/rows)-1;
           if ismember(axisInd,[1,3,5]);  cornerY = yIndent;
           else cornerY = yIndent + axisHeight; end
             
  ax(axisInd)= axes      ('Parent',mainFig,...
                          'units','centimeters',...
                          'ColorOrder',col,...  
                          'Fontsize',24,...
                          'XTick',0:30:90,...
                          'XTickLabel',[],...
                          'XLim',[0 90],...
                          'YTick',0:50:300,...
                          'YTickLabel',[],...
                          'YGrid','on',...
                          'YLim',[-50 350],...
                          'Position',...
                          [cornerX cornerY axisWidth axisHeight]);
                hold(ax(axisInd),'on')           
                axis(ax(axisInd),'manual')  
          
% Process Data
        if ismember(axisInd,[1 2]); results = data(1); 
        elseif ismember(axisInd,[3 4]); results = data(2); 
        else results = data(3); 
        end 

        rghtSpks = (results.SpkSumPerBin + results.SpkDiffPerBin)./2;
        lftSpks = results.SpkSumPerBin - rghtSpks;
  
        rghtSpks = rghtSpks - repmat(rghtSpks(1,:),size(rghtSpks,1),1);
        lftSpks = lftSpks - repmat(lftSpks(1,:),size(rghtSpks,1),1);

                    
% Plot Data
        if ismember(axisInd,[1 3 5])
        plot(ax(axisInd),azis,results.SpkDiffPerBin,...
                  'Linewidth',4);                   
                    
        else 
        plot(ax(axisInd),azis,rghtSpks,...
                  'Linewidth',4,...
                  'Linestyle','--')
                     
        plot(ax(axisInd),azis,lftSpks,...
                  'Linewidth',4);
        end
    
%Labels, legend, grid
      %Grid                   
      plot(ax(axisInd),azis,zeros(numel(azis)),...
                  'Color','black','Linewidth',2,'Linestyle','--');
            
      plot(ax(axisInd),azis,zeros(numel(azis)),...
                  'Color','black','Linewidth',3)

                    
      if ismember(axisInd,[2 4 6])        
                  %axes(ax(axisInd));   
                  text(3,-300,colLabel{axisInd/2},...
                        'Fontsize',24,...
                        'EdgeColor','black',...
                        'BackgroundColor','white');
      end                    
         
      if ismember(axisInd,2:2:6) 
                    ylim(ax(axisInd),[-300 120])
                    set(ax(axisInd),'YTick',-250:50:100);
      end
            
      if ismember(axisInd,2)
                    set(ax(axisInd),'YTickLabel',...
                        [{''},-200,{''},-100,{''}, 0,{''}, 100 ]);
      end
               
      if ismember(axisInd,1)
                    set(ax(axisInd),'YTickLabel',...
                        [0,{''},100,{''},200,{''},300]);
      end
               
      if ismember(axisInd,[1 3 5])
                    set(ax(axisInd),'XTickLabel',azisLabels);
      end
               
      if axisInd == 3
                    h = xlabel(ax(axisInd),xAxisLab);
                    p = get(h,'Position');%p(1) = p(1) - 50;
                    set(h,'position',p);
      end
      
      
      if axisInd == 2
                    h = ylabel(ax(axisInd),yAxisLab2);
                    p = get(h,'Position'); p(1) = p(1) - 5;
                    p(2) = p(2) - 250; set(h,'position',p);
      end

                
      if axisInd == 6
                    %axes(ax(axisInd));
                    text(100,-40,yAxisLab3,'Fontsize',24,...
                             'Rotation',90,...
                             'VerticalAlignment','middle',...
                             'HorizontalAlignment','center');
      end
                    
      if axisInd == 5
                    %axes(ax(axisInd));
                    text(102,150,yAxisLab1,'Fontsize',24,...
                             'Rotation',90,...
                             'VerticalAlignment','middle',...
                             'HorizontalAlignment','center');
      end
                     
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Plot Figure 7 from Kelvasa and Dietz 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig7

%Set plot layout for consecutive plots
rows =  1;
columns = 2; 
axisHeight = 5.5;
axisWidth = 11;
xIndent = 4.7;
yIndent = 1.9;

%Needed parameters
azis = 0:5:90;
binPos  = [ 10, 13, 16, 19, 22];

%Set plot labels
xAxisLab = 'Azimuthal angle';
yAxisLab2 = sprintf('AN response\n rate differences\n(spikes/sec)');  
colLabel{1} = sprintf('SSN\n65dB');
colLabel{2} = sprintf('PN\n65dB');
azisLabels = [num2str((0:30:90)'),repmat('^\circ',numel(0:30:90),1)];
%label = [num2str(binPos'),repmat(('mm'),5,1)];                         

%Set colors for each bin plot
col = hsv(numel(binPos));    
col(2,:) = col(2,:) + [0.2 -0.5 0];

%Generate Main Figure
mainFig = figure('units','centimeters','position',[0 0 30 8],...
    'PaperUnits','centimeters','Paperposition',[0 0 30 8]);

for axisInd = 1 : rows*columns
%Generate axis
           if axisInd == 1; cornerX = xIndent;
           else cornerX = 1.5*xIndent + axisWidth;end
           cornerY = yIndent; 
             
  ax(axisInd)= axes      ('Parent',mainFig,...
                          'units','centimeters',...
                          'ColorOrder',col,...  
                          'Fontsize',24,...
                          'XTick',0:30:90,...
                          'XTickLabel',[],...
                          'XLim',[0 90],...
                          'YLim',[-120 280],...
                          'YTick',-100:50:250,...
                          'YTickLabel',[],...
                          'YGrid','on',...
                          'Position',...
                          [cornerX cornerY axisWidth axisHeight]);
                hold(ax(axisInd),'on')           
                axis(ax(axisInd),'manual')  
          
% Process Data
        results = data(axisInd); 
                    
% Plot Data
        plot(ax(axisInd),azis,results.SpkDiffPerBin,...
                  'Linewidth',4);                   
                    
%Labels, legend, grid
      %Grid                   
      plot(ax(axisInd),azis,zeros(numel(azis)),...
                  'Color','black','Linewidth',2,'Linestyle','--');
            
      plot(ax(axisInd),azis,zeros(numel(azis)),...
                  'Color','black','Linewidth',3)

                            
     axes(ax(axisInd));   
                  text(3,150,colLabel{axisInd},...
                        'Fontsize',24,...
                        'EdgeColor','black',...
                        'BackgroundColor','white');
               
      if ismember(axisInd,1)
                    set(ax(axisInd),'YTickLabel',...
                        [-100,{''},0,{''},100,{''},200,{''}]);
      end
               
      if ismember(axisInd,[1 2])
                    set(ax(axisInd),'XTickLabel',azisLabels);
      end
               
      if axisInd == 1
                    h = xlabel(ax(axisInd),xAxisLab);
                    p = get(h,'Position');p(1) = p(1) + 50;
                    set(h,'position',p);
      end
      
      
      if axisInd == 1
                    h = ylabel(ax(axisInd),yAxisLab2);
      end

               
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Plot Figure 8a, 8b, 9a from Kelvasa and Dietz 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig8a || flags.do_fig8b || flags.do_fig9a

%Set plot layout for consecutive plots
axisHeight = 3.8;
axisWidth1 = 3.8;
axisWidth2 = 1; 
xIndent = 2.3;
yIndent = 2;

%Set plot labels
xAxisLab =  'Target angle';  
yAxisLab = 'Predicted angle';  
axisLabel = {'45dB','22mm','16mm','10mm','65dB','55dB','19mm','13mm'};
azisLabels = {'0^\circ','30^\circ','60^\circ','90^\circ'};

%Set colors for each bin plot
col = hsv(5);    
col(2,:) = col(2,:) + [0.2 -0.5 0];

%Generate Main Figure
mainFig = figure('units','centimeters','position',[1 1 12 28.5],...
    'PaperUnits','centimeters','PaperPosition',[1 1 12 28.5]);
% mainFig = figure('units','normalized',...
%                         'position',[0.05 0.05 0.5 0.9]);
                    
binPlotInd = [3,5,7,13,15]; binPlotMap(binPlotInd) = [5,3,1,4,2];

for axisInd = 1 : 2: 16
 
%Generate axis
        if ismember(axisInd,1:8); cornerX = xIndent; end 
        if ismember(axisInd,9:16); cornerX = xIndent + axisWidth1 +...
                                            axisWidth2; end 
        if ismember(axisInd,[1 9]); cornerY = yIndent; end 
        if ismember(axisInd,[3 11]); cornerY = yIndent + axisHeight; end 
        if ismember(axisInd,[5 13]); cornerY = yIndent + 2*axisHeight; end 
        if ismember(axisInd,[7 15]); cornerY = yIndent + 3*axisHeight; end 
 

ax(axisInd)=axes('Parent',mainFig,...
                'units','centimeters',...
                'fontsize', 24,...
                'xlim',[-10 100],...
                'ylim',[-10 100],...
                'XTick',0:30:90,...
                'XTickLabel',[],...
                'YTick',0:30:90,...
                'YTickLabel',[],...
                'Position',...
                [cornerX cornerY axisWidth1 axisHeight]);
                hold(ax(axisInd),'on')           
                axis(ax(axisInd),'manual')
                grid(ax(axisInd),'on')
             

          
ax(axisInd+1)=axes('Parent',mainFig,...
                 'units','centimeters',...
                 'xlim',[-10 100],...
                 'ylim',[0 20],...
                  'View',[90 -90],...
                 'Color',get(gcf,'color'),...
                 'xcolor',get(gcf,'color'),...
                 'ycolor',get(gcf,'color'),...
                 'color',get(gcf,'color'),...
                 'ytick',[],...
                 'xtick',[],...
                 'position',[cornerX+axisWidth1 cornerY ...
                    axisWidth2 axisHeight]);
                hold(ax(axisInd+1),'on')
                axis(ax(axisInd+1),'manual')
         
         
%Process data 
        
        if axisInd == 1; 
            predictions = data(1); color = 'black'; 
            predictions = predictions.groupedWeightedPrediction;
        end
        if axisInd == 9; 
            predictions = data(3); color = 'black'; 
            predictions = predictions.groupedWeightedPrediction; 
        end
        if axisInd == 11; 
            predictions= data(2);color = 'black'; 
            predictions = predictions.groupedWeightedPrediction; 
        end
        if ismember(axisInd,binPlotInd);  
            predictions = data(2);
            predictions = ...
               predictions.groupedBinPredictions{binPlotMap(axisInd)}; 
            color = col(binPlotMap(axisInd),:); 
        end
               
      rmsError = sqrt(sum((predictions(1,:) - ...
                                        predictions(3,:)).^2.*...
                                            predictions(2,:))./...
                                            sum(predictions(2,:)));

%Plot data
        scatter(predictions(3,:),predictions(1,:),...
            predictions(2,:).*3,color,'Linewidth',2,'Parent',ax(axisInd)) 
             
        hist(ax(axisInd+1),predictions(1,:),-10:5:100)
        h = findobj(ax(axisInd+1),'Type','patch');
        set(h,'FaceColor',color);
        set(h,'EdgeColor','black');
   
%Labels, legend, grid
         plot(ax(axisInd),-10:100,-10:100,'linestyle','-',...
                        'color','black')
                    
         axes(ax(axisInd))
         text(45,5,[num2str(round(rmsError)),'^\circ'],...
                        'Fontsize',24)

         text(-8,85,axisLabel{ceil(axisInd/2)},...
                        'Fontsize',18,...
                        'Edgecolor','black',...
                        'BackgroundColor','white')  
            
         if ismember(axisInd,1:2:8)
                        set(ax(axisInd),...
                        'YTickLabel',azisLabels,...
                        'Fontsize',20);
         end
               
         if ismember(axisInd,[1,9])
                       set(ax(axisInd),...
                       'XTickLabel',azisLabels,...
                       'Fontsize',20);
         end
               
         if axisInd == 9
                    h = xlabel(ax(axisInd),xAxisLab);
                    p = get(h,'Position');p(2) = p(2)-5;
                    p(1) = p(1)-60; set(h,'position',p);
         end
                
         if axisInd == 5
                    h = ylabel(ax(axisInd),yAxisLab);
                    p = get(h,'Position');p = p + [0 -65 0];
                    set(h,'position',p);
         end                 
end       
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Plot Figure 10 from Kelvasa and Dietz 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig10


%Set plot layout for consecutive plots
axisHeight = 6;
axisWidth1 = 6;
axisWidth2 = 1; 
xIndent = 2.4;
yIndent = 1.7;

%Set plot labels
xAxisLab =  'Target angle';  
yAxisLab = 'Predicted angle';  
azisLabels = {'0^\circ','30^\circ','60^\circ','90^\circ'};
axisLabel = {{'\fontsize{16}Continuous','Speech', '10s'},...
             {'\fontsize{16}Sp 55dB','200 ms'},...
             {'\fontsize{16}Continuous','Speech', '10s'},...
             {'\fontsize{16}Sp 55dB','200 ms'}};

%Generate Main Figure
mainFig = figure('units','centimeters','position',[0 0 17 16],...
    'PaperUnits','centimeters','PaperPosition',[0 0 17 16]);
                    
for axisInd = 1 : 2: 8
 
%Generate axis
        if ismember(axisInd,1:4); cornerX = xIndent; end 
        if ismember(axisInd,5:8); cornerX = xIndent + axisWidth1 +...
                                            axisWidth2; end 
        if ismember(axisInd,[1,5]); cornerY = yIndent; end  
        if ismember(axisInd,[3,7]); cornerY = yIndent + axisHeight; end  


ax(axisInd)=axes('Parent',mainFig,...
                'units','centimeters',...
                'fontsize', 24,...
                'xlim',[-10 100],...
                'ylim',[-10 100],...
                'XTick',0:30:90,...
                'XTickLabel',[],...
                'YTick',0:30:90,...
                'YTickLabel',[],...
                'Position',...
                [cornerX cornerY axisWidth1 axisHeight]);
                hold(ax(axisInd),'on')           
                axis(ax(axisInd),'manual')
                grid(ax(axisInd),'on')
             

          
ax(axisInd+1)=axes('Parent',mainFig,...
                 'units','centimeters',...
                 'xlim',[-10 100],...
                 'ylim',[0 20],...
                  'View',[90 -90],...
                 'Color',get(gcf,'color'),...
                 'xcolor',get(gcf,'color'),...
                 'ycolor',get(gcf,'color'),...
                 'color',get(gcf,'color'),...
                 'ytick',[],...
                 'xtick',[],...
                 'position',[cornerX+axisWidth1 cornerY ...
                    axisWidth2 axisHeight]);
                hold(ax(axisInd+1),'on')
                axis(ax(axisInd+1),'manual')
         
         
%Process data 
         
            predictions = data((axisInd+1)/2);
            predictions = ...
               predictions.groupedWeightedPrediction; 
            color = 'black'; 
               
            rmsError = sqrt(sum((predictions(1,:) - ...
                                        predictions(3,:)).^2.*...
                                            predictions(2,:))./...
                                            sum(predictions(2,:)));

%Plot data
        scatter(predictions(3,:),predictions(1,:),...
            predictions(2,:).*3,color,'Linewidth',2,'Parent',ax(axisInd)) 
             
        hist(ax(axisInd+1),predictions(1,:),-10:5:100)
        h = findobj(ax(axisInd+1),'Type','patch');
        set(h,'FaceColor',color);
        set(h,'EdgeColor','black');
   
%Labels, legend, grid
         plot(ax(axisInd),-10:100,-10:100,'linestyle','-',...
                        'color','black')
                    
         axes(ax(axisInd))
         text(45,5,[num2str(round(rmsError)),'^\circ'],...
                        'Fontsize',24)

         text(-8,80,axisLabel{ceil(axisInd/2)},...
                        'Fontsize',18,...
                        'Edgecolor','black',...
                        'BackgroundColor','white')  
            
         if ismember(axisInd,[1,3])
                        set(ax(axisInd),...
                        'YTickLabel',azisLabels,...
                        'Fontsize',20);
         end
               
         if ismember(axisInd,[1,5])
                       set(ax(axisInd),...
                       'XTickLabel',azisLabels,...
                       'Fontsize',20);
         end
               
         if axisInd == 1
                    h = xlabel(ax(axisInd),xAxisLab);
                    p = get(h,'Position');
                    p(1) = p(1)+55; set(h,'position',p);
         end
                
         if axisInd == 1
                    h = ylabel(ax(axisInd),yAxisLab);
                    p = get(h,'Position');
                    p(2) = p(2)+55; set(h,'position',p);
         end     
         
         if axisInd == 3
            title({'AN Response','difference m90'})
         end     
         if axisInd == 7
            title({'MLE'})
         end 
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Plot Figure 12 from Kelvasa and Dietz 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig12
 

%Set plot layout for consecutive plots
axisHeight = 6;
axisWidth1 = 6;
axisWidth2 = 1; 
xIndent = 2.4;
yIndent = 1.7;

%Set plot labels
xAxisLab =  ('Target angle');  
yAxisLab = ('Predicted angle');  
axisLabel = {'SSN',{'SSN','5dB SNR'}};
azisLabels = {'0^\circ','30^\circ','60^\circ','90^\circ'};

%Set colors for each bin plot
col = hsv(5);    
col(2,:) = col(2,:) + [0.2 -0.5 0];

%Generate Main Figure
mainFig = figure('units','centimeters','position',[0 0 17 8],...
    'PaperUnits','centimeters','PaperPosition',[0 0 17 8]);
                    
for axisInd = 1 : 2: 4
 
%Generate axis
        if ismember(axisInd,1:2); cornerX = xIndent; end 
        if ismember(axisInd,3:4); cornerX = xIndent + axisWidth1 +...
                                            axisWidth2; end 
        if ismember(axisInd,1:4); cornerY = yIndent; end  


ax(axisInd)=axes('Parent',mainFig,...
                'units','centimeters',...
                'fontsize', 24,...
                'xlim',[-10 100],...
                'ylim',[-10 100],...
                'XTick',0:30:90,...
                'XTickLabel',[],...
                'YTick',0:30:90,...
                'YTickLabel',[],...
                'Position',...
                [cornerX cornerY axisWidth1 axisHeight]);
                hold(ax(axisInd),'on')           
                axis(ax(axisInd),'manual')
                grid(ax(axisInd),'on')
             

          
ax(axisInd+1)=axes('Parent',mainFig,...
                 'units','centimeters',...
                 'xlim',[-10 100],...
                 'ylim',[0 20],...
                  'View',[90 -90],...
                 'Color',get(gcf,'color'),...
                 'xcolor',get(gcf,'color'),...
                 'ycolor',get(gcf,'color'),...
                 'color',get(gcf,'color'),...
                 'ytick',[],...
                 'xtick',[],...
                 'position',[cornerX+axisWidth1 cornerY ...
                    axisWidth2 axisHeight]);
                hold(ax(axisInd+1),'on')
                axis(ax(axisInd+1),'manual')
         
         
%Process data 
         
            predictions = data((axisInd+1)/2);
            predictions = ...
               predictions.groupedWeightedPrediction; 
            color = 'black'; 
               
            rmsError = sqrt(sum((predictions(1,:) - ...
                                        predictions(3,:)).^2.*...
                                            predictions(2,:))./...
                                            sum(predictions(2,:)));

%Plot data
        scatter(predictions(3,:),predictions(1,:),...
            predictions(2,:).*3,color,'Linewidth',2,'Parent',ax(axisInd)) 
             
        hist(ax(axisInd+1),predictions(1,:),-10:5:100)
        h = findobj(ax(axisInd+1),'Type','patch');
        set(h,'FaceColor',color);
        set(h,'EdgeColor','black');
   
%Labels, legend, grid
         plot(ax(axisInd),-10:100,-10:100,'linestyle','-',...
                        'color','black')
                    
         axes(ax(axisInd))
         text(45,5,[num2str(round(rmsError)),'^\circ'],...
                        'Fontsize',24)

         text(-8,85,axisLabel{ceil(axisInd/2)},...
                        'Fontsize',18,...
                        'Edgecolor','black',...
                        'BackgroundColor','white')  
            
         if ismember(axisInd,1)
                        set(ax(axisInd),...
                        'YTickLabel',azisLabels,...
                        'Fontsize',20);
         end
               
         if ismember(axisInd,[1,3])
                       set(ax(axisInd),...
                       'XTickLabel',azisLabels,...
                       'Fontsize',20);
         end
               
         if axisInd == 1
                    h = xlabel(ax(axisInd),xAxisLab);
                    p = get(h,'Position');
                    p(1) = p(1)+55; set(h,'position',p);
         end
                
         if axisInd == 1
                    h = ylabel(ax(axisInd),yAxisLab);
         end                 
end              

end
end

