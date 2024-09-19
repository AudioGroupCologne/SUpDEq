% AKplotMixingTime(tmp50_model_based, tmp95_model_based, tmp50_data_based, 
%                  tmp95_data_based, tmp50_interchannel_mean_data_based, tmp95_interchannel_mean_data_based,
%                  echo_dens, fs, do_print)
% Shows the results of the perceptual mixing time prediction from model
% and data based predictors. If data based prediction was selected it plots
% the curve of the echo density of each channel and marks the corresponding
% perceptual mixing times.
%
% see AKperceptualMixingTimeDemo for examples
%
% I N P U T:
% tmp50_model_based                  - average perceptual mixing time, model based prediction
% tmp95_model_based                  - 95%-point perceptual mixing time, model based prediction
% tmp50_data_based                   - average  mixing time, data based prediction
% tmp95_data_based                   - 95%-point perceptual mixing time, data based prediction
% tmp50_interchannel_mean_data_based - interchannel average of average perceptual mixing time, data based prediction, average from all channels
% tmp95_interchannel_mean_data_based - interchannel average of 95%-point perceptual mixing time, data based prediction, average from all channels
% echo_dens                          - echo density (cf. Abel & Huang (2006))
% fs                                 - sampling frequency
% do_print                           - save plot, [1/0]
%
% A. Lindau, L. Kosanke, 2011
% alexander.lindau@tu-berlin.de
% audio communication group
% Technical University of Berlin

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
%-------------------------------------------------------------------------%
function AKplotMixingTime(tmp50_model_based, tmp95_model_based, tmp50_data_based, tmp95_data_based, tmp50_interchannel_mean_data_based, tmp95_interchannel_mean_data_based, echo_dens, fs, do_print)

% figure properties
AKf

if isempty(tmp50_data_based)
    plot([4 4],[0.5 5.5],'k',[0 10],[4.6 4.6],'k')
    text(0.1,5,'MODEL-BASED','FontSize',12)
    text(0.1,4,['t_{mp50}= ',num2str(round(tmp50_model_based)),' ms'],'FontSize',10)
    text(0.1,3,['t_{mp95}= ',num2str(round(tmp95_model_based)),' ms'],'FontSize',10)
    text(4.5,5,'DATA-BASED','FontSize',12)
    text(4.5,3,'No results available.','FontSize',10)
    axis([0 12 0 6]),axis off
    
else
    % time vector for echo density plot
    t = linspace(0,length(echo_dens)/fs*1000,length(echo_dens));
    
    subplot(2,1,1)
    hold on
    
    % plot echo density
    colors = AKcolors;
    colors = colors.rgb([1 3:11],:);
    for nn = 1:size(echo_dens, 2)
        
        % plot non-NaN values
        plot(t,echo_dens(:,nn),'color', colors(mod(nn,10)+1,:), 'lineWidth', 1, 'handleVisibility', 'off')
        
        % check for NaN values in the echo density (occurs if the
        % calcualtion was speed up)
        id = isnan( echo_dens(:,nn) );
        echo_dens(id,nn) = 1;
        plot(t(id), echo_dens(id,nn), '--k', 'color', colors(mod(nn,10)+1,:), 'lineWidth', 1, 'handleVisibility', 'off')
        
    end
    xlabel('t in [ms]','FontSize',8)
    ylabel('echo density','FontSize',8), title('Data-based prediction using echo-densitiy-approach from Abel & Huang (2006, criterion I)','FontSize',8)
    box on
    
    % plot estimated tmp-values
    for n = 1:length(tmp50_data_based)
        
        if round(tmp50_data_based(n)/1000*fs) > length(echo_dens)
            plot(tmp50_data_based(n),1,'or','LineWidth',4)
        else
            plot(t(round(tmp50_data_based(n)/1000*fs)),echo_dens(round(tmp50_data_based(n)/1000*fs),n),'or','LineWidth',4)
        end
        
        if round(tmp95_data_based(n)/1000*fs) > length(echo_dens)
            plot(tmp95_data_based(n),1,'og','LineWidth',4)
        else
            plot(t(round(tmp95_data_based(n)/1000*fs)),echo_dens(round(tmp95_data_based(n)/1000*fs),n),'og','LineWidth',4)
        end
        
    end
    legend('t_{mp50}','t_{mp95}','Location','SouthEast')
    
    hold off
    set(gca,'XScale','log')
    axis([1 max(max(t),max(max(tmp50_data_based),max(tmp95_data_based))) 0 1.2]),grid on
    
    subplot(2,1,2)
    plot([4 4],[0.5 5.5],'k',[0 13.5],[4.6 4.6],'k')
    % write results of model based predicition into figure
    if isempty(tmp50_model_based)
        text(0.1,5,'MODEL-BASED','FontSize',12)
        text(0.1,3,'No results available.','FontSize',10)
    else
        text(0.1,5,'MODEL-BASED','FontSize',12)
        text(0.1,4,['t_{mp50}= ',num2str(round(tmp50_model_based)),' ms'],'FontSize',10)
        text(0.1,3,['t_{mp95}= ',num2str(round(tmp95_model_based)),' ms'],'FontSize',10)
    end
    
    % write results of data based prediction into figure
    
    % results of prediction for more than one channel
    if length(tmp50_data_based)>1
        for n = 1:length(tmp50_data_based)
            text(4.5,5,'DATA-BASED','FontSize',12)
            text(4.5+(n-1)*4,4,['t_{mp50} (ch. ',num2str(n)',') = ',num2str(round(tmp50_data_based(n))),' ms'],'FontSize',10)
            text(4.5+(n-1)*4,3,['t_{mp95} (ch. ',num2str(n)',') = ',num2str(round(tmp95_data_based(n))),' ms'],'FontSize',10)
            text(4.5,2,['t_{mp50} (interchannel mean) = ',num2str(round(tmp50_interchannel_mean_data_based)),' ms'],'FontSize',10)
            text(4.5,1,['t_{mp95} (interchannel mean) = ',num2str(round(tmp95_interchannel_mean_data_based)),' ms'],'FontSize',10)
        end
    else
        % results of prediction for one channel
        text(4.5,5,'DATA-BASED','FontSize',12)
        text(4.5,4,['t_{mp50} = ',num2str(round(tmp50_data_based)),' ms'],'FontSize',10)
        text(4.5,3,['t_{mp95} = ',num2str(round(tmp95_data_based)),' ms'],'FontSize',10)
    end
end

axis([0 12 0 6]),axis off

% print
if do_print == 1
    print('-dtiff','-r150','tmp_prediction_results')
end