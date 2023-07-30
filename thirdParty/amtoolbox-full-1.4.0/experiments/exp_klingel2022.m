function exp_klingel2022(varargin)
%EXP_KLINGEL2022 Experiments of Klingel & Laback (2022)
%   Usage: [] = exp_klingel2022(flag) 
%
%   EXP_KLINGEL2022(flag) reproduces figures of the study from  
%   Klingel & Laback (2022).
%
%
%   The following flags can be specified
%
%     'fig3'    Reproduce Fig.3:
%               Pretest ILD weights for each frequency band, averaged 
%               across azimuths. Blue circles show the results of 
%               experiment 1, red triangles of experiment 2. ILD weights 
%               gradually increase for single-band conditions from low 
%               (794-1260 Hz), to mid-low (1414-1782 Hz), to mid-high 
%               (2245-2828 Hz), to high (3564-4490 Hz). The multiband 
%               condition shows ILD weights similar to the mid-low band.
%
%     'fig4'    Reproduce Fig.4:
%               ILD weights as a function of azimuth for each frequency 
%               band (columns) and group (rows; top row shows ITD groups,
%               bottom row shows ILD groups). ILD weights gradually 
%               increase from low- to high-frequency stimuli and on average
%               are lower for lateral azimuths. Significant weight changes
%               from pre- to posttest were observed for the ILD group in 
%               experiment 1 for the high-frequency band as well as in 
%               experiment 2 for the mid-high-frequency band. Additionally, 
%               a significant weight change from pre- to posttest was 
%               observed for the ITD group in experiment 1 for the mid-low-
%               frequency band at 9 deg azimuth as well as in experiment 2 at 
%               39 deg azimuth.
%
%
%   Additional info about the data:
%
%     'wILD_exp1'   ILD weights of experiment 1
%
%                   - first dimension: azimuth (3, 9, 15, 21 deg)
%                   - second dimension: frequency band / testing time (multiband
%                     pretest, low pretest, mid-low pretest,
%                     mid-high pretest, high pretest, multi
%                     posttest, low posttest, mid-low posttest,
%                     mid-high posttest, high posttest)
%                   - third dimension: subject (odd numbers belong to ITD group,
%                     even numbers belong to ILD group)
%
%     'wILD_exp2'   ILD weights of experiment 2
%
%                   - first dimension: azimuth (3, 9, 15, 21, 27, 33, 39 deg)
%                   - second dimension: frequency band / testing time (mid-low 
%                     pretest, mid-high pretest, mid-low 
%                     posttest, mid-high posttest)
%                   - third dimension: subject (odd numbers belong to ITD group,
%                     even numbers belong to ILD group)
%
%
%   Examples:
%   ---------
%
%   To display Fig.3 use :
%
%     exp_klingel2022('fig3');
%
%   To display Fig.4 use :
%
%     exp_klingel2022('fig4');
%
%   See also: data_klingel2022
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_klingel2022.php


%   #Author: Maike Klingel (2022)
%   #Author: Clara Hollomey (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.flags.type = {'missingflag','fig3','fig4'};
[flags,~]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

wILD1 = amt_load('klingel2022', 'wILD_exp1.mat', 'wILD_exp1');
wILD_exp1 = wILD1.wILD_exp1;
wILD2 = amt_load('klingel2022', 'wILD_exp2.mat', 'wILD_exp2');
wILD_exp2 = wILD2.wILD_exp2;

if flags.do_fig3
    figure;
    hold on
    errorbar((1:4)-.02,nanmean(permute(nanmean(wILD_exp1(:,2:5,:)),[3 2 1])),(nanstd(permute(nanmean(wILD_exp1(:,2:5,:)),[3 2 1])))/sqrt(19),'bo-','MarkerFaceColor','b','LineWidth',1.5,'CapSize',8,'MarkerSize',7);
    errorbar((2:3)+.02,nanmean(permute(nanmean(wILD_exp2(1:4,1:2,:)),[3 2 1])),(nanstd(permute(nanmean(wILD_exp2(1:4,1:2,:)),[3 2 1])))/sqrt(20),'r^-','MarkerFaceColor','r','LineWidth',1.5,'CapSize',8,'MarkerSize',7);
    errorbar(2.5,nanmean(permute(nanmean(wILD_exp1(:,1,:)),[3 2 1])),(nanstd(permute(nanmean(wILD_exp1(:,1,:)),[3 2 1])))/sqrt(19),'bo','MarkerFaceColor','b','LineWidth',1.5,'CapSize',8,'MarkerSize',7);
    title('Mean ILD Weights Pretest','FontSize',18.7);
    ylabel('ILD Weight','FontSize',18.7);
    ylim([0 1]);
    xlabel('Band','FontSize',18.7);
    xlim([.8 4.2]);
    xticks([1 2 2.5 3 4]);
    xticklabels({'Low','Mid-Low','Multi','Mid-High','High'});
    %legend('Experiment 1','Experiment 2','location','northwest','FontSize',15.3,'EdgeColor','none');
    legend('Experiment 1','Experiment 2','location','northwest');
end

if flags.do_fig4
    bandname={'Low','Mid-Low','Mid-High','High'};
    figure;
    for group=1:2
        for band=1:5
            if band==1
                subplot(2,5,group*5-1+band);
            else
                subplot(2,5,group*5-6+band);
            end
            hold on
            errorbar((3:6:21)+.5,nanmean(wILD_exp1(:,band,group:2:end),3),permute(nanstd(permute(wILD_exp1(:,band,group:2:end),[3 2 1])),[3 2 1])/sqrt(sum(~isnan(wILD_exp1(1,band,group:2:end)))),'bo-','MarkerFaceColor','b','LineWidth',1,'CapSize',6,'MarkerSize',5)
            errorbar((3:6:21)-.25,nanmean(wILD_exp1(:,band+5,group:2:end),3),permute(nanstd(permute(wILD_exp1(:,band+5,group:2:end),[3 2 1])),[3 2 1])/sqrt(sum(~isnan(wILD_exp1(1,band,group:2:end)))),'bo--','MarkerFaceColor','w','LineWidth',1,'CapSize',6,'MarkerSize',5)
            xlim([0 45]);
            xlabel('Azimuth (deg)','FontSize',9.9);
            ylim([0 1]);
            ylabel('ILD Weight','FontSize',9.9);
            if group==1
                if band==1
                    title('Multiband','FontSize',9.9)
                else
                    title(bandname{band-1},'FontSize',9.9)
                end
            end
        end
        for band=1:2
            subplot(2,5,group*5-4+band);
            hold on
            errorbar((3:6:39)+.25,nanmean(wILD_exp2(:,band,group:2:end),3),permute(nanstd(permute(wILD_exp2(:,band,group:2:end),[3 2 1])),[3 2 1])/sqrt(sum(~isnan(wILD_exp2(1,band,group:2:end)))),'r^-','MarkerFaceColor','r','LineWidth',1,'CapSize',6,'MarkerSize',5)
            errorbar((3:6:39)-.75,nanmean(wILD_exp2(:,band+2,group:2:end),3),permute(nanstd(permute(wILD_exp2(:,band+2,group:2:end),[3 2 1])),[3 2 1])/sqrt(sum(~isnan(wILD_exp2(1,band+2,group:2:end)))),'r^--','MarkerFaceColor','w','LineWidth',1,'CapSize',6,'MarkerSize',5)
        end
    end
    for panel=1:10
        subplot(2,5,panel);
        hold on
        plot([0 45],[.5 .5],'k--')
    end
    subplot(2,5,7);
    %legend('Experiment 1 Pretest','Experiment 1 Posttest','Experiment 2 Pretest','Experiment 2 Posttest','location','north','EdgeColor','none','FontSize',8.1);
    legend('Experiment 1 Pretest','Experiment 1 Posttest','Experiment 2 Pretest','Experiment 2 Posttest','location','north','EdgeColor','none');
end


