function exp_ziegelwanger2014(varargin)
%EXP_ZIEGELWANGER2014   Figures from Ziegelwanger and Majdak (2014)
%   Usage: data = exp_ziegelwanger2014(flag)
%
%   EXP_ZIEGELWANGER2014(flags) reproduces figures of the paper from
%   Ziegelwanger and Majdak (2014).
%
%   The following flags can be specified:
%
%     'redo'    Recalculate results.
%     'cached'  Use cached results. Default. 
%
%     'fig3'    Reproduce Fig. 3:
%               
%               TOAs resulting from TOA estimators in the horizontal plane 
%               applied on calculated HRTFs of the objects Sphere, SAT,
%               STP, as well on measured HRTFs of an exemplary listener
%               (NH89, ARI).
%
%     'fig5'    Reproduce Fig. 5:
%
%               Estimated on-axis model parameters from HRTFs calculated
%               for the centered object Sphere. Condition: combination of
%               actual parameters used in HRTF simulation. Circle and
%               cross: Parameters estimated for the left and right ears,
%               respectively. Gray lines: Actual parameters.
%
%     'fig6'    Reproduce Fig. 6:
%
%               Estimated on-axis model parameters from left-ear (top row)
%               and right-ear (bottom row) HRTFs of human listeners. Lines:
%               normal distribution fitted to the data. phi_e: Positive and
%               negative values corresponding to the left and right ears,
%               respectively.
%
%     'fig7'    Reproduce Fig. 7:
%
%               Estimated on-axis model IRDs from acoustically measured
%               HRTFs of listeners. Line: normal distribution fitted to the
%               data.
%
%     'fig8'    Reproduce Fig. 8:
%
%               Relative TOAs (top row) and on-axis model fit residuals
%               (bottom row) from HRTFs of STP (left column) and of an
%               exemplary listener (NH89, ARI; right column). Black points:
%               data classified as outliers by the ESD test with an upper
%               bound of outlier rate of 1% (see ziegelwanger2014). The
%               reference for the relative TOAs is the smallest TOA in each
%               HRTF set. Horizontal lines: pm1 sampling interval.
%
%     'fig9'    Reproduce Fig. 9:
%
%               Relative TOAs of an exemplary listener (NH89, ARI) in the
%               interaural horizontal plane for the left (black) and right
%               (gray) ears as results from the MCM estimator (symbols) and
%               the on-axis model (lines). The reference for the relative
%               TOAs is the smallest TOA in HRTF sets of both ears.
%
%     'fig10'    Reproduce Fig. 10:
%
%               Estimated on-axis model parameters from HRTFs calculated
%               for the non-centered object Sphere. Other details as in
%               Fig. 5.
%
%     'fig12'   Reproduce Fig. 12:
%
%               Estimated on-axis model parameters from HRTFs calculated
%               for the non-centered object Sphere. Other details as in
%               Fig. 5. phi_e: Positive and negative values correspond to
%               the left and right ears, respectively.
%
%     'tab1'    Reproduce Tab. 1:
%
%               ANRs and parameter errors (average pm1 standard deviation)
%               resulting from fitting the on-axis model to TOAs estimated
%               from HRTFs of the Sphere. Parameter errors: Differences
%               between the estimated and actual parameters.
%
%     'tab2'    Reproduce Tab. 2:
%
%               Parameters and ANRs resulting from fitting the on-axis
%               model to TOAs estimated from HRTFs of the objects Sphere,
%               SAT, STP. Actual parameters: r=87,5 mm, phi_e=90deg, and
%               theta_e=0deg.
%
%     'tab3'    Reproduce Tab. 3:
%
%               Parameters (average pm1 standard deviation) and ANRs
%               (median) resulting from fitting the on-axis model to TOAs
%               estimated from acoustically measured HRTFs of human
%               listeners. L: Left ear. R: Right ear. All: Results for all
%               listeners. NH89: Results for a single listener (NH89, ARI).
%
%     'tab5'    Reproduce Tab. 5:
%
%               Parameters and ANRs (average pm1 standard deviation)
%               resulting from fitting the off-axis model to TOAs estimated
%               from HRTFs of SAT, STP, and all listeners (All). Full: Fits
%               to full TOA sets. O-A: Fits to the outlier-adjusted TOA
%               sets. L: Left ear. R: Right ear.
%
%     'tab6'    Reproduce Tab. 6:
%
%               ANRs and parameter errors (average pm1 standard deviation)
%               resulting from fitting the off-axis model to TOAs estimated
%               from HRTFs of the object Sphere. Centered: Conditions in
%               which r and vec(e) varied from M=0. Non-centered:
%               Conditions in which M varied. Other defails as in Tab. 5.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Optimization Toolbox for Matlab
%
%   3) Optional: Data in hrtf/ziegelwanger2014. Will be downloaded on demand.
%
%   Examples:
%   ---------
%
%   To display Fig. 3, use :
%
%     exp_ziegelwanger2014('fig3');
%
%   To display Fig. 5, use :
%
%     exp_ziegelwanger2014('fig5');
%
%   To display Fig. 6, use :
%
%     exp_ziegelwanger2014('fig6');
%
%   To display Fig. 7, use :
%
%     exp_ziegelwanger2014('fig7');
%
%   To display Fig. 8, use :
%
%     exp_ziegelwanger2014('fig8');
%
%   To display Fig. 9, use :
%
%     exp_ziegelwanger2014('fig9');
%
%   To display Fig. 10, use :
%
%     exp_ziegelwanger2014('fig10');
%
%   To display Fig. 12, use :
%
%     exp_ziegelwanger2014('fig12');
%
%   To display Tab. 1, use :
%
%     exp_ziegelwanger2014('tab1');
%
%   To display Tab. 2, use :
%
%     exp_ziegelwanger2014('tab2');
%
%   To display Tab. 3, use :
%
%     exp_ziegelwanger2014('tab3');
%
%   To display Tab. 5, use :
%
%     exp_ziegelwanger2014('tab5');
%
%   To display Tab. 6, use :
%
%     exp_ziegelwanger2014('tab6');
%
%   See also: ziegelwanger2014, ziegelwanger2014_onaxis,
%   ziegelwanger2014_offaxis, data_ziegelwanger2014
%
%   References:
%     H. Ziegelwanger and P. Majdak. Modeling the direction-continuous
%     time-of-arrival in head-related transfer functions. J. Acoust. Soc.
%     Am., 135:1278--1293, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_ziegelwanger2014.php


%   #Author: Harald Ziegelwanger (2014)
%   #Author: Clara Hollomey (2020)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% ------ Check input options --------------------------------------------

definput.flags.type = {'missingflag',...
'fig3','fig5','fig6','fig7','fig8','fig9',...
'fig10','fig12','tab1','tab2','tab3','tab5',...
'tab6'};
definput.import={'amt_cache'}; % get the flags of amt_cache

% Parse input options
[flags,~]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}), ...
        sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% Figure 3
if flags.do_fig3

%load data
    Obj{1}=data_ziegelwanger2014('Sphere',flags.cachemode);
    Obj{2}=data_ziegelwanger2014('SAT',flags.cachemode);
    Obj{3}=data_ziegelwanger2014('STP',flags.cachemode);
    Obj{4}=data_ziegelwanger2014('NH89',flags.cachemode);
    
%plot figure
    figure('Position',[ 520   221   732   577]);
    
    methodLabel=['MAX';'CTD';'AGD';'MCM'];
    sty=[': ';'-.';'- ';'--'];
    clr=[0 0 0;...
        0 0 0;...
        0 0 0;...
        0 0 0];
    lw=[3 1 1 1];

    for method=1:4
        subplot(2,2,method)
        for hrtf=1:4
            if hrtf==4
                Obj{4}.Data.toaEst{method}=Obj{4}.Data.toaEst{method}+110;
            end
            h=plot_ziegelwanger2014(Obj{hrtf},Obj{hrtf}.Data.toaEst{method},4,'k',0,1,1,sty(hrtf,:),lw(hrtf));
            set(h,'color',clr(method,:));
            hold on
        end
        ylabel('TOA (ms)','Fontname','Arial','Fontsize',14);
        xlabel('');
        grid off
        ylim([2.9,4.3])
        xlim([-10,370])
        set(gca,'Fontname','Arial','Fontsize',10)
        title('')
        switch(method)
            case 1
                xlabel('')
            case 2
                ylabel('')
                xlabel('')
            case 3
                l=legend('Sphere','SAT','STP','NH89','Location','NorthWest');
                set(l,'Fontsize',9,'Fontname','Arial')
            case 4
                ylabel('')
        end

        switch(method)
            case 1
                tmp=get(gca,'Position');
                set(gca,'Position',tmp+[0 0.08 0.05 -0.08]);
                set(gca,'xticklabel',[]);
            case 2
                tmp=get(gca,'Position');
                set(gca,'Position',tmp+[-0.05 0.08 0.05 -0.08]);
                set(gca,'yticklabel',[]);
                set(gca,'xticklabel',[]);
            case 3
                tmp=get(gca,'Position');
                set(gca,'Position',tmp+[0 0.285 0.05 -0.08]);
                set(gca,'xticklabel',{'0' '90' '180' '270' '360 '});
            case 4
                tmp=get(gca,'Position');
                set(gca,'Position',tmp+[-0.05 0.285 0.05 -0.08]);
                set(gca,'yticklabel',[]);
                set(gca,'xticklabel',{'0' '90' '180' '270' '360 '});
        end

        text(300,4.1,methodLabel(method,:),'Fontname','Arial','Fontsize',14)
    end
    
end

%% Figure 5
if flags.do_fig5

%load data
    data=data_ziegelwanger2014('SPHERE_ROT',flags.cachemode);
    for ii=1:length(data.results)
        p_onaxis{1}(:,:,ii)=data.results(ii).MAX{1}.p_onaxis;
        p_onaxis{2}(:,:,ii)=data.results(ii).CTD{1}.p_onaxis;
        p_onaxis{3}(:,:,ii)=data.results(ii).AGD{1}.p_onaxis;
        p_onaxis{4}(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
    end

%plot figure
    figure
    tmp=get(gcf,'Position');
    set(gcf,'Position',tmp.*[1 1 1 2]);
    sym='ox';%plot symbols
    ms=8;%markersize
    lw1=2;%linewidth1
    lw2=2;%linewidth2
    fs=18;%fontsize
    ls='k-';%linestyle
    lc=[0.5 0.5 0.5];
    h=[];
    
    % radii
    h(end+1)=subplot(411);
    var=[squeeze(p_onaxis{4}(1,1,:))*1000 squeeze(p_onaxis{4}(1,2,:))*1000 data.radius(:)];
    for ch=1:size(p_onaxis{4},2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
        hold on
    end
    plot(1:14,var(1:14,3),ls,'Linewidth',lw2,'color',lc)
    plot(15:28,var(15:28,3),ls,'Linewidth',lw2,'color',lc)
    plot(29:42,var(29:42,3),ls,'Linewidth',lw2,'color',lc)
    for ch=1:size(p_onaxis{4},2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
        hold on
    end
    clear var;
    ylabel('r (mm)','Fontname','Arial','Fontsize',fs)
    ylim([72,105])
    xlim([-1,44])
    set(gca,'xtick',1:42)
    set(gca,'ytick',[80 90 100])
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0 0 0]);

    %phi
    h(end+1)=subplot(412);
    var=[squeeze(p_onaxis{4}(2,1,:))/pi*180 squeeze(p_onaxis{4}(2,2,:))/pi*180 data.phi+ones(length(data.phi),1)*90 data.phi-ones(length(data.phi),1)*90];
    for ch=1:size(p_onaxis{4},2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
        hold on
    end
    for ch=1:size(p_onaxis{4},2)
        plot(1:9,var(1:9,2+ch),ls,'Linewidth',lw2,'color',lc)
        plot(10:14,var(10:14,2+ch),ls,'Linewidth',lw2,'color',lc)
        plot(15:23,var(15:23,2+ch),ls,'Linewidth',lw2,'color',lc)
        plot(24:28,var(24:28,2+ch),ls,'Linewidth',lw2,'color',lc)
        plot(29:37,var(29:37,2+ch),ls,'Linewidth',lw2,'color',lc)
        plot(38:42,var(38:42,2+ch),ls,'Linewidth',lw2,'color',lc)
    end
    for ch=1:size(p_onaxis{4},2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
        hold on
    end
    clear var;
    set(gca,'YTick',[-90 90])
    ylabel('   _e (deg)','Fontname','Arial','Fontsize',fs)
    xlim([-1,44])
    set(gca,'xtick',1:42)
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0.055 0 0]);

    %theta
    h(end+1)=subplot(413);
    var=[squeeze(p_onaxis{4}(3,1,:))/pi*180 squeeze(p_onaxis{4}(3,2,:))/pi*180 data.theta -data.theta];
    for ch=1:size(p_onaxis{4},2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
        hold on
    end
    for ch=1:size(p_onaxis{4},2)
        plot(1:9,var(1:9,2+ch),ls,'Linewidth',lw2,'color',lc)
        plot(10:14,var(10:14,2+ch),ls,'Linewidth',lw2,'color',lc)
        plot(15:23,var(15:23,2+ch),ls,'Linewidth',lw2,'color',lc)
        plot(24:28,var(24:28,2+ch),ls,'Linewidth',lw2,'color',lc)
        plot(29:37,var(29:37,2+ch),ls,'Linewidth',lw2,'color',lc)
        plot(38:42,var(38:42,2+ch),ls,'Linewidth',lw2,'color',lc)
    end
    for ch=1:size(p_onaxis{4},2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
        hold on
    end
    clear var;
    ylabel('\theta_e (deg)','Fontname','Arial','Fontsize',fs)
    xlabel('Condition','Fontname','Arial','Fontsize',fs)
    ylim([-15,15])
    xlim([-1,44])
    set(gca,'xtick',1:42)
    set(gca,'xticklabel',[])
    set(gca,'ytick',[-10 0 10])
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0.11 0 0]);

    for ii=1:length(h)
        set(h(ii),'fontsize',14)
        set(h(ii),'Linewidth',2)
        tmp=findobj(h(ii),'Type','patch');
        set(tmp,'EdgeColor','k');
    end
    set(gcf,'color',[1 1 1])
    
end

%% Figure 6
if flags.do_fig6

%load data
    hrtf={'ARI','CIPIC','LISTEN'};
    for kk=1:length(hrtf)
        data=data_ziegelwanger2014(hrtf{kk},flags.cachemode);
        if kk==3
            data.results=data.results([1:27 29:end]);
        end
        for ii=1:length(data.results)
            temp1(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
            temp2(:,1:size(data.results(ii).MCM{1}.p_offaxis,2),ii)=data.results(ii).MCM{1}.p_offaxis;
            temp3(ii)=mean([data.results(ii).MCM{1}.performance.on_axis{1}.resnormS ...
                data.results(ii).MCM{1}.performance.on_axis{2}.resnormS]);
            temp4(ii)=mean([data.results(ii).MCM{1}.performance.off_axis{1}.resnormS ...
                data.results(ii).MCM{1}.performance.off_axis{2}.resnormS]);
            temp5(ii)=mean([data.results(ii).MCM{1}.performance.on_axis{1}.resnormP ...
                data.results(ii).MCM{1}.performance.on_axis{2}.resnormP]);
            temp6(ii)=mean([data.results(ii).MCM{1}.performance.off_axis{1}.resnormP ...
                data.results(ii).MCM{1}.performance.off_axis{2}.resnormP]);
            temp8(ii)=data.results(ii).MCM{1}.performance.on_axis{1}.resnormS;
            temp9(ii)=data.results(ii).MCM{1}.performance.on_axis{2}.resnormS;
            temp10(ii)=data.results(ii).MCM{1}.performance.off_axis{1}.resnormS;
            temp11(ii)=data.results(ii).MCM{1}.performance.off_axis{2}.resnormS;
        end
        p_onaxis{kk}=temp1;
        p_offaxis{kk}=temp2;
        resnormS_onaxis{kk}=temp3;
        resnormS_offaxis{kk}=temp4;
        resnormP_onaxis{kk}=temp5;
        resnormP_offaxis{kk}=temp6;
        resnormS_onaxis_left{kk}=temp8;
        resnormS_onaxis_right{kk}=temp9;
        resnormS_offaxis_left{kk}=temp10;
        resnormS_offaxis_right{kk}=temp11;
        clear temp1 temp2 temp3 temp4 temp5 temp6 temp7 temp8 temp9 temp10 temp11
    end

%plot figure
    figure('PaperUnits','centimeters','PaperType','A4','Paperposition',[0, 0, 21, 29.7],'Units','centimeters','Position',[0 0 21 29.7],'Resize','off')
    fs=14;%fontsize
    lw=1.1;%linewidth
    h=[];

    %radii
    temp=1;
    for kk=1:length(hrtf)
        var(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            squeeze(mean(p_onaxis{kk}(1,:,:)*1000,2));
        varl(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            squeeze(p_onaxis{kk}(1,1,:)*1000);
        varr(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            squeeze(p_onaxis{kk}(1,2,:)*1000);
        temp=size(var,1)+1;
    end

    h(1)=subplot(631);
    colormap('gray');
    binranges = 45:5:135;
    bincounts=histc(varl,binranges);
    bar(binranges,bincounts/sum(bincounts)*100,'Facecolor',[0.6 0.6 0.6],'edgecolor',[0 0 0]);
    hold on
    box on
    if ~isoctave
      [mu,sigma]=normfit(1:length(varl), varl);
      y=normpdf(binranges,mu,sigma);
      plot(binranges,y*100*5,'k','linewidth',lw)
      plot([-500 500],[0 0],'k');
      xlim([45,135])
      ylim([-1,25])
      ylabel('Left ear ','Fontname','Arial','Fontsize',fs)
      set(gca,'xtick',[60 80 100 120])
      set(gca,'xticklabel',[]);
      tmp=get(gca,'Position');
      set(gca,'Position',tmp+[0 0 0.04 0]);
    else
      xlabel('Sorry: normfit not available in Octave.')
    end
    
    h(2)=subplot(634);
    binranges = 45:5:135;
    bincounts=histc(varr,binranges);
    bar(binranges,bincounts/sum(bincounts)*100,'Facecolor',[0.6 0.6 0.6],'edgecolor',[0 0 0]);
    hold on
    box on
    if ~isoctave
      [mu,sigma]=normfit(varr);  
      y=normpdf(binranges,mu,sigma);
      plot(binranges,y*100*5,'k','linewidth',lw)
      plot([-500 500],[0 0],'k');
      xlim([45,135])
      ylim([-1,25])
      set(gca,'xtick',[60 80 100 120])
      xlabel('r (mm)','Fontname','Arial','Fontsize',fs)
      ylabel('Right ear ','Fontname','Arial','Fontsize',fs)
      tmp=get(gca,'Position');
      set(gca,'Position',tmp+[0 0.03 0.04 0]);
    else
      xlabel('Sorry: normfit not available in Octave')
    end
    
    clear var varl varr;

    % phi_e
    temp=1;
    for kk=1:length(hrtf)
        varl(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            squeeze(p_onaxis{kk}(2,1,:))*180/pi;
        varr(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            mod(squeeze(p_onaxis{kk}(2,2,:))*180/pi,360);
        temp=size(varl,1)+1;
    end

    h(3)=subplot(632);
    binranges = 55:2:125;
    bincounts=histc(varl,binranges);
    bar(binranges,bincounts/sum(bincounts)*100,'Facecolor',[0.6 0.6 0.6],'edgecolor',[0 0 0]);
    hold on
    box on
    if ~isoctave
      [mu,sigma]=normfit(varl);
      y=normpdf(binranges,mu,sigma);
      plot(binranges,y*100*2,'k','linewidth',lw)
      plot([-500 500],[0 0],'k');
      xlim([55,125])
      ylim([-1,25])
      set(gca,'xticklabel',[]);
      set(gca,'yticklabel',[]);
      tmp=get(gca,'Position');
      set(gca,'Position',tmp+[-0.02 0 0.04 0]);
    else
      xlabel('Sorry: normfit not available in Octave.')
    end
    
    clear y

    varr=abs(varr-360);
    h(4)=subplot(635);
    binranges = 55:2:125;
    bincounts=histc(varr,binranges);
    bar(binranges,bincounts/sum(bincounts)*100,'Facecolor',[0.6 0.6 0.6],'edgecolor',[0 0 0]);
    hold on
    box on
    if ~isoctave
      [mu,sigma]=normfit(varr);
      y=normpdf(binranges,mu,sigma);
      plot(binranges,y*100*2,'k','linewidth',lw)
      plot([-500 500],[0 0],'k');
      xlim([55,125])
      ylim([-1,25])
      xlabel('   _e (deg)','Fontname','Arial','Fontsize',fs)
      set(gca,'xtick',[60 80 100 120]);
      set(gca,'xticklabel',{'\pm60','\pm80','\pm100','\pm120'});
      set(gca,'yticklabel',[]);
      tmp=get(gca,'Position');
      set(gca,'Position',tmp+[-0.02 0.03 0.04 0]);
    else
      xlabel('Sorry: normfit not available in Octave')
    end
    
    clear varl varr a

    % theta_e
    temp=1;
    for kk=1:length(hrtf)
        varl(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            squeeze(p_onaxis{kk}(3,1,:))*180/pi;
        varr(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            squeeze(p_onaxis{kk}(3,2,:))*180/pi;
        temp=size(varl,1)+1;
    end

    h(5)=subplot(633);
    binranges = -25:2:15;
    bincounts=histc(varl,binranges);
    bar(binranges,bincounts/sum(bincounts)*100,'Facecolor',[0.6 0.6 0.6],'edgecolor',[0 0 0]);
    hold on
    box on
    if ~isoctave
      [mu,sigma]=normfit(varl);
      y=normpdf(binranges,mu,sigma);
      plot(binranges,y*100*2,'k','linewidth',lw)
      plot([-500 500],[0 0],'k');
      xlim([-25,15])
      ylim([-1,25])
      set(gca,'xticklabel',[]);
      set(gca,'yticklabel',[]);
      tmp=get(gca,'Position');
      set(gca,'Position',tmp+[-0.04 0 0.04 0]);
    else
      xlabel('Sorry: normfit not available in Octave')
    end
    
    h(6)=subplot(636);
    binranges = -25:2:15;
    bincounts=histc(varr,binranges);
    bar(binranges,bincounts/sum(bincounts)*100,'Facecolor',[0.6 0.6 0.6],'edgecolor',[0 0 0]);
    hold on
    box on
    if ~isoctave
      [mu,sigma]=normfit(varr);
      y=normpdf(binranges,mu,sigma);
      plot(binranges,y*100*2,'k','linewidth',lw)
      plot([-500 500],[0 0],'k');
      xlim([-25,15])
      ylim([-1,25])
      xlabel('\theta_e (deg)','Fontname','Arial','Fontsize',fs)
      set(gca,'yticklabel',[]);
      tmp=get(gca,'Position');
      set(gca,'Position',tmp+[-0.04 0.03 0.04 0]);
    else
      xlabel('Sorry: normfit not available in Octave.')
    end
    
    clear varl varr

    for ii=1:length(h)
        set(h(ii),'linewidth',lw)
        set(h(ii),'TickLength',[0.015 0.015])
        set(h(ii),'Fontname','Arial','Fontsize',12)
        tmp=findobj(h(ii),'Type','patch');
        set(tmp,'EdgeColor','k');
    end
    set(gcf,'color',[1 1 1])
end

%% Figure 7
if flags.do_fig7
    
%load data
    hrtf={'ARI','CIPIC','LISTEN'};
    for kk=1:length(hrtf)
        data=data_ziegelwanger2014(hrtf{kk},flags.cachemode);
        if kk==3
            data.results=data.results([1:27 29:end]);
        end
        for ii=1:length(data.results)
            temp1(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
            temp2(:,:,ii)=data.results(ii).MCM{1}.p_offaxis;
            temp3(:,:,ii)=data.results(ii).MCM{2}.p_offaxis;
        end
        p_onaxis{kk}=temp1;
        p_offaxis{kk,1}=temp2;
        p_offaxis{kk,2}=temp3;
        clear data temp1 temp2 temp3
    end
    data=data_ziegelwanger2014('SPHERE_ROT',flags.cachemode);
    for ii=1:length(data.results)
        temp(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
    end
    p_onaxis{4}=temp;
    clear data temp
    Obj=data_ziegelwanger2014('Sphere',flags.cachemode);
    [~,temp]=ziegelwanger2014(Obj,Obj.Data.toaEst{4},0,1);
    p_onaxis{4}(:,:,end+1)=temp.p_onaxis;
    clear data temp Obj
    Obj=data_ziegelwanger2014('SAT',flags.cachemode);
    [~,temp]=ziegelwanger2014(Obj,Obj.Data.toaEst{4},0,1);
    p_onaxis{4}(:,:,end+1)=temp.p_onaxis;
    clear data temp Obj
    Obj=data_ziegelwanger2014('STP',flags.cachemode);
    [~,temp]=ziegelwanger2014(Obj,Obj.Data.toaEst{4},0,1);
    p_onaxis{4}(:,:,end+1)=temp.p_onaxis;
    clear data temp Obj
    data=data_ziegelwanger2014('SPHERE_DIS',flags.cachemode);
    for ii=1:length(data.results)
        temp1(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
        temp2(:,:,ii)=data.results(ii).MCM{1}.p_offaxis;
    end
    p_onaxis{5}=temp1;
    p_offaxis{5}=temp2;
    clear data temp1 temp2

%plot figure
    figure('PaperUnits','centimeters','PaperType','A4','Paperposition',[0, 0, 21, 29.7],'Units','centimeters','Position',[0 0 21 29.7],'Resize','off')
    fs=14;%fontsize

    temp=1;
    for kk=1:3
        varl(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            squeeze(p_onaxis{kk}(1,1,:)*1000);
        varr(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            squeeze(p_onaxis{kk}(1,2,:)*1000);
        temp=size(varl,1)+1;
    end

    h(1)=subplot(621);
    binranges = -105:5:105;
    bincounts=histc(varl(:,1)-varr(:,1),binranges);
    bar(binranges,bincounts/sum(bincounts)*100,'Facecolor',[0.6 0.6 0.6]);
    hold on
    if ~isoctave
      [mu,sigma]=normfit(varl(:,1)-varr(:,1));
      y=normpdf(binranges,mu,sigma);
      plot(binranges,y*100*5,'k','linewidth',1.1)
      ylabel('Rel. Freq.','Fontname','Arial','Fontsize',fs)
      xlabel('IRD (mm)','Fontname','Arial','Fontsize',fs)
      xlim([-80,85])
      ylim([-1,16])
      set(gca,'xtick',[-60 -30 0 30 60])
    else
      xlabel('Sorry: normfit not available in Octave')
    end
    
    clear varl varr y

    clear varl varr temp
    temp=1;
    for kk=1:3
        varl(temp:temp+size(p_offaxis{kk,2},3)-1,:)= ...
            squeeze(p_offaxis{kk,2}(1,1,:)*1000);
        varr(temp:temp+size(p_offaxis{kk,2},3)-1,:)= ...
            squeeze(p_offaxis{kk,2}(1,2,:)*1000);
        temp=size(varl,1)+1;
    end

    h(2)=subplot(623);
    binranges = -105:1:105;
    bincounts=histc(varl(:,1)-varr(:,1),binranges);
    bar(binranges,bincounts/sum(bincounts)*100,'Facecolor',[0 0 0]);
    hold on
    if ~isoctave
      [mu,sigma]=normfit(varl(:,1)-varr(:,1));
      y=normpdf(binranges,mu,sigma);
      % plot(binranges,y*100*0.5,'k','linewidth',1.1)
      xlabel('IRD (mm)')
      ylabel('Rel. Freq.','Fontname','Arial','Fontsize',fs)
      xlabel('IRD (mm)','Fontname','Arial','Fontsize',fs)
      xlim([-80,85])
      ylim([-3,48])
      set(gca,'xtick',[-60 -30 0 30 60])
      set(gca,'ytick',[0 20 40])
    else
      xlabel('Sorry: normfit not available in Octave')
    end
    
    clear varl varr y

    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0.03 0 0]);

    for ii=1:length(h)
        set(h(ii),'Linewidth',1.1)
        set(h(ii),'TickLength',[0.015 0.015])
        tmp=findobj(h(ii),'Type','patch');
        set(tmp,'EdgeColor','k');
    end
    set(gcf,'color',[1 1 1])
    
end

%% Figure 8
if flags.do_fig8
    
    siglevel=0.05;
    outlierrate=0.01; % upper bound for the outlier rate

    Obj=data_ziegelwanger2014('STP',flags.cachemode);
    idx=Obj.Data.toaEst{4}(:,1)>-100000;
    stpEst=Obj.Data.toaEst{4}(idx,1);
    [~,results]=ziegelwanger2014(Obj,Obj.Data.toaEst{4},0,2);
    stpMod=results.toa(idx,1);
    stpRes=stpEst-stpMod;

    Obj=data_ziegelwanger2014('NH89',flags.cachemode);
    fs=Obj.Data.SamplingRate*1e-6;
    idx=Obj.Data.toaEst{4}(:,1)>-100000;
    sp=Obj.SourcePosition(idx,:);
    [lat,~]=local_geo2horpolar(sp(:,1),sp(:,2));
    nhEst=Obj.Data.toaEst{4}(idx,1);
    [~,results]=ziegelwanger2014(Obj,Obj.Data.toaEst{4},0,2);
    nhMod=results.toa(idx,1);
    nhRes=nhEst-nhMod;

    res=stpRes;
    est=stpEst;

    [~,idx]=deleteoutliers(res,siglevel*outlierrate*length(stpEst));

    figure('Position',[620   229   560*2   497],'resize','off');
    h(1)=subplot(2,2,1);
    box on; hold on;
    plot(lat,est/fs-min(est/fs),'.','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',8);
    plot(lat(idx),est(idx)/fs-min(est/fs),'.k','MarkerSize',24);
    axis([-99 99 -50 995]);
    set(gca,'XTickLabel',[],'FontName','Arial');
    ylabel('Relative TOA (\mus)','FontName','Arial','FontSize',18);
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0 0.06 0]);

    h(2)=subplot(2,2,3);
    box on; hold on;
    plot(lat,res/fs,'.','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',8)
    plot(lat(idx),res(idx)/fs,'.k','MarkerSize',24);
    line([-99 99],1./[-fs -fs],'LineStyle','--','Color',[0 0 0]);
    line([-99 99],1./[fs fs],'LineStyle','--','Color',[0 0 0]);
    axis([-99 99 -80 195])
    xlabel('\Phi_h (deg)','FontName','Arial','FontSize',18);
    ylabel('Residual error (\mus)','FontName','Arial','FontSize',18);
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0.11 0.06 0]);

    res=nhRes;
    est=nhEst;

    [~,idx]=deleteoutliers(res,siglevel*outlierrate*length(nhEst));

    h(3)=subplot(2,2,2);
    box on; hold on;
    plot(lat,est/fs-min(est/fs),'.','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',8);
    plot(lat(idx),est(idx)/fs-min(est/fs),'.k','MarkerSize',24);
    axis([-99 99 -50 995]);
    set(gca,'XTickLabel',[],'FontName','Arial');
    set(gca,'yticklabel',[]);
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[-0.03 0 0.06 0]);

    h(4)=subplot(2,2,4);
    box on; hold on;
    plot(lat,res/fs,'.','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',8)
    plot(lat(idx),res(idx)/fs,'.k','MarkerSize',24);
    line([-99 99],1./[-fs -fs],'LineStyle','--','Color',[0 0 0]);
    line([-99 99],1./[fs fs],'LineStyle','--','Color',[0 0 0]);
    axis([-99 99 -80 195])
    xlabel('\Phi_h (deg)','FontName','Arial','FontSize',18);
    set(gca,'yticklabel',[]);
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[-0.03 0.11 0.06 0]);

    for ii=1:length(h)
        set(h(ii),'linewidth',1.1)
        set(h(ii),'Fontname','Arial','Fontsize',14)
    end
    set(gcf,'color',[1 1 1])
    
end

%% Figure 9
if flags.do_fig9
    
%load data
    Obj=data_ziegelwanger2014('NH89',flags.cachemode);
    toa1=Obj.Data.toaEst{4};
    [~,tmp]=ziegelwanger2014(Obj,Obj.Data.toaEst{4},0,1);
    toa2=tmp.toa;

%plot figure
    h=subplot(1,1,1);
    plot_ziegelwanger2014(Obj,toa2-min(min(toa1))-1,4,'k',0,1,1,'-',3);
    h2=plot_ziegelwanger2014(Obj,toa1-min(min(toa1))-1,4,'k',0,1,1,'o',1);
    set(h2,'MarkerFaceColor',[0 0 0]);
    h3=plot_ziegelwanger2014(Obj,toa2-min(min(toa1))-1,4,'k',0,2,1,'-',3);
    set(h3,'Color',[0.5 0.5 0.5]);
    h4=plot_ziegelwanger2014(Obj,toa1-min(min(toa1))-1,4,'k',0,2,1,'o',1);
    set(h4,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor',[.5 .5 .5]);

    xlim([-9 369])
    ylim([-0.04 0.98])
    grid off
    xlabel(' ');
    ylabel('Relative TOA (ms)','Fontname','Arial','Fontsize',22)
    title('')
    set(gca,'xticklabel',{'0' '90' '180' '270' '360   '});
    set(gca,'ytick',[0 0.2 0.4 0.6 0.8])
    set(gca,'Fontname','Arial','fontsize',18)
    set(h,'linewidth',2)
    
end

%% Figure 10
if flags.do_fig10

%load data
    data=data_ziegelwanger2014('SPHERE_DIS',flags.cachemode);
    for ii=1:length(data.results)
        p_onaxis(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
    end
    p1=p_onaxis(:,:,[1:3 length(data.xM)/3+1:length(data.xM)/3+3 length(data.xM)/3*2+1:length(data.xM)/3*2+3]);
    r1=data.radius([1:3 length(data.xM)/3+1:length(data.xM)/3+3 length(data.xM)/3*2+1:length(data.xM)/3*2+3]);
    yM1=data.yM([1:3 length(data.xM)/3+1:length(data.xM)/3+3 length(data.xM)/3*2+1:length(data.xM)/3*2+3]);

%plot figure
    figure
    set(gcf,'Position',[0 0 560 285]);
    sym='ox';%plot symbols
    ms=8;%markersize
    lw1=2;%linewidth
    lw2=2;%linewidth
    fs=18;%fontsize
    ls='k-';%linestyle
    lc=[0.5 0.5 0.5];
    h=[];
    
    %radii
    h(end+1)=subplot(211);
    var=[squeeze(p1(1,1,:))*1000 squeeze(p1(1,2,:))*1000 r1];
    for ch=1:size(p1,2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
        hold on
    end
    plot(1:3,r1(1:3),ls,'Linewidth',lw2,'color',lc)
    hold on
    plot(4:6,r1(4:6),ls,'Linewidth',lw2,'color',lc)
    plot(7:9,r1(7:9),ls,'Linewidth',lw2,'color',lc)
    set(gca,'xtick',1:9)
    set(gca,'xticklabel',[])
    clear var;
    ylabel('r (mm)','Fontname','Arial','Fontsize',fs)
    ylim([51,129])
    xlim([0.5,9.5])

    %yM
    h(end+1)=subplot(212);
    plot(1:3,-yM1(1:3)*1000,ls,'Linewidth',lw2,'color',lc)
    hold on
    plot(4:6,-yM1(4:6)*1000,ls,'Linewidth',lw2,'color',lc)
    plot(7:9,-yM1(7:9)*1000,ls,'Linewidth',lw2,'color',lc)
    set(gca,'xtick',1:9)
    set(gca,'xticklabel',[])
    clear var;
    xlabel('Condition','Fontname','Arial','Fontsize',fs)
    ylabel('y_M (mm) ','Fontname','Arial','Fontsize',fs)
    ylim([-5 25])
    xlim([0.5,9.5])
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0.02 0 0]);

    for ii=1:length(h)
        set(h(ii),'fontsize',14)
        set(h(ii),'Linewidth',2)
        tmp=findobj(h(ii),'Type','patch');
        set(tmp,'EdgeColor','k');
    end
    set(gcf,'color',[1 1 1])
    
end

%% Figure 12
if flags.do_fig12

%load data
    data=data_ziegelwanger2014('SPHERE_DIS',flags.cachemode);
    for ii=1:length(data.results)
        p_offaxis(:,:,ii)=data.results(ii).MCM{1}.p_offaxis;
    end

    center=[data.xM(1:length(data.xM)/3) data.yM(1:length(data.yM)/3) data.zM(1:length(data.zM)/3)];
    [~,idx3]=sort(squeeze(center(:,3)));
    [~,idx2]=sort(squeeze(center(:,1)));
    [~,idx1]=sort(squeeze(center(:,2)));
    idx=idx3(idx2(idx1));
    idx=[idx; idx+length(data.xM)/3; idx+length(data.xM)/3*2];
    data.radius=data.radius(idx);
    p_offaxis=p_offaxis(:,:,idx);
    data.xM=data.xM(idx);
    data.yM=data.yM(idx);
    data.zM=data.zM(idx);

%plot figure
    figure
    tmp=get(gcf,'Position');
    set(gcf,'Position',tmp.*[1 1 1.5 2]);
    sym='ox';%plot symbols
    ms=8;%markersize
    lw1=2;%linewidth
    lw2=2;%linewidth
    fs=18;%fontsize
    ls='k-';%linestyle
    lc=[0.5 0.5 0.5];
    h=[];

    %radii
    h(end+1)=subplot(611);
    var=[squeeze(p_offaxis(1,1,:))*1000 squeeze(p_offaxis(1,2,:))*1000 data.radius];
    plot(var(:,3),ls,'Linewidth',lw2,'color',lc)
    hold on
    for ch=1:size(p_offaxis,2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
    end
    ylabel('r (mm)','Fontname','Arial','Fontsize',fs)
    xlim([-1 size(var,1)+2])
    ylim([72 108])
    set(gca,'xtick',1:size(var,1))
    set(gca,'xticklabel',[])
    clear var;

    %xM
    h(end+1)=subplot(612);
    var=[squeeze(p_offaxis(2,1,:))*1000 squeeze(p_offaxis(2,2,:))*1000 -data.xM*1000];
    plot(var(:,3),ls,'Linewidth',lw2,'color',lc)
    hold on
    for ch=1:size(p_offaxis,2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
    end  
    ylabel('x_M (mm) ','Fontname','Arial','Fontsize',fs)
    xlim([-1 size(var,1)+2])
    ylim([-28 8])
    set(gca,'xtick',1:size(var,1))
    set(gca,'xticklabel',[])
    clear var;
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0.02 0 0]);

    %yM
    h(end+1)=subplot(613);
    var=[squeeze(p_offaxis(3,1,:))*1000 squeeze(p_offaxis(3,2,:))*1000 -data.yM*1000];
    plot(var(:,3),ls,'Linewidth',lw2,'color',lc)
    hold on
    for ch=1:size(p_offaxis,2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
    end
    ylabel('y_M (mm) ','Fontname','Arial','Fontsize',fs)
    xlim([-1 size(var,1)+2])
    ylim([-8 28])
    set(gca,'xtick',1:size(var,1))
    set(gca,'xticklabel',[])
    clear var;
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0.04 0 0]);

    %zM
    h(end+1)=subplot(614);
    var=[squeeze(p_offaxis(4,1,:))*1000 squeeze(p_offaxis(4,2,:))*1000 -data.zM*1000];
    plot([1 size(var,1)],var([1 3],3),ls,'Linewidth',lw2,'color',lc)
    hold on
    plot([1 size(var,1)],var([2 4],3),ls,'Linewidth',lw2,'color',lc)
    for ch=1:size(p_offaxis,2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
    end
    % set(gca,'YTick',[-10 0])
    ylabel('z_M (mm)','Fontname','Arial','Fontsize',fs)
    xlim([-1 size(var,1)+2])
    ylim([-14 4])
    set(gca,'xtick',1:size(var,1))
    set(gca,'xticklabel',[])
    clear var;
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0.06 0 0]);

    %phi
    h(end+1)=subplot(615);
    var=[squeeze(p_offaxis(6,1,:))/pi*180 squeeze(p_offaxis(6,2,:))/pi*180 ones(length(data.zM),1)*90 -ones(length(data.zM),1)*90];
    for ch=1:size(p_offaxis,2)
        plot(abs(var(:,2+ch)),ls,'Linewidth',lw2,'color',lc)
        hold on
    end
    for ch=1:size(p_offaxis,2)
        plot(abs(var(:,ch)),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
    end
    set(gca,'YTick',[-90 90])
    ylabel('   _e (deg)','Fontname','Arial','Fontsize',fs)
    xlim([-1 size(var,1)+2])
    ylim([82,98])
    set(gca,'xtick',1:size(var,1))
    set(gca,'xticklabel',[])
    set(gca,'ytick',[85 90 95])
    set(gca,'yticklabel',{'\pm85','\pm90','\pm95'})
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0.08 0 0]);
    clear var;

    %theta
    h(end+1)=subplot(616);
    var=[squeeze(p_offaxis(7,1,:))/pi*180 squeeze(p_offaxis(7,2,:))/pi*180 zeros(length(data.zM),1) zeros(length(data.zM),1)];
    for ch=1:size(p_offaxis,2)
        plot(var(:,2+ch),ls,'Linewidth',lw2,'color',lc)
        hold on
    end
    for ch=1:size(p_offaxis,2)
        plot(var(:,ch),['k' sym(ch)],'Markersize',ms,'Linewidth',lw1);
    end
    ylabel('\theta_e (deg)','Fontname','Arial','Fontsize',fs)
    xlabel('Condition','Fontname','Arial','Fontsize',fs)
    ylim([-7,7])
    xlim([-1 size(var,1)+2])
    set(gca,'xtick',1:size(var,1))
    set(gca,'xticklabel',[])
    set(gca,'ytick',[-5 0 5])
    tmp=get(gca,'Position');
    set(gca,'Position',tmp+[0 0.10 0 0]);
    clear var;

    for ii=1:length(h)
        set(h(ii),'fontsize',14)
        set(h(ii),'Linewidth',2)
        tmp=findobj(h(ii),'Type','patch');
        set(tmp,'EdgeColor','k');
    end
    set(gcf,'color',[1 1 1])
    
end

%% Table 1
if flags.do_tab1
    
%load data
    data=data_ziegelwanger2014('SPHERE_ROT',flags.cachemode);
    for ii=1:length(data.results)
        p_onaxis{1}(:,:,ii)=data.results(ii).MAX{1}.p_onaxis;
        p_onaxis{2}(:,:,ii)=data.results(ii).CTD{1}.p_onaxis;
        p_onaxis{3}(:,:,ii)=data.results(ii).AGD{1}.p_onaxis;
        p_onaxis{4}(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
    end
    resnormS=zeros(length(data.results),4);
    for ii=1:length(data.results)
        resnormS(ii,1)=mean([data.results(ii).MAX{1}.performance.on_axis{1}.resnormS ...
            data.results(ii).MAX{1}.performance.on_axis{2}.resnormS]);
        resnormS(ii,2)=mean([data.results(ii).CTD{1}.performance.on_axis{1}.resnormS ...
            data.results(ii).CTD{1}.performance.on_axis{2}.resnormS]);
        resnormS(ii,3)=mean([data.results(ii).AGD{1}.performance.on_axis{1}.resnormS ...
            data.results(ii).AGD{1}.performance.on_axis{2}.resnormS]);
        resnormS(ii,4)=mean([data.results(ii).MCM{1}.performance.on_axis{1}.resnormS ...
            data.results(ii).MCM{1}.performance.on_axis{2}.resnormS]);
    end
    
    out=zeros(4,8);
    out(1,7)=mean(resnormS(:,1))*1e6;
    out(2,7)=mean(resnormS(:,2))*1e6;
    out(3,7)=mean(resnormS(:,3))*1e6;
    out(4,7)=mean(resnormS(:,4))*1e6;
    out(1,8)=std(resnormS(:,1))*1e6;
    out(2,8)=std(resnormS(:,2))*1e6;
    out(3,8)=std(resnormS(:,3))*1e6;
    out(4,8)=std(resnormS(:,4))*1e6;

    idx=1:size(p_onaxis{1},3);
    % radii
    varl=[[squeeze(p_onaxis{1}(1,1,idx)) ...
        squeeze(p_onaxis{2}(1,1,idx))...
        squeeze(p_onaxis{3}(1,1,idx)) ...
        squeeze(p_onaxis{4}(1,1,idx))]*1000 ...
        data.radius(idx)];
    varr=[[squeeze(p_onaxis{1}(1,2,idx)) ...
        squeeze(p_onaxis{2}(1,2,idx))...
        squeeze(p_onaxis{3}(1,2,idx)) ...
        squeeze(p_onaxis{4}(1,2,idx))]*1000 ...
        data.radius(idx)];
    err=abs([varl(:,1:4); varr(:,1:4)]-[varl(:,5); varr(:,5)]*[1 1 1 1]);
    out(1,1)=mean(err(:,1));
    out(2,1)=mean(err(:,2));
    out(3,1)=mean(err(:,3));
    out(4,1)=mean(err(:,4));
    out(1,2)=std(err(:,1));
    out(2,2)=std(err(:,2));
    out(3,2)=std(err(:,3));
    out(4,2)=std(err(:,4));

    %phi
    varl=[[squeeze(p_onaxis{1}(2,1,idx)) ...
        squeeze(p_onaxis{2}(2,1,idx))...
        squeeze(p_onaxis{3}(2,1,idx)) ...
        squeeze(p_onaxis{4}(2,1,idx))]/pi*180 ...
        data.phi(idx)+ones(length(data.phi(idx)),1)*90];
    varr=[mod([squeeze(p_onaxis{1}(2,2,idx)) ...
        squeeze(p_onaxis{2}(2,2,idx))...
        squeeze(p_onaxis{3}(2,2,idx)) ...
        squeeze(p_onaxis{4}(2,2,idx))]/pi*180,360) ...
        mod(data.phi(idx)-ones(length(data.phi(idx)),1)*90,360)];
    err=abs([varl(:,1:4); varr(:,1:4)]-[varl(:,5); varr(:,5)]*[1 1 1 1]);
    out(1,3)=mean(err(:,1));
    out(2,3)=mean(err(:,2));
    out(3,3)=mean(err(:,3));
    out(4,3)=mean(err(:,4));
    out(1,4)=std(err(:,1));
    out(2,4)=std(err(:,2));
    out(3,4)=std(err(:,3));
    out(4,4)=std(err(:,4));

    %theta
    varl=[[squeeze(p_onaxis{1}(3,1,idx)) ...
        squeeze(p_onaxis{2}(3,1,idx))...
        squeeze(p_onaxis{3}(3,1,idx)) ...
        squeeze(p_onaxis{4}(3,1,idx))]/pi*180 ...
        data.theta(idx)];
    varr=[[squeeze(p_onaxis{1}(3,2,idx)) ...
        squeeze(p_onaxis{2}(3,2,idx))...
        squeeze(p_onaxis{3}(3,2,idx)) ...
        squeeze(p_onaxis{4}(3,2,idx))]/pi*180 ...
        -data.theta(idx)];
    err=abs([varl(:,1:4); varr(:,1:4)]-[varl(:,5); varr(:,5)]*[1 1 1 1]);
    out(1,5)=mean(err(:,1));
    out(2,5)=mean(err(:,2));
    out(3,5)=mean(err(:,3));
    out(4,5)=mean(err(:,4));
    out(1,6)=std(err(:,1));
    out(2,6)=std(err(:,2));
    out(3,6)=std(err(:,3));
    out(4,6)=std(err(:,4));

%display data
    rows=['MAX ||';...
        'CTD ||';...
        'AGD ||';...
        'MCM ||'];...
    fprintf('\nTab. I.:\n')
    fprintf('----------------------------------------------------------------------\n')
    fprintf('EST || r error (mm) | phi_e error (deg) | theta_e errror |  ANR (\mus)  \n')
    fprintf('----------------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------------\n')
    for ii=1:4
        fprintf(rows(ii,:))
        fprintf('   %4.1f \pm %3.1f |        %4.1f \pm %3.1f |     %4.1f \pm %3.1f | %4.1f \pm %3.1f\n',out(ii,:)');
    end
    fprintf('----------------------------------------------------------------------\n\n')

end

%% Table 2
if flags.do_tab2

%load data
    out=zeros(3,16);
    Obj1=data_ziegelwanger2014('Sphere',flags.cachemode);
    for method=1:4
        [Obj1,results]=ziegelwanger2014(Obj1,Obj1.Data.toaEst{method},0,1,[[0.08; pi/18*8; 0; 0] [0.08; -pi/18*10; 0; 0]]);
        out(1,method+0)=results.p_onaxis(1,1)*1000;
        out(1,method+4)=results.p_onaxis(2,1)*180/pi;
        out(1,method+8)=results.p_onaxis(3,1)*180/pi-1;
        out(1,method+12)=results.performance.on_axis{1}.resnormS*1e06;
    end
    Obj2=data_ziegelwanger2014('SAT',flags.cachemode);
    for method=1:4
        [Obj2,results]=ziegelwanger2014(Obj2,Obj2.Data.toaEst{method},0,1,[[0.08; pi/18*8; 0; 0] [0.08; -pi/18*10; 0; 0]]);
        out(2,method+0)=results.p_onaxis(1,1)*1000;
        out(2,method+4)=results.p_onaxis(2,1)*180/pi;
        out(2,method+8)=results.p_onaxis(3,1)*180/pi-1;
        out(2,method+12)=results.performance.on_axis{1}.resnormS*1e06;
    end
    Obj3=data_ziegelwanger2014('STP',flags.cachemode);
    for method=1:4
        [Obj3,results]=ziegelwanger2014(Obj3,Obj3.Data.toaEst{method},0,1,[[0.08; pi/18*8; 0; 0] [0.08; -pi/18*10; 0; 0]]);
        out(3,method+0)=results.p_onaxis(1,1)*1000;
        out(3,method+4)=results.p_onaxis(2,1)*180/pi;
        out(3,method+8)=results.p_onaxis(3,1)*180/pi-1;
        out(3,method+12)=results.performance.on_axis{1}.resnormS*1e06;
    end

%display data
    rows=['Sphere ||';...
        'SAT    ||';...
        'STP    ||'];
    fprintf('\nTab. II.:\n')
    fprintf('----------------------------------------------------------------------------------------------------------\n')
    fprintf('       ||     r error (mm)       |   phi_e error (deg)   |    theta_e errror     |       ANR (\mus)        \n')
    fprintf('       || MAX | CTD | AGD  | MCM | MAX | CTD | AGD | MCM | MAX | CTD | AGD | MCM | MAX | CTD | AGD | MCM \n')
    fprintf('----------------------------------------------------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------------------------------------------------\n')
    for ii=1:3
        fprintf(rows(ii,:))
        fprintf(' %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f| %4.1f\n',out(ii,:)');
    end
    fprintf('----------------------------------------------------------------------------------------------------------\n\n')

end

%% Table 3
if flags.do_tab3

%load data
    out=zeros(4,8);
    hrtf={'ARI','CIPIC','LISTEN'};
    for kk=1:length(hrtf)
        data=data_ziegelwanger2014(hrtf{kk},flags.cachemode);
        if kk==3
            data.results=data.results([1:27 29:end]);
        end
        for ii=1:length(data.results)
            temp1(:,:,ii)=data.results(ii).MCM{1}.p_onaxis;
            temp3(ii)=mean([data.results(ii).MCM{1}.performance.on_axis{1}.resnormS ...
                data.results(ii).MCM{1}.performance.on_axis{2}.resnormS]);
            temp5(ii)=mean([data.results(ii).MCM{1}.performance.on_axis{1}.resnormP ...
                data.results(ii).MCM{1}.performance.on_axis{2}.resnormP]);
            temp8(ii)=data.results(ii).MCM{1}.performance.on_axis{1}.resnormS;
            temp9(ii)=data.results(ii).MCM{1}.performance.on_axis{2}.resnormS;
        end
        p_onaxis{kk}=temp1;
        resnormS_onaxis{kk}=temp3;
        resnormP_onaxis{kk}=temp5;
        resnormS_onaxis_left{kk}=temp8;
        resnormS_onaxis_right{kk}=temp9;
        clear temp1 temp3 temp5 temp8 temp9
    end

    % radii
    out(1,1)=mean([mean(squeeze(p_onaxis{1}(1,1,:)*1000)) ...
        mean(squeeze(p_onaxis{2}(1,1,:)*1000)) ...
        mean(squeeze(p_onaxis{3}(1,1,:)*1000))]);
    out(1,2)=std([reshape(p_onaxis{1}(1,1,:)*1000,1,numel(p_onaxis{1}(1,1,:))) ...
        reshape(p_onaxis{2}(1,1,:)*1000,1,numel(p_onaxis{2}(1,1,:))) ...
        reshape(p_onaxis{3}(1,1,:)*1000,1,numel(p_onaxis{3}(1,1,:)))]);
    out(2,1)=mean([mean(squeeze(p_onaxis{1}(1,2,:)*1000)) ...
        mean(squeeze(p_onaxis{2}(1,2,:)*1000)) ...
        mean(squeeze(p_onaxis{3}(1,2,:)*1000))]);
    out(2,2)=std([reshape(p_onaxis{1}(1,2,:)*1000,1,numel(p_onaxis{1}(1,2,:))) ...
        reshape(p_onaxis{2}(1,2,:)*1000,1,numel(p_onaxis{2}(1,2,:))) ...
        reshape(p_onaxis{3}(1,2,:)*1000,1,numel(p_onaxis{3}(1,2,:)))]);

    % phi_e
    out(1,3)=mean([mean(squeeze(p_onaxis{1}(2,1,:)*180/pi)) ...
        mean(squeeze(p_onaxis{2}(2,1,:)*180/pi)) ...
        mean(squeeze(p_onaxis{3}(2,1,:)*180/pi))]);
    out(1,4)=std([reshape(p_onaxis{1}(2,1,:)*180/pi,1,numel(p_onaxis{1}(2,1,:))) ...
        reshape(p_onaxis{2}(2,1,:)*180/pi,1,numel(p_onaxis{2}(2,1,:))) ...
        reshape(p_onaxis{3}(2,1,:)*180/pi,1,numel(p_onaxis{3}(2,1,:)))]);
    out(2,3)=mean([mean(squeeze(p_onaxis{1}(2,2,:)*180/pi)) ...
        mean(squeeze(p_onaxis{2}(2,2,:)*180/pi)) ...
        mean(squeeze(p_onaxis{3}(2,2,:)*180/pi))]);
    out(2,4)=std([reshape(p_onaxis{1}(2,2,:)*180/pi,1,numel(p_onaxis{1}(2,2,:))) ...
        reshape(p_onaxis{2}(2,2,:)*180/pi,1,numel(p_onaxis{2}(2,2,:))) ...
        reshape(p_onaxis{3}(2,2,:)*180/pi,1,numel(p_onaxis{3}(2,2,:)))]);

    % theta_e
    out(1,5)=mean([mean(squeeze(p_onaxis{1}(3,1,:)*180/pi)) ...
        mean(squeeze(p_onaxis{2}(3,1,:)*180/pi)) ...
        mean(squeeze(p_onaxis{3}(3,1,:)*180/pi))]);
    out(1,6)=std([reshape(p_onaxis{1}(3,1,:)*180/pi,1,numel(p_onaxis{1}(3,1,:))) ...
        reshape(p_onaxis{2}(3,1,:)*180/pi,1,numel(p_onaxis{2}(3,1,:))) ...
        reshape(p_onaxis{3}(3,1,:)*180/pi,1,numel(p_onaxis{3}(3,1,:)))]);
    out(2,5)=mean([mean(squeeze(p_onaxis{1}(3,2,:)*180/pi)) ...
        mean(squeeze(p_onaxis{2}(3,2,:)*180/pi)) ...
        mean(squeeze(p_onaxis{3}(3,2,:)*180/pi))]);
    out(2,6)=std([reshape(p_onaxis{1}(3,2,:)*180/pi,1,numel(p_onaxis{1}(3,2,:))) ...
        reshape(p_onaxis{2}(3,2,:)*180/pi,1,numel(p_onaxis{2}(3,2,:))) ...
        reshape(p_onaxis{3}(3,2,:)*180/pi,1,numel(p_onaxis{3}(3,2,:)))]);

    out(1,7)=mean([resnormS_onaxis_left{1} resnormS_onaxis_left{2} resnormS_onaxis_left{3}])*1e6;
    out(1,8)=std([resnormS_onaxis_left{1} resnormS_onaxis_left{2} resnormS_onaxis_left{3}])*1e6;
    out(2,7)=mean([resnormS_onaxis_right{1} resnormS_onaxis_right{2} resnormS_onaxis_right{3}])*1e6;
    out(2,8)=std([resnormS_onaxis_right{1} resnormS_onaxis_right{2} resnormS_onaxis_right{3}])*1e6;

    Obj=data_ziegelwanger2014('NH89',flags.cachemode);
    [~,results]=ziegelwanger2014(Obj,4,0,1);

    out(3:4,[1 3 5])=results.p_onaxis(1:3,:)'.*([1;1]*[1000 180/pi 180/pi]);
    out(3,7)=results.performance.on_axis{1}.resnormS*1e6;
    out(4,7)=results.performance.on_axis{2}.resnormS*1e6;

%display data
    rows=['All  |  L  ||';...
        '     |  R  ||';...
        'NH89 |  L  ||';...
        '     |  R  ||'];...
    fprintf('\nTab. III.:\n')
    fprintf('-----------------------------------------------------------------------\n')
    fprintf('     | Ear || r (mm)       | phi_e (deg)  | theta_e (deg) |  ANR (\mus)  \n')
    fprintf('-----------------------------------------------------------------------\n')
    fprintf('-----------------------------------------------------------------------\n')
    for ii=1:4
        fprintf(rows(ii,:))
        fprintf(' %5.1f \pm %4.1f | %5.1f \pm %4.1f |   %5.1f \pm %3.1f | %4.1f \pm %4.1f\n',out(ii,:)');
    end
    fprintf('-----------------------------------------------------------------------\n\n')

end

%% Table 5
if flags.do_tab5
    
%load data
    out=zeros(12,14);
    Obj=data_ziegelwanger2014('SAT',flags.cachemode);
    [~,results]=ziegelwanger2014(Obj,Obj.Data.toaEst{4},0,1);
    out(1:2,[1 7 9 11 3 5])=results.p_offaxis([1:4 6:7],:)'.*([1;1]*[1000 1000 1000 1000 180/pi 180/pi]);
    out(1,13)=results.performance.off_axis{1}.resnormS*1e6;
    out(2,13)=results.performance.off_axis{2}.resnormS*1e6;

    [~,results]=ziegelwanger2014(Obj,Obj.Data.toaEst{4},[0.05 0.01],1);
    out(3:4,[1 7 9 11 3 5])=results.p_offaxis([1:4 6:7],:)'.*([1;1]*[1000 1000 1000 1000 180/pi 180/pi]);
    out(3,13)=results.performance.off_axis{1}.resnormS*1e6;
    out(4,13)=results.performance.off_axis{2}.resnormS*1e6;

    Obj=data_ziegelwanger2014('STP',flags.cachemode);
    [~,results]=ziegelwanger2014(Obj,Obj.Data.toaEst{4},0,1);
    out(5:6,[1 7 9 11 3 5])=results.p_offaxis([1:4 6:7],:)'.*([1;1]*[1000 1000 1000 1000 180/pi 180/pi]);
    out(5,13)=results.performance.off_axis{1}.resnormS*1e6;
    out(6,13)=results.performance.off_axis{2}.resnormS*1e6;

    [~,results]=ziegelwanger2014(Obj,Obj.Data.toaEst{4},[0.05 0.01],1);
    out(7:8,[1 7 9 11 3 5])=results.p_offaxis([1:4 6:7],:)'.*([1;1]*[1000 1000 1000 1000 180/pi 180/pi]);
    out(7,13)=results.performance.off_axis{1}.resnormS*1e6;
    out(8,13)=results.performance.off_axis{2}.resnormS*1e6;

    out(1:8,5)=out(1:8,5)-1;

    hrtf={'ARI','CIPIC','LISTEN'};
    for kk=1:length(hrtf)
        data=data_ziegelwanger2014(hrtf{kk},flags.cachemode);
        if kk==3
            data.results=data.results([1:27 29:end]);
        end
        for jj=1:2
            for ii=1:length(data.results)
                temp1(:,:,ii)=data.results(ii).MCM{jj}.p_offaxis;
                temp3(ii)=mean([data.results(ii).MCM{jj}.performance.off_axis{1}.resnormS ...
                    data.results(ii).MCM{jj}.performance.off_axis{2}.resnormS]);
                temp5(ii)=mean([data.results(ii).MCM{jj}.performance.off_axis{1}.resnormP ...
                    data.results(ii).MCM{jj}.performance.off_axis{2}.resnormP]);
                temp8(ii)=data.results(ii).MCM{jj}.performance.off_axis{1}.resnormS;
                temp9(ii)=data.results(ii).MCM{jj}.performance.off_axis{2}.resnormS;
            end
            p_offaxis{kk,jj}=temp1;
            resnormS_offaxis{kk,jj}=temp3;
            resnormP_offaxis{kk,jj}=temp5;
            resnormS_offaxis_left{kk,jj}=temp8;
            resnormS_offaxis_right{kk,jj}=temp9;
            clear temp1 temp3 temp5 temp8 temp9
        end
    end

    % radii
    for jj=1:2
        out(9+(jj-1)*2,1)=mean([mean(squeeze(p_offaxis{1,jj}(1,1,:)*1000)) ...
            mean(squeeze(p_offaxis{2,jj}(1,1,:)*1000)) ...
            mean(squeeze(p_offaxis{3,jj}(1,1,:)*1000))]);
        out(9+(jj-1)*2,2)=std([reshape(p_offaxis{1,jj}(1,1,:)*1000,1,numel(p_offaxis{1,jj}(1,1,:))) ...
            reshape(p_offaxis{2,jj}(1,1,:)*1000,1,numel(p_offaxis{2,jj}(1,1,:))) ...
            reshape(p_offaxis{3,jj}(1,1,:)*1000,1,numel(p_offaxis{3,jj}(1,1,:)))]);
        out(10+(jj-1)*2,1)=mean([mean(squeeze(p_offaxis{1,jj}(1,2,:)*1000)) ...
            mean(squeeze(p_offaxis{2,jj}(1,2,:)*1000)) ...
            mean(squeeze(p_offaxis{3,jj}(1,2,:)*1000))]);
        out(10+(jj-1)*2,2)=std([reshape(p_offaxis{1,jj}(1,2,:)*1000,1,numel(p_offaxis{1,jj}(1,2,:))) ...
            reshape(p_offaxis{2,jj}(1,2,:)*1000,1,numel(p_offaxis{2,jj}(1,2,:))) ...
            reshape(p_offaxis{3,jj}(1,2,:)*1000,1,numel(p_offaxis{3,jj}(1,2,:)))]);
    end

    % phi_e
    for jj=1:2
    out(9+(jj-1)*2,3)=mean([mean(squeeze(p_offaxis{1,jj}(6,1,:)*180/pi)) ...
        mean(squeeze(p_offaxis{2,jj}(6,1,:)*180/pi)) ...
        mean(squeeze(p_offaxis{3,jj}(6,1,:)*180/pi))]);
    out(9+(jj-1)*2,4)=std([reshape(p_offaxis{1,jj}(6,1,:)*180/pi,1,numel(p_offaxis{1,jj}(6,1,:))) ...
        reshape(p_offaxis{2,jj}(6,1,:)*180/pi,1,numel(p_offaxis{2,jj}(6,1,:))) ...
        reshape(p_offaxis{3,jj}(6,1,:)*180/pi,1,numel(p_offaxis{3,jj}(6,1,:)))]);
    out(10+(jj-1)*2,3)=mean([mean(squeeze(p_offaxis{1,jj}(6,2,:)*180/pi)) ...
        mean(squeeze(p_offaxis{2,jj}(6,2,:)*180/pi)) ...
        mean(squeeze(p_offaxis{3,jj}(6,2,:)*180/pi))]);
    out(10+(jj-1)*2,4)=std([reshape(p_offaxis{1,jj}(6,2,:)*180/pi,1,numel(p_offaxis{1,jj}(6,2,:))) ...
        reshape(p_offaxis{2,jj}(6,2,:)*180/pi,1,numel(p_offaxis{2,jj}(6,2,:))) ...
        reshape(p_offaxis{3,jj}(6,2,:)*180/pi,1,numel(p_offaxis{3,jj}(6,2,:)))]);
    end

    % theta_e
    for jj=1:2
    out(9+(jj-1)*2,5)=mean([mean(squeeze(p_offaxis{1,jj}(7,1,:)*180/pi)) ...
        mean(squeeze(p_offaxis{2,jj}(7,1,:)*180/pi)) ...
        mean(squeeze(p_offaxis{3,jj}(7,1,:)*180/pi))]);
    out(9+(jj-1)*2,6)=std([reshape(p_offaxis{1,jj}(7,1,:)*180/pi,1,numel(p_offaxis{1,jj}(7,1,:))) ...
        reshape(p_offaxis{2,jj}(7,1,:)*180/pi,1,numel(p_offaxis{2,jj}(7,1,:))) ...
        reshape(p_offaxis{3,jj}(7,1,:)*180/pi,1,numel(p_offaxis{3,jj}(7,1,:)))]);
    out(10+(jj-1)*2,5)=mean([mean(squeeze(p_offaxis{1,jj}(7,2,:)*180/pi)) ...
        mean(squeeze(p_offaxis{2,jj}(7,2,:)*180/pi)) ...
        mean(squeeze(p_offaxis{3,jj}(7,2,:)*180/pi))]);
    out(10+(jj-1)*2,6)=std([reshape(p_offaxis{1,jj}(7,2,:)*180/pi,1,numel(p_offaxis{1,jj}(7,2,:))) ...
        reshape(p_offaxis{2,jj}(7,2,:)*180/pi,1,numel(p_offaxis{2,jj}(7,2,:))) ...
        reshape(p_offaxis{3,jj}(7,2,:)*180/pi,1,numel(p_offaxis{3,jj}(7,2,:)))]);
    end

    % xM
    for jj=1:2
        out(9+(jj-1)*2,7)=mean([mean(squeeze(p_offaxis{1,jj}(2,1,:)*1000)) ...
            mean(squeeze(p_offaxis{2,jj}(2,1,:)*1000)) ...
            mean(squeeze(p_offaxis{3,jj}(2,1,:)*1000))]);
        out(9+(jj-1)*2,8)=std([reshape(p_offaxis{1,jj}(2,1,:)*1000,1,numel(p_offaxis{1,jj}(2,1,:))) ...
            reshape(p_offaxis{2,jj}(2,1,:)*1000,1,numel(p_offaxis{2,jj}(2,1,:))) ...
            reshape(p_offaxis{3,jj}(2,1,:)*1000,1,numel(p_offaxis{3,jj}(2,1,:)))]);
        out(10+(jj-1)*2,7)=mean([mean(squeeze(p_offaxis{1,jj}(2,2,:)*1000)) ...
            mean(squeeze(p_offaxis{2,jj}(2,2,:)*1000)) ...
            mean(squeeze(p_offaxis{3,jj}(2,2,:)*1000))]);
        out(10+(jj-1)*2,8)=std([reshape(p_offaxis{1,jj}(2,2,:)*1000,1,numel(p_offaxis{1,jj}(2,2,:))) ...
            reshape(p_offaxis{2,jj}(2,2,:)*1000,1,numel(p_offaxis{2,jj}(2,2,:))) ...
            reshape(p_offaxis{3,jj}(2,2,:)*1000,1,numel(p_offaxis{3,jj}(2,2,:)))]);
    end

    % yM
    for jj=1:2
        out(9+(jj-1)*2,9)=mean([mean(squeeze(p_offaxis{1,jj}(3,1,:)*1000)) ...
            mean(squeeze(p_offaxis{2,jj}(3,1,:)*1000)) ...
            mean(squeeze(p_offaxis{3,jj}(3,1,:)*1000))]);
        out(9+(jj-1)*2,10)=std([reshape(p_offaxis{1,jj}(3,1,:)*1000,1,numel(p_offaxis{1,jj}(3,1,:))) ...
            reshape(p_offaxis{2,jj}(3,1,:)*1000,1,numel(p_offaxis{2,jj}(3,1,:))) ...
            reshape(p_offaxis{3,jj}(3,1,:)*1000,1,numel(p_offaxis{3,jj}(3,1,:)))]);
        out(10+(jj-1)*2,9)=mean([mean(squeeze(p_offaxis{1,jj}(3,2,:)*1000)) ...
            mean(squeeze(p_offaxis{2,jj}(3,2,:)*1000)) ...
            mean(squeeze(p_offaxis{3,jj}(3,2,:)*1000))]);
        out(10+(jj-1)*2,10)=std([reshape(p_offaxis{1,jj}(3,2,:)*1000,1,numel(p_offaxis{1,jj}(3,2,:))) ...
            reshape(p_offaxis{2,jj}(3,2,:)*1000,1,numel(p_offaxis{2,jj}(3,2,:))) ...
            reshape(p_offaxis{3,jj}(3,2,:)*1000,1,numel(p_offaxis{3,jj}(3,2,:)))]);
    end

    % zM
    for jj=1:2
        out(9+(jj-1)*2,11)=mean([mean(squeeze(p_offaxis{1,jj}(4,1,:)*1000)) ...
            mean(squeeze(p_offaxis{2,jj}(4,1,:)*1000)) ...
            mean(squeeze(p_offaxis{3,jj}(4,1,:)*1000))]);
        out(9+(jj-1)*2,12)=std([reshape(p_offaxis{1,jj}(4,1,:)*1000,1,numel(p_offaxis{1,jj}(4,1,:))) ...
            reshape(p_offaxis{2,jj}(4,1,:)*1000,1,numel(p_offaxis{2,jj}(4,1,:))) ...
            reshape(p_offaxis{3,jj}(4,1,:)*1000,1,numel(p_offaxis{3,jj}(4,1,:)))]);
        out(10+(jj-1)*2,11)=mean([mean(squeeze(p_offaxis{1,jj}(4,2,:)*1000)) ...
            mean(squeeze(p_offaxis{2,jj}(4,2,:)*1000)) ...
            mean(squeeze(p_offaxis{3,jj}(4,2,:)*1000))]);
        out(10+(jj-1)*2,12)=std([reshape(p_offaxis{1,jj}(4,2,:)*1000,1,numel(p_offaxis{1,jj}(4,2,:))) ...
            reshape(p_offaxis{2,jj}(4,2,:)*1000,1,numel(p_offaxis{2,jj}(4,2,:))) ...
            reshape(p_offaxis{3,jj}(4,2,:)*1000,1,numel(p_offaxis{3,jj}(4,2,:)))]);
    end

    for jj=1:2
        out(9+(jj-1)*2,13)=mean([resnormS_offaxis_left{1,jj} resnormS_offaxis_left{2,jj} resnormS_offaxis_left{3,jj}])*1e6;
        out(9+(jj-1)*2,14)=std([resnormS_offaxis_left{1,jj} resnormS_offaxis_left{2,jj} resnormS_offaxis_left{3,jj}])*1e6;
        out(10+(jj-1)*2,13)=mean([resnormS_offaxis_right{1,jj} resnormS_offaxis_right{2,jj} resnormS_offaxis_right{3,jj}])*1e6;
        out(10+(jj-1)*2,14)=std([resnormS_offaxis_right{1,jj} resnormS_offaxis_right{2,jj} resnormS_offaxis_right{3,jj}])*1e6;
    end

%display data

    rows=['SAT |  Full  |  L  ||';...
        '    |  Full  |  R  ||';...
        '    |  O-A   |  L  ||';...
        '    |  O-A   |  R  ||';...
        'STP |  Full  |  L  ||';...
        '    |  Full  |  R  ||';...
        '    |  O-A   |  L  ||';...
        '    |  O-A   |  R  ||';...
        'ALL |  Full  |  L  ||';...
        '    |  Full  |  R  ||';...
        '    |  O-A   |  L  ||';...
        '    |  O-A   |  R  ||'];
    fprintf('\nTab. V.:\n')
    fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------\n')
    fprintf('    | TOAset | Ear || r (mm)     | phi_e (deg) | theta_e (deg)| x_M (mm)    | y_M (mm)    | z_M (mm)    |  ANR (\mus)  \n')
    fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------\n')
    fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------\n')
    for ii=1:12
        fprintf(rows(ii,:))
        fprintf(' %4.1f \pm %3.1f | %5.1f \pm %3.1f |  %5.1f \pm %4.1f | %4.1f \pm %4.1f | %4.1f \pm %4.1f | %4.1f \pm %4.1f | %4.1f \pm %4.1f\n',out(ii,:)');
    end
    fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------\n\n')

end

%% Table 6
if flags.do_tab6
    
%load data
    out=zeros(8,14);

    hrtf={'SPHERE_ROT','SPHERE_DIS'};
    for kk=1:length(hrtf)
        data=data_ziegelwanger2014(hrtf{kk},flags.cachemode);
        for jj=1:2
            for ii=1:length(data.results)
                temp1(:,:,ii)=data.results(ii).MCM{jj}.p_offaxis;
                temp3(ii)=mean([data.results(ii).MCM{jj}.performance.off_axis{1}.resnormS ...
                    data.results(ii).MCM{jj}.performance.off_axis{2}.resnormS]);
                temp8(ii)=data.results(ii).MCM{jj}.performance.off_axis{1}.resnormS;
                temp9(ii)=data.results(ii).MCM{jj}.performance.off_axis{2}.resnormS;
            end
            p_offaxis{kk,jj}=temp1;
            resnormS_offaxis{kk,jj}=temp3;
            resnormS_offaxis_left{kk,jj}=temp8;
            resnormS_offaxis_right{kk,jj}=temp9;
            if kk==1
                radius{kk,jj}=[data.radius data.radius];
                phi{kk,jj}=[90+data.phi 270+data.phi];
                theta{kk,jj}=[data.theta -data.theta];
                xM{kk,jj}=zeros(length(data.radius),2);
                yM{kk,jj}=zeros(length(data.radius),2);
                zM{kk,jj}=zeros(length(data.radius),2);
            else
                radius{kk,jj}=[data.radius data.radius];
                phi{kk,jj}=ones(length(data.radius),1)*[90 270];
                theta{kk,jj}=zeros(length(data.radius),2);
                xM{kk,jj}=-[data.xM data.xM];
                yM{kk,jj}=-[data.yM data.yM];
                zM{kk,jj}=-[data.zM data.zM];
            end
        end
        clear data temp1 temp3 temp8 temp9
    end

    for kk=1:2
        for jj=1:2
            for ch=1:2
                err=squeeze(p_offaxis{kk,jj}(1,ch,:))*1000-radius{kk,jj}(:,ch);
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),1)=mean(err);
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),2)=std(err);

                err=mod(squeeze(p_offaxis{kk,jj}(6,ch,:))*180/pi,360)-phi{kk,jj}(:,ch);
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),3)=mean(err);
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),4)=std(err);

                err=squeeze(p_offaxis{kk,jj}(7,ch,:))*180/pi-theta{kk,jj}(:,ch);
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),5)=mean(err);
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),6)=std(err);

                err=(squeeze(p_offaxis{kk,jj}(2,ch,:))-xM{kk,jj}(:,ch))*1000;
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),7)=mean(err);
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),8)=std(err);

                err=(squeeze(p_offaxis{kk,jj}(3,ch,:))-yM{kk,jj}(:,ch))*1000;
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),9)=mean(err);
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),10)=std(err);

                err=(squeeze(p_offaxis{kk,jj}(4,ch,:))-zM{kk,jj}(:,ch))*1000;
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),11)=mean(err);
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),12)=std(err);

                out(1+(kk-1)*4+(jj-1)*2+(ch-1),13)=mean(resnormS_offaxis_left{kk,jj})*1e6;
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),14)=std(resnormS_offaxis_left{kk,jj})*1e6;
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),13)=mean(resnormS_offaxis_right{kk,jj})*1e6;
                out(1+(kk-1)*4+(jj-1)*2+(ch-1),14)=std(resnormS_offaxis_right{kk,jj})*1e6;
            end
        end
    end

%display data
    rows=['Centered     |  Full  |  L  ||';...
        '             |  Full  |  R  ||';...
        '             |  O-A   |  L  ||';...
        '             |  O-A   |  R  ||';...
        'Non-centered |  Full  |  L  ||';...
        '             |  Full  |  R  ||';...
        '             |  O-A   |  L  ||';...
        '             |  O-A   |  R  ||'];
    fprintf('\nTab. V.:\n')
    fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------\n')
    fprintf('Condition    | TOAset | Ear || r error (mm) | phi_e error (deg) | theta_e errror | x_M error (mm) | y_M error (mm) | z_M error (mm) |  ANR (\mus)  \n')
    fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------\n')
    fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------\n')
    for ii=1:8
        fprintf(rows(ii,:))
        fprintf('   %4.1f \pm %3.1f |        %4.1f \pm %3.1f |     %4.1f \pm %3.1f |      %4.1f \pm %3.1f |     %4.1f \pm %3.1f |     %4.1f \pm %3.1f | %4.1f \pm %3.1f\n',out(ii,:)');
    end
    fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------\n\n')

end

end

%function idx=ARI_FindPosition(data,azimuth,elevation)
%    psi=sin(elevation/180*pi).*sin(data.APV(:,2)/180*pi) + ...
%        cos(elevation/180*pi).*cos(data.APV(:,2)/180*pi).*...
%        cos(azimuth/180*pi-data.APV(:,1)/180*pi);
%    [~,idx]=min(acos(psi));
%end

function [lat,pol]=local_geo2horpolar(azi,ele)
    warning('off');

    azi=mod(azi+360,360);
    ele=mod(ele+360,360);

    razi = deg2rad(azi);
    rele = deg2rad(ele);
    rlat=asin(sin(razi).*cos(rele));
    rpol=asin(sin(rele)./cos(rlat));
    idx=find(cos(rlat)==0);
    rpol(idx)=0;
    pol = rad2deg(rpol);
    lat = rad2deg(rlat);

    idx = find(razi>pi/2 & razi < 3*pi/2 & (rele < pi/2 | rele > 3*pi/2));
    pol(idx)=180-pol(idx);
    idx = find(~(razi>pi/2 & razi < 3*pi/2) & rele > pi/2 & rele < 3*pi/2);
    pol(idx)=180-pol(idx);
end

