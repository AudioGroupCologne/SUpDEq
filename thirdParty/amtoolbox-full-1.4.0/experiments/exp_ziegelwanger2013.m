function varargout=exp_ziegelwanger2013(varargin)
%EXP_ZIEGELWANGER2013   Figures from Ziegelwanger and Majdak (2013)
%   Usage: data = exp_ziegelwanger2013(flag)
%
%   EXP_ZIEGELWANGER2013(flags) reproduces figures of the paper from
%   Ziegelwanger and Majdak (2013).
%
%   Optional fields of output data structure:
%
%   The following flags can be specified:
%
%     'plot'       Plot the output of the experiment. This is the default.
%
%     'no_plot'    Don't plot, print results in the console.
%
%     'cached'     Reload previously calculated results. This is the default.
%
%     'redo'       Recalculate results.
%
%
%     'fig1a'    Reproduce Fig. 1 from ziegelwanger2013:
%                Model parameters (sphere radius, sphere center) resulting
%                from fitting the off-axis model to HRTFs of human
%                listeners.
%
%     'fig1b'    Reproduce Fig. 1 from majdak2013:
%               
%                Left panel: 
%                Normalized HRIRs of NH89 (ARI database). Sound source was
%                placed 45 deg left in the horizontal plane.
%                Solid line: for the left ear
%                Dashed line: for the right ear
%                Vertical lines: Estimated TOAs
%                Arrows: Resulting ITDs from the corresponding TOAs
%                Circles, Red: Minimum-Phase Cross-Correlation Method
%                Triangle, Green: Time-Position of the HRIR-Maximum
%                Diamonds, Blue: Centroid of the HRIR
%                Squares, Magenta: Average Group Delay (1-5 kHz)
%               
%                Right panel: 
%                Estimated TOAs of NH89 (ARI database) in the horizontal
%                interaural plane.
%                Black: Minimum-Phase Cross-Correlation Method
%                Blue: Time-Position of the HRIR-Maximum
%                Green: Centroid of the HRIR
%                Red: Average Group Delay (1-5 kHz)
%
%     'fig2b'    Reproduce Fig. 2 from majdak2013:
%               
%                Right panel: 
%                Estimated TOAs, detected outliers and outlier adjusted set
%                of TOAs for NH89 (ARI database) in the horizontal plane.
%                Line: Estimated TOAs
%                Triangles (down), Blue: Detected outliers for the azimuthal
%                slope criterion
%                Triangles (up), Blue: Detected outliers for the sagittal
%                variance criterion
%                Dots, Red: Outlier-adjusted set
%
%     'fig3b'    Reproduce Fig. 3 from majdak2013:
%
%                Left panel:
%                Model parameters (sphere radius, ear position) resulting
%                from fitting the on-axis model to HRTFs of a rigid sphere.
%                Squares: for the left ear
%                Diamonds: for the right ear
%                Dashed lines: set values used in the numerical HRTF
%                calculation
%
%                Right panel:
%                Model parameters (sphere radius) resulting from fitting the
%                on-axis model to HRTFs of an off-axis placed rigid sphere.
%
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Optimization Toolbox for Matlab
%
%   3) Data in hrtf/ziegelwanger2013
%
%   Examples:
%   ---------
%
%   To display Fig. 1a, use :
%
%     exp_ziegelwanger2013('fig1a');
%
%   To display Fig. 1b, use :
%
%     exp_ziegelwanger2013('fig1b');
%
%   To display Fig. 2b, use :
%
%     exp_ziegelwanger2013('fig2b');
%
%   To display Fig. 3b, use :
%
%     exp_ziegelwanger2013('fig3b');
%
%   See also: ziegelwanger2013, ziegelwanger2013_onaxis,
%   ziegelwanger2013_offaxis, data_ziegelwanger2013
%
%   References:
%     P. Majdak and H. Ziegelwanger. Continuous-direction model of the
%     broadband time-of-arrival in the head-related transfer functions. In
%     ICA 2013 Montreal, volume 19, page 050016, Montreal, Canada, 2013. ASA.
%     
%     H. Ziegelwanger and P. Majdak. Modeling the broadband time-of-arrival
%     of the head-related transfer functions for binaural audio. In
%     Proceedings of the 134th Convention of the Audio Engineering Society,
%     page 7, Rome, 2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_ziegelwanger2013.php


%   #Author: Harald Ziegelwanger (2013)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% ------ Check input options --------------------------------------------

    definput.flags.type = {'missingflag',...
    'fig1a','fig1b','fig2b','fig3b'};
    definput.flags.plot = {'plot','no_plot'};
    definput.import={'amt_cache'}; % get the flags of amt_cache
    
    % Parse input options
    [flags,kv]  = ltfatarghelper({},definput,varargin);

    if flags.do_missingflag
        flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}), ...
            sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
        error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
    end;

%% Figure 1b
if flags.do_fig1b
    
    data=data_ziegelwanger2013('NH89',flags.cachemode);

    subplot(122)
    %---------------------------Threshold---------------------------
    [~,tmp]=ziegelwanger2013(data,1,0);
    MAX=tmp.toa;

    %---------------------------Centroid----------------------------
    [~,tmp]=ziegelwanger2013(data,2,0);
    CTD=tmp.toa;

    %---------------------------Groupdelay--------------------------
    [~,tmp]=ziegelwanger2013(data,3,0);
    AGD=tmp.toa;

    %---------------------------Minimal-Phase-----------------------
    [~,tmp]=ziegelwanger2013(data,4,0);
    MCM=tmp.toa;
    clear tmp

    plot_ziegelwanger2013(data,MCM(:,1),4,[0 0 0]/255,0,1,1,'-',1);
    hold on
    plot_ziegelwanger2013(data,MAX(:,1),4,[0 0 80]/255,0,1,1,'--',1);
    plot_ziegelwanger2013(data,CTD(:,1),4,[50 220 50]/255,0,1,1,'-',1);
    plot_ziegelwanger2013(data,AGD(:,1),4,[250 80 80]/255,0,1,1,'--',1);
    xlim([-10 370])
    ylim([0.65 2.05])
    grid off
    legend(' MCM',' MAX',' CTD',' AGD','Location','NorthWest');
    legend boxoff
    xlabel('Azimuth (deg) ')
    ylabel('TOA (ms) ')
%     set(gca,'FontSize',ticklabelsize,'LineWidth',bw,'LineWidth',bw);%,'YTick',[-0.8 -0.4 0 0.4 0.8]
    title('');
    
    subplot(121)
    time=(0:data.API.N-1)/data.Data.SamplingRate*1000;
    MAX=round(MAX);
    CTD=round(CTD);
    AGD=round(AGD);
    MCM=round(MCM);
    idx=local_ARI_FindPosition(data,85,0);

    fprintf(['MAX: ITD is ' num2str(time(diff(MAX(idx,:),1,2))) ' ms\n'])
    fprintf(['CTD: ITD is ' num2str(time(diff(CTD(idx,:),1,2))) ' ms\n'])
    fprintf(['AGD: ITD is ' num2str(time(diff(AGD(idx,:),1,2))) ' ms\n'])
    fprintf(['MCM: ITD is ' num2str(time(diff(MCM(idx,:),1,2))) ' ms\n'])

    plot(time,squeeze(data.Data.IR(idx,1,:))/max(abs(data.Data.IR(idx,1,:))),'k-');
    hold on
    plot(time,squeeze(data.Data.IR(idx,2,:))/max(abs(data.Data.IR(idx,2,:))),'k--')
    h=stem([time(MCM(idx,1)) time(MCM(idx,1))],[-2 data.Data.IR(idx,1,MCM(idx,1))/max(abs(data.Data.IR(idx,1,:)))],'r-','BaseValue',-1);
    stem([time(CTD(idx,1)) time(CTD(idx,1))],[-2 data.Data.IR(idx,1,CTD(idx,1))/max(abs(data.Data.IR(idx,1,:)))],'g-','BaseValue',-1)
    stem([time(MAX(idx,1)) time(MAX(idx,1))],[-2 data.Data.IR(idx,1,MAX(idx,1))/max(abs(data.Data.IR(idx,1,:)))],'b-','BaseValue',-1)
    stem([time(AGD(idx,1)) time(AGD(idx,1))],[-2 data.Data.IR(idx,1,AGD(idx,1))/max(abs(data.Data.IR(idx,1,:)))],'m-','BaseValue',-1)
    stem([time(MCM(idx,2)) time(MCM(idx,2))],[-2 data.Data.IR(idx,2,MCM(idx,2))/max(abs(data.Data.IR(idx,2,:)))],'r-','BaseValue',-1)
    stem([time(CTD(idx,2)) time(CTD(idx,2))],[-2 data.Data.IR(idx,2,CTD(idx,2))/max(abs(data.Data.IR(idx,2,:)))],'g-','BaseValue',-1)
    stem([time(AGD(idx,2)) time(AGD(idx,2))],[-2 data.Data.IR(idx,2,AGD(idx,2))/max(abs(data.Data.IR(idx,2,:)))],'m-','BaseValue',-1)
    stem([time(MAX(idx,2)) time(MAX(idx,2))],[-2 data.Data.IR(idx,2,MAX(idx,2))/max(abs(data.Data.IR(idx,2,:)))],'b-','BaseValue',-1)
    plot(time(MAX(idx,1)),data.Data.IR(idx,1,MAX(idx,1))/max(abs(data.Data.IR(idx,1,:))),'b^','MarkerFaceColor','b')
    plot(time(MCM(idx,1)),data.Data.IR(idx,1,MCM(idx,1))/max(abs(data.Data.IR(idx,1,:))),'ro','MarkerFaceColor','r')
    plot(time(CTD(idx,1)),data.Data.IR(idx,1,CTD(idx,1))/max(abs(data.Data.IR(idx,1,:))),'gd','MarkerFaceColor','g')
    plot(time(AGD(idx,1)),data.Data.IR(idx,1,AGD(idx,1))/max(abs(data.Data.IR(idx,1,:))),'ms','MarkerFaceColor','m')
    plot(time(MAX(idx,2)),data.Data.IR(idx,2,MAX(idx,2))/max(abs(data.Data.IR(idx,2,:))),'b^','MarkerFaceColor','b')
    plot(time(MCM(idx,2)),data.Data.IR(idx,2,MCM(idx,2))/max(abs(data.Data.IR(idx,2,:))),'ro','MarkerFaceColor','r')
    plot(time(CTD(idx,2)),data.Data.IR(idx,2,CTD(idx,2))/max(abs(data.Data.IR(idx,2,:))),'gd','MarkerFaceColor','g')
    plot(time(AGD(idx,2)),data.Data.IR(idx,2,AGD(idx,2))/max(abs(data.Data.IR(idx,2,:))),'ms','MarkerFaceColor','m')
    xlim([0.4 2.1])
    ylim([-1.1 1.1])
    xlabel('Time (ms) ')
    ylabel('Amplitude ')
%     set(gca,'FontSize',ticklabelsize,'LineWidth',bw,'LineWidth',bw,'XTick',[2.5 3 3.5 4]);
    set(get(h,'Baseline'),'Visible','off')
end

%% Figure 2b
if flags.do_fig2b
        
    data=data_ziegelwanger2013('NH89',flags.cachemode);
    ch=1;

    p0_onaxis=[[0.0875; pi/2; 0; 0.0001] [0.0875; -pi/2; 0; 0.0001]];
    p0_onaxis=transpose(p0_onaxis);
    p_onaxis=zeros(size(p0_onaxis));
    p0_offaxis=zeros(2,7);
    p_offaxis=p0_offaxis;

    indicator=zeros(data.API.M,data.API.R);
    indicator_hor=indicator;
    indicator_sag=indicator;
    pos=zeros(data.API.M,8);
    pos(:,1:2)=data.SourcePosition(:,1:2);
    [pos(:,6),pos(:,7)]=sph2hor(data.SourcePosition(:,1),data.SourcePosition(:,2));
    pos(:,8)=cumsum(ones(data.API.M,1));
    [~,tmp]=ziegelwanger2013(data,4,0);
    toaEst=tmp.toa;

    % Outlier detection: smooth TOA in horizontal planes
    epsilon=5;
    slope=zeros(data.API.M,1);
    for ele=min(pos(:,2)):epsilon:max(pos(:,2)) %calculate slope for each elevation along azimuth
        idx=find(pos(:,2)>ele-epsilon/2 & pos(:,2)<=ele+epsilon/2);
        if numel(idx)>1
            idx(length(idx)+1)=idx(1);
            slope(idx(1:end-1),1)=diff(toaEst(idx,ch))./abs(diff(pos(idx,1)));
        end
    end
    sloperms=sqrt(sum(slope.^2)/length(slope));
    if sloperms<30/(length(find(pos(:,2)==0))/2)
        sloperms=30/(length(find(pos(:,2)==0))/2);
    end
    for ele=min(pos(:,2)):epsilon:max(pos(:,2))
        idx=find(pos(:,2)>ele-epsilon/2 & pos(:,2)<=ele+epsilon/2);
        for ii=1:length(idx)-1
            if abs(slope(idx(ii)))>sloperms
                for jj=0:1
                    if ii+jj==0 || ii+jj==length(idx)
                        indicator_hor(idx(end),ch)=1;
                    else
                        indicator_hor(idx(mod(ii+jj,length(idx))),ch)=1;
                    end
                end
            end
        end
        clear idx
    end

    % Outlier detection: constant TOA in sagittal planes
    epsilon=2;
    for ii=1:20
        sag_dev=zeros(data.API.M,1);
        for lat=-90:epsilon:90
            idx=find(pos(:,6)>lat-epsilon/2 & pos(:,6)<=lat+epsilon/2); 
            idx2=find(pos(:,6)>lat-epsilon/2 & pos(:,6)<=lat+epsilon/2 & indicator_hor(:,ch)==0 & indicator(:,ch)==0);
            if length(idx2)>2
                sag_dev(idx,1)=toaEst(idx,ch)-mean(toaEst(idx2,ch));
            end
        end
        sag_var=sqrt(sum(sag_dev.^2)/length(sag_dev));
        if sag_var<2
            sag_var=2;
        end
        indicator(:,ch)=zeros(size(indicator,1),1);
        indicator_sag(:,ch)=zeros(size(indicator_sag,1),1);
        indicator_sag(abs(sag_dev)>sag_var,ch)=ones(length(find(abs(sag_dev)>sag_var)),1);
        indicator(abs(sag_dev)>sag_var | indicator_hor(:,ch)==1,ch)=ones(length(find(abs(sag_dev)>sag_var | indicator_hor(:,ch)==1)),1);
    end

    clear sag_dev; clear sag_var;

    plot_ziegelwanger2013(data,toaEst,4,'k',0,1,1,{'-'},1);
    hold on
    h=plot_ziegelwanger2013(data,indicator_sag.*toaEst,3,'w',0,1,1,{'^'},4);
    set(h,'MarkerFaceColor','w','MarkerEdgeColor','w');
    h=plot_ziegelwanger2013(data,indicator_hor.*toaEst,3,'w',0,1,1,{'v'},4);
    set(h,'MarkerFaceColor','w','MarkerEdgeColor','w');
    h=plot_ziegelwanger2013(data,indicator_sag.*toaEst,3,'b',0,1,1,{'^'},4);
    set(h,'LineWidth',2);
    h=plot_ziegelwanger2013(data,indicator_hor.*toaEst,3,'b',0,1,1,{'v'},4);
    set(h,'LineWidth',2);
    h=plot_ziegelwanger2013(data,(-indicator+1).*toaEst,3,'r',0,1,1,{'o'},4);
    set(h,'MarkerFaceColor','r','MarkerEdgeColor','r');
    ylabel('TOA (ms)')
    xlabel('Azimuth (deg)')
    grid off
    xlim([-10 370])
    ylim([0.65 1.65])
    title('')
    set(gca,'YTick',[0.7 1 1.3 1.6])
end

%% Figure 3b
if flags.do_fig3b

    figure
    sym='sdo'; %plot symbols
    clr=[0,0,255; 255,0,0; 255,255,67]/255; %plot colors
    meclr=[0,0,255; 255,0,0; 255,255,67]/255; %marker edge colors
    
    data=data_ziegelwanger2013('SPHERE_ROT',flags.cachemode);
    
    % radii
    subplot(311)
    var=[squeeze(data.results.p_onaxis(1,1,:))*100 squeeze(data.results.p_onaxis(1,2,:))*100 data.radius(:)/10];
    for ch=1:size(data.results.p_onaxis,2)
        plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
        hold on
    end
    plot(var(:,3),'k--')
    clear var;
    ylabel('r in cm')

    %phi
    subplot(312)
    var=[squeeze(data.results.p_onaxis(2,1,:))/pi*180 squeeze(data.results.p_onaxis(2,2,:))/pi*180 -data.phi+ones(length(data.phi),1)*90 -data.phi-ones(length(data.phi),1)*90];
    for ch=1:size(data.results.p_onaxis,2)
        plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
        hold on
    end
    for ch=1:size(data.results.p_onaxis,2)
        plot(1:9,var(1:9,2+ch),'k--')
        plot(10:14,var(10:14,2+ch),'k--')
        plot(15:23,var(15:23,2+ch),'k--')
        plot(24:28,var(24:28,2+ch),'k--')
        plot(29:37,var(29:37,2+ch),'k--')
        plot(38:42,var(38:42,2+ch),'k--')
    end
    clear var;
    set(gca,'YTick',[-90 90])
    ylabel('\phi_e in deg')

    %theta
    subplot(313)
    var=[squeeze(data.results.p_onaxis(3,1,:))/pi*180 squeeze(data.results.p_onaxis(3,2,:))/pi*180 data.theta -data.theta];
    for ch=1:size(data.results.p_onaxis,2)
        plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
        hold on
    end
    for ch=1:size(data.results.p_onaxis,2)
        plot(1:9,var(1:9,2+ch),'k--')
        plot(10:14,var(10:14,2+ch),'k--')
        plot(15:23,var(15:23,2+ch),'k--')
        plot(24:28,var(24:28,2+ch),'k--')
        plot(29:37,var(29:37,2+ch),'k--')
        plot(38:42,var(38:42,2+ch),'k--')
    end
    clear var;
    ylabel('\theta_e in deg')
    xlabel('Condition')
    
    figure
    sym='sdo'; %plot symbols
    clr=[0,0,255; 255,0,0; 255,255,67]/255; %plot colors
    meclr=[0,0,255; 255,0,0; 255,255,67]/255; %marker edge colors
    
    data=data_ziegelwanger2013('SPHERE_DIS',flags.cachemode);
    
    p1=data.results.p_onaxis(:,:,[1:3 length(data.xM)/3+1:length(data.xM)/3+3 length(data.xM)/3*2+1:length(data.xM)/3*2+3]);
    r1=data.radius([1:3 length(data.xM)/3+1:length(data.xM)/3+3 length(data.xM)/3*2+1:length(data.xM)/3*2+3]);
    yM1=data.yM([1:3 length(data.xM)/3+1:length(data.xM)/3+3 length(data.xM)/3*2+1:length(data.xM)/3*2+3]);

    %radii
    subplot(211)
%         ymax2=round(max(max(squeeze(p1(1,:,:)*100))))+1;
    var=[squeeze(p1(1,1,:))*100 squeeze(p1(1,2,:))*100 r1/10];
    for ch=1:size(p1,2)
        plot(var(:,ch),sym(ch),'MarkerEdgeColor',meclr(ch,:),'MarkerFaceColor',clr(ch,:));
        hold on
    end
    plot(r1/10,'k--')
    clear var;
    ylabel('r in cm')

    %yM
    subplot(212)
    plot(-yM1*100,'k--')
    clear var;
    xlabel('Condition')
    ylabel('y_M in cm ')
end

%% Figure 1a
if flags.do_fig1a
    
    hrtf={'ARI','CIPIC','LISTEN'};
    sym='ods'; %plot symbols

    %-------------------------------Load Data----------------------------------
    for kk=1:length(hrtf)
        data=data_ziegelwanger2013(hrtf{kk},flags.cachemode);
        if kk==3
            data.results=data.results([1:27 29:end]);
        end
        temp1=zeros(size(data.results(1).meta.p_onaxis,1),size(data.results(1).meta.p_onaxis,2),length(data.results));
        temp3=zeros(size(data.results(1).meta.p_offaxis,1),size(data.results(1).meta.p_offaxis,2),length(data.results));
        for ii=1:length(data.results)
            temp1(:,:,ii)=data.results(ii).meta.p_onaxis;
            temp3(:,1:size(data.results(ii).meta.p_offaxis,2),ii)=data.results(ii).meta.p_offaxis;
        end
        p_onaxis{kk}=temp1;
        p_offaxis{kk}=temp3;
    end
    
    %radii
    temp=1;
    for kk=1:length(hrtf)
        var(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            [squeeze(mean(p_offaxis{kk}(1,:,:)*100,2)) kk* ...
             ones(size(p_onaxis{kk},3),1) transpose(1: ...
                                                    size(p_onaxis{kk},3))];
        varl(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            [squeeze(p_offaxis{kk}(1,1,:)*100) kk* ...
             ones(size(p_onaxis{kk},3),1) transpose(1: ...
                                                    size(p_onaxis{kk},3))];
        varr(temp:temp+size(p_onaxis{kk},3)-1,:)= ...
            [squeeze(p_offaxis{kk}(1,2,:)*100) kk* ...
             ones(size(p_onaxis{kk},3),1) transpose(1: ...
                                                    size(p_onaxis{kk},3))];
        temp=size(var,1)+1;
    end
    [~,idx]=sort(var(:,1));
    var=var(idx,:);
    varl=varl(idx,:);
    varr=varr(idx,:);
    var(:,3)=transpose(1:size(var,1));
    varl(:,3)=transpose(1:size(var,1));
    varr(:,3)=transpose(1:size(var,1));
    subplot(411)
    for kk=1:length(hrtf)
        h{kk}=plot(var(var(:,2)==kk,3),var(var(:,2)==kk,1),sym(kk),'MarkerEdgeColor','b','MarkerFaceColor','b');
        hold on
    end
    ylabel('Average r (cm)')
    xlim([-1.5 temp+1.5])
    ylim([7.5 10.57])
%     legend([h{1},h{2},h{3}],'ARI','CIPIC','LISTEN','Location','NorthWest');
    legend boxoff

    subplot(412)
    for kk=1:length(hrtf)
        plot(var(var(:,2)==kk,3),abs(varl(varl(:,2)==kk,1)-varr(varr(:,2)==kk,1)),sym(kk),'MarkerEdgeColor','b','MarkerFaceColor','b');
        hold on
    end
    ylabel('Binaural difference of r (cm)')
    xlim([-1.5 temp+1.5])
    ylim([-0.5 3.5])
    clear var varl varr;
    
    % yM
    temp=1;
    for kk=1:length(hrtf)
        var(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(mean(p_offaxis{kk}(3,:,:)*100,2)) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
        varl(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_offaxis{kk}(3,1,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
        varr(temp:temp+size(p_onaxis{kk},3)-1,:)=[squeeze(p_offaxis{kk}(3,2,:)*100) kk*ones(size(p_onaxis{kk},3),1) transpose(1:size(p_onaxis{kk},3))];
        temp=size(var,1)+1;
    end
    var=var(idx,:);
    varl=varl(idx,:);
    varr=varr(idx,:);
    var(:,3)=transpose(1:size(var,1));
    varl(:,3)=transpose(1:size(var,1));
    varr(:,3)=transpose(1:size(var,1));
    subplot(413)
    for kk=1:length(hrtf)
        plot(var(var(:,2)==kk,3),var(var(:,2)==kk,1),sym(kk),'MarkerEdgeColor','b','MarkerFaceColor','b');
        hold on
    end
    ylabel('y_M in cm')
    xlim([-1.5 temp+1.5])
    ylim([-3 4])

    subplot(414)
    for kk=1:length(hrtf)
        plot(var(var(:,2)==kk,3),abs(varl(varl(:,2)==kk,1)-varr(varr(:,2)==kk,1)),sym(kk),'MarkerEdgeColor','b','MarkerFaceColor','b');
        hold on
    end
    ylabel('Binaural difference of r (cm)')
    xlim([-1.5 temp+1.5])
    ylim([-0.5 3.5])
    clear var varl varr;

    xlabel('Subjects')

end



end

function idx=local_ARI_FindPosition(data,azimuth,elevation)
    psi=sin(elevation/180*pi).*sin(data.SourcePosition(:,2)/180*pi) + ...
        cos(elevation/180*pi).*cos(data.SourcePosition(:,2)/180*pi).*...
        cos(azimuth/180*pi-data.SourcePosition(:,1)/180*pi);
    [~,idx]=min(acos(psi));
end


