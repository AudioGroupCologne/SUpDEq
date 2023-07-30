function [Obj,results]=ziegelwanger2013(Obj,method,model,p0_onaxis)
%ZIEGELWANGER2013 Direction-continuous model of time-of-arrival (TOA) in HRTFs (simple)
%   Usage: [Obj,results]=ziegelwanger2013(Obj,method,model,p0_onaxis) 
%
%   Input parameters:
%     Obj                 : SOFA object
%
%     method              : (Optional) Select one of the estimation methods
%
%     model               : (Optional) Correct estimated toa, using geometrical
%                           TOA-Model. If model=0 use TOA estimated, 
%                           if model=1 use TOA modeled (default).
%
%     p0_onaxis           : (optional) Starting values for lsqcurvefit
% 
%   Output parameters:
%     Obj                 : SOFA Object 
%
%                           - results.toa   data matrix with time of arrival (TOA) for each impulse response (IR):
%                           - results.p_onaxis   estimated on-axis model-parameters
%                           - results.p_offaxis   estimated off-axis model-parameters
%
%   Estimates the Time-of-Arrival for each measurement in Obj (SOFA) and
%   corrects the results with a geometrical model of the head.
%
%   The value of method is an integer choosing one of the following
%   methods. XXX Explain for each method a little about how they work:
%
%   1) Threshold-Detection
%
%   2) Centroid of squared IR
%
%   3) Mean Groupdelay
%
%   4) Minimal-Phase Cross-Correlation (Max) (default)
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Optimization Toolbox for Matlab
%
%
%   Examples:
%   ---------
% 
%   To calculate the model parameters for the on-axis time-of-arrival model
%   (p_onaxis) and for the off-axis time-of-arrival model (p_offaxis) for a
%   given HRTF set (SOFA object, 'Obj') with the minimum-phase
%   cross-correlation method, use:
%
%       [Obj,results]=ziegelwanger2013(Obj,4,1);
%
%   See also: ziegelwanger2013_onaxis, ziegelwanger2013_offaxis,
%   data_ziegelwanger2013, exp_ziegelwanger2013 plot_ziegelwanger2013
%   ziegelwanger2014
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/ziegelwanger2013.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA M-Signal M-Optimization
%   #Author: Harald Ziegelwanger (2013), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% ----------------------------convert to SOFA-----------------------------
if ~isfield(Obj,'GLOBAL_Version')
    Obj=SOFAconvertARI2SOFA(Obj.hM,Obj.meta,Obj.stimPar);
end

%% ----------------------------check variables-----------------------------

if ~exist('method','var')
    method=4;
else if isempty(method)
        method=4;
    end
end

if ~exist('model','var')
    model=1;
else if isempty(model)
        model=1;
    end
end

if ~exist('p0_onaxis','var')
    p0_onaxis=[[0.0875; pi/2; 0; 0] [0.0875; -pi/2; 0; 0]];
end

%% -------------------------initialize variables---------------------------
p0_onaxis=transpose(p0_onaxis);
p_onaxis=zeros(size(p0_onaxis));
p0_offaxis=zeros(2,7);
p_offaxis=p0_offaxis;

toa=zeros(Obj.API.M,Obj.API.R);
toaEst=zeros(Obj.API.M,Obj.API.R);
indicator=zeros(Obj.API.M,Obj.API.R);
indicator_hor=indicator;
indicator_sag=indicator;
pos=zeros(Obj.API.M,8);
pos(:,1:2)=Obj.SourcePosition(:,1:2);
[pos(:,6),pos(:,7)]=sph2hor(Obj.SourcePosition(:,1),Obj.SourcePosition(:,2));
pos(:,8)=cumsum(ones(Obj.API.M,1));
% for ii=1:Obj.API.M
%     for jj=1:Obj.API.R
%         if isnan(Obj.Data.IR_min(1,ii,jj))
%             hM_min(:,ii,jj)=ARI_MinimalPhase(Obj.hM(:,ii,jj));
%         end
%     end
% end

%% -----------------------estimate time-of-arrival-------------------------
switch method
    case 1 %---------------------------Threshold---------------------------
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                toaEst(ii,jj)=find(abs(Obj.Data.IR(ii,jj,:))==max(abs(Obj.Data.IR(ii,jj,:))),1);
            end
        end
    case 2 %---------------------------Centroid----------------------------
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                toaEst(ii,jj)=find(cumsum(Obj.Data.IR(ii,jj,:).^2)>(sum(Obj.Data.IR(ii,jj,:).^2)/2),1);
            end
        end
    case 3 %---------------------------Groupdelay--------------------------
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                [Gd,F]=grpdelay(transpose(double(squeeze(Obj.Data.IR(ii,jj,:)))),1,Obj.API.N,Obj.Data.SamplingRate);
                toaEst(ii,jj)=median(Gd(find(F>500):find(F>2000)));
            end
        end
    case 4 %---------------------------Minimal-Phase-----------------------
        Obj=ARI_MinimalPhase(Obj);
        corrcoeff=zeros(Obj.API.M,Obj.API.R);
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                [c,lag]=xcorr(squeeze(Obj.Data.IR(ii,jj,:)),squeeze(Obj.Data.IR_min(ii,jj,:)),Obj.API.N-1,'none');
                [corrcoeff(ii,jj),idx]=max(abs(c));
                corrcoeff(ii,jj)=corrcoeff(ii,jj)/sum(Obj.Data.IR(ii,jj,:).^2);
                toaEst(ii,jj)=lag(idx);
            end
        end
end

%% ----------------------Fit-Models-to-estimated-TOA-----------------------
for ch=1:Obj.API.R

    % Outlier detection: smooth TOA in horizontal planes
    epsilon=5;
    slope=zeros(Obj.API.M,1);
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
        sag_dev=zeros(Obj.API.M,1);
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
        indicator_sag(:,ch)=abs(sag_dev)>sag_var;
        indicator(:,ch)=(abs(sag_dev)>sag_var | indicator_hor(:,ch));
    end
    clear sag_dev; clear sag_var;
end

performance.indicator=indicator;
performance.outliers=sum(sum(indicator))/Obj.API.M/2*100;
performance.outliersl=sum(indicator(:,1))/Obj.API.M*100;
performance.outliersr=sum(indicator(:,2))/Obj.API.M*100;

if model
    % Fit on-axis model to outlier adjusted set of estimated TOAs
    for ch=1:Obj.API.R
        p0_onaxis(ch,4)=min(toaEst(indicator(:,ch)==0,ch))/Obj.Data.SamplingRate;
        p0offset_onaxis=[0.06 pi/4 pi/4 0.001];

        idx=find(indicator(:,ch)==0);
        x=pos(idx,1:2)*pi/180;
        y=toaEst(idx,ch)/Obj.Data.SamplingRate;
        if isoctave
            [~,p_onaxis(ch,:)]=leasqr(x,y,p0_onaxis(ch,:),@ziegelwanger2013_onaxis);
        else
            p_onaxis(ch,:)=lsqcurvefit(@ziegelwanger2013_onaxis,p0_onaxis(ch,:),x,y,p0_onaxis(ch,:)-p0offset_onaxis,p0_onaxis(ch,:)+p0offset_onaxis,optimset('Display','off','TolFun',1e-6));
        end
        toa(:,ch)=ziegelwanger2013_onaxis(p_onaxis(ch,:),pos(:,1:2)*pi/180)*Obj.Data.SamplingRate;
    end

    % Fit off-axis model to outlier adjusted set of estimated TOAs
    TolFun=[1e-5; 1e-6];
    for ii=1:size(TolFun,1)
        for ch=1:Obj.API.R
            idx=find(indicator(:,ch)==0);
            x=pos(idx,1:2)*pi/180;
            y=toaEst(idx,ch)/Obj.Data.SamplingRate;
            p0_offaxis(ch,:)=[p0_onaxis(ch,1) 0 0 0 p0_onaxis(ch,4) p0_onaxis(ch,2) p0_onaxis(ch,3)];
            p0offset_offaxis=[0.05 0.05 0.05 0.05 0.001 pi pi];
            if isoctave
                [~,p_offaxis(ch,:)]=leasqr(x,y,p0_offaxis(ch,:),@ziegelwanger2013_offaxis);
            else
                p_offaxis(ch,:)=lsqcurvefit(@ziegelwanger2013_offaxis,p0_offaxis(ch,:),x,y,p0_offaxis(ch,:)-p0offset_offaxis,p0_offaxis(ch,:)+p0offset_offaxis,optimset('Display','off','TolFun',TolFun(ii,1)));
            end
            toa(:,ch)=ziegelwanger2013_offaxis(p_offaxis(ch,:),pos(:,1:2)*pi/180)*Obj.Data.SamplingRate;
        end
        if abs(diff(p_offaxis(:,1)))>0.003 || abs(diff(p_offaxis(:,3)))>0.003
            p_offaxis(:,[1 3])=p_offaxis([2 1],[1 3]);
            for ch=1:Obj.API.R
                idx=find(indicator(:,ch)==0);
                x=pos(idx,1:2)*pi/180;
                y=toaEst(idx,ch)/Obj.Data.SamplingRate;
                p0_offaxis(ch,:)=[p_offaxis(ch,1) mean(p_offaxis(:,2)) p_offaxis(ch,3) mean(p_offaxis(:,4)) mean(p_offaxis(:,5)) p_offaxis(ch,6) p_offaxis(ch,7)];
                p0offset_offaxis=[0.05 0.05 0.05 0.05 0.001 pi/2 pi/2];
                if isoctave
                    [~,p_offaxis(ch,:)]=leasqr(x,y,p0_offaxis(ch,:),@ziegelwanger2013_offaxis);
                else
                    p_offaxis(ch,:)=lsqcurvefit(@ziegelwanger2013_offaxis,p0_offaxis(ch,:),x,y,p0_offaxis(ch,:)-p0offset_offaxis,p0_offaxis(ch,:)+p0offset_offaxis,optimset('Display','off','TolFun',TolFun(ii,1)));
                end
                toa(:,ch)=ziegelwanger2013_offaxis(p_offaxis(ch,:),pos(:,1:2)*pi/180)*Obj.Data.SamplingRate;
            end
        end
        if abs(diff(p_offaxis(:,1)))<0.003 && abs(diff(p_offaxis(:,2)))<0.003 && abs(diff(p_offaxis(:,3)))<0.003 && abs(diff(p_offaxis(:,4)))<0.003
            break
        end
    end
else
    toa=toaEst;
    p_offaxis=p0_offaxis;
end

Obj.Data.Delay=toa;
Obj.Data.p_onaxis=transpose(p_onaxis);
Obj.Data.p_offaxis=transpose(p_offaxis);
Obj.Data.performance=performance;

results.toa=toa;
results.p_onaxis=transpose(p_onaxis);
results.p_offaxis=transpose(p_offaxis);
results.performance=performance;

end %of function

function Obj=ARI_MinimalPhase(Obj)
    Obj.Data.IR_min=zeros(size(Obj.Data.IR));

    for jj=1:Obj.API.R
        for ii=1:Obj.API.M
%             h=squeeze(Obj.Data.IR(ii,jj,:));
            h=[squeeze(Obj.Data.IR(ii,jj,:)); zeros(4096-Obj.API.N,1)];
            % decompose signal
            amp1=abs(fft(h));

            % transform
            amp2=amp1;
            an2u=-imag(hilbert(log(amp1))); % minimal phase

            % reconstruct signal from amp2 and an2u
            % build a symmetrical phase 
            an2u=an2u(1:floor(length(h)/2)+1);
            an2u=[an2u; -flipud(an2u(2:end+mod(length(h),2)-1))];
            an2=an2u-round(an2u/2/pi)*2*pi;  % wrap around +/-pi: wrap(x)=x-round(x/2/pi)*2*pi
            % amplitude
            amp2=amp2(1:floor(length(h)/2)+1);
            amp2=[amp2; flipud(amp2(2:end+mod(length(h),2)-1))];
            % back to time domain
            h2=real(ifft(amp2.*exp(1i*an2)));
            Obj.Data.IR_min(ii,jj,:)=h2(1:Obj.API.N);
        end
    end
end



