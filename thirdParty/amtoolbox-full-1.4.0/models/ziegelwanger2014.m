function [Obj,results]=ziegelwanger2014(Obj,varargin)
%ZIEGELWANGER2014 Time of arrival estimates
%   Usage: [Obj,results]=ziegelwanger2014(data,estimation,outlierDetection,model,p0_onaxis) 
%
%
%   Input parameters:
%
%      Obj        : SOFA object
%
%      estimation : (optional) TOA estimation method:
%
%                   [TOAest]	pre-estimated TOAs
%
%                   1 		Built-in maximum detection
%
%                   2 	        Built-in centroid of squared IR,
%
%                   3 	        Built-in mean group delay,
%
%                   4 	        Built-in maximum of the minimum-phase cross-correlation (default)
%
%                   5          Built-in maxium detection after 3k-low-pass filtering (similar to Andreopoulou and Katz, 2017)
%
%                   String	(Experimental) Estimator from the 'itdestimator' function, requires further parameters
%
%      outlierDetection : (optional) detect outliers in estimated TOAs
%
%                         0           off
%
%                         1           on (default values: [0.05;0.01])
%
%                         [alpha r]   reject outliers using the extreme Studentized deviance test with the significance level of ALPHA and upper bound of outlier rate R
%
%      model: (optional) correct estimated toa, using geometrical TOA-Model
%
%             0  no model, just estimate the TOAs
%
%             1  (default) off-axis TOA modeled
%
%             2  on-axis TOA modeled
%
%      p0_onaxis: (optional) startvalues for lsqcurvefit
%
%      lowpass:   (optional) bandwidth setting when used with the estimator from 'itdestimator' ('lp' for lowpass, 'bb' for broadband)
%
%      upper_cutfreq: (optional) lowpass cutoff (Hz) when used with the estimator from 'itdestimator'
%
%      threshlvl: (optional) threshold level (dB) when used with the 'Threshold' estimator from 'itdestimator'
%   
%   Output parameters:
%
%      Obj: SOFA Object
%
%      results: struct with fields
%
%               .toa: data matrix with time of arrival (TOA) for each impulse response (IR)
%
%               .p_onaxis: estimated on-axis model-parameters
%
%               .p_offaxis: estimated off-axis model-parameters
%
%   Estimates the Time-of-Arrival for each measurement in Obj (SOFA) and
%   corrects the results with a geometrical model of the head.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Optimization Toolbox for Matlab
%
%   Examples:
%   ---------
% 
%   To calculate the model parameters for the on-axis time-of-arrival model
%   (p_onaxis) and for the off-axis time-of-arrival model (p_offaxis) for a
%   given HRTF set (SOFA object, 'Obj') with the minimum-phase
%   cross-correlation estimation, use:
%
%       [Obj,results]=ziegelwanger2014(Obj,4,1);
%
%   See also: ziegelwanger2014_onaxis, ziegelwanger2014_offaxis,
%             data_ziegelwanger2014, exp_ziegelwanger2014, itdestimator
%
%   References:
%     H. Ziegelwanger and P. Majdak. Modeling the direction-continuous
%     time-of-arrival in head-related transfer functions. J. Acoust. Soc.
%     Am., 135:1278--1293, 2014.
%     
%
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/ziegelwanger2014.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA M-Signal M-Optimization
%   #Author: Harald Ziegelwanger (2014), Acoustics Research Institute, Vienna, Austria
%   #Author: Robert Baumgartner (2018)
%   #Author: Laurin Steidle (2018)

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

definput.flags.lowpass = {'lp','bb'};
definput.flags.peak = {'hp','fp'};
definput.flags.disp = {'no_debug','debug'};
    
definput.keyvals.estimation = 4;
definput.keyvals.outlierDetection = [0.05;0.01];
definput.keyvals.model=1e-6;
definput.keyvals.p0_onaxis=[[0.0875; pi/2; 0; 0] [0.0875; -pi/2; 0; 0]];

definput.keyvals.threshlvl = -10;
definput.keyvals.upper_cutfreq = 3000;
definput.keyvals.butterpoly = 10;

[flags,kv]=ltfatarghelper({'estimation','outlierDetection','model','p0_onaxis'},definput,varargin);



estimation = kv.estimation;
outlierDetection = prod(kv.outlierDetection);

%% -------------------------initialize variables---------------------------
p0_onaxis=transpose(kv.p0_onaxis);
p_onaxis=zeros(size(p0_onaxis));
p0_offaxis=zeros(2,7);
p_offaxis=p0_offaxis;

toa=zeros(Obj.API.M,Obj.API.R);
toa_onaxis=toa;
toa_offaxis=toa;
indicator=zeros(Obj.API.M,Obj.API.R);

pos(:,1:2)=Obj.SourcePosition(:,1:2);

%% -----------------------estimate time-of-arrival-------------------------
if isscalar(estimation)
    hM=Obj.Data.IR;
    
    toaEst=zeros(Obj.API.M,Obj.API.R);
    switch estimation
        case 1 %---------------------------Threshold---------------------------
            for ii=1:Obj.API.M
                for jj=1:Obj.API.R
                    toaEst(ii,jj)=find(abs(hM(ii,jj,:))==max(abs(hM(ii,jj,:))),1);
                end
            end
        case 2 %---------------------------Centroid----------------------------
            for ii=1:Obj.API.M
                for jj=1:Obj.API.R
                    toaEst(ii,jj)=find(cumsum(hM(ii,jj,:).^2)>(sum(hM(ii,jj,:).^2)/2),1);
                end
            end
        case 3 %---------------------------Groupdelay--------------------------
            for ii=1:Obj.API.M
                for jj=1:Obj.API.R
                    [Gd,F]=grpdelay(transpose(double(squeeze(hM(ii,jj,:)))),...
                        1,Obj.API.N*4,Obj.Data.SamplingRate*4);
                    toaEst(ii,jj)=mean(Gd(find(F>1000):find(F>5000)));
                end
            end
        case 4 %---------------------------Minimal-Phase-----------------------
            hMmin=ARI_MinimalPhase(Obj);
            corrcoeff=zeros(Obj.API.M,Obj.API.R);
            for ii=1:Obj.API.M
                for jj=1:Obj.API.R
                    [c,lag]=xcorr(squeeze(hM(ii,jj,:)),squeeze(hMmin(ii,jj,:)),...
                        Obj.API.N*4-1,'none');
                    [corrcoeff(ii,jj),idx]=max(abs(c));
                    corrcoeff(ii,jj)=corrcoeff(ii,jj)/sum(hM(ii,jj,:).^2);
                    toaEst(ii,jj)=lag(idx);
                end
            end
      case 5 % ------------ Threshold on low-passed IRs

        amt_disp('  Applying Butterworth low pass',flags.disp)
        amt_disp(strcat('  Polynomial order of Butterworth filter: ',num2str(kv.butterpoly)),flags.disp)
        amt_disp(strcat('  Cut-off frequency is: ',num2str(kv.upper_cutfreq),' Hz'),flags.disp)
        cut_off_freq_norm = kv.upper_cutfreq/(Obj.Data.SamplingRate/2);
        [lp_a,lp_b] = butter(kv.butterpoly,cut_off_freq_norm);
%         f_ir = zeros(pos,ear,Ns);
        for ii=1:Obj.API.M
            for jj=1:Obj.API.R
                sir = squeeze( hM(ii,jj,:) );
                f_sir = filter(lp_a,lp_b,sir);
%                 f_ir(ii,jj,:) = f_sir;
                toaEst(ii,jj)=(find(abs(f_sir)==max(abs(f_sir)),1))-round(0.4e-3*Obj.Data.SamplingRate); % compensate for the delay of the filter, estimated for default parameters
            end
        end
    end
elseif ischar(estimation)
    [~,toaEst] = itdestimator(Obj,estimation,flags.lowpass,flags.peak,'guesstoa');
    toaEst = toaEst*Obj.Data.SamplingRate;
    
else
    toaEst=estimation;
end

%% --------------------Detect-Outliers-in-estimated-TOA--------------------
if outlierDetection>0
    for ch=1:Obj.API.R
        p0_onaxis(ch,4) = min(toaEst(indicator(:,ch)==0,ch))/Obj.Data.SamplingRate;
        p0offset_onaxis = [0.06 pi pi/2 0.001];
        x               = pos(:,1:2)*pi/180;
        y               = toaEst(:,ch)/Obj.Data.SamplingRate;

            tmp=lsqcurvefit(@ziegelwanger2014_onaxis,p0_onaxis(ch,:),...
                x,y,p0_onaxis(ch,:)-p0offset_onaxis,...
                p0_onaxis(ch,:)+p0offset_onaxis,...
                optimset('Display','off','TolFun',1e-6));

        outliers = toaEst(:,ch)-ziegelwanger2014_onaxis(tmp,pos(:,1:2)*pi/180)*Obj.Data.SamplingRate;
        [~,idx]=deleteoutliers(outliers, outlierDetection*Obj.API.M);
        indicator(idx,ch)=ones(length(idx),1);
    end
end
    
%% ----------------------Fit-Models-to-estimated-TOA-----------------------
if kv.model>0
    % Fit on-axis model to outlier adjusted set of estimated TOAs
    for ch=1:Obj.API.R
        p0_onaxis(ch,4) = min(toaEst(indicator(:,ch)==0,ch))/Obj.Data.SamplingRate;
        p0offset_onaxis = [0.06 pi pi/2 0.001];
        idx             = find(indicator(:,ch)==0);
        x               = pos(idx,1:2)*pi/180;
        y               = toaEst(idx,ch)/Obj.Data.SamplingRate;

            [p_onaxis(ch,:),performance.on_axis{ch}.resnormS,...
             performance.on_axis{ch}.residualS,...
             performance.on_axis{ch}.exitflag,...
             performance.on_axis{ch}.output] =...
                lsqcurvefit(@ziegelwanger2014_onaxis,p0_onaxis(ch,:),x,y,...
                    p0_onaxis(ch,:)-p0offset_onaxis,...
                    p0_onaxis(ch,:)+p0offset_onaxis,...
                    optimset('Display','off','TolFun',1e-6));
            toa(:,ch)=ziegelwanger2014_onaxis(p_onaxis(ch,:),pos(:,1:2)*pi/180)*Obj.Data.SamplingRate;

        performance.on_axis{ch}.resnormS = ...
            sqrt(performance.on_axis{ch}.resnormS/(Obj.API.M-sum(indicator(:,ch))));
        performance.on_axis{ch}.resnormP = ...
            norm((toaEst(:,ch)-toa(:,ch))/Obj.Data.SamplingRate)/sqrt(Obj.API.M);
    end
    toa_onaxis=toa;

    % Fit off-axis model to outlier adjusted set of estimated TOAs
    if kv.model~=2
        for ch=1:Obj.API.R
            idx = find(indicator(:,ch)==0);
            x   = pos(idx,1:2)*pi/180;
            y   = toaEst(idx,ch)/Obj.Data.SamplingRate;
            p0_offaxis(ch,:) = [mean(p_onaxis(:,1)) 0.001 ...
                                -diff(p_onaxis(:,1))/2 0.001 ...
                                mean(p_onaxis(:,4)) ...
                                p_onaxis(ch,2) ...
                                p_onaxis(ch,3)];
            p0offset_offaxis = [abs(diff(p_onaxis(:,1))/4) 0.1 0.1 0.1 0.001 pi/4 pi/4];
 
                [p_offaxis(ch,:),performance.off_axis{ch}.resnormS,...
                 performance.off_axis{ch}.residualS,...
                 performance.off_axis{ch}.exitflag,...
                 performance.off_axis{ch}.output] = ...
                    lsqcurvefit(@ziegelwanger2014_offaxis,...
                        p0_offaxis(ch,:),x,y,...
                        p0_offaxis(ch,:)-p0offset_offaxis,...
                        p0_offaxis(ch,:)+p0offset_offaxis,...
                        optimset('Display','off','TolFun',kv.model));
                toa(:,ch)=ziegelwanger2014_offaxis(p_offaxis(ch,:),pos(:,1:2)*pi/180)*Obj.Data.SamplingRate;

            performance.off_axis{ch}.resnormS = ...
                sqrt(performance.off_axis{ch}.resnormS/(Obj.API.M-sum(indicator(:,ch))));
            performance.off_axis{ch}.resnormP = ...
                norm((toaEst(:,ch)-toa(:,ch))/Obj.Data.SamplingRate)/sqrt(Obj.API.M);
        end
        toa_offaxis=toa;
    end
else
    toa=toaEst;
    p_offaxis=p0_offaxis;
end

%Save to output variables
performance.outliers=indicator;
for ii=1:size(indicator,2)
  performance.outlierRate(ii)=sum(indicator(:,ii))/Obj.API.M*100;
end
results.toa=toa;
results.toaEst=toaEst;
results.toa_onaxis=toa_onaxis;
results.toa_offaxis=toa_offaxis;
results.p_onaxis=transpose(p_onaxis);
results.p_offaxis=transpose(p_offaxis);
results.performance=performance;
if exist('corrcoeff','var')
    results.performance.corrcoeff=corrcoeff;
end

end %of function

function hMmin=ARI_MinimalPhase(Obj)
    hM=Obj.Data.IR;
    hMmin=hM;

    for jj=1:Obj.API.R
        for ii=1:Obj.API.M
            h=[squeeze(hM(ii,jj,:)); zeros(4096*4-size(hM,3),1)];
            amp1=abs(fft(h));
            amp2=amp1;
            an2u=-imag(hilbert(log(amp1)));
            an2u=an2u(1:floor(length(h)/2)+1);
            an3u=[an2u; -flipud(an2u(2:end+mod(length(h),2)-1))];
            an3=an3u-round(an3u/2/pi)*2*pi;
            amp2=amp2(1:floor(length(h)/2)+1);
            amp3=[amp2; flipud(amp2(2:end+mod(length(h),2)-1))];
            h2=real(ifft(amp3.*exp(1i*an3)));
            hMmin(ii,jj,:)=h2(1:Obj.API.N);
        end
    end
end


