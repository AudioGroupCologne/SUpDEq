function [E] = li2020(target_direct,template_direct,target_brir,template_brir,varargin)
%LI2020 Sound externalization in reverberant spaces
%
%   Usage:    E = li2020(target_direct,template_direct,target_brir,template_brir)
% 
%   Input parameters:
%     target        : binaural impulse response or head-related impulse response
%                     of the target sound. Matrix dimensions: time x receiver.
%     template      : binaural impulse response or head-related impulse response
%                     of the template sound. Matrix dimensions: time x receiver.
%                     Note that the dimensions of target and template must be
%                     identical.
%     target_brir   : target binaural room impulse response
%     template_brir : template binaural room impulse response
% 
%   Output parameters: 
%     E         : predicted degree of externalization
%
%
%   Examples:
%
%   Anechoic signal (rendered with HRTFs or direct sound part of BRIRs):
%
%     E = Li2020(tar_HRIR,tem_HRIR,[],[],'stim',stimulus,'fs',44100,'fsstim',44100,'flow',200, 'fhigh',16000, 'space',1);  
%
%   Reverberant signal (rendered with BRIRs with unkonwn direct sound part):
%
%     E = Li2020([],[],tar_BRIR,tem_BRIR,'stim',stimulus,'fs',44100,'fsstim',44100,'flow',200, 'fhigh',16000, 'space',1);  
%
%   Reverberant signal (rendered with BRIRs with konwn direct part (HRTFs)):
%
%     E = Li2020(tar_BRIR_dir,tem_BRIR_dir,tar_BRIR,tem_BRIR,'stim',stimulus,'fs',44100,'fsstim',44100,'flow',200, 'fhigh',16000, 'space',1);  
%
%   The implementation is based on baumgartner2021.m
%   See also: exp_li2020 baumgartner2021 exp_baumgartner2021 data_li2020 sig_li2020
%
%
%   References:
%     S. Li, R. Baumgartner, and J. Peissig. Modeling perceived
%     externalization of a static, lateral sound image. Acta Acustica, 4(5),
%     2021.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/li2020.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: M-Signal M-Image
%   #Author: Song Li (2020), Institute of Communications Technology, Leibniz University of Hannover, Germany

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Check input
definput.import={'baumgartner2014','baumgartner2014_pmv2ppp','localizationerror','amt_cache'};
definput.keyvals.JND = 1;
definput.keyvals.cueWeights = [2.5,1.7,0.9,2.8,0.6]; % weightings for the acoustic cues:
%[b_ILD, b_SG, b_w, b_ILD_TSD, b_lamda] 
definput.keyvals.ILD_sd_ref = 2.1; % ILD_TSD_ref = ([2.2 2.0 2.2 2.0 2.2]);
definput.keyvals.ILD_sd_temp = 0.67;% ILD_TSD_temp = ([0.7254 0.6011 0.7415 0.6635 0.6432]); ILD TSDs in anechoic chamber

[flags,kv]=ltfatarghelper(...
  {'fs','stim','space','do','flow','fhigh','fsstim','bwcoef'},definput,varargin);

  
  if(isempty(target_direct)||isempty(template_direct))
  % direct sound part is truncated from BRIRs      
    win_ramp_down = hann(round(0.001*kv.fs),'symmetric'); 
    win_echo=[ones(round(0.0025*kv.fs),1);win_ramp_down(round(length(win_ramp_down)/2):end);zeros(length(target_brir)-round(0.0025*kv.fs)-round(length(win_ramp_down)/2-1)-2,1)];
       
    target = permute(target_brir.*win_echo,[1,3,2]);
    template = permute(template_brir.*win_echo,[1,3,2]); 

    
  else
        % direct sound part and BRIRs are given 
        if(~isempty(target_brir)||~isempty(template_brir))

            win_ramp_down = hann(round(0.01*kv.fs)-1,'symmetric'); 
            win_ramp_up = hann(round(0.005*kv.fs)-1,'symmetric'); 
            win_echo=[ones(round(0.0025*kv.fs),1); win_ramp_up(round(0.0025*kv.fs):end); zeros(round(0.01*kv.fs)-1-round(0.0025*kv.fs),1);...
            win_ramp_down(1:round(0.005*kv.fs)*1);ones(length(target_brir)-round(0.01*kv.fs)-round(0.005*kv.fs)-round(0.0025*kv.fs),1)];
            
            if (length(target_brir)>length(target_direct))
                target_direct=padarray(target_direct,length(target_brir)-length(target_direct),0,'post');
            else
                target_brir=padarray(target_brir,length(target_direct)-length(target_brir),0,'post');
            end
            % remove additional echo
            target = permute(target_direct.*win_echo,[1,3,2]);
            template = permute(template_direct.*win_echo,[1,3,2]);
        else
            % only HRTFs are given, no BRIRs
            target = permute(target_direct,[1,3,2]);
            template = permute(template_direct,[1,3,2]);
        end
  end
  


%% DTF filtering

dimtar = size(target); % for lconv dim check
Nang = size(template,2);

if not(isempty(kv.stim))
  target = lconv(target,kv.stim);
  template = lconv(template,kv.stim);  
  
  target(length(kv.stim)+1:end,:) =[];
  template(length(kv.stim)+1:end,:) =[]; 
end

if size(target,2) ~= dimtar(2)
 target = reshape(target,[size(target,1),dimtar(2:end)]);
 template = reshape(template,[size(template,1),dimtar(2:end)]); 
end
%% Spectral Analysis

[tar.mp,fc] = baumgartner2014_spectralanalysis(target(1:length(target),:,:),'argimport',flags,kv);
[tem.mp,fc] = baumgartner2014_spectralanalysis(template(1:length(template),:,:),'argimport',flags,kv);

 %
%% ILD cues
    tar.ild = -diff(tar.mp,1,3); % ILD = left - right
    tem.ild = -diff(tem.mp,1,3);

    % target-template comparison -> ILD deviation
    dILD = abs(tem.ild-repmat(tar.ild,1,Nang));
    dILD(dILD < kv.JND) = 0; % limit minimum ILD difference according to JND  
    % overall normalized ILD deviation
    ILD_deviation = mean(dILD./abs(tem.ild));

%% SG cue

    % Spectral gradient extraction
    nrep.tem = baumgartner2014_gradientextraction(tem.mp,fc,'both');
    nrep.tar = baumgartner2014_gradientextraction(tar.mp,fc,'both');


    dNrep = (repmat(nrep.tar,1,Nang)-nrep.tem);  
    dNrep(abs(dNrep) < kv.JND) = 0;
        
    SG_deviation = mean(abs(dNrep))./mean((abs(nrep.tem)));
  
     
     
%% ILD TSD cue
% echo suppression  
    
  if(~isempty(target_brir)||~isempty(template_brir))

    win_ramp_down = hann(round(0.01*kv.fs)-1,'symmetric'); 
    win_ramp_up = hann(round(0.005*kv.fs)-1,'symmetric'); 
    win_echo=[ones(round(0.0025*kv.fs),1); win_ramp_up(round(0.0025*kv.fs):end); zeros(round(0.01*kv.fs)-1-round(0.0025*kv.fs),1);...
             win_ramp_down(1:round(0.005*kv.fs)*1);ones(length(target_brir)-round(0.01*kv.fs)-round(0.005*kv.fs)-round(0.0025*kv.fs),1)];
    
    target = permute(target_brir.*win_echo,[1,3,2]);
    template = permute(template_brir.*win_echo,[1,3,2]);
     
  else
      
    win_ramp_down = hann(round(0.001*kv.fs)); 
    win_echo=[ones(round(0.0025*kv.fs),1);win_ramp_down(round(0.0005*kv.fs):end);zeros(length(target_direct)-round(0.0025*kv.fs)-round(0.0005*kv.fs)-1,1)];

    target = permute(target_direct.*win_echo,[1,3,2]);
    template = permute(template_direct.*win_echo,[1,3,2]); 
    
  end
  
    dimtar = size(target); % for lconv dim check

    if not(isempty(kv.stim))
      target = lconv(target,kv.stim);
      template = lconv(template,kv.stim);  

      target(length(kv.stim)+1:end,:) =[];
      template(length(kv.stim)+1:end,:) =[];

    end

    % check that lconv preserved matrix dimensions
    if size(target,2) ~= dimtar(2)
     target = reshape(target,[size(target,1),dimtar(2:end)]);
     template = reshape(template,[size(template,1),dimtar(2:end)]); 
    end

    buf_frame = buffer(squeeze(target(:,1,1))', 0.02*kv.fs, 0.01*kv.fs, 'nodelay');
    dim_len = size(buf_frame);
    frameLength = dim_len(1);
    Nframes  = dim_len(2);

    if Nframes == 0
      Nframes = 1;
      frameLength = size(target,1);
      target = cat(1,target,zeros(frameLength-size(target,1),size(target,2),size(target,3)));
    end
    % window overlap  
      win_overlap = hann(round(0.02*kv.fs));
      win_overlap = permute(win_overlap,[1,3,2]);
      Nang = size(template,2);
    for iframe = 1:Nframes-1
            idt = (1:frameLength) + (iframe-1)*frameLength/2;

        %% Spectral Analysis
        tem =[];
        tar =[];
        [tar.mp,fc] = baumgartner2014_spectralanalysis(target(idt,:,:).*win_overlap,'argimport',flags,kv);
        [tem.mp,fc] = baumgartner2014_spectralanalysis(template(idt,:,:).*win_overlap,'argimport',flags,kv);

        %% ILDs
          tar.ild = diff(tar.mp,1,3);
          tem.ild = diff(tem.mp,1,3);

        %% target-template comparison -> ILD deviation
         ILD_temp(:,iframe) = tem.ild;
         ILD_tar(:,iframe) = repmat(tar.ild,1,Nang);

    end
    
    range=size(ILD_temp); 
    ILD_low_num=1;  ILD_high_num=range(2); 
    ILD_temp_std=(std(ILD_temp(:,ILD_low_num:ILD_high_num)'))';
    
    ILD_tar_std=(std(ILD_tar(:,ILD_low_num:ILD_high_num)'))';
    
    if (~isempty(target_brir)||~isempty(template_brir)) % reverberant condition
        norm_ILD_std_mean = abs(mean(ILD_temp_std - ILD_tar_std)) / (mean (ILD_temp_std));
        ILD_sd_temp_mean  = mean (ILD_temp_std);
        
    else
        ILD_TSD_ref = ([2.2 2.0 2.2 2.0 2.2]); % ILD TSD reference
        ILD_TSD_remaind = 0.0777; % Mal-adapted ILD TSD for anechoic conditions, which is used to compensate for the offset of the max. E-rating                            
        norm_ILD_std_mean = abs(mean(ILD_temp_std - ILD_tar_std + ones(34,1)*ILD_TSD_remaind)) / (mean (ILD_TSD_ref));                        
        ILD_sd_temp_mean =mean(ILD_temp_std + ones(34,1)*ILD_TSD_remaind); 
    end

 
%% Calculation of externalization ratings
    b_ILD = kv.cueWeights(1);
    b_SG = kv.cueWeights(2);
    w = kv.cueWeights(3);
    b_ILDTSD = kv.cueWeights(4);
    b_lamda = kv.cueWeights(5);

    % 
    ILD_sd_reference =kv.ILD_sd_ref; 
    if (~isempty(target_brir)||~isempty(template_brir)) % reverberant condition
        ILD_sd_template =kv.ILD_sd_ref;
    else
        ILD_sd_template =kv.ILD_sd_temp;
    end
    gamma = 1- b_lamda * ILD_sd_template/ILD_sd_reference;  
    delta_xi = w * SG_deviation(:,:,1) + (1-w) * SG_deviation(:,:,2);
    delta_m = gamma * ( b_ILD * ILD_deviation + b_SG * delta_xi) + b_ILDTSD * norm_ILD_std_mean;

    E = 2 * exp(-delta_m) + 1;
    
end


