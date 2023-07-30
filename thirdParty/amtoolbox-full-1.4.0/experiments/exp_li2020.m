function exp_li2020(varargin)
%EXP_LI2020 experimental results for externalization model from Li et al. (2020)
%
%   EXP_LI2020(flag) demonstrates how to apply the model to estimate the degree of externalization 
%   by comparing template and target HRTFs
%
%   Examples:
%   ---------
%   To display measured and predicted results of experiment A use :
%     exp_li2020('fig2');
%
%   To display measured and predicted results of experiment B use :
%     exp_li2020('fig3');
%
%   To display measured and predicted results of experiment C use :
%     exp_li2020('fig5');
%
%   To display measured and predicted results of experiment D use :
%     exp_li2020('fig6');
%
%   To display measured and predicted results of experiment E use :
%     exp_li2020('fig7');
%
%   Requirements:
%   -------------
%
%   1) Data in hrtf/li2020 and auxdata/li2020 (noise_input,..., ExpResults_LI2020, modified_BRIR_HRTF_dataset)
%
%   References:
%     S. Li, R. Baumgartner, and J. Peissig. Modeling perceived
%     externalization of a static, lateral sound image. Acta Acustica, 4(5),
%     2021.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_li2020.php


%   #Author: Song Li, Institute of Communications Technology, Leibniz
%   University of Hannover, Germany

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   To make the scripts clear, the general weighting factors are used for
%   calculating the fitted and predicted reults. The fitted results are
%   slightly different than the results in the paper, since the fitted results
%   shown in the paper are calcualted based on the individual weighting
%   factors.
 
%% check input    
definput.keyvals.Obj = [];
definput.flags.experiment = {'fig2','fig3','fig5','fig6','fig7'};
[flags,kv]=ltfatarghelper({'Obj'},definput,varargin);

%% Get listener's data
    %stimulus = load('AMT Implementation\data\noise_input.mat');
    stimulus = amt_load('li2020', 'noise_input.mat');
    data = data_li2020('hrtf');   % load HRTFs of listener pool
    Nsubj = 5;                    % 5 listeners
    fs = 44100;

    %% Run model to get externalization scores for every conditions/experiment
if flags.do_fig2
    %% Experiment 1
    ILD_expansion = 5;
    E_predict_exp1=zeros(5,5,3); % 5 subjects x 5 ILD expansion levels x 3 conditions
    for ids = 1:Nsubj
        % condition 1: broadband BB
        for k =1:ILD_expansion
            template = squeeze(data.hrir_bb_modified_Exp1(ids,:,1,:));
            target = squeeze(data.hrir_bb_modified_Exp1(ids,:,k,:));
            [E_predict_exp1(ids,k,1)] = li2020(target,template,[],[],'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
        end

        % condition 2: Low frequency LO
        for k =1:ILD_expansion
            template = squeeze(data.hrir_low_modified_Exp1(ids,:,1,:));
            target = squeeze(data.hrir_low_modified_Exp1(ids,:,k,:));
            [E_predict_exp1(ids,k,2)] = li2020(target,template,[],[],'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
        end

        % condition 3: High frequency HI    
        for k =1:ILD_expansion
            template = squeeze(data.hrir_high_modified_Exp1(ids,:,1,:));
            target = squeeze(data.hrir_high_modified_Exp1(ids,:,k,:));
            [E_predict_exp1(ids,k,3)] = li2020(target,template,[],[],'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
        end  

    end

     CI95 = [ bsxfun(@times, std(E_predict_exp1(:,:,1))/sqrt(length(E_predict_exp1(:,:,1))),...
              tinv([0.025 0.975], length(E_predict_exp1(:,:,1))-1)');
              bsxfun(@times, std(E_predict_exp1(:,:,2))/sqrt(length(E_predict_exp1(:,:,2))),...
              tinv([0.025 0.975], length(E_predict_exp1(:,:,2))-1)');
              bsxfun(@times, std(E_predict_exp1(:,:,3))/sqrt(length(E_predict_exp1(:,:,3))),...
              tinv([0.025 0.975], length(E_predict_exp1(:,:,3))-1)')];    
    CI95 = CI95(2:2:6,:);
    
    data_li2020('exp1','plot');
    x_conditions = 1:5;
    legend off
    
    errorbar(x_conditions-0.05, mean(E_predict_exp1(:,:,1)), CI95(1,:),'r','LineStyle','none','LineWidth',1.5); hold on;
    errorbar(x_conditions, mean(E_predict_exp1(:,:,2)), CI95(2,:),'r','LineStyle','none','LineWidth',1.5); 
    errorbar(x_conditions+0.05, mean(E_predict_exp1(:,:,3)), CI95(3,:),'r','LineStyle','none','LineWidth',1.5); 
    
    plot(x_conditions-0.05, median(E_predict_exp1(:,:,1)),'ro--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r'); 
    plot(x_conditions,median(E_predict_exp1(:,:,2)),'rs--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r'); 
    plot(x_conditions+0.05,median(E_predict_exp1(:,:,3)),'rd--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r'); 

   


end

if flags.do_fig3
    %% Experiment 2
    bandwidth = 5;
    E_predict_exp2=zeros(5,5); % 5 subjects x 5 conditions
    for ids = 1:Nsubj

        for k =1:bandwidth
            template = squeeze(data.hrir_modified_Exp2(ids,:,1,:));
            target = squeeze(data.hrir_modified_Exp2(ids,:,k,:));
            [E_predict_exp2(ids,k)] = li2020(target,template,[],[],'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
        end

    end
  
                 
    CI95 = bsxfun(@times, std(E_predict_exp2)/sqrt(length(E_predict_exp2)),...
        tinv([0.025 0.975], length(E_predict_exp2)-1)');  
    
    data_li2020('exp2','plot');
    x_conditions = 1:5;
    legend off
    
    errorbar(x_conditions, mean(E_predict_exp2), CI95(1,:),'r','LineStyle','none','LineWidth',1.5); hold on;
    plot(x_conditions,median(E_predict_exp2),'ro--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r'); 
    
 
end

if flags.do_fig5
    %% Experiment 3
    compressed_ILD = 5;
    E_predict_exp3=zeros(5,5,2); % 5 subjects x 5 ILD compression levels x 2 conditions
    for ids = 1:Nsubj
        % condition 1: ipsi
        for k =1:compressed_ILD
            template = squeeze(data.hrir_ipsi_modified_Exp3(ids,:,1,:));
            target = squeeze(data.hrir_ipsi_modified_Exp3(ids,:,k,:));
            [E_predict_exp3(ids,k,1)] = li2020(target,template,[],[],'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
        end

        % condition 2: contra
        for k =1:compressed_ILD
            template = squeeze(data.hrir_contra_modified_Exp3(ids,:,1,:));
            target = squeeze(data.hrir_contra_modified_Exp3(ids,:,k,:));
            [E_predict_exp3(ids,k,2)] = li2020(target,template,[],[],'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
        end

    end


       CI95 = [ bsxfun(@times, std(E_predict_exp3(:,:,1))/sqrt(length(E_predict_exp3(:,:,1))),...
              tinv([0.025 0.975], length(E_predict_exp3(:,:,1))-1)');
              bsxfun(@times, std(E_predict_exp3(:,:,2))/sqrt(length(E_predict_exp3(:,:,2))),...
              tinv([0.025 0.975], length(E_predict_exp3(:,:,2))-1)')];    
        CI95 = CI95(2:2:4,:);
    
    
    
    data_li2020('exp3','plot');
    x_conditions = 1:5;
    legend off
    
    errorbar(x_conditions-0.05, mean(E_predict_exp3(:,:,1)), CI95(1,:),'r','LineStyle','none','LineWidth',1.5); hold on;
    errorbar(x_conditions+0.05, mean(E_predict_exp3(:,:,2)), CI95(2,:),'r','LineStyle','none','LineWidth',1.5); 
    plot(x_conditions-0.05,median(E_predict_exp3(:,:,1)),'ro--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r'); 
    plot(x_conditions+0.05,median(E_predict_exp3(:,:,2)),'rs--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r'); 

end

if flags.do_fig6
    %% Experiment 4
    bandwidth = 5;
    E_predict_exp4=zeros(5,5,3); % 5 subjects x 5 smoothing levels x 3 conditions
    for ids = 1:Nsubj
        % condition 1: bilateral
        for k =1:bandwidth
            template = squeeze(data.hrir_both_modified_Exp4(ids,:,1,:));
            target = squeeze(data.hrir_both_modified_Exp4(ids,:,k,:));
            [E_predict_exp4(ids,k,1)] = li2020(target,template,[],[],'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
        end

        % condition 2: ipsilateral
        for k =1:bandwidth
            template = squeeze(data.hrir_ipsi_modified_Exp4(ids,:,1,:));
            target = squeeze(data.hrir_ipsi_modified_Exp4(ids,:,k,:));
            [E_predict_exp4(ids,k,2)] = li2020(target,template,[],[],'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
        end

        % condition 3: contralateral  
        for k =1:bandwidth
            template = squeeze(data.hrir_contra_modified_Exp4(ids,:,1,:));
            target = squeeze(data.hrir_contra_modified_Exp4(ids,:,k,:));
            [E_predict_exp4(ids,k,3)] = li2020(target,template,[],[],'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
        end  

    end

    CI95 = [ bsxfun(@times, std(E_predict_exp4(:,:,1))/sqrt(length(E_predict_exp4(:,:,1))),...
              tinv([0.025 0.975], length(E_predict_exp4(:,:,1))-1)');
              bsxfun(@times, std(E_predict_exp4(:,:,2))/sqrt(length(E_predict_exp4(:,:,2))),...
              tinv([0.025 0.975], length(E_predict_exp4(:,:,2))-1)');
              bsxfun(@times, std(E_predict_exp4(:,:,3))/sqrt(length(E_predict_exp4(:,:,3))),...
              tinv([0.025 0.975], length(E_predict_exp4(:,:,3))-1)')];    
    CI95 = CI95(2:2:6,:);
    
    
    data_li2020('exp4','plot');
    x_conditions = 1:5;
    legend off
   
    
    errorbar(x_conditions-0.05, mean(E_predict_exp4(:,:,1)), CI95(1,:),'r','LineStyle','none','LineWidth',1.5); hold on
    errorbar(x_conditions, mean(E_predict_exp4(:,:,2)), CI95(2,:),'r','LineStyle','none','LineWidth',1.5);
    errorbar(x_conditions+0.05, mean(E_predict_exp4(:,:,3)), CI95(3,:),'r','LineStyle','none','LineWidth',1.5); 
    
    plot(x_conditions-0.05,median(E_predict_exp4(:,:,1)),'ro--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r'); 
    plot(x_conditions,median(E_predict_exp4(:,:,2)),'rs--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r'); 
    plot(x_conditions+0.05,median(E_predict_exp4(:,:,3)),'rd--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r'); 


end

if flags.do_fig7

    %% Experiment 5
    reverbLevel = 5;
    % reverbLevel =1: original BRIRs
    % reverbLevel =5: Anechoic condition/HRTF
    E_predict_exp5=zeros(5,5,5); % 5 subjects x 5 bandwidth factors x 5 reverb levels 
    for ids = 1:Nsubj
        % condition 1: B = 0
        for k =1:reverbLevel
            template = squeeze(data.brir_both_modified_Exp5(ids,:,1,1,:));  % with given BRIRs
            target = squeeze(data.brir_both_modified_Exp5(ids,:,1,k,:));
            
             template_dir = squeeze(data.brir_both_modified_Exp5(ids,:,1,reverbLevel,:));  % with given direct sound part
             target_dir = squeeze(data.brir_both_modified_Exp5(ids,:,1,reverbLevel,:));   
            
            [E_predict_exp5(ids,1,k)] = li2020(target_dir,template_dir,target,template,'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);
            
        %    [E_predict_exp5(ids,1,k)] = Li2020([],[],target,template,'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
        %        'flow',200, 'fhigh',16000, 'space',1);  % option
            
            
        end

        % condition 2: B = 1
        for k =1:reverbLevel
            template = squeeze(data.brir_both_modified_Exp5(ids,:,1,1,:)); % with given BRIRs
            target = squeeze(data.brir_both_modified_Exp5(ids,:,2,k,:));
            
             template_dir = squeeze(data.brir_both_modified_Exp5(ids,:,1,reverbLevel,:)); % with given direct sound part
             target_dir = squeeze(data.brir_both_modified_Exp5(ids,:,2,reverbLevel,:));

            [E_predict_exp5(ids,2,k)] = li2020(target_dir,template_dir,target,template,'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);  
            
        %    [E_predict_exp5(ids,2,k)] = Li2020([],[],target,template,'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
        %        'flow',200, 'fhigh',16000, 'space',1);  % option
            
        end

        % condition 3: B = 4
        for k =1:reverbLevel
            template = squeeze(data.brir_both_modified_Exp5(ids,:,1,1,:)); % with given BRIRs
            target = squeeze(data.brir_both_modified_Exp5(ids,:,3,k,:));
            
            template_dir = squeeze(data.brir_both_modified_Exp5(ids,:,1,reverbLevel,:)); % with given direct sound part
            target_dir = squeeze(data.brir_both_modified_Exp5(ids,:,3,reverbLevel,:));
            
            [E_predict_exp5(ids,3,k)] = li2020(target_dir,template_dir,target,template,'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
            
        %    [E_predict_exp5(ids,3,k)] = Li2020([],[],target,template,'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
        %        'flow',200, 'fhigh',16000, 'space',1);  % option
            
            
        end

        % condition 4: B = 16
        for k =1:reverbLevel 
            template = squeeze(data.brir_both_modified_Exp5(ids,:,1,1,:)); % with given BRIRs
            target = squeeze(data.brir_both_modified_Exp5(ids,:,4,k,:));
          
             template_dir = squeeze(data.brir_both_modified_Exp5(ids,:,1,reverbLevel,:)); % with given direct sound part
             target_dir = squeeze(data.brir_both_modified_Exp5(ids,:,4,reverbLevel,:));
            
             [E_predict_exp5(ids,4,k)] = li2020(target_dir,template_dir,target,template,'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                 'flow',200, 'fhigh',16000, 'space',1);   
          
        %    [E_predict_exp5(ids,4,k)] = Li2020([],[],target,template,'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
        %        'flow',200, 'fhigh',16000, 'space',1);  % option
             

        end

        % condition 5: B = 64
        for k =1:reverbLevel
            template = squeeze(data.brir_both_modified_Exp5(ids,:,1,1,:)); % with given BRIRs
            target = squeeze(data.brir_both_modified_Exp5(ids,:,5,k,:));
            
             template_dir = squeeze(data.brir_both_modified_Exp5(ids,:,1,reverbLevel,:)); % with given direct sound part
             target_dir = squeeze(data.brir_both_modified_Exp5(ids,:,5,reverbLevel,:));
            
            [E_predict_exp5(ids,5,k)] = li2020(target_dir,template_dir,target,template,'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
                'flow',200, 'fhigh',16000, 'space',1);   
       
        %    [E_predict_exp5(ids,5,k)] = Li2020([],[],target,template,'stim',stimulus.noise_input,'fs',fs,'fsstim',fs, ...
        %        'flow',200, 'fhigh',16000, 'space',1);  % option
            
        end

    end

   CI95 =[bsxfun(@times, std(squeeze(E_predict_exp5(:,1,:)))/sqrt(length(E_predict_exp5(:,1,:))),...
          tinv([0.025 0.975], length(E_predict_exp5(:,1,:))-1)');
          bsxfun(@times, std(squeeze(E_predict_exp5(:,2,:)))/sqrt(length(E_predict_exp5(:,2,:))),...
          tinv([0.025 0.975], length(E_predict_exp5(:,2,:))-1)');
          bsxfun(@times, std(squeeze(E_predict_exp5(:,3,:)))/sqrt(length(E_predict_exp5(:,3,:))),...
          tinv([0.025 0.975], length(E_predict_exp5(:,3,:))-1)');
          bsxfun(@times, std(squeeze(E_predict_exp5(:,4,:)))/sqrt(length(E_predict_exp5(:,4,:))),...
          tinv([0.025 0.975], length(E_predict_exp5(:,4,:))-1)');
          bsxfun(@times, std(squeeze(E_predict_exp5(:,5,:)))/sqrt(length(E_predict_exp5(:,5,:))),...
          tinv([0.025 0.975], length(E_predict_exp5(:,5,:))-1)')];    
    CI95 = CI95(2:2:10,:);
    
    
    data_li2020('exp5','plot');
    x_conditions = 1:5;
    legend off
    
    errorbar(x_conditions-0.1,  mean(squeeze(E_predict_exp5(:,1,:))), CI95(1,:),'r', 'LineStyle','none','LineWidth',1.5); hold on,
    errorbar(x_conditions-0.05, mean(squeeze(E_predict_exp5(:,2,:))), CI95(2,:),'r', 'LineStyle','none','LineWidth',1.5);
    errorbar(x_conditions,      mean(squeeze(E_predict_exp5(:,3,:))), CI95(3,:),'r', 'LineStyle','none','LineWidth',1.5);
    errorbar(x_conditions+0.05, mean(squeeze(E_predict_exp5(:,4,:))), CI95(4,:),'r', 'LineStyle','none','LineWidth',1.5);
    errorbar(x_conditions+0.1,  mean(squeeze(E_predict_exp5(:,5,:))), CI95(5,:),'r', 'LineStyle','none','LineWidth',1.5);
    
    plot(x_conditions-0.1, median(squeeze(E_predict_exp5(:,1,:))),'ro--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r'); hold on,
    plot(x_conditions-0.05, median(squeeze(E_predict_exp5(:,2,:))),'rs--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r');
    plot(x_conditions, median(squeeze(E_predict_exp5(:,3,:))),'rd--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r');
    plot(x_conditions+0.05, median(squeeze(E_predict_exp5(:,4,:))),'rv--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r');
    plot(x_conditions+0.1, median(squeeze(E_predict_exp5(:,5,:))),'r^--', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','r');




end

end



