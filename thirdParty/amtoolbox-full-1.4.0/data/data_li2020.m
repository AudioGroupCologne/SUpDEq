function varargout = data_li2020(varargin)
%data_li2020 Load experimental results from Li et al. (2020)
%
%   Usage: data = data_li2020(dataFlag,measure)
%          data = data_li2020(fig)
%
%   Output parameters:
%     data    : structure containing either HRTFs/BRIRs or externalization results
%
%   DATA_LI2020 provides individually measured/modified HRTFs/BRIRs and
%   experimental results from Li et al. (2020). Use the fig flag 
%   to obtain data shown in figures from Li et al. (2020).
%
%   The dataFlag flag may be used to choose between HRTFs and various
%   experimental results:
%
%     'hrtf'   HRTFs used in all experiments.
%
%     'exp1'    Experimental results of Exp. A. Default.
%
%     'exp2'    Experimental results of Exp. B.
%
%     'exp3'    Experimental results of Exp. C.
%
%     'exp4'    Experimental results of Exp. D.
%
%     'exp5'    Experimental results of Exp. E.
%
%
%   Additional flags may be:
%
%     'plot'    Plot results as published.
%
%     'no_plot'  No plots. Default.
%
%
%   For experiments 1 to 4, data contains HRTFs/BRIRs (exp1 - exp4 (size of out matrix): 5x256x5x2 <---> Nr. subjects x HRIR length x conditions x left/right, for experiment 5 (size of out matrix): 5x16384x5x5x2 <---> Nr. subjects x BRIR length x smoothing condition x compression condition x left/right ).
%   The externalization results comprise (mean, median, and 95% CI).      
%
%
%   Examples:
%   ---------
%
%   To display results of experiment A:
%
%     data_li2020('exp1','plot');
%
%   To display results of experiment B:
%
%     data_li2020('exp2','plot');
%
%   To display results of experiment C:
%
%     data_li2020('exp3','plot');
%
%   To display results of experiment D:
%
%     data_li2020('exp4','plot');
%
%   To display results of experiment E:
%
%     data_li2020('exp5','plot');
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_li2020.php


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

%   References:
%     S. Li, R. Baumgartner, and J. Peissig. 
%     Modeling perceived externalization of a static, lateral sound image.
%     Acta Acust.,4(5) (2020) 
%     
%     References
%     1. https://doi.org/10.1051/aacus/2020020 
    


% definput.import={'amt_cache'};
definput.flags.expResults = {'exp1','exp2','exp3','exp4','exp5','hrtf'};
definput.flags.plot = {'no_plot','plot'};

[flags]=ltfatarghelper({},definput,lower(varargin));

%% HRTFs
if flags.do_hrtf
    %%
    % load original HRTFs and BRIRs and + processing stage
    %%
    varargout{1} = amt_load('li2020', 'modified_BRIR_HRTF_dataset.mat');
    return
  end

  %% Experimental Results
  if flags.do_exp1 || flags.do_exp2 || flags.do_exp3 || flags.do_exp4 || flags.do_exp5
      
      x_conditions =1:5;

      Results = amt_load('li2020', 'ExpResults_LI2020.mat');   
      if flags.do_exp1 
        varargout{1} = Results.E_result_Exp1;
        %% Plot
        if flags.do_plot
          out.fig = figure;
            
            % median
            e_exp1_bb_median = plot(x_conditions-0.05,  Results.E_result_Exp1.E_modified_Exp1_bb.median','ko-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k'); hold on   
            e_exp1_low_median = plot(x_conditions,      Results.E_result_Exp1.E_modified_Exp1_low.median','ks-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k'); 
            e_exp1_high_median = plot(x_conditions+0.05,Results.E_result_Exp1.E_modified_Exp1_high.median','kd-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k'); 

            % CI 95
            e_exp1_bb=errorbar(x_conditions-0.05, Results.E_result_Exp1.E_modified_Exp1_bb.mean, Results.E_result_Exp1.E_modified_Exp1_bb.CI95,'k','LineStyle','none','LineWidth',1.5); 
            e_exp1_low=errorbar(x_conditions, Results.E_result_Exp1.E_modified_Exp1_low.mean, Results.E_result_Exp1.E_modified_Exp1_low.CI95,'k','LineStyle','none','LineWidth',1.5); 
            e_exp1_high=errorbar(x_conditions+0.05, Results.E_result_Exp1.E_modified_Exp1_high.mean, Results.E_result_Exp1.E_modified_Exp1_high.CI95,'k','LineStyle','none','LineWidth',1.5); 


            duration=[{'0'}; {'5'}; {'10'}; {'15'}; {'20'}];
            xticks([1 2 3 4 5]);
            set(gca,'XTickLabel',duration); 

            yticks([0 1 2 3]);
            set(gca,'yTickLabel',[{'0'}; {'1'}; {'2'}; {'3'}]); 

            xlim([0.5 5.5]), ylim([-0.1 3.1]),
            xlabel('ILD expansion [dB]');
            ylabel('Externalization rating')

            legend('BB','LO','HI')  
            set(gca, 'FontSize',26, 'FontName', 'Times New Roman')

        end
      end

      
    if flags.do_exp2 
        varargout{1} = Results.E_result_Exp2;
        %% Plot
        if flags.do_plot
          out.fig = figure;

          
            e_exp2_bb_median = plot(x_conditions, Results.E_result_Exp2.median','ko-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k'); hold on
            e_exp2_bb=errorbar(x_conditions, Results.E_result_Exp2.mean, Results.E_result_Exp2.CI95, 'k','LineStyle','none','LineWidth',1.5); hold on,


            duration=[{'0'}; {'1'}; {'4'}; {'16'}; {'64'}];
            xticks([1 2 3 4 5]);
            set(gca,'XTickLabel',duration); 
            yticks([0 1 2 3]);
            set(gca,'yTickLabel',[{'0'}; {'1'}; {'2'}; {'3'}]); 

            xlim([0.5 5.5]), ylim([-0.1 3.1]),

            xlabel('Spectral smoothing [ERB]');
            ylabel('Externalization rating')

            legend('measured')  
            set(gca, 'FontSize',26, 'FontName', 'Times New Roman')

        end
     end     
      
    if flags.do_exp3 
        varargout{1} = Results.E_result_Exp3;
        %% Plot
        if flags.do_plot
          out.fig = figure;

            % median values
            e_exp3_bb_ipsi_median = plot(x_conditions-0.05,Results.E_result_Exp3.E_modified_ipsi_Exp3.median','ko-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k'); hold on
            e_exp3_bb_contra_median = plot(x_conditions+0.05,Results.E_result_Exp3.E_modified_contra_Exp3.median','ks-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k'); 

            % CI 95
            e_exp3_bb_ipsi=errorbar(x_conditions-0.05, Results.E_result_Exp3.E_modified_ipsi_Exp3.mean, Results.E_result_Exp3.E_modified_ipsi_Exp3.CI95,'k','LineStyle','none','LineWidth',1.5); 
            e_exp3_bb_contra=errorbar(x_conditions+0.05, Results.E_result_Exp3.E_modified_contra_Exp3.mean, Results.E_result_Exp3.E_modified_contra_Exp3.CI95,'k','LineStyle','none','LineWidth',1.5);

            duration=[{'0'}; {'25'}; {'50'}; {'75'}; {'100'}];
            xticks([1 2 3 4 5]);
            set(gca,'XTickLabel',duration); 
            yticks([0 1 2 3]);
            set(gca,'yTickLabel',[{'0'}; {'1'}; {'2'}; {'3'}]); 
            xlim([0.5 5.5]), ylim([-0.1 3.1]),


            xlabel('Compressed ILD contrast [%]');
            ylabel('Externalization rating')


            legend('ipsi','contra')  
            set(gca, 'FontSize',26, 'FontName', 'Times New Roman')

        end
    end           
      
  if flags.do_exp4 
        varargout{1} = Results.E_result_Exp4;
        %% Plot
        if flags.do_plot
          out.fig = figure;

            e_exp4_bb_median = plot(x_conditions-0.05, Results.E_result_Exp4.E_modified_Exp4_bb.median','ko-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k'); hold on
            e_exp4_ipsi_median = plot(x_conditions, Results.E_result_Exp4.E_modified_Exp4_ipsi.median','ks-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k');
            e_exp4_contra_median = plot(x_conditions+0.05, Results.E_result_Exp4.E_modified_Exp4_contra.median','kd-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k'); 


            e_exp4_bb=errorbar(x_conditions-0.05, Results.E_result_Exp4.E_modified_Exp4_bb.mean, Results.E_result_Exp4.E_modified_Exp4_bb.CI95,'k','LineStyle','none','LineWidth',1.5);
            e_exp4_ipsi=errorbar(x_conditions, Results.E_result_Exp4.E_modified_Exp4_ipsi.mean, Results.E_result_Exp4.E_modified_Exp4_ipsi.CI95,'k','LineStyle','none','LineWidth',1.5); 
            e_exp4_contra=errorbar(x_conditions+0.05, Results.E_result_Exp4.E_modified_Exp4_contra.mean, Results.E_result_Exp4.E_modified_Exp4_contra.CI95,'k','LineStyle','none','LineWidth',1.5); 


         
            legend('bi','ipsi','contra');

            duration=[{'0'}; {'1'}; {'4'}; {'16'}; {'64'}];
            xticks([1 2 3 4 5]);
            set(gca,'XTickLabel',duration); 
            yticks([0 1 2 3]);
            set(gca,'yTickLabel',[{'0'}; {'1'}; {'2'}; {'3'}]); 
            xlim([0.5 5.5]), ylim([-0.1 3.1]),
            set(gca, 'FontSize',26, 'FontName', 'Times New Roman')

            xlabel('Spectral smoothing [ERB]');
            ylabel('Externalization rating')

        end
  end                
      
 
  if flags.do_exp5 
        varargout{1} = Results.E_result_Exp5;
        %% Plot
        if flags.do_plot
          out.fig = figure;

            e_exp5_bb1_median=plot(x_conditions-0.1,  Results.E_result_Exp5.E_modified_bb1_Exp5.median',  'ko-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k'); hold on,
            e_exp5_bb2_median=plot(x_conditions-0.05, Results.E_result_Exp5.E_modified_bb2_Exp5.median', 'ks-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k');
            e_exp5_bb3_median=plot(x_conditions,      Results.E_result_Exp5.E_modified_bb3_Exp5.median', 'kd-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k');
            e_exp5_bb4_median=plot(x_conditions+0.05, Results.E_result_Exp5.E_modified_bb4_Exp5.median', 'kv-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k');
            e_exp5_bb5_median=plot(x_conditions+0.1,  Results.E_result_Exp5.E_modified_bb5_Exp5.median', 'k^-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerFaceColor','k');

            
            e_exp5_bb1=errorbar(x_conditions-0.1,  Results.E_result_Exp5.E_modified_bb1_Exp5.mean, Results.E_result_Exp5.E_modified_bb1_Exp5.CI95,'k','LineStyle','none','LineWidth',1.5); 
            e_exp5_bb2=errorbar(x_conditions-0.05, Results.E_result_Exp5.E_modified_bb2_Exp5.mean, Results.E_result_Exp5.E_modified_bb2_Exp5.CI95,'k','LineStyle','none','LineWidth',1.5);
            e_exp5_bb3=errorbar(x_conditions,      Results.E_result_Exp5.E_modified_bb3_Exp5.mean, Results.E_result_Exp5.E_modified_bb3_Exp5.CI95,'k','LineStyle','none','LineWidth',1.5);
            e_exp5_bb4=errorbar(x_conditions+0.05, Results.E_result_Exp5.E_modified_bb4_Exp5.mean, Results.E_result_Exp5.E_modified_bb4_Exp5.CI95,'k','LineStyle','none','LineWidth',1.5);
            e_exp5_bb5=errorbar(x_conditions+0.1,  Results.E_result_Exp5.E_modified_bb5_Exp5.mean, Results.E_result_Exp5.E_modified_bb5_Exp5.CI95,'k','LineStyle','none','LineWidth',1.5);


         
            legend('B = 0','B = 1', 'B = 4','B = 16', 'B = 64');

            duration=[{'0'}; {'25'}; {'50'}; {'75'}; {'100'}];
            xticks([1 2 3 4 5]);
            set(gca,'XTickLabel',duration); 
            yticks([0 1 2 3]);
            set(gca,'yTickLabel',[{'0'}; {'1'}; {'2'}; {'3'}]); 

            xlim([0.5 5.5]), ylim([-0.1 3.1]),

            xlabel('Reverberation reduction [%]');
            ylabel('Externalization rating')
            set(gca, 'FontSize',26, 'FontName', 'Times New Roman')

        end
   end                
     
     
 
end
  




