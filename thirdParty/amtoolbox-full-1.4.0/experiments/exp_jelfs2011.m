function output = exp_jelfs2011(varargin)
%EXP_JELFS2011  Figures from Jelfs et al. (2011)
%   Usage: output = exp_jelfs2011(flag);
%
%   EXP_JELFS2011 reproduces result from the paper Jelfs et al. (2011).
%
%   The following flags can be specified;
%
%     'plot'     plot the output of the experiment. This is the default.
%
%     'no_plot'   Don't plot, only return data.
%
%     'fig4'     Reproduce Fig. 4. The generated data is compared against
%                data from Hawley et at. (2004). 
%
%   Examples:
%   ---------
%
%   To display Figure 4 use :
%
%     exp_jelfs2011('fig4');
%
%   See also: jelfs2011, culling2004
%
%   References:
%     S. Jelfs, J. Culling, and M. Lavandier. Revision and validation of a
%     binaural model for speech intelligibility in noise. Hearing Research,
%     2011.
%     
%     M. Hawley, R. Litovsky, and J. Culling. The benefit of binaural hearing
%     in a cocktail party: Effect of location and type of interferer. J.
%     Acoust. Soc. Am., 115:833--843, 2004.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_jelfs2011.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: MATLAB
%   #Author: Peter Sondergaard (2013)
%   #Author: Matthieu Lavandier (2021)
%   #Author: Clara Hollomey (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% TODO: Extract the human data from this function and put it into a
% dedicated function, either data_jelfs2011 or data_hawley2004

definput.flags.type={'fig4'};
definput.flags.plot={'plot','no_plot'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_fig4
  bin_data = [-3.5 -10 -12.5 -10.5 -1.25 -2.25 -8 -7.25 1.5 0 -4.5 -5.5];
  mon_data = [-3.4688 -0.2813 -8.9063 -6.9375 -0.9063 0.9375 -5.2813 -4.1563 1.1875 1.2813 -2.7813 -1.9688 ];
  offset = [0 0 0 0 3.01 3.01 3.01 3.01 4.77 4.77 4.77 4.77]';
  bSNR = zeros(12,3);
  database = 'kemar';
  
  mode = 'both';
  
  bSNR(1,:)  = jelfs2011({0,database},{[0],database},mode);
  bSNR(2,:)  = jelfs2011({0,database},{[330],database},mode);
  bSNR(3,:)  = jelfs2011({0,database},{[60],database},mode);
  bSNR(4,:)  = jelfs2011({0,database},{[90],database},mode);
  bSNR(5,:)  = jelfs2011({0,database},{[0 0],database},mode);
  bSNR(6,:)  = jelfs2011({0,database},{[330 90],database},mode);
  bSNR(7,:)  = jelfs2011({0,database},{[60 90],database},mode);
  bSNR(8,:)  = jelfs2011({0,database},{[90 90],database},mode);
  bSNR(9,:)  = jelfs2011({0,database},{[0 0 0],database},mode);
  bSNR(10,:) = jelfs2011({0,database},{[330 60 90],database},mode);
  bSNR(11,:) = jelfs2011({0,database},{[30 60 90],database},mode);
  bSNR(12,:) = jelfs2011({0,database},{[90 90 90],database},mode);
  predictions1 = -bSNR(:,1)+offset;
    
  mode = 'left';
  bSNR(1,:)  = jelfs2011({0,database},{[0],database},mode);
  bSNR(2,:)  = jelfs2011({0,database},{[330],database},mode);
  bSNR(3,:)  = jelfs2011({0,database},{[60],database},mode);
  bSNR(4,:)  = jelfs2011({0,database},{[90],database},mode);
  bSNR(5,:)  = jelfs2011({0,database},{[0 0],database},mode);
  bSNR(6,:)  = jelfs2011({0,database},{[330 90],database},mode);
  bSNR(7,:)  = jelfs2011({0,database},{[60 90],database},mode);
  bSNR(8,:)  = jelfs2011({0,database},{[90 90],database},mode);
  bSNR(9,:)  = jelfs2011({0,database},{[0 0 0],database},mode);
  bSNR(10,:) = jelfs2011({0,database},{[330 60 90],database},mode);
  bSNR(11,:) = jelfs2011({0,database},{[30 60 90],database},mode);
  bSNR(12,:) = jelfs2011({0,database},{[90 90 90],database},mode);
  predictions2 = -bSNR(:,1)+offset;
  
  if flags.do_plot
    figure;
    hold all;
    scatter(predictions1(1:4),bin_data(1:4),'or');
    scatter(predictions2(1:4),mon_data(1:4),'ob');
    scatter(predictions1(5:8),bin_data(5:8),'^r');
    scatter(predictions2(5:8),mon_data(5:8),'^b');
    scatter(predictions1(9:12),bin_data(9:12),'hr');
    scatter(predictions2(9:12),mon_data(9:12),'hb');
    
    line([-14,5],[-14,6],'LineStyle','--');
    
    legend('1 int., Binaural',...
           '1 int., Monaural',...
           '2 int., Binaural',...
           '2 int., Monaural',...
           '3 int., Binaural',...
           '3 int., Monaural',...
           'Location','NorthWest');
    
    xlabel('Observed SRT (dB)');
    ylabel('Predicted SRT (dB)');
    
  end;
  output=[predictions1,predictions2];
end;

