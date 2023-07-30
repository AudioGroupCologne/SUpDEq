function data = data_lindemann1986(varargin)
%DATA_LINDEMANN1986 Data points from the Lindemann (1986a) paper
%   Usage: data = data_lindemann1986_data(flag)
%
%   Output parameters:
%       data    : the data points from the given figure as a matrix with x and y
%                 data or (in cases where no clear x values were given such as
%                 T/2) as a vector containing only the y data
%
%   DATA_LINDEMANN1986(flag) returns data points from the Lindemann
%   1986 paper. The flag may be one of:
%
%     'no_plot'        Don't plot, only return data. This is the default.
%
%     'plot'          Plot the data.
%
%     'fig11_yost'    Return data from Fig. 11. with condition "yost".
%                     The data is ILD vs. lateral displacement (10 is max
%                     Displacement).
%
%     'fig11_sayers'  Return data from Fig. 11. with condition "sayers". The
%                     data is ILD vs. lateral displacement (10 is max
%                     displacement).
%
%     'fig12_400'     Return data from Fig. 12. for the 400 Hz pure tone. The
%                     data is ILD vs. ITD.
%
%     'fig12_600'     Return data from Fig. 12. for the 600 Hz pure tone. The
%                     data is ILD vs. ITD. 
%
%     'fig13'         Return data from Fig. 13. The output data format is
%                     x-axis, -3dB, 3dB, 9dB, 15dB, 25dB. The data is ILD vs.
%                     ITD.
%
%     'fig16'         Return data from Fig. 16. The data is ILD vs. ITD.
%
%     'fig17'         Return data from Fig. 17. The output data format is
%                     x-axis, 0ms, 0.09ms, 0.18ms, 0.27ms. The data is ITD vs.
%                     ILD.
%
%   If no flag is given, the function will print the list of valid flags.
%
%   Examples:
%   ---------
%
%   Figure 11 with the "yost" condition can be displayed using :
%
%     data_lindemann1986('fig11_yost','plot');
%
%   Figure 11 with the "sayers" condition can be displayed using :
%
%     data_lindemann1986('fig11_sayers','plot');
%
%   Figure 12 for a 400 Hz pure tone can be displayed using :
%
%     data_lindemann1986('fig12_400','plot');
%
%   Figure 12 for a 600 Hz pure tone can be displayed using :
%
%     data_lindemann1986('fig12_600','plot');
%
%   Figure 13 can be displayed by using :
%
%     data_lindemann1986('fig13','plot');
%
%   Figure 16 can be displayed using :
%
%     data_lindemann1986('fig16','plot');
%
%   Figure 17 can be displayed using :
%
%     data_lindemann1986('fig17','plot');
% 
%   References:
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. I. Simulation of lateralization for
%     stationary signals. J. Acoust. Soc. Am., 80:1608--1622, 1986.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_lindemann1986.php


%   #Author: Hagen Wierstorf

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% TODO: explain Data in description;


%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.type={'missingflag',...
    'fig11_yost','fig11_sayers',...
    'fig12_400','fig12_600',...
    'fig13',...
    'fig16',...
    'fig17',...
    };
definput.flags.plot = {'no_plot','plot'};

% Parse input options
[flags,keyvals]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% ------ Data points from the paper ------------------------------------
% The following data points are guessed after the plots in the paper from
% Lindemann 1986.
%
% Data for the given figure

% --- 'fig11_yost' ---
if flags.do_fig11_yost
  data = [ 0  0.0
           3  2.1
           6  4.4
           9  6.4
           12  8.6
           15  9.2
           18 10.0 ];
  % Plotting
  if flags.do_plot
    plot(data(:,1),data(:,2)*0.058,'+');
    axis([0 25 0 0.8]);
    legend('500 Hz, Yost');
    xlabel('interaural level difference (dB)');
    ylabel('displacement of the centroid d');
  end;
end;

% --- 'fig11_sayers' ---
if flags.do_fig11_sayers
  data = [ 6  3.9
           9  4.9
           12  8.9 ];
  % Plotting
  if flags.do_plot
    plot(data(:,1),data(:,2)*0.058,'*');
    axis([0 25 0 0.8]);
    legend('600 Hz, Sayers');
    xlabel('interaural level difference (dB)');
    ylabel('displacement of the centroid d');
  end;
end;

% --- 'fig12_400' ---
if flags.do_fig12_400
  data = [ -1.00  0.4
           -0.83  3.3
           -0.67  6.4
           -0.50  8.3
           -0.33  6.7
           -0.17  3.3
           0.00  0.1 ];
  % Plotting
  if flags.do_plot
    plot(data(:,1),data(:,2),'+');
    axis([-1 0 0 10]);
    legend('Trading experiment of Young, 400 Hz');
    xlabel('interaural time difference (ms)');
    ylabel('interaural level difference (dB)');
  end;
end;

% --- 'fig12_600' ---
if flags.do_fig12_600
  data = [ -1.00  0.6
           -0.83  2.6
           -0.67  5.8
           -0.50  7.2
           -0.33  7.0
           -0.17  2.6
           0.00  0.4 ];
  % Plotting
  if flags.do_plot
    plot(data(:,1),data(:,2),'*');
    axis([-1 0 0 10]);
    legend('Trading experiment of Young, 600 Hz');
    xlabel('interaural time difference (ms)');
    ylabel('interaural level difference (dB)');
  end;
end

% --- 'fig13' ---
if flags.do_fig13
  % data format: x-axis, -3dB, 3dB, 9dB, 15dB, 25dB
  data = [ -1.0 -10  18  29  37  38
           -0.9 -10   9  25  31  38
           -0.7 -10  -2  11  27  33
           -0.5 -10  -9   4  15  30
           -0.3 -10  -7   3  13  29
           -0.1  -7  -1   6  14  29
           0.1   1   8  13  19  30
           0.3  10  14  19  26  32
           0.5  14  20  26  31  34
           0.7  15  24  30  35  36
           0.9  11  25  31  34  39
           1.0   2  19  30  37  38 ];
  % Plotting
  if flags.do_plot
    plot(data(:,1),data(:,2),'x-r', ...   % -3dB
         data(:,1),data(:,3),'x-b', ...   %  3dB
         data(:,1),data(:,4),'x-g', ...   %  9dB
         data(:,1),data(:,5),'x-b', ...   % 15dB
         data(:,1),data(:,6),'x-r');      % 25dB
    axis([-1 1 -10 40]);
    set(gca,'XTick',-1:0.4:1);
    legend('25dB','15dB','9dB','3dB','-3dB');
    xlabel('interaural time difference (ms)');
    ylabel('interaural level difference (dB)');
  end;
end

% --- 'fig16' ---
if flags.do_fig16
  data = [ -1.000  0.0
           -0.875  3.3
           -0.750  4.9
           -0.625  6.6
           -0.500  7.2
           -0.375  7.5
           -0.250  9.0
           -0.125  9.5
           0.000  11.1 ];
  % Plotting
  if flags.do_plot
    plot(data(:,1),data(:,2),'+');
    axis([-1 0 0 15]);
    legend('f = 500 Hz');
    xlabel('interaural time difference (ms)');
    ylabel('interaural level difference (dB)');
  end;
end

% --- 'fig17' ---
if flags.do_fig17
  % data format: x-axis, 0ms, 0.09ms, 0.18ms, 0.27ms
  data = [ -9 -0.18 -0.09  0.00  0.09
           -6 -0.14 -0.05  0.05  0.14
           -3 -0.09  0.00  0.10  0.18
           0  0.00  0.09  0.17  0.27
           3  0.08  0.17  0.23  0.35
           6  0.10  0.23  0.28  0.37
           9  0.17  0.25  0.34  0.42 ];
  % Plotting
  if flags.do_plot
    plot(data(:,1),data(:,5),'x-r', ...   %  0 ms
         data(:,1),data(:,4),'x-b', ...   %  0.09 ms
         data(:,1),data(:,3),'x-g', ...   %  0.18 ms
         data(:,1),data(:,2),'x-b');      %  0.27 ms
    axis([-9 9 -0.36 0.72]);
    set(gca,'XTick',-9:3:9);
    legend('0.27 ms','0.18 ms','0.09 ms','0 ms');
    xlabel('interaural level difference (dB)');
    ylabel('interaural time difference (ms)');
  end;
end;


