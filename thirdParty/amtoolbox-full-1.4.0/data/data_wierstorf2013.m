function [data,description] = data_wierstorf2013(varargin)
%DATA_WIERSTORF2013 Data points from the Wierstorf (2013) book chapter
%   Usage: data = data_wierstorf2013_data(flag)
%
%   Output parameters:
%       data        : the data points from the given figure as columns of a
%                     matrix. The x axis is alway the first column, followed by
%                     data and for example confident intervals. Which kind of
%                     data is stored in which column is indicated by the
%                     description string.
%       description : string that describes the type and unit of the data
%                     represented in every column of the data matrix
%
%   DATA_WIERSTORF2013(flag) returns data points from the Wierstorf 2013 book
%   chapter. The flag may be one of:
%
%     'no_plot'        Don't plot, only return data. This is the default.
%
%     'plot'          Plot the data.
%
%     'itd2angle_lookuptable' Return the data for the ITD-to-angle look-up.
%
%     'fig6'          Return data from Fig. 6. The data describes the difference
%                     between the azimuth angle of the auditory event and the
%                     azimuth angle of the sound event for different incidence
%                     angles of the sound event. The data was collected with
%                     real loudspeakers as sound events and for via binaural
%                     synthesis simulated loudspeakers as sound event.
%
%     'fig7'          Return data from Fig. 7. The data descirbes the result of
%                     a listening test measuring the localization performance
%                     for Wave Field Synthesis setups at 16 different positions
%                     in the listening area. The data comes for three different
%                     loudspeaker arrays consisting of  3, 8, or 15 loudspeakers.
%                     Note: these are the same data as from Fig. 10.
%
%     'fig10'         Return data from Fig. 10. The data descirbes the result of
%                     a listening test measuring the localization performance
%                     for Wave Field Synthesis setups at 16 different positions
%                     in the listening area. The data comes for three different
%                     loudspeaker arrays consisting of  3, 8, or 15 loudspeakers.
%                     Note: these are the same data as from Fig. 7.
%
%   If no flag is given, the function will print the list of valid flags.
%
%   Examples:
%   ---------
%
%   Figure 6 can be displayed using :
%
%     data_wierstorf2013('fig6','plot');
%
%   Figure 7 can be displayed using :
%
%     data_wierstorf2013('fig7','plot');
%
%   Figure 10 can be displayed using :
%
%     data_wierstorf2013('fig10','plot');
%
%   References:
%     H. Wierstorf, A. Raake, and S. Spors. Binaural assessment of
%     multi-channel reproduction. In J. Blauert, editor, The technology of
%     binaural listening, chapter 10. Springer, Berlin--Heidelberg--New York
%     NY, 2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_wierstorf2013.php


%   #Author: Hagen Wierstorf

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.type={'missingflag','fig6','fig7','fig10','itd2angle_lookuptable'};
definput.flags.plot = {'no_plot','plot'};

% Parse input options
[flags,keyvals]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;
%% Load ITD-to-angle look-up

if flags.do_itd2angle_lookuptable
  data=amt_load('wierstorf2013','itd2angle_lookuptable.mat');
end

%% ------ Data points from the paper ------------------------------------
% The following data points are the original data points from the plots in the book
% chapter from Wierstorf 2013.
%
% Data for the given figure
% --- 'fig6' ---
if flags.do_fig6
    description = { ...
        'all','phi_{sound event}/deg'; ...
        'loudspeaker','phi_{auditory event}-phi_{sound event}/deg'; ...
        'loudspeaker','confidence interval'; ...
        'binaural synthesis','phi_{auditory event}-phi_{sound event}/deg'; ...
        'binaural synthesis','confidence interval';
    };
    data = [ 41.60 -4.60 0.82 -2.21 2.81
             34.80 -2.84 1.13 -0.62 1.36
             21.80 -2.01 1.30  0.04 1.27
             11.40 -0.98 0.81 -0.85 1.47
              5.70 -1.46 0.97 -1.31 0.84
              0.00 -0.31 0.77 -0.67 0.76
            -11.40  2.57 0.97 -0.56 0.84
            -16.80  2.54 0.90 -1.01 0.99
            -31.20  1.66 1.42  0.77 1.48
            -38.90  2.91 1.01  2.57 1.45
            -42.20  2.51 1.33  3.86 2.23
    ];
    % Plotting
    if flags.do_plot
        figure;
        subplot(2,1,1);
        errorbar(data(:,1),data(:,2),data(:,3),'o');
        axis([-45 45 -7 7]);
        legend(description{2,1});
        xlabel(description{1,2});
        ylabel(description{2,2});
        subplot(2,1,2);
        errorbar(data(:,1),data(:,4),data(:,5),'^');
        axis([-45 45 -7 7]);
        legend(description{4,1});
        xlabel(description{1,2});
        ylabel(description{4,2});
    end;
end;

% --- fig7 ---
if flags.do_fig7
    description = { ...
         'all','x/m'; ...
         'all','y/m'; ...
         '1.43m','dx/m'; ...
         '1.43m','dy/m'; ...
         '0.41m','dx/m'; ...
         '0.41m','dy/m'; ...
         '0.20m','dx/m'; ...
         '0.20m','dy/m';
    };
    data = [ ...
      -0.00 -1.50  0.03  2.50  0.01  2.50  0.01  2.50
      -0.25 -1.50  0.38  2.50  0.38  2.50  0.24  2.50
      -0.50 -1.50  0.82  2.50  0.40  2.50  0.50  2.50
      -0.75 -1.50  1.21  2.50  0.67  2.50  0.67  2.50
      -1.00 -1.50  1.61  2.50  1.05  2.50  0.98  2.50
      -1.25 -1.50  2.02  2.50  1.52  2.50  1.20  2.50
      -1.50 -1.50  2.18  2.50  1.51  2.50  1.50  2.50
      -1.75 -1.50  1.76  2.50  1.85  2.50  1.73  2.50
      -0.00 -2.00  0.03  3.00  0.02  3.00 -0.01  3.00
      -0.25 -2.00  0.35  3.00  0.32  3.00  0.23  3.00      
      -0.50 -2.00  0.72  3.00  0.53  3.00  0.45  3.00      
      -0.75 -2.00  1.05  3.00  0.71  3.00  0.81  3.00      
      -1.00 -2.00  1.44  3.00  1.05  3.00  1.06  3.00      
      -1.25 -2.00  1.88  3.00  1.42  3.00  1.28  3.00      
      -1.50 -2.00  2.26  3.00  1.77  3.00  1.51  3.00      
      -1.75 -2.00  2.42  3.00  1.87  3.00  1.72  3.00      
    ];
    if flags.do_plot
        figure;
        subplot(1,3,1);
        quiver(data(:,1),data(:,2),data(:,3),data(:,4),20,'.b');
        title(description{3,1});
        axis([-2.13 1.63 -2.2 1.2]);
        xlabel('x/m');
        ylabel('y/m');
        subplot(1,3,2);
        quiver(data(:,1),data(:,2),data(:,5),data(:,6),20,'.b');
        title(description{5,1});
        axis([-2.13 1.63 -2.2 1.2]);
        xlabel('x/m');
        ylabel('y/m');
        subplot(1,3,3);
        quiver(data(:,1),data(:,2),data(:,7),data(:,8),20,'.b');
        title(description{7,1});
        axis([-2.13 1.63 -2.2 1.2]);
        xlabel('x/m');
        ylabel('y/m');
    end
end

% --- fig10 ---
if flags.do_fig10
    description = { ...
        'all','all','x/m'; ...
        'Y=1.5m','1.43m','phi_{auditory event}-phi_{sound event}/deg'; ...
        'Y=1.5m','1.43m','confidence interval'; ...
        'Y=1.5m','0.41m','phi_{auditory event}-phi_{sound event}/deg'; ...
        'Y=1.5m','0.41m','confidence interval'; ...
        'Y=1.5m','0.20m','phi_{auditory event}-phi_{sound event}/deg'; ...
        'Y=1.5m','0.20m','confidence interval'; ...
        'Y=2.0m','1.43m','phi_{auditory event}-phi_{sound event}/deg'; ...
        'Y=2.0m','1.43m','confidence interval'; ...
        'Y=2.0m','0.41m','phi_{auditory event}-phi_{sound event}/deg'; ...
        'Y=2.0m','0.41m','confidence interval'; ...
        'Y=2.0m','0.20m','phi_{auditory event}-phi_{sound event}/deg'; ...
        'Y=2.0m','0.20m','confidence interval'; ...
    };
    data = [
        0.00  -0.80 1.57  -0.29  2.49 -0.30  2.04  -0.48  0.75 -0.45  1.84  0.13  1.22
       -0.25  -3.03 1.52  -3.02  1.62  0.18  1.62  -1.94  0.80 -1.27  1.32  0.33  0.97
       -0.50  -6.94 1.43   2.28  2.01  0.02  1.67  -4.00  1.41 -0.56  1.24  0.93  1.01
       -0.75  -9.22 1.29   1.62  2.97  1.72  1.40  -5.32  2.04  0.74  1.95 -1.02  1.34
       -1.00 -10.99 2.14  -1.07  1.47  0.45  1.16  -7.18  2.09 -0.80  1.23 -1.08  1.53
       -1.25 -12.40 2.53  -4.78  1.34  0.95  1.34  -9.51  2.10 -2.73  1.54 -0.50  1.20
       -1.50 -10.07 3.24  -0.25  1.16 -0.01  1.19 -10.38  2.07 -3.94  1.29 -0.21  1.30
       -1.75  -0.11 2.56  -1.56  1.28  0.28  1.91  -8.61  2.05 -1.68  1.27  0.38  1.91
    ];
    if flags.do_plot
        figure
        subplot(3,1,1)
        h = errorbar(data(:,1)-0.025,data(:,2),data(:,3),'o');
        set(h,'MarkerEdgeColor','none','MarkerFaceColor','b')
        hold on;
        errorbar(data(:,1)+0.025,data(:,8),data(:,9),'o');
        axis([-1.85 0.125 -16 7]);
        legend(description{2,1},description{8,1});
        title(description{2,2});
        xlabel(description{1,3});
        ylabel(description{2,3});
        subplot(3,1,2)
        h = errorbar(data(:,1)-0.025,data(:,4),data(:,5),'o');
        set(h,'MarkerEdgeColor','none','MarkerFaceColor','b')
        hold on;
        errorbar(data(:,1)+0.025,data(:,10),data(:,11),'o');
        axis([-1.85 0.125 -16 7]);
        legend(description{4,1},description{10,1});
        title(description{4,2});
        xlabel(description{1,3});
        ylabel(description{2,3});
        subplot(3,1,3)
        h = errorbar(data(:,1)-0.025,data(:,6),data(:,7),'o');
        set(h,'MarkerEdgeColor','none','MarkerFaceColor','b')
        hold on;
        errorbar(data(:,1)+0.025,data(:,12),data(:,13),'o');
        axis([-1.85 0.125 -16 7]);
        legend(description{6,1},description{12,1});
        title(description{6,2});
        xlabel(description{1,3});
        ylabel(description{2,3});
    end
end


