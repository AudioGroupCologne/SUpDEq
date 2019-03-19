function [delay,data_mean,data_std]  = data_elberling2010(varargin)
%DATA_ELBERLING2010 ABR wave V data as functon of level and sweeping rate
%   Usage: data = data_elberling2010(flag)
%
%   Output parameters:
%     delay      : "x-axis" - sweeping rate delay between 710Hz and 5700Hz
%     data_mean  : Mean of data 
%     data_std   : Standard deviation of data
%
%   DATA_ELBERLING2010(flag) returns data points from the Elberling et
%   al. (2010)
%
%   The flag may be one of:
%
%     'noplot'  Don't plot, only return data. This is the default.
%
%     'plot'    Plot the data.
%  
%     'fig4'    Data from Fig. 4, Amplitude of wave V.
%
%     'fig5'    Data from Fig. 5, Latency of wave V
%
%     'stim'    Return the stimulus and the sampling frequency. XXX Describe
%               the data. XXX change 'stim' to figure no., if applicable.
%
%   Examples:
%   ---------
%
%   Figure 4 can be displayed using:
%
%     data_elberling2010('fig4','plot');
%
%   Figure 5 can be displayed using:
%
%     data_elberling2010('fig5','plot');
%
%   References:
%     C. Elberling, J. Calloe, and M. Don. Evaluating auditory brainstem
%     responses to different chirp stimuli at three levels of stimulation. J.
%     Acoust. Soc. Am., 128(1):215-223, 2010.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_elberling2010.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% TODO: explain Data in description;
    
% Define input flags
definput.flags.type={'missingflag','fig4','fig5','stim'};
definput.flags.plot = {'noplot','plot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
    flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
        sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
    error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end

% Font size
ftsz=12;
%color
col=[0.7,0.7,0.7];

fs=30000;

delay = [0 1.86 2.56 3.32 4.14 5.04];
delay2 = [delay;delay;delay];

if flags.do_fig4
  data_mean = [256 311 339 349 382 391; ...
               368 560 613 645 631 581; ...
               407 624 654 596 516 389]';

  data_std  = [61   86  86  72  83  91;...
               92  135 119  79  79 100;...
               91  162 118 118 96  120]';
  
  if flags.do_plot
    set(gca,'fontsize',ftsz);
    
    figure;
    hold on;
    errorbar(delay2',data_mean,data_std/sqrt(20),'-', 'linewidth',1.5, ...
             'color',col);
    %,ylim([0 800])
    axis([-1.2 5.5 0 800]);
    xlabel('Change of delay [ms]');
    ylabel('ABR amplitude [nv]')
    text(-.7,data_mean(1,1), '20','fontsize',ftsz)
    text(-.7,data_mean(1,2), '40','fontsize',ftsz)
    text(-.7,data_mean(1,3), '60','fontsize',ftsz)
    text(-.9,data_mean(1,3)+70, 'dB nHL','fontsize',ftsz)
    text(-.2,50, 'Click','fontsize',ftsz)
    text(delay(2),75, '1','fontsize',ftsz)
    text(delay(3),75, '2','fontsize',ftsz)
    text(delay(4),75, '3','fontsize',ftsz)
    text(delay(5),75, '4','fontsize',ftsz)
    text(delay(6),75, '5','fontsize',ftsz)
    text(3,40, 'Chirps','fontsize',ftsz)
  end;
  
end;

if flags.do_fig5
  
  data_mean    = [...
      7.60 7.23 7.12 6.97 6.86 6.66; ...
      6.59 6.15 5.96 5.76 5.48 4.98; ...
      5.88 4.97 4.57 4.12 3.50 2.87]';
  data_std     = [...
      0.79 0.36 0.32 0.34 0.36 0.36; ...
      0.74 0.29 0.23 0.27 0.34 0.59; ...
      0.25 0.28  0.35 0.51 0.51 0.38]';
  
  if flags.do_plot
    figure; 
    
    set(gca,'fontsize',ftsz);
    axis([-1.2 6.5 0 10]);
    xlabel('Change of delay [ms]');
    ylabel('ABR latency [ms]')
    text(-.7,5.88, '60','fontsize',ftsz,'color',col)
    text(-.7,6.59, '40','fontsize',ftsz,'color',col)
    text(-.7,7.6, '20','fontsize',ftsz,'color', col)
    text(-.9,8.7, 'dB nHL','fontsize',ftsz,'color',col)
    text(5.5,data_mean(6,1)/fs*1000-15, '20','fontsize',ftsz)
    text(5.5,data_mean(6,2)/fs*1000-15, '40','fontsize',ftsz)
    text(5.5,data_mean(6,3)/fs*1000-15, '60','fontsize',ftsz)
    text(5.3,data_mean(6,1)/fs*1000-15+.6, 'dB nHL','fontsize',ftsz)
    text(-.2,.5, 'Click','fontsize',ftsz)
    text(delay(2),1, '1','fontsize',ftsz)
    text(delay(3),1, '2','fontsize',ftsz)
    text(delay(4),1, '3','fontsize',ftsz)
    text(delay(5),1, '4','fontsize',ftsz)
    text(delay(6),1, '5','fontsize',ftsz)
    text(3,.5, 'Chirps','fontsize',ftsz)
    box on;
    hold on;
    set(gca,'fontsize',ftsz);
    errorbar(delay2',data_mean,data_std,'-',...
             'linewidth',1.5,'color',col)
    %,ylim([0 800])
  end;
end;

if flags.do_stim
  % delay is the first output parameter, use it to return the stimulus.  
  delay=amt_load('elberling2010','stim.mat');
  
  % fs is stored in data_mean
  data_mean=fs;
  
  if flags.do_plot
    amt_disp('XXX Plot is missing.');
  end;
end;
