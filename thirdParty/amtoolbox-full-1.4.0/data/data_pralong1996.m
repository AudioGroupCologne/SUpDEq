function data = data_pralong1996(varargin)
%DATA_PRALONG1996 Head phone data from Pralong & Carlile (1996)
%   Usage: data = data_pralong1996(flag)
%
%   DATA_PRALONG1996(flag) returns data points from the Pralong & Carlile
%   (1996) paper.
%
%   The flag may be one of:
%
%     'no_plot'  Don't plot, only return data. This is the default.
%
%     'plot'    Plot the data.
%  
%     'fig1e'   Data from Fig. 1(e), Gain of Sennheiser 250 Linear
%               circumaural headphones. This is the default
%
%   Examples:
%   ---------
%
%   Figure 1e can be displayed using:
%
%     data_pralong1996('plot');
%
%   References:
%     D. Pralong and S. Carlile. The role of individualized headphone
%     calibration for the generation of high fidelity virtual auditory space.
%     J. Acoust. Soc. Am., 100:3785--3793, 1996.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_pralong1996.php


%   #Author: Peter L. Soendergaard (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
  
% Define input flags
definput.flags.type={'fig1e'};
definput.flags.plot = {'no_plot','plot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_fig1e
  data = [...
      125,	1; ...
      250,	1; ...
      500,	1; ...
      1000,	0.994850557; ...
      1237.384651,	0.994850557; ...
      1531.120775,	0.994850557; ...
      1894.585346,	1.114513162; ...
      2002.467159,	1.235743262; ...
      2344.330828,	1.867671314; ...
      2721.273584,	2.822751493; ...
      3001.403462,	2.180544843; ...
      3589.453635,	1.442755787; ...
      4001.342781,	1.173563859; ...
      4441.534834,	1.37016005; ...
      5004.212211,	1.599690164; ...
      5495.887031,	1.37016005; ...
      5997.423738,	1.114513162; ...
      6800.526258,	0.648125625; ...
      6946.931144,	0.631609176; ...
      7995.508928,	0.276505667; ...
      8414.866811,	0.084335217; ...
      9008.422743,	0.084335217; ...
                 ];
  
  
  if flags.do_plot
    figure;
    x=data(:,1);
    freqresp=20*log10(data(:,2));
    semilogx(x,freqresp);
    xlabel('Frequency (Hz)');
    ylabel('Gain (dB)');
  end;
  
end;


