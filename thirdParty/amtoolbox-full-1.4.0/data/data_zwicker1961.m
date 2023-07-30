function data = data_zwicker1961(varargin)
%DATA_ZWICKER1961  Data for the Bark scale
%   Usage: data = data_zwicker1961;
%
%   DATA_ZWICKER1961 returns the data points defining the notion of
%   critical bands. The output data consists of the frequency limit of
%   each band in Hz.
%
%   To get the bandwidth of each channel, simply use:
%
%     bw = diff(data_zwicker1961);
%
%   The first entry has been modified from the original paper: It was
%   originally 20 Hz, but in Zwicker and Fastl 1999 this was changed to 0
%   Hz.
%
%   The data can be plotted using :
%
%     data_zwicker1961('plot');
%  
%   References:
%     E. Zwicker. Subdivision of the audible frequency range into critical
%     bands (frequenzgruppen). J. Acoust. Soc. Am., 33(2):248--248, 1961.
%     
%     E. Zwicker and H. Fastl. Psychoacoustics: Facts and models, volume 254.
%     Springer Berlin, 1999.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_zwicker1961.php


%   #Author: Peter L. Soendergaard (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% TODO: describe Data in description
  
definput.flags.plot = {'no_plot','plot'};

[flags,kv]  = ltfatarghelper({},definput,varargin);

data = [0,100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,...
        2320,2700,3150,3700,4400,5300,6400,7700,9500,12000,15500].';

if flags.do_plot
  dl=length(data);
  semilogx(data(2:end),1:dl-1);
  xlim([100, 16000]);
  xlabel('Frequency in Hz.');
  ylabel('Critical-band function (Tonheit) in Bark.');
end;

