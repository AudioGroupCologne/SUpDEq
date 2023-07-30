function [hartestim,fs]  = data_harte2009()
%DATA_HARTE2009 Tone burst stimuli from Harte et al. (2009)
%
%   Usage: [hartestim,fs] = data_harte2009;
%
%   DATA_HARTE2009 returns the tone burst stimuli from
%   Harte et al. (2009) and the sampling frequency, fs=48000.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_harte2009.php


%   #Author: The AMT Team (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%   References:
%   J. Harte, G. Pigasse, and T. Dau. Comparison of cochlear delay
%   estimates using otoacoustic emissions and auditory brainstem responses.
%   J. Acoust. Soc. Am., 126(3):1291--1301, 2009.



% TODO: explain stimuli in description;

hartestim = amt_load('harte2009','stim.mat');
  
fs = 48e3;
  


