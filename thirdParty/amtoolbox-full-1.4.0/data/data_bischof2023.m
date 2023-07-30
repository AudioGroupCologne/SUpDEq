function [bischofdata,fs]  = data_bischof2023()
%DATA_BISCHOF2023 data from Bischof et al. (2023)
%   Usage: [bischofdata,fs]  = data_bischof2023;
%
%   DATA_BISCHOF2023
%   returns the reverberant harmonic complex tone target and noise interferer
%   signals stimuli from Bischof et al. (2023) and the used sampling 
%   frequency, fs=44100.
%
%   References:
%     N. Bischof, P. Aublin, and B. Seeber. Fast processing models effects of
%     reflections on binaural unmasking. Acta Acustica, 2023.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_bischof2023.php


%   #Author: Norbert F. Bischof (2023)
%   #Author: Pierre G. Aublin
%   #Author: Bernhard Seeber (2023)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.

bischofdata = amt_load('bischof2023','expdata_bischof2023.mat');
fs = 44100;

