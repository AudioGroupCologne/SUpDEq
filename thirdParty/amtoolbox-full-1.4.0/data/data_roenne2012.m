function [ur,fs]  = data_roenne2012(varargin)
%DATA_ROENNE2012 Unitary response
%   Usage: ur = data_roenne2012()
%
%   [ur,fs]=DATA_ROENNE2012 returns the unitary response from Roenne
%   (2012) and its sampling frequency, fs=30000.
%
%   Examples:
%   ---------
%
%   The first plot shows the unitary response in the time-domain:
%
%     [ur,fs]  = data_roenne2012;
%     plot((0:length(ur)-1)/fs,ur);
%     xlabel('Time / seconds');
%     ylabel('Amplitude');
%
%   The second plot shows the magnitude response of the unitary
%   response, normalized so the highest peak reaches 0-dB:
%
%     [ur,fs]  = data_roenne2012;
%     magresp(ur,fs,90,'fir','1');
%
%   References:
%     F. M. Rønne, T. Dau, J. Harte, and C. Elberling. Modeling auditory
%     evoked brainstem responses to transient stimuli. The Journal of the
%     Acoustical Society of America, 131(5):3903--3913, 2012. [1]http ]
%     
%     References
%     
%     1. http://scitation.aip.org/content/asa/journal/jasa/131/5/10.1121/1.3699171
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_roenne2012.php


%   #Author: Peter L. Soendergaard (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% TODO: explain "unitary response"
  
s=amt_load('roenne2012','ur.mat');
ur=s.ur;

fs=30000;




