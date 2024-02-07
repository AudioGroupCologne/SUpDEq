function data = data_boyd2012
%DATA_BOYD2012 data from Boyd et al. (JASA-EL, 2012)
%
%   Usage: data = data_boyd2012
%
%   Output parameters:
%     data    : structure
%
%   The 'data' struct has the following fields:
%
%     'ID'               subject ID
%
%     'Resp'             externalization responses 
%
%     'BRIR'             binaural room impulse responses for 4 positions
%
%     'Target'           target stimuli (single and four talker conditions)
%
%     'Reference_1T'     reference stimuli (single talker)
%
%     'Reference_4T'     reference stimuli (four talkers)
%
%     'fs'               sampling rate in Hz
%
%   Mean externalization scores of NH listeners extracted from top panels 
%   (1 talker condition) of Fig. 1 
%
%   References:
%     A. W. Boyd, W. M. Whitmer, J. J. Soraghan, and M. A. Akeroyd. Auditory
%     externalization in hearing-impaired listeners: The effect of pinna cues
%     and number of talkers. J. Acoust. Soc. Am., 131(3):EL268--EL274, 2012.
%     [1]www: ]
%     
%     References
%     
%     1. http://dx.doi.org/10.1121/1.3687015
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_boyd2012.php


%   #Author: Robert Baumgartner

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

tmp = amt_load('boyd2012','boyd2012.mat');
data = tmp.data;

end


