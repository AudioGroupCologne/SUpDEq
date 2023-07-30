function [processed, mpar] = eurich2022(mRef, mTest,mpar)
%EURICH2022  Binaural detection model based on interaural coherence
%   Usage: [processed, mpar] = eurich2022(mRef, mTest,mpar)
%
%   Input parameters:
%     mRef     : Binaural reference signal
%     mTest    : Binaural test signal
%     mpar     : Model parameters
%
%   Output parameters:
%     processed  : Processed model outputs (reference and test, concatenated)
%     mpar       : Modified model parameters
%
%   EURICH2022 is a binaural detection model based 
%   on interaural coherence
%
%
%   See also: exp_eurich2022 eurich2022_processing
%
%   References:
%     B. Eurich, J. Encke, S. D. Ewert, and M. Dietz. Lower interaural
%     coherence in off-signal bands impairs binaural detection. The Journal
%     of the Acoustical Society of America, 151(6):3927--3936, 06 2022.
%     [1]arXiv | [2]http ]
%     
%     References
%     
%     1. http://arxiv.org/abs/https://pubs.aip.org/asa/jasa/article-pdf/151/6/3927/16528275/3927\_1\_online.pdf
%     2. https://doi.org/10.1121/10.0011673
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/eurich2022.php


%   #StatusDoc: Submitted
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Author: Bernhard Eurich (2022)
%   #Author: Piotr Majdak (2023): integration in the AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 



% process Reference and Test signals
[processed_mRef,mpar]  = eurich2022_processing(mRef,mpar);
[processed_mTest,mpar] = eurich2022_processing(mTest,mpar);

processed = cat(1,processed_mRef, processed_mTest);



