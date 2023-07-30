function out = moore2016_agcnextframe( dLastFrame, dThisInput, aA, aR )
%MOORE2016_AGCNEXTFRAME adjusts successive short term loudness frames
%
%   Input parameters:
%     dLastFrame : last idx
%     dThisInput : current idx
%     aA : attack coefficient
%     aR : release coefficient
%
%   Output parameters:
%     out : adjusted loudness frame
%
%   This code corresponds to the AGC stage in the binaural loudness model moore2016
%   in the version for TVL 2016 based on ANSI S3.4-2007 and Moore & Glasberg (2007).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/moore2016_agcnextframe.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author: Josef Schlittenlacher (2018): original code
%   #Author: Clara Hollomey (2021): integration in the AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if ( dThisInput > dLastFrame )                            % attack
    out = aA * dThisInput + ( 1 - aA ) * dLastFrame;
else                                                      %release
    out = aR * dThisInput + ( 1 - aR ) * dLastFrame;
end


