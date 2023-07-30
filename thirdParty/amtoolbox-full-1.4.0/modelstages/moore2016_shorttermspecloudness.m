function [ShortTermSpecificLoudness, ShortTermLoudness] = moore2016_shorttermspecloudness( InstantaneousSpecificLoudness )
%MOORE2016_SHORTTERMSPECLOUDNESS calculates the short-term specific loudness
%   
%
%   Input parameters:
%     InstantaneousSpecificLoudness : monaural instantaneous specific loudness
%
%   Output parameters:
%     ShortTermSpecificLoudness : short term specific loudness
%     ShortTermLoudness         : short term loudness
%
%   This code calculates the short term specific loudness for the binaural 
%   loudness model moore2016 in the version for TVL 2016 based on ANSI S3.4-2007 and Moore & Glasberg (2007).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/moore2016_shorttermspecloudness.php


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


aA = 0.045;
aR = 0.033;

ShortTermSpecificLoudness = zeros( size(InstantaneousSpecificLoudness) );

ShortTermSpecificLoudness(1,:) = moore2016_agcnextframeofvector( zeros(1,150) , InstantaneousSpecificLoudness(1,:), aA, aR );
for i=2:size( InstantaneousSpecificLoudness, 1 )
    ShortTermSpecificLoudness(i,:) = moore2016_agcnextframeofvector( ShortTermSpecificLoudness(i-1,:), InstantaneousSpecificLoudness(i,:), aA, aR );
end

ShortTermLoudness = sum( ShortTermSpecificLoudness, 2 ) / 4;

end

function out = moore2016_agcnextframeofvector( vLastFrame, vThisInput, aA, aR )

outThisIsLarger = aA * vThisInput + ( 1 - aA ) * vLastFrame;   % attack
outLastIsLarger = aR * vThisInput + ( 1 - aR ) * vLastFrame;   % release

out = ( vThisInput > vLastFrame ) .* outThisIsLarger + ( vThisInput <= vLastFrame ) .* outLastIsLarger;

end


