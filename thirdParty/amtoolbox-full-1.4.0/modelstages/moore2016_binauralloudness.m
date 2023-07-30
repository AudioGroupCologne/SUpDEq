function [Loudness, LoudnessLeft, LoudnessRight] = moore2016_binauralloudness( SpecificLoudnessLeftMon, SpecificLoudnessRightMon )
%MOORE2016_BINAURALLOUDNESS Calculate the binaural loudness out of monaural loudness at 0.25-ERB steps
%
%   Input parameters:
%     SpecificLoudnessLeftMon : monaural specific loudness left ear
%     SpecificLoudnessRightMon : monaural specific loudness right ear
%
%   Output parameters:
%     Loudness      : binaural loudness
%     LoudnessLeft  : loudness at left ear
%     LoudnessRight : loudness at right ear
%
%   This code weights the respective monaural loudness taking account binaural
%   inhibition for calculating the binaural loudness in moore2016
%   in the version for TVL 2016 based on ANSI S3.4-2007 and Moore & Glasberg (2007).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/moore2016_binauralloudness.php


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


g = -18:0.25:18; 
W = exp( -(0.08 .* g) .^2 );

SpecificLoudnessLeftSmoothed  = conv(W,SpecificLoudnessLeftMon)/sum(W);
SpecificLoudnessLeftSmoothed  = SpecificLoudnessLeftSmoothed(73:222);   
SpecificLoudnessRightSmoothed = conv(W,SpecificLoudnessRightMon)/sum(W); 
SpecificLoudnessRightSmoothed = SpecificLoudnessRightSmoothed(73:222);

SpecificLoudnessLeftSmoothed  = SpecificLoudnessLeftSmoothed + 10^-13;
SpecificLoudnessRightSmoothed = SpecificLoudnessRightSmoothed + 10^-13;

p = 1.5978;

InhibLeft = 2 ./ ( 1 + ( sech( SpecificLoudnessRightSmoothed ./ SpecificLoudnessLeftSmoothed ) ) .^ p ); 
InhibRight = 2 ./ ( 1 + ( sech( SpecificLoudnessLeftSmoothed ./ SpecificLoudnessRightSmoothed ) ) .^ p );

SpecificLoudnessLeft = SpecificLoudnessLeftMon ./ InhibLeft; 
SpecificLoudnessRight = SpecificLoudnessRightMon ./ InhibRight;            

LoudnessLeft = sum(SpecificLoudnessLeft) / 4;   
LoudnessRight = sum(SpecificLoudnessRight) / 4;

Loudness = LoudnessLeft + LoudnessRight;  


