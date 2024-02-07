%DEMO_MOORE2016 calculates binaural loudness
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_moore2016.php


%   #Author: Clara Hollomey (2021): for testing

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

t = 1:32000;
s = sind(1000 * t).';
s = scaletodbspl(s, 100, 100);
%be careful when comparing to sone, scaletodbspl scales to rms
%so a 40 dBSPL sine set here will not lead exactly to 1 sone
Fs = 32000;
dBMax = 100;


filenameFilter = 'ff_32000.mat';

%calculate the outer/middleear filtering
s = moore2016_cochlea(s, filenameFilter);

%calculate the short term loudness
[InstantaneousSpecificLoudnessLeft, InstantaneousSpecificLoudnessRight] = moore2016_monauralinstspecloudness( s, Fs, dBMax );

%remove NAs (only necessary because the first value tends to be one, at least in Octave,
%and then, by successive summing, everything gets messed up)
InstantaneousSpecificLoudnessLeft(find(isnan(InstantaneousSpecificLoudnessLeft ))) = 0;
InstantaneousSpecificLoudnessRight(find(isnan(InstantaneousSpecificLoudnessRight ))) = 0;


ShortTermSpecificLoudnessLeft  = moore2016_shorttermspecloudness( InstantaneousSpecificLoudnessLeft );
ShortTermSpecificLoudnessRight = moore2016_shorttermspecloudness( InstantaneousSpecificLoudnessRight );

ShortTermLoudnessLeft  = zeros( size( ShortTermSpecificLoudnessLeft, 1 ), 1 );
ShortTermLoudnessRight = zeros( size( ShortTermSpecificLoudnessRight, 1 ), 1 );

%calculate the binaural loudness
for i = 1:size( ShortTermSpecificLoudnessLeft, 1 )
    [monauralShortTermLoudness(i), ShortTermLoudnessLeft(i), ShortTermLoudnessRight(i)] = moore2016_binauralloudness( ShortTermSpecificLoudnessLeft(i,:), ShortTermSpecificLoudnessRight(i,:) );
end

%calculate the long term loudness
LongTermLoudnessLeft   = moore2016_longtermloudness( ShortTermLoudnessLeft );
LongTermLoudnessRight  = moore2016_longtermloudness( ShortTermLoudnessRight );

LongTermLoudness = LongTermLoudnessLeft + LongTermLoudnessRight;

Loudness = max(LongTermLoudness);

ShortTermLoudness = ShortTermLoudnessLeft + ShortTermLoudnessRight;

% plot the results
out = plot_moore2016(ShortTermLoudness, LongTermLoudness);

