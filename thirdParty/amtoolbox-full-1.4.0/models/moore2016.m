function [shorttermloudness, longtermloudness, loudness] = moore2016(inputsignal)
%MOORE2016 Binaural loudness model
%   Usage: [shorttermloudness, longtermloudness, loudness] = moore2016(inputsignal)
%
%   Input parameters:
%     inputsignal : a vector or 2 dimensional matrix containing the input signal
%                   sampled to 32 kHz
%
%   Output parameters:
%     shorttermloudness : shortterm loudness [phon]
%     longtermloudness  : longterm loudness [phon]
%     loudness          : maximum loudness [phon]
%
%   For each ear, the model includes: an outer and middle ear filter; short-term 
%   spectral analysis; calculation of an excitation pattern, a compressive nonlinearity, 
%   and smoothing over time.
%   The short-term loudness is calculated as the sum of the short-term loudness 
%   values for the two ears. The long-term loudness for each ear is obtained from
%   the short-term loudness. The overall loudness impression is calculated as the 
%   sum of the long-term loudness of both ears. 
%   The Matlab code provided calculates loudness according to the model described by Moore et 
%   al. (2016), but with the modified time constants described by Moore et al. (2018). It was 
%   developed from C code for the same model, and Matlab code written for ANSI S3.4-2007, 
%   based on Moore et al. (1997) and Glasberg and Moore (2006) and ISO 532-2 (2017), based 
%   on Moore and Glasberg (2007).
%   The code may be used with wav files (one or two channels). If a one-channel file is used, the 
%   program assumes diotic presentation. To calculate the loudness of a monaural signal, a 
%   second channel filled with zeros must be added. 
%
%   See also: f2erb plot_moore2016 moore2016_cochlea moore2016_monauralinstspecloudness
%             moore2016_agcnextframe moore2016_longtermloudness moore2016_binauralloudness
%             moore2016_spectrum moore2016_excitationpattern 
%             moore2016_shorttermspecloudness
%
%   References:
%     B. R. Glasberg and B. C. J. Moore. A Model of Loudness Applicable to
%     Time-Varying Sounds. J. Audio Eng. Soc, 50(5):331--342, 2002.
%     
%     B. C. J. Moore, B. R. Glasberg, and T. Baer. A Model for the Prediction
%     of Thresholds, Loudness, and Partial Loudness. J. Audio Eng. Soc,
%     45(4):224--240, 1997.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/moore2016.php


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


%defaults
Fs = 32000;
dBMax = 100;

filenameFilter = 'ff_32000.mat';

%calculate the outer/middleear filtering
s = moore2016_cochlea(inputsignal, filenameFilter);

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

longtermloudness = LongTermLoudnessLeft + LongTermLoudnessRight;

loudness = max(longtermloudness);

shorttermloudness = ShortTermLoudnessLeft + ShortTermLoudnessRight;

