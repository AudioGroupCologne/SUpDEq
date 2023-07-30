function [InstantaneousSpecificLoudnessLeft, InstantaneousSpecificLoudnessRight] = moore2016_monauralinstspecloudness(signal, Fs, dBMax)
%MOORE2016_MONAURALINSTSPECLOUDNESS calculates instantaneous specific loudness over time
%
%   Input parameters:
%     signal : input signal
%     fs : sampling frequency
%     dbMax : the RMS SPL of a sinusoid with a peak amplitude of 1.
%
%   Output parameters:
%     InstantaneousSpecificLoudnessLeft  : instantaneous specific loudness left ear
%     InstantaneousSpecificLoudnessRight : instantaneous specific loudness right ear
%
%
%   This code calculates the monaural instantaneous specific loudness for the binaural 
%   loudness model moore2016 in the version for TVL 2016 based on ANSI S3.4-2007 and Moore & Glasberg (2007).
%
%   The left and the right input signal and sampling rate Fs must be passed, the signal must
%   contain two columns. Fs must be 32 kHz, the reading of the wav file and resampling must 
%   be done before calling this function. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/moore2016_monauralinstspecloudness.php


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

npts               = Fs / 1000 * 64;                     % points for FFT, 2048
nSegmentDuration   = 1;                                  % duration of segment in ms
nSamplesPerSegment = Fs / 1000 * nSegmentDuration;       % 32
nSegmentsInSignal  = floor( (length(signal) - npts) / Fs * 1000 / nSegmentDuration );

% Hann windows for 6 FFTs; 1st column 64 ms, 6th column 2 ms
wHann = zeros(npts,6);
for i = 1:6
    wHann(:,i) = [ zeros( (1 - 1/2^(i-1) ) / 2 * npts, 1 ); hann( npts / 2^(i-1) ); zeros( (1 - 1/2^(i-1) ) / 2 * npts, 1 ) ];
end

% indices which shall be used from the ffts, 20-80 Hz from the first...
vLimitingF = [20 80 500 1250 2540 4050 15000]; 
vLimitingIndices = zeros(1,7);
for i=1:7
    vLimitingIndices(i) = floor( vLimitingF(i) / (Fs/npts) ) + 1;
end

InstantaneousSpecificLoudnessLeft  = zeros(nSegmentsInSignal + 1, length(1.75:0.25:39) );
InstantaneousSpecificLoudnessRight = zeros(nSegmentsInSignal + 1, length(1.75:0.25:39) );

%% loop: window segment, apply stationary loudness model, move window
%% by 1 ms

% tic
for iSegment = 0:nSegmentsInSignal
%     
%     if ( mod(iSegment,50) == 0 )
%         disp([num2str(iSegment) ' ms of ' num2str(nSegmentsInSignal) ' ms analyzed'])
%     end
    
    [fLeftRelevant, LLeftRelevant, fRightRelevant, LRightRelevant] = moore2016_spectrum(signal( (iSegment*nSamplesPerSegment+1) : (iSegment*nSamplesPerSegment+npts), : ), Fs, dBMax, wHann, vLimitingIndices);
    
    if ( isempty( LLeftRelevant ) )
        SpecificLoudnessLeft = zeros(1,150);
    else
        ExcitationLevelsLeft = moore2016_excitationpattern(fLeftRelevant', LLeftRelevant');
        SpecificLoudnessLeft = moore2016_specloudnessbinaural( ExcitationLevelsLeft );
    end
    if ( isempty( LRightRelevant ) )
        SpecificLoudnessRight = zeros(1,150);
    else
        ExcitationLevelsRight = moore2016_excitationpattern(fRightRelevant', LRightRelevant');
        SpecificLoudnessRight = moore2016_specloudnessbinaural( ExcitationLevelsRight );
    end
%     [ ThisSpecificLoudness, ThisLoudness ] = MonauralSpecificLoudness2BinauralSpecificLoudness025( SpecificLoudnessLeft, SpecificLoudnessRight );
    InstantaneousSpecificLoudnessLeft(iSegment+1,:) = SpecificLoudnessLeft;
    InstantaneousSpecificLoudnessRight(iSegment+1,:) = SpecificLoudnessRight;
%     InstantaneousLoudness(iSegment+1)         = ThisLoudness;
end
% toc
end 

function out = moore2016_specloudnessbinaural( ExcitationLevels )

% version for TVL 2016 based on ANSI S3.4-2007 and Moore & Glasberg (2007)
% Calculate the specific loudness out of excitation patterns at 0.25-ERB
% steps, using the constant C of the model that accounts for binaural
% loudness

% ERB-scale
ERBc = 1.75:0.25:39;
fc = erb2fc(ERBc);

% Parameters
C = 0.06267;
C = 0.063026;
C = 0.0631;
G = get_G(fc);
Alpha = get_Alpha(fc);
A = get_A(fc);

% Excitation and threshold in linear power units
Excitation = 10 .^ ( ExcitationLevels ./ 10 );
Threshold = excitationthreshold(fc);
Threshold = 10 .^ ( Threshold ./ 10 );

% 3 ranges - calculate all for each case (computationally cheap)
N1 = C .* ( ( G .* Excitation + A ) .^ Alpha - A .^Alpha );
N2 = C .* ( 2 .* Excitation ./ ( Excitation + Threshold ) ) .^ 1.5 .* ( ( G .* Excitation + A ) .^ Alpha - A .^Alpha );
N3 = C .* ( Excitation ./ 1.0707 ) .^ 0.2;

out = ( Excitation > Threshold & Excitation < 10^10 ) .* N1 + ( Excitation <= Threshold ) .* N2 + (Excitation >= 10^10 ) .* N3;

% figure;
% hold on;
% semilogy(ExcitationLevels, out);
% xlim([0 110]);
% ylim([0.005 10]);
end

%end

function out = get_Alpha(f)

TableG2Alpha = [ -25.0   -20.0   -15.0   -10.0   -5.0    0.0;
                 0.26692 0.25016 0.23679 0.22228 0.21055 0.20000];

G = get_G(f);
G = 10 .* log(G) ./ log(10);

out = interpolation(TableG2Alpha(1,:), TableG2Alpha(2,:), G,'linear');

end

function out = get_A(f)

TableG2A = [ -24.54531 -23.78397 -22.78169 -21.76854 -20.74442 -19.78305 -18.90431 -18.01605 -17.11816 -16.21055 -15.32375 -14.59341 ...
             -13.91727 -13.29726 -12.73537 -12.23364 -11.75255 -11.23866 -10.75136 -10.29164  -9.86051  -9.45902  -9.08823  -8.72191 ...
              -8.35715  -8.01199  -7.68715  -7.38338  -7.10145  -6.84213  -6.60623  -6.39458  -6.14589  -5.89392  -5.65071  -5.41661 ...
              -5.19198  -4.97718  -4.77258  -4.57857  -4.39555  -4.20148  -4.00538  -3.81442  -3.62882  -3.27454  -3.10633  -2.94438 ...
              -2.78894  -2.64027  -2.50042  -2.37015  -2.24820  -2.13487  -2.03046  -1.93531  -1.84973  -1.77405  -1.70863  -1.65382 ...
              -1.60997  -1.57745  -1.51786  -1.44522  -1.37466  -1.30624  -1.24006  -1.17621  -1.11478  -1.05587  -0.99956  -0.94596 ...
              -0.89518  -0.84731  -0.80246  -0.75663  -0.69834  -0.64029  -0.58251  -0.52501  -0.46783  -0.41099  -0.35451  -0.29842 ...
              -0.24274  -0.18750  -0.13273  -0.07845  -0.02470   0.00000;
               8.85200   8.63150   8.35840   8.10120   7.85850   7.65258   7.49124   7.32853   7.16458   6.99954   6.83957   6.70559 ...
               6.58115   6.46869   6.36813   6.27970   6.19583   6.10719   6.02404   5.94640   5.87876   5.82490   5.77551   5.72719 ...
               5.67939   5.63443   5.59236   5.55322   5.51708   5.48399   5.45402   5.42723   5.39588   5.36425   5.33386   5.30472 ...
               5.27688   5.25064   5.22800   5.20667   5.18660   5.16539   5.14401   5.12326   5.10314   5.06490   5.04681   5.02944 ...
               5.01280   4.99693   4.98203   4.96818   4.95524   4.94324   4.93219   4.92215   4.91312   4.90515   4.89827   4.89251 ...
               4.88790   4.88449   4.87824   4.87063   4.86324   4.85609   4.84918   4.84252   4.83611   4.82998   4.82412   4.81854 ...
               4.81327   4.80830   4.80365   4.79890   4.79286   4.78685   4.78088   4.77494   4.76904   4.76318   4.75736   4.75159 ...
               4.74586   4.74019   4.73456   4.72899   4.72349   4.72096 ];

G = get_G(f);
G = 10 .* log(G) ./ log(10);

out = interpolation(TableG2A(1,:), TableG2A(2,:), G);

end

function out = get_G(f)
%addpath(fullfile (pwd, 'functions'));
% low level gain of the cochlear amplifier relative to the gain at 500 Hz

LinearThreshold = 10.^( excitationthreshold(f) ./ 10 );
out = 10^(3.63/10) ./ LinearThreshold;
end

function out = excitationthreshold( f )
% version for TVL 2016 based on ANSI S3.4-2007 and Moore & Glasberg (2007)         
Above500 = 3.63;
Threshold = [ 50    63    80    100   125   160   200  250  315  400  500  630  750  800  1000;
              28.18 23.90 19.20 15.68 12.67 10.09 8.08 6.30 5.30 4.50 3.63 3.63 3.63 3.63 3.63 ];

Below500 = interpolation(Threshold(1,:), Threshold(2,:), f, 'linear');
          
out = ( f < 500 ) .* Below500 + ( f >= 500 ) .* Above500;
end


