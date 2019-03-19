% AMT - General functions used within the AMT
%
%  Signal levels and thresholds
%     DBSPL            - SPL of signal measured in dB.
%     SETDBSPL         - Specify SPL of signal.
%     ABSOLUTETHRESHOLD    - Absolute threshold of hearing
%
%  Emulation of experiments
%     EMUEXP         - Emulate psychoacoustic experiments
%
%  Plotting
%     AUDSPECGRAM      - Auditory spectrogram.
%     MODSPECGRAM      - Temporal modulation spectrogram
%     STMODSPECGRAM    - Spectro-temporal modulation spectrogram
%
%  Filters
%     GAMMATONE        - Calculate Gammatone filter coefficients
%     AUDITORYFILTERBANK - Linear auditory filterbank.
%     IHCENVELOPE        - Inner hair cell envelope extration.
%     ADAPTLOOP          - Adaptation loops.
%     MODFILTERBANK      - Modulation filter bank.
%     HEADPHONEFILTER    - FIR filter to model headphones 
%     MIDDLEEARFILTER    - FIR filter to model the middle ear
%     UFILTERBANKZ     - Apply multiple filters
%     FILTERBANKZ      - Apply multiple filters with non-equidistant downsampling
%
%  Sound localization
%     LOCALIZATIONERROR - Calculates various localization errors from localization responses
%     ITD2ANGLE        - Convert ITD to an angle using a lookup table
%     itd2angle_lookuptable - Create the lookup table
%
%  Speech
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/Contents.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%     SIIWEIGHTINGS        - Speech intelligibility weighted by frequency
