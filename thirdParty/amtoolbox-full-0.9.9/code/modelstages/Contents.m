% AMT - Various stages of auditory models
%
%  HRTF processing
%     ziegelwanger2014_onaxis      - On-axis version of Ziegelwanger et al. (2014)
%     ziegelwanger2014_offaxis     - Off-axis version of Ziegelwanger et al. (2014)
%     ziegelwanger2013_onaxis      - On-axis version of Ziegelwanger et al. (2013) (legacy only)
%     ziegelwanger2013_offaxis     - Off-axis version of Ziegelwanger et al. (2013) (legacy only)
%
%  Peripheral processing stages
%     ZILANY2014_FFGN   - Fast (exact) fractional Gaussian noise and Brownian motion generator.
%
%  Binaural processing stages
%     langendijk2002_comp     - Comparison process from Langendijk 20002.
%     lindemann1986_bincorr - Running cross-correlation between two signals.
%     culling2005_bmld    - BMLD calculation from Culling et al. (2005).
%     dietz2011_unwrapitd - IPD to ITD transformation for the Dietz model
%     dietz2011_filterbank - filterbank of Dietz 2011 binaural model
%     dietz2011_interauralfunctions - Calculate interaural parameters for Dietz 2011 model
%     wierstorf2013_estimateazimuth - Azimuth position estimation based on dietz2011 or lindemann1986
%     ziegelwanger2014_offaxis - synthesis of TOAs for the off-axis model
%     ziegelwanger2014_onaxis  - synthesis of TOAs for the on-axis model
%     langendijk2002_likelihood - Likelihood estimation
%
%  
%   Modulation detection and masking stages
%     dau1996_preproc - signal detection model from Dau et al. (1996)
%     dau1997_preproc - modulation-detection and modulation-masking processing from Dau et al. (1997)
%     jepsen2008_preproc - intensity discrimination, tone-in-noise detection, spectral masking, forward masking, and amplitude-modulation detection from Jepsen et al. (2008)
%
%     lindemann1986_centroid - Centroid of the cross-correlation activation
%
%     roenne2012_chirp          - Simulate chirp evoked ABRs
%     roenne2012_click          - Simulate ABR respone to click
%     roenne2012_tonebursts     - Simulate tone burst evoked ABR wave V latencies
%
%   Gammatone filterbank framework from Hohmann (2002)
%     hohmann2002_delay            - Create a delay object
%     hohmann2002_mixer            - Create a mixer object
%     hohmann2002_synth            - Create a synthesis object 
%     hohmann2002_filter           - Create a single Gammatone filter object
%     hohmann2002_process          - Process the input signals by the corresponding filterbank object
%     hohmann2002_clearstate       - Clears the state of the filterbank object
%     hohmann2002_freqz            - Calculates frequency response of a filterbank object
%     hohmann2002_clearstate       - Clears the state of the filterbank objects
%     hohmann2002_freqz            - Calculates frequency response of a filterbank object
%
%  The Takanen 2013 model   - Takanen 2013 binaural auditory model
%     takanen2013_contracomparison - Enhance contrast between hemispheres
%     takanen2013_cueconsistency - Check consistency before cue combination
%     takanen2013_directionmapping - Map the directional cues to directions
%     takanen2013_formbinauralactivitymap - Steer cues on a topographic map
%     takanen2013_lso        - Model of the lateral superior olive
%     takanen2013_mso        - Model of the medial superior olive
%     takanen2013_onsetenhancement - Emphasize onsets on direction analysis
%     takanen2013_periphery  - Process input through the model of periphery
%     takanen2013_wbmso      - Wideband medial superior olive model
%     takanen2013_weightedaveragefilter - Part of the takanen2013 model
%
%  Binaural masking level differences from Breebaart et al. (2001)
%     breebaart2001_preproc  - Preprocessor of a binaural signal with EI-cell output
%     breebaart2001_eicell - Excitation-inhibition cell from Breebaart et al. (2001).
%     breebaart2001_centralproc - Central processor taking decision in an experiment
%     breebaart2001_outmiddlefilter - Outer and middle-ear filter used by the model Breebaart et al. (2001)
% 
%  Sagittal-plane localization models
%     langendijk2002_spectralanalysis - FFT-based filter bank with constant relative bandwidth
%     baumgartner2013_pmv2ppp - Calculate performance predictions from PMVs for baumgartner2013
%     baumgartner2014_spectralanalysis - Approximation of spectral analysis by auditory periphery
%     baumgartner2014_gradientextraction - Extraction of positive spectral gradients
%     baumgartner2014_comparisonprocess - Comparison with direction-specific templates
%     baumgartner2014_similarityestimation - Similarity estimation with listener-specific sensitivity
%     baumgartner2014_binauralweighting - Binaural combination of monaural similarity estimates
%     baumgartner2014_sensorimotormapping - Response scatter induced by localization task
%     baumgartner2014_calibration - Calibration of the model (linear periphery)
%     baumgartner2014_virtualexp - Performs a virtual sound-localization experiment
%     baumgartner2014_pmv2ppp - Performance predictions from PMVs of baumgartner2014
%     baumgartner2014_likelistat - Likelihood statistics for evaluation of model performance
%     baumgartner2016_calibration - Calibration of the model (nonlinear periphery)
%     baumgartner2016_spectralanalysis - Approximation of spectral analysis by auditory periphery
%     baumgartner2016_gradientextraction - Extraction of positive spectral gradients
%     baumgartner2016_comparisonprocess - Comparison with direction-specific templates
%
%  Cochlear-implant localization model
%     kelvasa2015_anbinning - AN and time binning from Kelvasa and Dietz 2015 binaural model
%     kelvasa2015_anprocessing -  AN model used in Kelvasa and Dietz 2015 binaural model
%     kelvasa2015_calibratemapping - Produces necessary mappings for localization model
%     kelvasa2015_ciprocessing - CI ACE processing strategy used in Kelvasa and Dietz 2015 binaural model
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/Contents.php

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
%     kelvasa2015_localize - uses calibration data to map bilateral spike rate differences to an azimuthal angle
