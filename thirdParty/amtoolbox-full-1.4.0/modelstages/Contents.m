% AMT - Various stages of auditory models
%
%   Basilar membrane velocity
%     hohmann2002_clearstate               - Clears the state of the filterbank object         
%     hohmann2002_delay                    - Create a delay object    
%     hohmann2002_filter                   - Create a single Gammatone filter object     
%     hohmann2002_freqz                    - Calculates frequency response of a filterbank object    
%     hohmann2002_mixer                    - Create a mixer object
%     hohmann2002_process                  - Process the input signals by the corresponding filterbank object
%     hohmann2002_synth                    - Create a synthesis object
%     lyon2011_agcstep                     - Active gain control update step
%     lyon2011_carstep                     - One sample-time update step for the filter part of the model
%     lyon2011_closeagcloop                - Active gain control loop
%     lyon2011_crosscouple                 - Adjust the intensity of the ear signals
%     lyon2011_design                      - Computes all the filter coefficients needed to run the CARFAC model
%     lyon2011_detect                      - Calculates conductance using a sigmoidal detection nonlinearity
%     lyon2011_ihcstep                     - Update step of inner-hair-cell (IHC) model
%     lyon2011_init                        - Allocates the state vector storage in the CARFAC model
%     lyon2011_spatialsmooth               - Spatial smoothing of FIR coefficients
%     lyon2011_stageg                      - Return the stage gain g needed to get unity gain at DC
%
%   Auditory nerve
%     bruce2018_ffgn                       - Fast (exact) fractional Gaussian noise and Brownian motion generator.
%     bruce2018_fitaudiogram              - Fit audiogram to threshold shifts      
%     bruce2018_generateanpopulation       - Generate auditory nerve fibre population
%     bruce2018_innerhaircells             - Calculation of inner haircell potential
%     bruce2018_synapse                    - Calculation of synapse output
%     zilany2014_ffgn                      - Fast (exact) fractional Gaussian noise and Brownian motion generator.
%     zilany2014_innerhaircells            - Calculation of inner haircell potential
%     zilany2014_synapse                   - Calculation of synapse output
%
%   Temporal modulation sensitivity
%     carney2015_fitaudiogram             - Fit audiogram to threshold shifts     
%     carney2015_generateneurogram         - Generate a neurogram
%     carney2015_getalphanorm              - Returns filter coefficients for a normalized alpha function
%     king2019_modfilterbank               - Modulation filterbank
%     relanoiborra2019_drnl                - Dual-resonance nonlinear filterbank
%     relanoiborra2019_mfbtd               - Modulation filterbank
%     relanoiborra2019_featureextraction             - Preprocessing stage
%     roenne2012_chirp                     - Simulate chirp evoked ABRs
%     roenne2012_click                     - Simulate ABR respone to click
%     roenne2012_tonebursts                - Simulate tone burst evoked ABR wave V latencies
%     verhulst2018_auditorynerve           - Auditory nerve models used by Verhulst et al. 2018 and 2015
%     verhulst2015_cn                      - Cochlear nucleus model
%     verhulst2015_ic                      - Inferior colliculus model
%     verhulst2018_ihctransduction         - IHC transduction in Verhulst et al. 2018
%
%   Binaural processing
%     bischof2023_filterbank               - Filterbank of Bischof et al. model
%     breebaart2001_centralproc            - Central processor taking decision in an experiment
%     breebaart2001_eicell                 - Excitation-inhibition cell from Breebaart et al. (2001)
%     breebaart2001_outmiddlefilter        - Outer and middle-ear filter used by the model Breebaart et al. (2001)
%     dietz2011_filterbank                 - Filterbank of Dietz 2011 binaural model  
%     dietz2011_interauralfunctions        - Calculate interaural parameters for Dietz 2011 model
%     dietz2011_unwrapitd                  - IPD to ITD transformation for the Dietz model
%     eurich2022_processing                - Binaural processing stage of the Eurich et al. (2022) model
%     eurich2022_decision                  - Decision stage of the Eurich et al. (2022) model 
%     lindemann1986_bincorr                - Running cross-correlation between two signals
%     lindemann1986_centroid               - Centroid of the cross-correlation activation
%     tabuchi2016_estimatethreshold        - Threshold estimation
%     takanen2013_contracomparison         - Enhance contrast between hemispheres
%     takanen2013_cueconsistency           - Check consistency before cue combination
%     takanen2013_directionmapping         - Map the directional cues to directions
%     takanen2013_formbinauralactivitymap  - Steer cues on a topographic map
%     takanen2013_lso                      - Model of the lateral superior olive
%     takanen2013_mso                      - Model of the medial superior olive
%     takanen2013_onsetenhancement         - Emphasize onsets on direction analysis
%     takanen2013_periphery                - Process input through the model of periphery
%     takanen2013_wbmso                    - Wideband medial superior olive model
%     takanen2013_weightedaveragefilter    - Part of the takanen2013 model
%     vicente2020_betterearsnrframe        - Better ear SNR from audiogram
%     vicente2020_buadvantage              - Better ear SNR from audiogram
%     vicente2020_internalnoise            - Internal noise calculation
%
%   Loudness
%     moore2016_agcnextframe               - Adjusts successive short term loudness frames
%     moore2016_binauralloudness           - Calculate the binaural loudness
%     moore2016_cochlea                    - Outer and middle ear filtering
%     moore2016_excitationpattern          - Calculate the excitation patterns
%     moore2016_longtermloudness           - Calculate the long term loudness
%     moore2016_monauralinstspecloudness   - Calculate instantaneous specific loudness over time
%     moore2016_shorttermspecloudness      - Calculate the short-term specific loudness
%     moore2016_spectrum                   - Calculate the spectrum for an audio segment of 2048x2 samples
%
%   Monaural speech perception
%     joergensen2011_combineinformation    - Combines the SNRenv across modulation and audio filters 
%     joergensen2011_multchansnrenv        - Calculates the SNRenv
%     joergensen2011_overlapadd3           - Overlap-add calculation for an FFT matrix
%     joergensen2011_pctodsrt              - Calculates the SRT and change in SRT from the simulated percent correct
%     joergensen2011_sim                   - Simulate the experiments of Jørgensen and Dau (2011)
%     joergensen2011_sepsub               - Calculates an estimate of the clean signal using spectral subtraction
%     joergensen2013_sim                   - Simulate the experiments of Jørgensen, Ewert and Dau (2013)
%
%   Binaural speech perception
%     hauth2020_ecprocess4optsigs          - Applies the Equalization-Cancellation process in the frequency domain
%     hauth2020_fftcon                     - Calculates the interaural delay of the dominant source
%     hauth2020_sii                        - Calculates the SII according to ANSI S3.5-1997
%     hauth2020_srmr                       - Computes the speech-to-reverberation modulation energy ratio
%
%   Spatial perception
%     barumerli2023_metrics                - Extract localization metrics
%     barumerli2023_featureextraction                - Extract HRTF using gammatone frequency bands and ITDs from SOFA object
%     baumgartner2013_calibration          - Calibration of the model (linear periphery)
%     baumgartner2013_pmv2ppp              - Calculate performance predictions from PMVs for baumgartner2013
%     baumgartner2014_binauralweighting    - Binaural combination of monaural similarity estimates
%     baumgartner2014_calibration          - Calibration of the model (linear periphery)
%     baumgartner2014_comparisonprocess    - Comparison with direction-specific templates
%     baumgartner2014_gradientextraction   - Extraction of positive spectral gradients
%     baumgartner2014_likelistat           - Likelihood statistics for evaluation of model performance
%     baumgartner2014_parametrization      - Joint optimization of model parameters
%     baumgartner2014_pmv2ppp              - Performance predictions from PMVs of baumgartner2014
%     baumgartner2014_sensorimotormapping  - Response scatter induced by localization task
%     baumgartner2014_similarityestimation - Similarity estimation with listener-specific sensitivity
%     baumgartner2014_spectralanalysis     - Approximation of spectral analysis by auditory periphery
%     baumgartner2014_virtualexp           - Performs a virtual sound-localization experiment
%     baumgartner2016_calibration          - Calibration of the model (nonlinear periphery)
%     baumgartner2016_comparisonprocess    - Comparison with direction-specific templates
%     baumgartner2016_gradientextraction   - Extraction of positive spectral gradients
%     baumgartner2016_parametrization      - Joint optimization of model parameters
%     baumgartner2016_spectralanalysis     - Approximation of spectral analysis by auditory periphery
%     baumgartner2021_mapping              - Mapping of externalization scores        
%     kelvasa2015_anbinning                - AN and time binning from Kelvasa and Dietz 2015 binaural model
%     kelvasa2015_anprocessing             - AN model used in Kelvasa and Dietz 2015 binaural model
%     kelvasa2015_calibratemapping         - Produces necessary mappings for localization model
%     kelvasa2015_ciprocessing             - CI ACE processing strategy used in Kelvasa and Dietz 2015 binaural model
%     kelvasa2015_localize                 - Uses calibration data to map bilateral spike rate differences to an azimuthal angle
%     langendijk2002_comp                  - Comparison process from Langendijk 20002.
%     langendijk2002_likelihood            - Likelihood estimation
%     langendijk2002_spectralanalysis      - FFT-based filter bank with constant relative bandwidth
%     llado2022_binauralfeats              - binaural auditory model based on the example 13.6.2 in: Pulkki and Karjalainen: Communication  acoustics.
%     llado2022_evaluatenn                 - evaluate the neural network
%     llado2022_extractirs                 - extract impulse responses
%     llado2022_trainnn                    - train the neural network
%     llado2022_weightsanalysis            - analyses the weights learnt by the neural network
%     may2011_cbarlabel                    - Sets the labels of a cbar plot
%     may2011_classifygmm                  - Gaussian mixture model
%     may2011_estazimuthgmm                - Estimate azimuth from Gaussian Mixture Model output
%     may2011_findlocalpeaks               - Finds peaks with optional quadratic interpolation
%     may2011_fireprint                    - Colormap that increases linearly in lightness
%     may2011_framedata                    - Frame data
%     may2011_gammatoneinit                - Initialize gammatone filterbank structure
%     may2011_gammatone                    - Gammatone filterbank
%     may2011_interpolateparabolic         - Multi-channel parabolic interpolation
%     may2011_neuraltransduction           - Calculates the auditory nerve response
%     may2011_xcorrnorm                    - Normalized time-domain cross-correlation function
%     mclachlan2021_metrics                - Extract localization metrics
%     mclachlan2021_preproc                - Preprocessing for Mclachlan et al. (2021)
%     mclachlan2021_rotatedirs             - Rotate a set of coordinates on a sphere
%     reijniers2014_metrics                - Extract localization metrics (errors) as in Reijniers et al. 2014
%     reijniers2014_featureextraction                - Extract HRTF using gammatone, frequency bands, and ITDs from a SOFA object
%     wierstorf2013_estimateazimuth        - Azimuth position estimation based on dietz2011 or lindemann1986
%     ziegelwanger2013_offaxis             - Off-axis version of Ziegelwanger et al. (2013) (legacy only)
%     ziegelwanger2013_onaxis              - On-axis version of Ziegelwanger et al. (2013) (legacy only)
%     ziegelwanger2014_offaxis             - Off-axis version of Ziegelwanger et al. (2014)
%     ziegelwanger2014_onaxis              - On-axis version of Ziegelwanger et al. (2014)
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/Contents.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


