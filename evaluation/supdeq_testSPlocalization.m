%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function results = supdeq_testSPlocalization( testHRTFdataset, referenceHRTFdataset, samplingGrid, S, lat, fs, printResults )
%
% This function models sagittal plane localization for the testHRTFdataset
% and the referenceHRTFdataset on a pre-defined or given sampling grid. 
% The function uses the model according to Baumgartner2014 (see references), 
% which is part of the Auditory Modeling Toolbox (AMT). The result is a 
% struct with quadrant errors (qe) in percent and polar errors (pe) in degree 
% for the test/reference HRTF dataset and the difference between both estimations.
%
% Output:
% results               - Output struct with qe (quadrant errors in %) and pe
%                         (local polar errors in degree) of the
%                         test/referenceHRTFdataset
%
% Input:        
% testHRTFdataset       - Struct with SH-coefficients of the test 
%                         HRTF dataset (for example the de-equalized HRTF set
%                         or any other HRTF set based on the referenceHRTFdataset) 
%                         for the left (Hl_nm) and right (Hr_nm) channel/ear, 
%                         absolute frequency scale f, transform order N, and FFToversize
% referenceHRTFdataset  - Struct with SH-coefficients of the reference 
%                         HRTF dataset for the left (Hl_nm) and right (Hr_nm) channel/ear, 
%                         absolute frequency scale f, transform order N, and FFToversize
% samplingGrid          - Spatial sampling grid. Must be a Qx2 matrix where 
%                         the first column holds the azimuth and the second 
%                         the elevation.
%                         If no samplingGrid (or []) is passed, the default
%                         sampling grid is used with equidistant values
%                         (steps of 1° in azimuth and elevation, azimuth
%                         range from 0-359°, elevation range from -30 - +90)
% S                     - Listener-specific sensitivity threshold 
%                         (threshold of the sigmoid link function representing 
%                         the psychometric link between transformation from the
%                         distance metric and similarity index) to S. 
%                         Default: 0.76 (according to Baumgartner2014)
% lat                   - Apparent lateral angle of the sound sources to be
%                         tested (azimuth)
%                         Default: 0 (0 degree - sound source in the front)   
% fs                    - Sampling Rate
%                         Default: 48000
% printResults          - Print results true/false
%                         Default:false
%
% Dependencies: Auditory Modeling Toolbox (AMT)
%
% References:
% R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source 
% localization in sagittal planes for human listeners. 
% The Journal of the Acoustical Society of America, 136(2):791-802, 2014. 
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function results = supdeq_testSPlocalization( testHRTFdataset, referenceHRTFdataset, samplingGrid, S, lat, fs, printResults )

% SamplingGrid according to Baumgartner2014
% (-30 to +60 elevation, steps of 1 degree in azimuth and elevation
if nargin < 3 || isempty(samplingGrid)
    [~, samplingGrid] = AKsubGrid(1, 'transverse', 90:-1:-30);
    %Transform elevation to required range
    samplingGrid(:,2) = 90-samplingGrid(:,2);
end

if nargin < 4 || isempty(S)
    S = 0.76;
end

if nargin < 5 || isempty(lat)
    lat = 0;
end

if nargin < 6 || isempty(fs)
    fs = 48000;
end

if nargin < 7
    printResults = false;
end

% Get values of sampling grid for defined lateral value (azimuth)
latIDFront  = find(samplingGrid(:,1) == lat);
latIDBack   = find(samplingGrid(:,1) == 180-lat); 
if isempty(latIDFront) || isempty(latIDBack)
    error('Could not find lateral target angle (lat, front or back) in sampling grid. Please re-define lateral target angle or sampling grid!');
else
    %Get also values for back (latID = 180)
    samplingGrid = [samplingGrid(latIDFront,:);samplingGrid(latIDBack,:)];
end

%% (1) - Get HRIRs according to sampling grid

fprintf('Extracting 2 x %d HRIRs. This may take some time...\n',size(samplingGrid,1));
[HRIRs_Test_L,HRIRs_Test_R] = supdeq_getArbHRIR(testHRTFdataset,samplingGrid);
[HRIRs_Ref_L,HRIRs_Ref_R]   = supdeq_getArbHRIR(referenceHRTFdataset,samplingGrid);
fprintf('2 x %d HRIRs extracted...\n',size(samplingGrid,1))
    
%% (2) - Save in SOFA format

testSOFA = supdeq_writeSOFAobj(HRIRs_Test_L',HRIRs_Test_R',samplingGrid, fs);
refSOFA  = supdeq_writeSOFAobj(HRIRs_Ref_L',HRIRs_Ref_R',samplingGrid, fs);

%Apply diffuse field compensation
testSOFA = SOFAhrtf2dtf(testSOFA);
refSOFA = SOFAhrtf2dtf(refSOFA);

%% (3) - Run model baumgartner2014

disp('Running model according to Baumgartner2014...');

qe_pe_test  = baumgartner2014(testSOFA, testSOFA, 'QE_PE_EB', 'S', S, 'lat', lat, ...
                                         'fs', fs, 'fsstim', fs);
qe_pe_ref  = baumgartner2014(refSOFA, refSOFA, 'QE_PE_EB', 'S', S, 'lat', lat, ...
                                         'fs', fs, 'fsstim', fs);
qe_pe_testVSref  = baumgartner2014(testSOFA, refSOFA, 'QE_PE_EB', 'S', S, 'lat', lat, ...
                                         'fs', fs, 'fsstim', fs);
                                     
qe_diff = qe_pe_ref.qe - qe_pe_testVSref.qe;
pe_diff = qe_pe_ref.pe - qe_pe_testVSref.pe;

%% (4) - Write results struct

results.qe_test         = qe_pe_test.qe;
results.pe_test         = qe_pe_test.pe;
results.qe_ref          = qe_pe_ref.qe;
results.pe_ref          = qe_pe_ref.pe;
results.qe_testVSref    = qe_pe_testVSref.qe;
results.pe_testVSref    = qe_pe_testVSref.pe;
results.qe_diff         = qe_diff;
results.pe_diff         = pe_diff;
results.samplingGrid    = samplingGrid;
results.S               = S;
results.lat             = [lat, 180-lat];

disp('Done with sagittal plane localization test...');

%% (5) - Print results

if printResults
    fprintf('\nRESULTS \n');
    disp('------------------------------------------------------');
    fprintf('Quadrant errors testHRTFdataset (%%) \t\t %4.1f \n',results.qe_test)
    fprintf('Local polar RMS error testHRTFdataset (deg) \t %4.1f \n',results.pe_test)
    disp('------------------------------------------------------');
    fprintf('Quadrant errors referenceHRTFdataset (%%) \t %4.1f \n',results.qe_ref)
    fprintf('Local polar RMS error referenceHRTFdataset (deg)  %4.1f \n',results.pe_ref)
    disp('------------------------------------------------------');
    fprintf('Quadrant errors testVSreference (%%) \t\t %4.1f \n',results.qe_testVSref)
    fprintf('Local polar RMS error testVSreference (deg) \t  %4.1f \n',results.pe_testVSref)
    disp('------------------------------------------------------');
    fprintf('Quadrant errors difference (%%) \t\t\t %4.1f \n',results.qe_diff)
    fprintf('Local polar RMS error difference (deg) \t\t %4.1f \n',results.pe_diff) 
    disp('------------------------------------------------------');
end

end

