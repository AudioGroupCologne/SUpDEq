%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function results = supdeq_testHPlocalization( testHRTFdataset, referenceHRTFdataset, samplingGrid, fs, printResults)
%
% This function models horizontal plane localization for the testHRTFdataset
% and the referenceHRTFdataset on a pre-defined or given sampling grid. 
% The function uses the GMM-based estimation of azimuth direction 
% according to May2011 (see references), which is part of the 
% Auditory Modeling Toolbox (AMT). The result is a struct with the azimuthal 
% errors in degree.
%
% Output:
% results               - Output struct with errors and mean absolute
%                         errors of the test/referenceHRTFdataset
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
%                         (steps of 5° in azimuth, elevation = 90° (sound source in the front), 
%                         azimuth range +-90°)
% fs                    - Sampling Rate
%                         Default: 48000
% printResults          - Print results true/false
%                         Default:false
%
% Dependencies: Auditory Modeling Toolbox (AMT)
%
% References:
% T. May, S. van de Par, and A. Kohlrausch. A probabilistic model for
% robust localization based on a binaural auditory front-end. IEEE Trans
% Audio Speech Lang Proc, 19:1-13, 2011.
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function results = supdeq_testHPlocalization( testHRTFdataset, referenceHRTFdataset, samplingGrid, fs, printResults)

%Sampling grid on horizontal plane with +-90 degree in steps of 5 degree
if nargin < 3 || isempty(samplingGrid)
    [~, samplingGrid] = AKsubGrid(5, 'transverse', 0);
    %Transform elevation to required range
    samplingGrid(:,2) = 90-samplingGrid(:,2);
    %Set az to range between -90 and 90 degree
    ids = [find(samplingGrid(:,1)<=90);find(samplingGrid(:,1)>=270)];
    samplingGrid = samplingGrid(ids,:);    
end

if nargin < 4 || isempty(fs)
    fs = 48000;
end

if nargin < 5
    printResults = false;
end

%% Run model may2011

%Generate noise with a length of 1 second
noise = randn(fs,1);
%Set to -6dB peak level
noise = noise*((10^(-6/20))/(max(abs(noise))));

% Get HRIRs according to sampling grid
fprintf('Extracting 2 x %d HRIRs. This may take some time...\n',size(samplingGrid,1));
[HRIRs_Test_L,HRIRs_Test_R] = supdeq_getArbHRIR(testHRTFdataset,samplingGrid,[],[],'ak');
[HRIRs_Ref_L,HRIRs_Ref_R]   = supdeq_getArbHRIR(referenceHRTFdataset,samplingGrid,[],[],'ak');
fprintf('2 x %d HRIRs extracted...\n',size(samplingGrid,1))

disp('Running model according to May2011...');
az_Test = zeros(1,size(samplingGrid,1));
az_Ref  = zeros(1,size(samplingGrid,1));
for i = 1:size(samplingGrid,1)

    %Create binaural noise signal
    binNoise_Test(:,1) = conv(noise,HRIRs_Test_L(:,i));
    binNoise_Test(:,2) = conv(noise,HRIRs_Test_R(:,i));
    binNoise_Ref(:,1)  = conv(noise,HRIRs_Ref_L(:,i));
    binNoise_Ref(:,2)  = conv(noise,HRIRs_Ref_R(:,i));
    
    %Test azimuthal position with may2011
    azStruct_Test   = may2011(binNoise_Test,fs);
    azStruct_Ref    = may2011(binNoise_Ref,fs);
    %Histogram analysis of frame-based localization estimates
    az_Test(:,i)    = may2011_estAzimuthGMM(azStruct_Test,'HIST',1,0);
    az_Ref(:,i)     = may2011_estAzimuthGMM(azStruct_Ref,'HIST',1,0);
    
    if ~mod(i,5)
        fprintf('%d of %d azimuthal positions evaluated...\n',i,size(samplingGrid,1))
    end
    
end
    
%% (2) - Prepare az values of samplingGrid for output format of may2011

azSamplingGrid                      = samplingGrid(:,1);
azSamplingGrid(azSamplingGrid<=90)  = -1*azSamplingGrid(azSamplingGrid<=90);
azSamplingGrid(azSamplingGrid>90)   = 360-azSamplingGrid(azSamplingGrid>90);

%% (3) - Calculate errors and write results struct

results.azError_Test                = az_Test' - azSamplingGrid;
results.azAbsError_Test             = abs(results.azError_Test );
results.azError_Ref                 = az_Ref'  - azSamplingGrid;
results.azAbsError_Ref              = abs(results.azError_Ref  );
results.azMeanAbsError_Test         = mean(results.azAbsError_Test);
results.azMeanAbsError_Ref          = mean(results.azAbsError_Ref);
results.diff_azMeanAbsError         = results.azMeanAbsError_Ref - results.azMeanAbsError_Test;
results.azNonZeroMeanAbsError_Test  = mean(nonzeros(results.azAbsError_Test));
results.azNonZeroMeanAbsError_Ref   = mean(nonzeros(results.azAbsError_Ref));
results.diff_azNonZeroMeanAbsError  = results.azNonZeroMeanAbsError_Ref - results.azNonZeroMeanAbsError_Test;
results.azSamplingGrid              = azSamplingGrid;
results.samplingGrid                = samplingGrid;

disp('Done with horizontal plane localization test...');

%% (4) - Print results

if printResults
    fprintf('\nRESULTS \n');
    disp('------------------------------------------------------');
    fprintf('Mean absolute azimuthal error testHRTFdataset (deg) \t\t %4.1f \n',results.azMeanAbsError_Test)
    fprintf('Mean absolute azimuthal error referenceHRTFdataset (deg) \t %4.1f \n',results.azMeanAbsError_Ref)
    fprintf('Mean absolute azimuthal error difference (deg) \t\t\t %4.1f \n',results.diff_azMeanAbsError)
    disp('------------------------------------------------------');
    fprintf('Non zero mean absolute azimuthal error testHRTFdataset (deg) \t %4.1f \n',results.azNonZeroMeanAbsError_Test)
    fprintf('Non zero mean absolute azimuthal error referenceHRTFdataset (deg) %4.1f \n',results.azNonZeroMeanAbsError_Ref)
    fprintf('Non zero mean absolute azimuthal error difference (deg) \t %4.1f \n',results.diff_azNonZeroMeanAbsError)
    disp('------------------------------------------------------');
end

end

