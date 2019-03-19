%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function results = supdeq_calcSpectralDifference( testHRTFdataset, referenceHRTFdataset, samplingGrid, fs, filterMode, plotResults)
%
% This function calculates the spectral differences between a test and a 
% reference HRTF dataset at given spatial sampling points
%
% Output:
% results               - Results struct with spectral difference averaged 
%                         over all spatial sampling points in dB (20*log10) 
%                         for left and right channel (specDifference) as well as
%                         mean of specDifference in dB for left and right channel
%                         --> Spectral difference averaged over all spatial 
%                         sampling points and frequencies
%                         (meanSpecDifference). Additionally, 
%                         meanSpecDifference_freqRange provides the mean
%                         spectral difference in a frequency range between
%                         100Hz and 12kHz. "f" is a frequency vector
%
%                         Optional: filterMode = true
%                         Results struct with spectral difference averaged
%                         over all spatial sampling points in dB
%                         (20*log10) for the left and the right channel 
%                         and for every of the 40 filter bands
%                         ("specDiffPerBand_L/R"). Additionally,
%                         "specDiffL_meanPerBand" returns the mean spectral
%                         difference per band (mean per band of
%                         specDiffPerBand_L/R), fc returns the center
%                         frequencies of the auditory filterbank, and f
%                         is a frequency vector.
%
% Input:        
% testHRTFdataset       - Struct with (SH-Coefficients) of the test 
%                         HRTF dataset (for example the de-equalized HRTF set
%                         or any other HRTF set based on the referenceHRTFdataset) 
%                         for the left (Hl_nm) and right (Hr_nm), channel/ear, 
%                         absolute frequency scale f, transform order N, and FFToversize
% referenceHRTFdataset  - Struct with (SH-Coefficients) of the reference 
%                         HRTF dataset for the left (Hl_nm) and right (Hr_nm), channel/ear, 
%                         absolute frequency scale f, transform order N, and FFToversize
%                         NOTE: testHRTFdataset and referenceHRTFdataset
%                         need to have the same FFT size and FFToversize!!
% samplingGrid          - Spatial sampling grid. Must be a Qx2 matrix where 
%                         the first column holds the azimuth and the second 
%                         the elevation.
%                         If no samplingGrid (or []) is passed, the default
%                         sampling grid is passed (Lebedev grid with 2702
%                         nodes)
% fs                    - Sampling Rate
%                         Default: 48000
% filterMode            - Filter mode true/false. Use auditory filterbank
%                         with 40 bands (fLow = 50 Hz, fHigh = 20kHz)
%                         Default: false
% plotResults          -  plotResults true/false. Simple plot...
%                         Default: false
%
% Dependencies: Auditory Modeling Toolbox (AMT), if filterMode = true
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function results = supdeq_calcSpectralDifference( testHRTFdataset, referenceHRTFdataset, samplingGrid, fs, filterMode, plotResults)

%Default sampling grid (Lebedev grid with 2702 nodes)
if nargin < 3 || isempty(samplingGrid)
    samplingGrid = supdeq_lebedev(2702);
end

if nargin < 4 || isempty(fs)
    fs = 48000;
end

if nargin < 5 || isempty(filterMode)
    filterMode = false;
end

if nargin < 6
    plotResults = false;
end


%% Filter mode off
if ~filterMode
    
    %Get HRTFs according to sampling grid
    fprintf('Extracting 2 x %d HRTFs. This may take some time...\n',size(samplingGrid,1));
    [HRTFs_Test_L,HRTFs_Test_R] = supdeq_getArbHRTF(testHRTFdataset,samplingGrid);
    [HRTFs_Ref_L,HRTFs_Ref_R]   = supdeq_getArbHRTF(referenceHRTFdataset,samplingGrid);
    fprintf('2 x %d HRTFs extracted...\n',size(samplingGrid,1))
    
    %Calculate spectral difference
    disp('Calculating spectral differences...');
    specDiffL = abs( 20*log10(abs(HRTFs_Ref_L.')) - 20*log10(abs(HRTFs_Test_L.')));
    specDiffR = abs( 20*log10(abs(HRTFs_Ref_R.')) - 20*log10(abs(HRTFs_Test_R.')));
    specDiffL = sum(specDiffL,2);
    specDiffR = sum(specDiffR,2);
    specDiffL = specDiffL/size(samplingGrid,1);
    specDiffR = specDiffR/size(samplingGrid,1);
    
    %Get values from 100Hz - 12 kHz
    lowBin  = round(size(referenceHRTFdataset.f,2) / referenceHRTFdataset.f(end) * 100) + 1;
    highBin = round(size(referenceHRTFdataset.f,2) / referenceHRTFdataset.f(end) * 12000) + 1; 

    %Write results struct
    results.specDifference                  = [specDiffL,specDiffR];
    results.meanSpecDifference              = mean(results.specDifference);
    results.meanSpecDifference_freqRange    = mean(results.specDifference(lowBin:highBin,:));
    results.f                               = referenceHRTFdataset.f;

    disp('Done with calculation of spectral differences...');
    
    if plotResults
        figure;
        semilogx(results.f,results.specDifference(:,1),'k','LineWidth',1.5);
        hold on;
        semilogx(results.f,results.specDifference(:,2),'r','LineWidth',1.5);
        grid on;
        xlim([20,20000]);
        xlabel('Frequency [Hz]');
        ylabel('Spectral difference [dB]');
        legend('Left','Right','Location','NorthWest');
    end
    
end

%% Filter mode on
if filterMode
    
    %Get HRIRs according to sampling grid
    HRIRs_Test = zeros((size(testHRTFdataset.f,2)*2-2) / testHRTFdataset.FFToversize,2,size(samplingGrid,1));
    HRIRs_Ref = zeros((size(referenceHRTFdataset.f,2)*2-2) / referenceHRTFdataset.FFToversize,2,size(samplingGrid,1));
    
    fprintf('Extracting 2 x %d HRIRs. This may take some time...\n',size(samplingGrid,1));
    [HRIRs_Test_L,HRIRs_Test_R] = supdeq_getArbHRIR(testHRTFdataset,samplingGrid);
    [HRIRs_Ref_L,HRIRs_Ref_R]   = supdeq_getArbHRIR(referenceHRTFdataset,samplingGrid);
    %Kind of a workaround, but 4-D arrays needed later....
    HRIRs_Test(:,1,:)   = reshape(HRIRs_Test_L,[size(HRIRs_Test,1),1,size(HRIRs_Test,3)]);
    HRIRs_Test(:,2,:)   = reshape(HRIRs_Test_R,[size(HRIRs_Test,1),1,size(HRIRs_Test,3)]);
    HRIRs_Ref(:,1,:)    = reshape(HRIRs_Ref_L,[size(HRIRs_Ref,1),1,size(HRIRs_Ref,3)]); 
    HRIRs_Ref(:,2,:)    = reshape(HRIRs_Ref_R,[size(HRIRs_Ref,1),1,size(HRIRs_Ref,3)]);
    fprintf('2 x %d HRIRs extracted...\n',size(samplingGrid,1))
    
    %Zero-pad to 512 samples (if length is smaller) to get higher frequency
    %resolution. This only works if HRIRs_Test and HRIRs_Ref have equal
    %length!
    newLength = 512;
    if size(HRIRs_Test,1) < newLength
        zeroPad = zeros(newLength-size(HRIRs_Test,1),2,size(samplingGrid,1));
        HRIRs_Test  = [HRIRs_Test;zeroPad];
        HRIRs_Ref   = [HRIRs_Ref;zeroPad];
        %Calc f needed later for results struct
        f = linspace(0,fs/2,newLength/2+1);
    else
        f = linspace(0,fs/2,size(HRIRs_Test,1)/2+1);
    end
       
    %Filter HRIRs with auditory filter bank (settings lead to 40 bands)
    %Output is a 4D array with [samples X filterBand x channel x
    %samplingPoint]
    disp('Applying auditory filter bank. This may take some time...');
    HRIRs_Test_filt = auditoryfilterbank(HRIRs_Test,fs,'flow',50,'fhigh',20000);
    [HRIRs_Ref_filt,fc]  = auditoryfilterbank(HRIRs_Ref,fs,'flow',50,'fhigh',20000);
    disp('Done with filtering...');
    
    %Perform FFT
    disp('Calculating spectral differences...');
    HRTFs_Test_filt = fft(HRIRs_Test_filt);
    HRTFs_Test_filt = HRTFs_Test_filt(1:size(HRTFs_Test_filt,1)/2+1,:,:,:);
    HRTFs_Ref_filt  = fft(HRIRs_Ref_filt);
    HRTFs_Ref_filt = HRTFs_Ref_filt(1:size(HRTFs_Ref_filt,1)/2+1,:,:,:);
    
    %Calculate spectral difference for each band
    specDiffL = squeeze(abs(20*log10(abs(HRTFs_Ref_filt(:,:,1,:))) - 20*log10(abs(HRTFs_Test_filt(:,:,1,:)))));
    specDiffR = squeeze(abs(20*log10(abs(HRTFs_Ref_filt(:,:,2,:))) - 20*log10(abs(HRTFs_Test_filt(:,:,2,:)))));
    %Average across spatial sampling points
    specDiffL = sum(specDiffL,3);
    specDiffR = sum(specDiffR,3);
    for i = 1:size(specDiffL,2)
        specDiffL(:,i) = specDiffL(:,i)/size(samplingGrid,1);
        specDiffR(:,i) = specDiffR(:,i)/size(samplingGrid,1);
    end
    %Get mean per band
    for i = 1:size(specDiffL,2)
        specDiffL_meanPerBand(i,:) = mean(specDiffL(:,i));
        specDiffR_meanPerBand(i,:) = mean(specDiffR(:,i));
    end
    
    %Write results struct
    results.specDiffPerBand_L       = specDiffL;
    results.specDiffPerBand_R       = specDiffR;
    results.specDiffL_meanPerBand   = specDiffL_meanPerBand;
    results.specDiffR_meanPerBand   = specDiffR_meanPerBand;
    results.fc                      = fc;
    results.f                       = f;

    disp('Done with calculation of spectral differences...');
    
    if plotResults
        figure;
        semilogx(results.fc,results.specDiffL_meanPerBand,'k','LineWidth',1.5);
        hold on;
        semilogx(results.fc,results.specDiffR_meanPerBand,'r','LineWidth',1.5);
        grid on;
        xlim([20,20000]);
        xlabel('Frequency [Hz]');
        ylabel('Spectral Difference [dB]');
        legend('Left','Right','Location','NorthWest');
    end
    
end


end

