%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function results = supdeq_calcSpectralDifference_HRIR( testHRIRdataset, referenceHRIRdataset, fs, FFToversize, filterMode, plotResults)
%
% This function calculates the spectral differences between a test and a 
% reference HRTF dataset at given spatial sampling points
%
% Output:
% results               - Results struct with spectral difference per
%                         frequency and sampling point in dB (20*log10)
%                         (specDifferencePerPoint), spectral 
%                         difference averaged over all spatial sampling 
%                         points in dB (specDifference) as well 
%                         as mean of specDifference in dB 
%                         --> Spectral difference averaged over all 
%                         spatial sampling points and frequencies
%                         (meanSpecDifference). Additionally, 
%                         meanSpecDifference_freqRange provides the mean
%                         spectral difference in a frequency range between
%                         100Hz and 12kHz. "f" is a frequency vector
%
%                         Optional: filterMode = true
%                         Results struct with spectral difference averaged
%                         over all spatial sampling points in dB
%                         (20*log10) and for every of the 40 filter bands
%                         ("specDiffPerBand"). Additionally,
%                         "specDiffL_meanPerBand" returns the mean spectral
%                         difference per band (mean per band of
%                         specDiffPerBand), fc returns the center
%                         frequencies of the auditory filterbank, and f
%                         is a frequency vector.
%
% Input:        
% testHRIRdataset       - [N x M] array with HRIRs of the test HRIR
%                         dataset, with N samples and M channels                       
% referenceHRIRdataset  - [N x M] array with HRIRs of the reference HRIR
%                         dataset, with N samples and M channels  
% fs                    - Sampling Rate
%                         Default: 48000
% FFToversize           - FFToversize rises the FFT Blocksize [default = 1]
%                         A FFT of the blocksize (FFToversize*NFFT) is applied
%                         to the time domain data,  where  NFFT is determinded
%                         as the next power of two of the signalSize  which is
%                         signalSize = (lastSample-firstSample).
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

function results = supdeq_calcSpectralDifference_HRIR( testHRIRdataset, referenceHRIRdataset, fs, FFToversize, filterMode, plotResults)

%Default sampling grid (Lebedev grid with 2702 nodes)
if nargin < 3 || isempty(fs)
    fs = 48000;
end

if nargin < 4 || isempty(FFToversize)
    FFToversize = 1;
end

if nargin < 5 || isempty(filterMode)
    filterMode = false;
end

if nargin < 6
    plotResults = false;
end

%% Get NFFT based on number of samples in HRIR arrays and FFToversize
NFFT = size(testHRIRdataset,1);
if size(referenceHRIRdataset,1) > NFFT
    NFFT = size(referenceHRIRdataset,1);
end

NFFT    = (2^nextpow2(NFFT));
if FFToversize > 1
    NFFT = NFFT*FFToversize;
end
%Get frequency vector
f = linspace(0,fs/2,NFFT/2+1);

if size(testHRIRdataset,2) ~= size(referenceHRIRdataset,2)
    error('Number of sampling points in both HRIR datasets must be equal!');
end

numPoints = size(testHRIRdataset,2);

%% Filter mode off
if ~filterMode
    
    %Transform HRIRs to frequency domain
    fprintf('Transforming 2 x %d HRIRs to frequency domain. This may take some time...\n',size(testHRIRdataset,2));
    HRTFs_Test  = abs(fft(testHRIRdataset,NFFT));
    HRTFs_Test  = HRTFs_Test(1:end/2+1,:); 
    HRTFs_Ref   = abs(fft(referenceHRIRdataset,NFFT)); 
    HRTFs_Ref   = HRTFs_Ref(1:end/2+1,:);
        
    %Calculate spectral difference
    disp('Calculating spectral differences...');
    specDiff = abs( 20*log10(abs(HRTFs_Ref)) - 20*log10(abs(HRTFs_Test)));
    specDiffPerPoint = specDiff;
    specDiff = sum(specDiff,2);
    specDiff = specDiff/numPoints;
    
    %Get values from 100Hz - 12 kHz
    lowBin  = round(size(f,2) / f(end) * 100) + 1;
    highBin = round(size(f,2) / f(end) * 12000) + 1; 

    %Write results struct
    results.specDifferencePerPoint          = specDiffPerPoint;        
    results.specDifference                  = specDiff;
    results.meanSpecDifference              = mean(results.specDifference);
    results.meanSpecDifference_freqRange    = mean(results.specDifference(lowBin:highBin,:));
    results.f                               = f;

    disp('Done with calculation of spectral differences...');
    
    if plotResults
        figure;
        semilogx(results.f,results.specDifference(:,1),'k','LineWidth',1.5);
        grid on;
        xlim([20,20000]);
        xlabel('Frequency [Hz]');
        ylabel('Spectral difference [dB]');
    end
    
end

%% Filter mode on
if filterMode
    
    %Transform HRIRs to 3D array
    HRIRs_Test  = zeros(size(testHRIRdataset,1), 1, numPoints);
    HRIRs_Ref   = zeros(size(referenceHRIRdataset,1), 1, numPoints);
    
    %Kind of a workaround, but 4-D arrays needed later....
    HRIRs_Test(:,1,:)   = reshape(testHRIRdataset,[size(HRIRs_Test,1),1,size(HRIRs_Test,3)]);
    HRIRs_Ref(:,1,:)    = reshape(referenceHRIRdataset,[size(HRIRs_Ref,1),1,size(HRIRs_Ref,3)]); 
    
    %Zero-pad to 512 samples (if length is smaller) to get higher frequency
    %resolution. This only works if HRIRs_Test and HRIRs_Ref have equal
    %length!
    if NFFT < 512
        newLength = 512;
    else
        newLength = NFFT;
    end
    zeroPad = zeros(newLength-size(HRIRs_Test,1),1,numPoints);
    HRIRs_Test  = [HRIRs_Test;zeroPad];
    HRIRs_Ref   = [HRIRs_Ref;zeroPad];
    %Calc f needed later for results struct
    f = linspace(0,fs/2,newLength/2+1);
       
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
    specDiff = squeeze(abs(20*log10(abs(HRTFs_Ref_filt(:,:,1,:))) - 20*log10(abs(HRTFs_Test_filt(:,:,1,:)))));
    specDiffPerPoint = specDiff;
    %Average across spatial sampling points
    specDiff = sum(specDiff,3);
    for i = 1:size(specDiff,2)
        specDiff(:,i) = specDiff(:,i)/numPoints;
    end
    %Get mean per band
    for i = 1:size(specDiff,2)
        specDiff_meanPerBand(i,:) = mean(specDiff(:,i));
    end
    
    %Write results struct
    results.specDiffPerPoint        = specDiffPerPoint; 
    results.specDiffPerBand         = specDiff;
    results.specDiff_meanPerBand    = specDiff_meanPerBand;
    results.fc                      = fc;
    results.f                       = f;

    disp('Done with calculation of spectral differences...');
    
    if plotResults
        figure;
        semilogx(results.fc,results.specDiff_meanPerBand,'k','LineWidth',1.5);
        grid on;
        xlim([20,20000]);
        xlabel('Frequency [Hz]');
        ylabel('Spectral Difference [dB]');
    end
    
end


end

