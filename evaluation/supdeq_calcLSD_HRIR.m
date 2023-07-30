%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function results = supdeq_calcLSD_HRIR( testHRIRdataset, referenceHRIRdataset, fs, FFToversize, freqRange, freqRange2)
%
% This function calculates the log-spectral differences (LSD) between a test and a 
% reference HRTF dataset at given spatial sampling points
%
% Output:
% results               - Results struct with LSD per sampling point in dB 
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
% freqRange             - Define frequency range [lowFrequency highFrequency] 
%                         in Hz to calculate LSD in a specific frequency
%                         range.
%                         Default: [20 20000]
% freqRange2            - Optional another frequency range
%                         Default: nan
%
% Dependencies: -
%
% (C) 2023 by JMA, Johannes M. Arend
%             TU Berlin
%             Audio Communication Group
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function results = supdeq_calcLSD_HRIR( testHRIRdataset, referenceHRIRdataset, fs, FFToversize, freqRange, freqRange2)

%Default sampling grid (Lebedev grid with 2702 nodes)
if nargin < 3 || isempty(fs)
    fs = 48000;
end

if nargin < 4 || isempty(FFToversize)
    FFToversize = 1;
end

if nargin < 5 || isempty(freqRange)
    freqRange = [20 20000];
end

if nargin < 6 || isempty(freqRange2)
    freqRange2 = nan;
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

%Transform HRIRs to frequency domain
fprintf('Transforming 2 x %d HRIRs to frequency domain. This may take some time...\n',size(testHRIRdataset,2));
HRTFs_Test  = abs(fft(testHRIRdataset,NFFT));
HRTFs_Test  = HRTFs_Test(1:floor(end/2)+1,:); 
HRTFs_Ref   = abs(fft(referenceHRIRdataset,NFFT)); 
HRTFs_Ref   = HRTFs_Ref(1:floor(end/2)+1,:);

%Get frequency bins for f_low and f_high
[~,lowBin]  = min(abs(freqRange(1)-f));
[~,highBin]  = min(abs(freqRange(2)-f));
    
%Calculate LSD
disp('Calculating log-spectral differences...');
lsd_freq_point = sqrt((20*log10(HRTFs_Ref./HRTFs_Test)).^2);
lsd_freq = mean(lsd_freq_point,2); %Sum of frequency outside of sqrt to get frequency-dependent value
lsd_point = sqrt(mean((20*log10(HRTFs_Ref./HRTFs_Test)).^2,1)).'; %Sum of frequency inside of sqrt
lsd_point_freqRange = sqrt(mean((20*log10(HRTFs_Ref(lowBin:highBin,:)./HRTFs_Test(lowBin:highBin,:))).^2,1)).'; %Sum of frequency inside of sqrt, for frequency range
lsd = mean(lsd_point); %Sum of point outside of sqrt
lsd_freqRange = mean(lsd_point_freqRange); %Sum of point outside of sqrt

if ~isnan(freqRange2)
    %Get frequency bins for f_low and f_high
    [~,lowBin2]  = min(abs(freqRange2(1)-f));
    [~,highBin2]  = min(abs(freqRange2(2)-f));
    
    lsd_point_freqRange2 = sqrt(mean((20*log10(HRTFs_Ref(lowBin2:highBin2,:)./HRTFs_Test(lowBin2:highBin2,:))).^2,1)).'; %Sum of frequency inside of sqrt, for frequency range
    lsd_freqRange2 = mean(lsd_point_freqRange2); %Sum of point outside of sqrt
end

%Write results struct    
results.lsd                 = lsd;
results.lsd_freq            = lsd_freq;
results.lsd_point           = lsd_point;
results.lsd_freq_point      = lsd_freq_point;
results.lsd_freqRange       = lsd_freqRange; 
results.lsd_point_freqRange = lsd_point_freqRange;
if ~isnan(freqRange2)
results.lsd_freqRange2      = lsd_freqRange2; 
results.lsd_point_freqRange2= lsd_point_freqRange2;
end
results.f                   = f;
results.freqRange           = freqRange;
if ~isnan(freqRange2)
results.freqRange2           = freqRange2;
end

disp('Done with calculation of log-spectral differences...');

end

