%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function hrtfDataset = supdeq_sofa2hrtf(SOFAobj, N, samplingGrid, FFToversize)
%
% This function transformes a 'SimpleFreeFieldHRIR' SOFA file to frequency 
% domain and writes a hrtf dataset struct suitable for the SUpDEq toolbox.
%
% Output:
% hrtfDataset   - Struct with HRTFsfor the left (HRTF_L) and right (HRTF_R) 
%                 channel/ear, absolute frequency scale f, 
%                 maximum transform order N, FFToversize, and samplingGrid
%
% Input:
% SOFAobj       - Measurement data in form of a SOFA object
% N             - Maximal transform order N
% samplingGrid  - Optional input [default = nan]
%                 If a spatial sampling grid is passed, it must be a Qx3 
%                 matrix where the first column holds the azimuth, the 
%                 second the elevation, and the third the sampling weights.
%                 In this case, the sampling grid with weights will be
%                 written to the struct. Otherwise, the sampling grid
%                 stored in the SOFA file, transformed to SH coordinates,
%                 will be written to the struct.
%                 Azimuth positions in degree. Can be scalar or vector of size(el)
%                 (0=front, 90=left, 180=back, 270=right)
%                 (0 points to positive x-axis, 90 to positive y-axis)
%                 Elevations in degree. Can be scalar or vector of size(az)
%                 (0=North Pole, 90=front, 180=South Pole)
%                 (0 points to positive z-axis, 180 to negative z-axis)
% FFToversize   - FFToversize rises the FFT Blocksize [default = 1]
%                 A FFT of the blocksize (FFToversize*NFFT) is applied
%                 to the time domain data,  where  NFFT is determinded
%                 as the next power of two of the signalSize  which is
%                 signalSize = (lastSample-firstSample).
%
% Dependencies: SOFA API
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function hrtfDataset = supdeq_sofa2hrtf(SOFAobj, N, samplingGrid, FFToversize)

if nargin < 3
    samplingGrid = [];
end

if nargin < 4
    FFToversize = 1;
end

%Check FFToversize parameter
if FFToversize < 1
    warning('FFToversize needs to >= 1. Set to default = 1.');
    FFToversize = 1;
end

%Check grid if passed
if ~isempty(samplingGrid)
    if size(samplingGrid,2) > size(samplingGrid,1)
        samplingGrid = samplingGrid';
    end
    
    if size(samplingGrid,2) ~= 3
        warning('No sampling weights passed. Using sampling grid defined in SOFA object.')
        samplingGrid = [];
    end
end

%Check SOFA object
if ~strcmp(SOFAobj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR')
    error('Function only valid for SOFA convention "SimpleFreeFieldHRIR"');
end

%% Transform data to frequency domain

%Get IRs
irL = squeeze(SOFAobj.Data.IR(:,1,:)); 
irR = squeeze(SOFAobj.Data.IR(:,2,:));

%Zero pad to get length according to FFToversize
if FFToversize > 1
    NFFT = 2^nextpow2(size(irL,2));
    NFFT = NFFT*FFToversize;
    irL = [irL, zeros(size(irL,1),NFFT-size(irL,2))];
    irR = [irR, zeros(size(irR,1),NFFT-size(irR,2))];
else
    if size(irL,2) == 2^nextpow2(size(irL,2))
        NFFT = size(irL,2);
    else
        NFFT = 2^nextpow2(size(irL,2));
        irL = [irL, zeros(size(irL,1),NFFT-size(irL,2))];
        irR = [irR, zeros(size(irR,1),NFFT-size(irR,2))];
    end
end

%Transform HRIRs to frequency domain (HRTFs)
HRTF_L = fft(irL,[],2);
HRTF_R = fft(irR,[],2);
%Cut at fs/2
HRTF_L = HRTF_L(:,1:end/2+1);
HRTF_R = HRTF_R(:,1:end/2+1);

%If no samplingGrid was passed, use samplingGrid from SOFA
if isempty(samplingGrid)
    %Get samplingGrid from SOFA object
    samplingGrid = SOFAobj.SourcePosition(:,1:2);

    %Transform samplingGrid to SH coordinate system
    samplingGrid(:,2) = 90-samplingGrid(:,2);
end

%Write f
f = linspace(0,SOFAobj.Data.SamplingRate/2,NFFT/2+1);

%Write output struct
hrtfDataset.HRTF_L          = HRTF_L;
hrtfDataset.HRTF_R          = HRTF_R;
hrtfDataset.f               = f;
hrtfDataset.Nmax            = N;
hrtfDataset.FFToversize     = FFToversize;
hrtfDataset.samplingGrid    = samplingGrid;
hrtfDataset.sourceDistance  = SOFAobj.SourcePosition(1,3); 

end

