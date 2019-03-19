%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function HRIRs_sfd = supdeq_sofa2sfd(SOFAobj, N, samplingGrid, FFToversize)
%
% This function transformes spherical HRIR measurements in SOFA format
% ('SimpleFreeFieldHRIR' convention) to the SH-domain by spherical Fourier 
% transform and writes a output struct suitable for the SUpDEq toolbox.
%
% Output:
% HRIRs_sfd     - Struct with Spherical Harmonic coefficients 
%                 (SH-coefficients) for the left (Hl_nm) and right (Hr_nm) 
%                 channel/ear, absolute frequency scale f, 
%                 transform order N, and FFToversize
%
% Input:
% SOFAobj       - Measurement data in form of a SOFA object
% N             - Transform order N
% samplingGrid  - Optional input [default = nan]
%                 If a spatial sampling grid is passed, it must be a Qx3 
%                 matrix where the first column holds the azimuth, the 
%                 second the elevation, and the third the sampling weights.
%                 In this case, the SOFiA functions "sofia_fdt" and 
%                 "sofia_stc" are used.
%                 If no samplingGrid (or []) is passed, the sampling grid
%                 according to the SOFA object is used for a spherical
%                 Fourier transform without sampling weights. In this case,
%                 the function AKsht is used.
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
% Dependencies: SOFiA toolbox, AKtools, SOFA API
%
% References:
% Benjamin Bernschütz: Microphone Arrays and Sound Field Decomposition 
% for Dynamic Binaural Recording. Ph.D. dissertation, Technical University
% Berlin (2016).
%
% Boaz Rafaely: Fundamentals of spherical array processing. In. Springer 
% topics in signal processing. Benesty, J.; Kellermann, W. (Eds.), 
% Springer, Heidelberg et al. (2015).
% 
% Franz Zotter: Analysis and synthesis of sound-radiation with spherical 
% arrays. Ph.D. dissertation, University of Music and Performing arts (2009). 
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function HRIRs_sfd = supdeq_sofa2sfd(SOFAobj, N, samplingGrid, FFToversize)

if nargin < 3
    samplingGrid = [];
end

if nargin < 4
    FFToversize = 1;
end

%Check FFToversize parameter
if FFToversize < 1
    warning('FFToversize needs to >= 1. Set to default = 1');
    FFToversize = 1;
end

%Check grid if passed
if ~isempty(samplingGrid)
    if size(samplingGrid,2) > size(samplingGrid,1)
        samplingGrid = samplingGrid';
    end
    
    if size(samplingGrid,2) ~= 3
        warning('No sampling weights passed. Using sampling grid defined in SOFA object and AKsht performing SHT with pseudo-inverse SH-matrix.')
        samplingGrid = [];
    end
end

%Check SOFA object
if ~strcmp(SOFAobj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR')
    error('Function only valid for SOFA convention "SimpleFreeFieldHRIR".');
end

%% If sampling grid passed - use sofia_stc

if ~isempty(samplingGrid)
        
    %Convert az and el to rad again
    samplingGrid(:,1:2) = samplingGrid(:,1:2) * pi / 180;
    
    %Write timeData struct needed for sofia_fdt
    timeDataL.impulseResponses  = squeeze(SOFAobj.Data.IR(:,1,:));
    timeDataL.FS                = SOFAobj.Data.SamplingRate;
    timeDataL.radius            = abs(SOFAobj.ReceiverPosition(3));
    timeDataL.averageAirTemp    = 29.3100; %Number according to miro object HRIR_L2702
    
    timeDataR.impulseResponses  = squeeze(SOFAobj.Data.IR(:,2,:));
    timeDataR.FS                = SOFAobj.Data.SamplingRate;
    timeDataR.radius            = abs(SOFAobj.ReceiverPosition(3));
    timeDataR.averageAirTemp    = 29.3100; %Number according to miro object HRIR_L2702
    
    %Transform to Fourier domain with regular FFT (sofia_fdt)
    fftDataL = sofia_fdt(timeDataL,FFToversize);
    [fftDataR, ~, f] = sofia_fdt(timeDataR,FFToversize);
    
    %Transform to SH-domain with spherical Fourier transform (sofia_stc)
    %Get SH-coefficients for left and right channel
    Hl_nm = sofia_stc(N,fftDataL,samplingGrid);
    Hr_nm = sofia_stc(N,fftDataR,samplingGrid);

    %Write output struct
    HRIRs_sfd.Hl_nm         = Hl_nm;
    HRIRs_sfd.Hr_nm         = Hr_nm;
    HRIRs_sfd.f             = f;
    HRIRs_sfd.N             = N;
    HRIRs_sfd.FFToversize   = FFToversize;
    
    disp('Transformation done with SOFiA toolbox');
end

%% If no sampling grid passed - Use AKsht
if isempty(samplingGrid)
    
    %Get samplingGrid from SOFA object
    samplingGrid = SOFAobj.SourcePosition(:,1:2);
    
    %Transform samplingGrid back to AKsht coordinate system
    samplingGrid(:,2) = 90-samplingGrid(:,2);

    %Get IRs needed for AKsht
    %Arrays need to be flipped because AKsht needs vectors with [M x N]
    %(columns = channels)
    irL = squeeze(SOFAobj.Data.IR(:,1,:))'; 
    irR = squeeze(SOFAobj.Data.IR(:,2,:))';
        
    %Zero pad to get length according to FFToversize
    if FFToversize > 1
        NFFT = 2^nextpow2(size(irL,1));
        NFFT = NFFT*FFToversize;
        irL = [irL; zeros(NFFT-size(irL,1),size(irL,2))];
        irR = [irR; zeros(NFFT-size(irR,1),size(irR,2))];
    else
        if size(irL,1) == 2^nextpow2(size(irL,1))
            NFFT = size(irL,1);
        else
            NFFT = 2^nextpow2(size(irL,1));
            irL = [irL; zeros(NFFT-size(irL,1),size(irL,2))];
            irR = [irR; zeros(NFFT-size(irR,1),size(irR,2))];
        end
    end
    
    %Transform to SH-domain with spherical Fourier transform (AKsht)
    %Get SH-coefficients for left and right channel
    Hl_nm = AKsht(irL,true,samplingGrid,N,'complex',SOFAobj.Data.SamplingRate);
    [Hr_nm, f] = AKsht(irR,true,samplingGrid,N,'complex',SOFAobj.Data.SamplingRate);
    
    %Write output struct
    HRIRs_sfd.Hl_nm         = Hl_nm;
    HRIRs_sfd.Hr_nm         = Hr_nm;
    HRIRs_sfd.f             = f';
    HRIRs_sfd.N             = N;
    HRIRs_sfd.FFToversize   = FFToversize;
    
    disp('Transformation done with AKTools');
end

end

