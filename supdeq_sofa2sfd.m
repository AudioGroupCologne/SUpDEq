%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function HRIRs_sfd = supdeq_sofa2sfd(SOFAobj, N, samplingGrid, FFToversize, transformCore, tikhEps)
%
% This function transformes spherical HRIR measurements in SOFA format
% ('SimpleFreeFieldHRIR' convention) to the SH-domain by spherical Fourier 
% transform and writes a output struct suitable for the SUpDEq toolbox.
%
% Output:
% HRIRs_sfd     - Struct with Spherical Harmonic coefficients 
%                 (SH-coefficients) for the left (Hl_nm) and right (Hr_nm) 
%                 channel/ear, absolute frequency scale f, 
%                 transform order N, FFToversize, and sourceDistance in m.
%
% Input:
% SOFAobj       - Measurement data in form of a SOFA object
% N             - Transform order N
% samplingGrid  - Optional input [default = nan]
%                 If a spatial sampling grid is passed, it should be a Qx3 
%                 matrix where the first column holds the azimuth, the 
%                 second the elevation, and the third the sampling weights.
%                 In this case, the SH transform is performed with weights.
%                 If no samplingGrid (or []) is passed, the sampling grid
%                 according to the SOFA object is used for a spherical
%                 Fourier transform without sampling weights. In this case,
%                 the transformCore 'ak' is used by default.
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
% transformCore - String to define method to be used for the spherical 
%                 Fourier transform. 
%                 'sofia - sofia_itc from SOFiA toolbox
%                 'ak'   - AKisht from AKtools 
%                 The results are exactly the same, but AKisht is faster 
%                 with big sampling grids
%                 Default: 'sofia'
% tikhEps       - Define epsilon of Tikhonov regularization if
%                 regularization should be applied
%                 Applying the Tikhonov regularization will always result 
%                 in a least-square fit solution for the SH transform. 
%                 Variable 'transformCore' is neglected when 'tikhEps' is 
%                 defined as the regularized least-square spherical Fourier 
%                 transform is applied directly without any third party 
%                 toolbox. Depending on the sampling grids, weights are
%                 applied or not.
%                 Default: 0 (no Tikhonov regularization)
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
% R. Duraiswami, D. N. Zotkin, and N. A. Gumerov, ?Interpolation and range 
% extrapolation of HRTFs,? in Proceedings of the IEEE International 
% Conference on Acoustics, Speech, and Signal Processing, 2004, pp. IV45?IV48.
%
% (C) 2018/2019 by JMA, Johannes M. Arend
%               TH Köln - University of Applied Sciences
%               Institute of Communications Engineering
%               Department of Acoustics and Audio Signal Processing

function HRIRs_sfd = supdeq_sofa2sfd(SOFAobj, N, samplingGrid, FFToversize, transformCore, tikhEps)

if nargin < 3 || isempty(samplingGrid)
    samplingGrid = [];
end

if nargin < 4 || isempty(FFToversize)
    FFToversize = 1;
end

if nargin < 5 || isempty(transformCore)
    transformCore = 'sofia';
end

if nargin < 6 || isempty(tikhEps)
    tikhEps = 0;
end

%Check FFToversize parameter
if FFToversize < 1
    warning('FFToversize needs to >= 1. Set to default = 1');
    FFToversize = 1;
end

%Check grid if passed
weightsPassed = true;
if ~isempty(samplingGrid)
    if size(samplingGrid,2) > size(samplingGrid,1)
        samplingGrid = samplingGrid';
        warning('Assuming grid was passed in wrong dimensions - Grid flipped');
    end
    
    if size(samplingGrid,2) ~= 3
        weightsPassed = false;
    end
end

%Set transformCore to 'ak' if no weights passed but sofia is choosen
if ~weightsPassed && strcmp(transformCore,'sofia')
    disp('Changed transformCore to ak as no sampling weights were passed!');
    transformCore = 'ak';
end
    
%Check SOFA object
if ~strcmp(SOFAobj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR')
    error('Function only valid for SOFA convention "SimpleFreeFieldHRIR".');
end

%% If sampling grid with weights passed and transformCore sofia

if tikhEps == 0 %Without Tikhonov Regularization
    
    if ~isempty(samplingGrid) && weightsPassed && strcmp(transformCore,'sofia')

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
        HRIRs_sfd.Hl_nm             = Hl_nm;
        HRIRs_sfd.Hr_nm             = Hr_nm;
        HRIRs_sfd.f                 = f;
        HRIRs_sfd.N                 = N;
        HRIRs_sfd.FFToversize       = FFToversize;
        HRIRs_sfd.sourceDistance    = SOFAobj.SourcePosition(1,3);

        disp('Transformation done with SOFiA toolbox - Sampling weights passed');
    end
end

%% If sampling grid passed with or without weights and transformCore ak

if tikhEps == 0 %Without Tikhonov Regularization

    if ~isempty(samplingGrid) && strcmp(transformCore,'ak')

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
        HRIRs_sfd.sourceDistance    = SOFAobj.SourcePosition(1,3);

        if weightsPassed
            disp('Transformation done with AKTools - Sampling weights passed');
        else
            disp('Transformation done with AKTools - No sampling weights passed');
        end
    end
end

%% If no sampling grid passed - Use AKsht

if tikhEps == 0 %Without Tikhonov Regularization

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
        HRIRs_sfd.sourceDistance    = SOFAobj.SourcePosition(1,3);

        disp('Transformation done with AKTools - Applied sampling grid from SOFA object without weights');
    end
end

%% If tikhEps is not zero - Calculate SH transform with Tikhonov regularization

if tikhEps ~= 0
    
    sofaGrid = false;
    if isempty(samplingGrid) %If no sampling grid passed
        %Get samplingGrid from SOFA object
        samplingGrid = SOFAobj.SourcePosition(:,1:2);
        %Transform samplingGrid back to SUpDEq coordinate system
        samplingGrid(:,2) = 90-samplingGrid(:,2);
        %Set weightsPassed
        weightsPassed = false;
        sofaGrid = true;
    end
        
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

    %Write f
    f = linspace(0,SOFAobj.Data.SamplingRate/2,NFFT/2+1);

    %Get SH functions
    [Ynm,n] = AKsh(N,[],samplingGrid(:,1),samplingGrid(:,2));
    n = n';
    nSH = (N+1).^2;

    %Create diagonal matrix according to Duraiswami2004
    I = eye(nSH);
    D = diag(1 + n.*(n+1)) .* I;

    if weightsPassed %% For Least-Square SH transform with Tikhonov regularization and weights
        W = diag(samplingGrid(:,3)); %Weights
        YnmInvTik = (Ynm' * W * Ynm + tikhEps*D)^-1 * Ynm' * W;
    else %% For Least-Square SH transform with Tikhonov regularization
        YnmInvTik = (Ynm' * Ynm + tikhEps*D)^-1 * Ynm';
    end
    
    %Get SH-coefficients for left and right channel
    Hl_nm = YnmInvTik*HRTF_L;
    Hr_nm = YnmInvTik*HRTF_R;
    
    %Write output struct
    HRIRs_sfd.Hl_nm         = Hl_nm;
    HRIRs_sfd.Hr_nm         = Hr_nm;
    HRIRs_sfd.f             = f;
    HRIRs_sfd.N             = N;
    HRIRs_sfd.FFToversize   = FFToversize;
    HRIRs_sfd.sourceDistance    = SOFAobj.SourcePosition(1,3);
    HRIRs_sfd.tikhEps = tikhEps;
    
    if weightsPassed 
        disp('Transformation done with least-squares method with Tikhonov regularization - Sampling weights passed');
    else %No weights
        if sofaGrid
            disp('Transformation done with least-squares method with Tikhonov regularization - Applied sampling grid from SOFA object without weights');
        else
            disp('Transformation done with least-squares method with Tikhonov regularization - No sampling weights passed');
        end
    end
    
end
end

