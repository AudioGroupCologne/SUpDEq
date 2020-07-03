%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function HRIRs_sfd = supdeq_hrtf2sfd(HRTF_L, HRTF_R, N, samplingGrid, fs, transformCore, tikhEps)
%
% This function transformes HRTFs to the SH-domain by spherical Fourier transform.
%
% Output:
% HRIRs_sfd     - Struct with Spherical Harmonic coefficients 
%                 (SH-Coefficients) for the left (Hl_nm) and right (Hr_nm) 
%                 channel/ear, absolute frequency scale f, 
%                 transform order N, and FFToversize
%
% Input:
% HRTF_L/R      - HRTF for left and right ear as [n x m] samples, with n =
%                 number of channels / sampling points and m = FFT bins
% N             - Transform order N
% samplingGrid  - Spatial sampling grid. Can be a Qx3 matrix where the first 
%                 column holds the azimuth, the second the elevation, and 
%                 the third the sampling weights. In this case, the SOFiA 
%                 functions "sofia_fdt" and "sofia_stc" as well as "AKsht"
%                 can be used
%                 If no weights are passed (Qx2 matrix), the function "AKsht" 
%                 is used providing a least-squares solution of the spherical
%                 Fourier transform.
%                 Azimuth positions in degree. Can be scalar or vector of size(el)
%                 (0=front, 90=left, 180=back, 270=right)
%                 (0 points to positive x-axis, 90 to positive y-axis)
%                 Elevations in degree. Can be scalar or vector of size(az)
%                 (0=North Pole, 90=front, 180=South Pole)
%                 (0 points to positive z-axis, 180 to negative z-axis)
% fs            - Sampling rate. Only needed to write frequency vector or
%                 if AKsht is used.  
%                 Default: 48000
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
% Dependencies: SOFiA toolbox, AKtools
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

function HRIRs_sfd = supdeq_hrtf2sfd(HRTF_L, HRTF_R, N, samplingGrid, fs, transformCore, tikhEps)

if nargin < 5 || isempty(fs)
    fs = 48000;
end

if nargin < 6 || isempty(transformCore)
    transformCore = 'sofia';
end

if nargin < 7 || isempty(tikhEps)
    tikhEps = 0;
end

%Check passed grid for weights  
weightsPassed = true;
if size(samplingGrid,2) ~= 3
    %warning('No sampling weights passed. Using AKsht to perform SHT with pseudo-inverse SH-matrix')
    weightsPassed = false;
end

%% If sampling grid with weights passed - use sofia_stc or AKsht

if tikhEps == 0 %Without Tikhonov Regularization

    if weightsPassed

        if strcmp(transformCore,'sofia')

            %Convert az and el to rad again
            samplingGrid(:,1:2) = samplingGrid(:,1:2) * pi / 180;

            %Transform to SH-domain with spherical Fourier transform (sofia_stc)
            %Get SH-coefficients for left and right channel
            Hl_nm = sofia_stc(N,HRTF_L,samplingGrid);
            Hr_nm = sofia_stc(N,HRTF_R,samplingGrid);

            %Write output struct
            HRIRs_sfd.Hl_nm         = Hl_nm;
            HRIRs_sfd.Hr_nm         = Hr_nm;
            HRIRs_sfd.f             = linspace(0,fs/2,size(HRTF_L,2));
            HRIRs_sfd.N             = N;
            HRIRs_sfd.FFToversize   = 1;

            disp('Transformation done with SOFiA toolbox - Sampling weights passed');

        elseif strcmp(transformCore,'ak')

            %Transform to SH-domain with spherical Fourier transform (AKsht)
            %Get SH-coefficients for left and right channel
            Hl_nm = AKsht(HRTF_L.',false,samplingGrid,N,'complex',fs);
            Hr_nm = AKsht(HRTF_R.',false,samplingGrid,N,'complex',fs);

            %Write output struct
            HRIRs_sfd.Hl_nm         = Hl_nm;
            HRIRs_sfd.Hr_nm         = Hr_nm;
            HRIRs_sfd.f             = linspace(0,fs/2,size(HRTF_L,2));
            HRIRs_sfd.N             = N;
            HRIRs_sfd.FFToversize   = 1;

            disp('Transformation done with AKTools - Sampling weights passed');

        end
    end
end

%% If no weights passed - Use always AKsht
if tikhEps == 0 %Without Tikhonov Regularization
    if ~weightsPassed 
    
        %Transform to SH-domain with spherical Fourier transform (AKsht)
        %Get SH-coefficients for left and right channel
        Hl_nm = AKsht(HRTF_L.',false,samplingGrid,N,'complex',fs);
        Hr_nm = AKsht(HRTF_R.',false,samplingGrid,N,'complex',fs);

        %Write output struct
        HRIRs_sfd.Hl_nm         = Hl_nm;
        HRIRs_sfd.Hr_nm         = Hr_nm;
        HRIRs_sfd.f             = linspace(0,fs/2,size(HRTF_L,2));
        HRIRs_sfd.N             = N;
        HRIRs_sfd.FFToversize   = 1;

        disp('Transformation done with AKTools - No sampling weights passed');
    end
end

%% If tikhEps is not zero - Calculate SH transform with Tikhonov regularization

if tikhEps ~= 0
    
    %Get SH functions, weights omitted
    [Ynm,n] = AKsh(N,[],samplingGrid(:,1),samplingGrid(:,2));
    n = n';
    nSH = (N+1).^2;

    %Create diagonal matrix according to Duraiswami2004
    I = eye(nSH);
    D = diag(1 + n.*(n+1)) .* I;
    
    % Inverse SH matrix for Least-Square SH transform with Tikhonov regularization
    YnmInvTik = (Ynm' * Ynm + tikhEps*D)^-1 * Ynm';
    
    %Get SH-coefficients for left and right channel
    Hl_nm = YnmInvTik*HRTF_L;
    Hr_nm = YnmInvTik*HRTF_R;
    
    %Write output struct
    HRIRs_sfd.Hl_nm         = Hl_nm;
    HRIRs_sfd.Hr_nm         = Hr_nm;
    HRIRs_sfd.f             = linspace(0,fs/2,size(HRTF_L,2));
    HRIRs_sfd.N             = N;
    HRIRs_sfd.FFToversize   = 1;
    HRIRs_sfd.tikhEps       = tikhEps;
    
    disp('Transformation done with least-squares method with Tikhonov regularization');
    
end

end

