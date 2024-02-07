%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function [interpHRTFset, HRTF_interp_L, HRTF_interp_R] = supdeq_interpHRTF(HRTFset, ipSamplingGrid, ppMethod, ipMethod, mc, headRadius, tikhEps, limitMC, mcKnee, mcMinPhase, limFade) 
%
% This function performs HRTF interpolation of the input HRTFset according
% to the passed interpolation sampling grid. Various interpolation and
% time-alignment (pre-processing) methods can be applied. Furthermore,
% post-interpolation magnitude correction for time-aligned interpolation
% (MCA) can be applied to further improve the interpolation results. See 
% also supdeq_demo_MCA for some examples of how to use the function.
%
% Output:
% interpHRTFset         - Struct with the interpolated HRTF for the 
%                         left (HRTF_L) and right (HRTF_R) ear, absolute 
%                         frequency scale f, various interpolation parameters
%                         according to the applied method, and FFToversize
% HRTF_L/R_ip           - Interpolated complex HRTFs (single-sided spectrum)
%
% Input:
% HRTFset               - Struct with (sparse) HRTF dataset in frequency
%                         domain (single-sided complex spectra). Struct needs to
%                         provide HRTF_L/HRTF_R, the samplingGrid, Nmax,
%                         absolute frequency scale f, and FFToversize
% ipSamplingGrid        - Spatial sampling grid (Q x 2 or Q x 3 matrix), 
%                         defining the interpolation points, where the first 
%                         column holds the azimuth, the second the elevation 
%                         (both in degree), and optionally the third the 
%                         sampling weights.
%                         Azimuth in degree
%                         (0=front, 90=left, 180=back, 270=right)
%                         (0 points to positive x-axis, 90 to positive y-axis)
%                         Elevation in degree
%                         (0=North Pole, 90=front, 180=South Pole)
%                         (0 points to positive z-axis, 180 to negative z-axis)
% ppMethod              - String defining the pre/post-processing applied
%                         to the HRTFs before/after interpolation. Besides
%                         'MagPhase', all of the pre-processing methods are
%                         time-alignment approaches.
%                         'SUpDEq'          - SUpDEq method [3][5] 
%                         'SUpDEq_Lim'      - Modified SUpDeq method where the eqDataset is limited to frequencies above the spatial aliasing frequency fA [5]
%                         'SUpDEq_AP'       - Modified SUpDEq method applying rigid sphere allpass filter (only the phase components of the rigid sphere transfer functions)
%                         'SUpDEq_Lim_AP'   - Modified SUpDeq method where the eqDataset is limited as above (SUpDEq_Lim) and only the phase components are used (SUpDEq_AP)
%                         'PC'              - Phase-Correction (open sphere transfer functions) [4][5]
%                         'OBTA'            - Onset-Based Time-Alignment [1][2][5]
%                         'MagPhase'        - Split HRTF to magnitude and unwrapped phase - No time-alignment
%                         'None'            - Perform interpolation without any pre/post-processing 
%                          Default: 'SUpDEq'
% ipMethod              - String defining the ipMethod
%                         'SH'      - Spherical Harmonics interpolation [1-5]
%                         'NN'      - Natural Neighbor interpolation with Voronoi weights [6]
%                         'Bary'    - Barycentric interpolation [7]
%                         'SARITA'  - Spherical Array Interpolation by Time Alignment [8]; employs Natural Neighbor interpolation with Voronoi weights
%                          Default: 'SH'
% mc                     - Define maximum boost in dB if magnitude correction according to [9] should be applied to interpolated HRTFs in a further postprocessing step
%                          Set to nan if magnitude correction should not be applied
%                          Set to inf if no limiting should be applied, i.e., unlimited boost
%                          Default: inf
% headRadius             - Head radius in m. Required for time-alignment methods SUpDEq [3][5] and PC [4][5]
%                          Default: 0.0875
% tikhEps                - Define epsilon of Tikhonov regularization if regularization should be applied
%                          Applying the Tikhonov regularization will always result in a least-square fit 
%                          solution for the SH transform. Only relevant if ipMethod is 'SH'.
%                          Default: 0 (no Tikhonov regularization)
% limitMC                - Boolean to set magnitude correction filters to 0 dB below spatial aliasing frequency fA, 
%                          assuming that (SH) interpolation below fA is physically correct. 
%                          Default: true (correction filter operate only above fA)
%                          Only applies if 'mc' is not nan
% mcKnee                 - Knee in dB of the soft-limiting applied to the correction filters. Soft-limiting is applied according to maximum boost 'mc' 
%                          Only applies if 'mc' is not nan or inf
%                          Default: 0 dB (no knee)
% mcMinPhase             - Boolean to design correction filters as minimum-phase filters. If set to false, correction filters are zero-phase filters
%                          Default: true (minimum-phase filters) as this provides a better performance in terms of ITD errors in combination with classical SUpDEq processing
%                          When applying other preprocessing methods, using zero-phase filters can be more appropriate. 
% limFade                - String to define whether to apply a linear fade upwards from aliasing frequency fA to fAt = fA + 1/3 Oct, or downwards, i.e., from fA to fAt = fA - 1/3 Oct
%                          'fadeUp' - upwards
%                          'fadeDown' - downwards
%                          Default: 'fadeDown'
% ILDComp                - Enable/disable the ILD error compensation for the magnitude correction. 
%                          Only applies if 'mc' is not nan.
%                          Default: true
%
% Dependencies: SUpDEq toolbox, AKtools, SOFiA toolbox, SFS toolbox, TriangleRayIntersection
%
% References: (Pre-Processing methods)
% Onset-Based Time-Alignment (OBTA)
% [1] M. J. Evans, J. A. S. Angus, and A. I. Tew, "Analyzing head-related transfer function measurements 
% using surface spherical harmonics,"
% J. Acoust. Soc. Am., vol. 104, no. 4, pp. 2400-2411, 1998,
% doi: https://doi.org/10.1121/1.423749.
% [2] F. Brinkmann and S. Weinzierl, "Comparison of head-related transfer functions pre-processing 
% techniques for spherical harmonics decomposition,"
% in Proceedings of the AES Conference on Audio for Virtual and Augmented Reality, 2018, pp. 1-10.
%
% Spatial Upsampling by Directional Equalization (SUpDEq)
% [3] C. Pörschmann, J. M. Arend, and F. Brinkmann, 
% "Directional Equalization of Sparse Head-Related Transfer Function Sets for Spatial Upsampling," 
% IEEE/ACM Trans. Audio, Speech, Lang. Process., vol. 27, no. 6, pp. 1060-1071, 2019,
% doi: https://doi.org/10.1109/TASLP.2019.2908057.
%
% Phase-Correction (PC)
% [4] Z. Ben-Hur, D. Lou Alon, R. Mehra, and B. Rafaely, 
% "Efficient Representation and Sparse Sampling of Head-Related Transfer Functions 
% Using Phase-Correction Based on Ear Alignment," 
% IEEE/ACM Trans. Audio, Speech, Lang. Process., vol. 27, no. 12, pp. 2249-2262, 2019,
% doi: https://doi.org/10.1109/TASLP.2019.2945479.
%
% Comparison of time-alignment methods & perceptual evaluation
% [5] J. M. Arend, F. Brinkmann, and C. Pörschmann, “Assessing Spherical 
% Harmonics Interpolation of Time-Aligned Head-Related Transfer Functions,” 
% J. Audio Eng. Soc., vol. 69, no. 1/2, pp. 104–117, 2021, 
% doi: https://doi.org/10.17743/jaes.2020.0070.
% 
% Spherical harmonics interpolation
% [1-5]
%
% Natural-Neighbor Interpolation
% [6] R. Sibson, “A brief description of natural neighbor interpolation,” 
% in Interpreting Multivariate Data, V. Barnett, Chichester, England: John Wiley & Sons, 1981, pp. 21–36.
% Good read about NN interpolation: https://github.com/gwlucastrig/Tinfour/wiki/Introduction-to-Natural-Neighbor-Interpolation
% 
% Barycentric Coordinates / Interpolation
% [7] P. Shirley and S. Marschner, Fundamentals of Computer Graphics, 3rd ed. Boca Raton, FL: Taylor & Francis, 2009.
% Good read about Barycentric interpolation: https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates
% 
% [8] Spherical Array Interpolation by Time Alignment (SARITA)
% T. Lübeck, J. M. Arend, and C. Pörschmann, “Spatial Upsampling of Sparse Spherical Microphone Array Signals,” 
% IEEE/ACM Trans. Audio Speech Lang. Proc., vol. 31, pp. 1163–1174, 2023,
% doi: https://doi.org/10.1109/TASLP.2023.3250841.
%
% Magnitude-Corrected and Time-Aligned Interpolation (MCA)
% [9] J. M. Arend, C. Pörschmann, S. Weinzierl, and F. Brinkmann, 
% "Magnitude-Corrected and Time-Aligned Interpolation of Head-Related Transfer Functions,” 
% IEEE/ACM Trans. Audio Speech Lang. Proc., vol. 31, pp. 3783–3799, 2023,
% doi: https://doi.org/10.1109/TASLP.2023.3313908.
%
%
% (C) 2020-2023 by JMA, Johannes M. Arend
%               TU Berlin, Audio Communication Group
%               TH Köln, Institute of Communications Engineering

function [interpHRTFset, HRTF_L_ip, HRTF_R_ip] = supdeq_interpHRTF(HRTFset, ipSamplingGrid, ppMethod, ipMethod, mc, headRadius, tikhEps, limitMC, mcKnee, mcMinPhase, limFade, ILDComp) 


if nargin < 3 || isempty(ppMethod)
    ppMethod = 'SUpDEq';
end

if nargin < 4 || isempty(ipMethod)
    ipMethod = 'SH';
end

if nargin < 5 || isempty(mc)
    mc = inf;
end

if nargin < 6 || isempty(headRadius)
    headRadius = 0.0875;
end

if nargin < 7 || isempty(tikhEps)
    tikhEps = 0;
end

if nargin < 8 || isempty(limitMC)
    limitMC = true;
end

if nargin < 9 || isempty(mcKnee)
    mcKnee = 0;
end

if nargin < 10 || isempty(mcMinPhase)
    mcMinPhase = true;
end

if nargin < 11 || isempty(limFade)
    limFade = 'fadeDown';
end

if nargin < 12 || isempty(ILDComp)
    ILDComp = true;
end
   
%% Get some required variables

fs = HRTFset.f(end)*2;
f = HRTFset.f;
NFFT = length(HRTFset.f)*2-2;
samplingGrid = HRTFset.samplingGrid;
HRTF_L = HRTFset.HRTF_L;
HRTF_R = HRTFset.HRTF_R;
c = 343;

%Get HRIRs - Required for threshold detection and mc
%Get mirror spectrum
HRIR_L = AKsingle2bothSidedSpectrum(HRTF_L.');
HRIR_R = AKsingle2bothSidedSpectrum(HRTF_R.');
%Transform back to time domain
HRIR_L = real(ifft(HRIR_L));
HRIR_R = real(ifft(HRIR_R));
hrirLength = size(HRIR_L,1)/HRTFset.FFToversize;

%% Perform preprocessing

switch ppMethod
    case 'SUpDEq'
        disp('Preprocessing: SUpDEq')
        %Apply SUpDEq - Equalize HRTFs

        %Get eqDataset at N = 44
        eqDataset = supdeq_getEqDataset(44,2*headRadius,NFFT,fs);

        %Equalization 
        [eqTF_L,eqTF_R] = supdeq_getEqTF(eqDataset,samplingGrid,'DEG',2,'ak',false);
        pHRTF_L = HRTF_L./eqTF_L;
        pHRTF_R = HRTF_R./eqTF_R;
        
    case 'SUpDEq_Lim'
        disp('Preprocessing: SUpDEq_Lim')
        %Apply SUpDEq - Equalize HRTFs

        %Get eqDataset at N = 44
        eqDataset = supdeq_getEqDataset(44,2*headRadius,NFFT,fs);
        
        %Limit eqDataset (Small improvement explained in [5] leading to better
        %results below the spatial aliasing frequency fA)
        eqDataset = supdeq_limitEqDataset(eqDataset,HRTFset.Nmax,headRadius);

        %Equalization 
        [eqTF_L,eqTF_R] = supdeq_getEqTF(eqDataset,samplingGrid,'DEG',2,'ak',false);
        pHRTF_L = HRTF_L./eqTF_L;
        pHRTF_R = HRTF_R./eqTF_R;

    case 'SUpDEq_AP'
        disp('Preprocessing: SUpDEq_AP')
        %Apply SUpDEq_AP - Equalize HRTFs with rigid sphere allpass functions

        %Get eqDataset at N = 44
        eqDataset = supdeq_getEqDataset(44,2*headRadius,NFFT,fs);

        %Equalization 
        [eqTF_L,eqTF_R] = supdeq_getEqTF(eqDataset,samplingGrid,'DEG',2,'ak',true);
        pHRTF_L = HRTF_L./eqTF_L;
        pHRTF_R = HRTF_R./eqTF_R;
        
    case 'SUpDEq_Lim_AP'
        disp('Preprocessing: SUpDEq_Lim_AP')
        %Apply SUpDEq_AP - Equalize HRTFs with limited rigid sphere allpass functions

        %Get eqDataset at N = 44
        eqDataset = supdeq_getEqDataset(44,2*headRadius,NFFT,fs);
        
        %Limit eqDataset (Small improvement explained in [5] leading to better
        %results below the spatial aliasing frequency fA)
        eqDataset = supdeq_limitEqDataset(eqDataset,HRTFset.Nmax,headRadius);

        %Equalization 
        [eqTF_L,eqTF_R] = supdeq_getEqTF(eqDataset,samplingGrid,'DEG',2,'ak',true);
        pHRTF_L = HRTF_L./eqTF_L;
        pHRTF_R = HRTF_R./eqTF_R;

    case 'PC'
        disp('Preprocessing: PC')
        %Apply Phase-Correction - Equalize HRTFs with open sphere allpass functions

        %Get wavenumber k
        k  = 2*pi*f/c; k = k.';

        %Transform sampling grid to radiant
        sg = samplingGrid * pi / 180;

        %Get phase correction term for left/right ear according to Eq. 13 in [4]
        cosThetaL = cos(sg(:,2)')*cos(pi/2) + sin(sg(:,2)')*sin(pi/2) .* cos(sg(:,1)'-pi/2); %Left ear with -pi/2
        phaseCorL = exp(-1j*headRadius * k .* cosThetaL);
        cosThetaR = cos(sg(:,2)')*cos(pi/2) + sin(sg(:,2)')*sin(pi/2) .* cos(sg(:,1)'+pi/2); %Right ear with +pi/2
        phaseCorR = exp(-1j*headRadius * k .* cosThetaR);

        %Apply to HRTFs
        pHRTF_L = HRTF_L .* phaseCorL.'; %Just transpose, no conjugate complex
        pHRTF_R = HRTF_R .* phaseCorR.';

    case 'OBTA' 
        disp('Preprocessing: OBTA')
        %Apply Onset-Based Time-Alignment with fractional delay filters

        %Not compatible with SARITA interpolation. Intercepted here.
        if strcmp(ipMethod,'SARITA')
            warning('OBTA not compatible with SARITA interpolation. Set to ''None'' for preprocessing')
            ppMethod = 'None';

            disp('Preprocessing: None')
            %Copy HRTF to pHRTF_L as in case 'None'
            pHRTF_L = HRTF_L;
            pHRTF_R = HRTF_R;
        
        else

            %Estimate TOA on low-passed and 10 times upsampled HRIRs according to [2]
            toa_L = AKonsetDetect(HRIR_L, 10, -20, 'rel', [3e3 fs]);
            toa_R = AKonsetDetect(HRIR_R, 10, -20, 'rel', [3e3 fs]);
    
            %Time-Align HRIRs (remove TOAs) with fractional delay
            pHRIR_L = AKfractionalDelayCyclic(HRIR_L, -toa_L);
            pHRIR_R = AKfractionalDelayCyclic(HRIR_R, -toa_R);
    
            %Transform to Fourier domain with same FFToversize as defined in struct
            pHRTF_L = AKboth2singleSidedSpectrum (fft(pHRIR_L)).';
            pHRTF_R = AKboth2singleSidedSpectrum (fft(pHRIR_R)).';
        
        end

    case 'MagPhase'
        disp('Preprocessing: MagPhase')

        %Not compatible with SARITA interpolation. Intercepted here.
        if strcmp(ipMethod,'SARITA')
            warning('MagPhase not compatible with SARITA interpolation. Set to ''None'' for preprocessing')
            ppMethod = 'None';

            disp('Preprocessing: None')
            %Copy HRTF to pHRTF_L as in case 'None'
            pHRTF_L = HRTF_L;
            pHRTF_R = HRTF_R;

        else
        
            pHRTF_L = abs(HRTF_L);
            pHRTF_R = abs(HRTF_R);
            pHRTF_L_ph = unwrap(angle(HRTF_L).').';
            pHRTF_R_ph = unwrap(angle(HRTF_R).').';

        end
        
    case 'None' 
        disp('Preprocessing: None')

        %Simply copy HRTF to pHRTF_L - Easier for further processing
        pHRTF_L = HRTF_L;
        pHRTF_R = HRTF_R;
end

%Save intermediate result (pre-processed HRTFs) in output struct
interpHRTFset.p.pHRTF_L = pHRTF_L;
interpHRTFset.p.pHRTF_R = pHRTF_R;
if strcmp(ppMethod,'OBTA')
    interpHRTFset.p.toa_L = toa_L;
    interpHRTFset.p.toa_R = toa_R;
end
if strcmp(ppMethod,'MagPhase')
    interpHRTFset.p.pHRTF_L_ph = pHRTF_L_ph;
    interpHRTFset.p.pHRTF_R_ph = pHRTF_R_ph;
end

%% Perform interpolation

switch ipMethod
    case 'SH'
        disp('HRTF Interpolation: SH')
        
        %Transform preprocessed HRTFs to SH domain
        %Tikhonov regularization can be applied if weights in sampling grid are erased and tikhEps ~= 0
        pHRTF_L_nm = AKsht(pHRTF_L.',false,samplingGrid,HRTFset.Nmax,'complex',fs,false,'complex',tikhEps);
        pHRTF_R_nm = AKsht(pHRTF_R.',false,samplingGrid,HRTFset.Nmax,'complex',fs,false,'complex',tikhEps);
        
        %Apply SH interpolation to ipSamplingGrid
        pHRTF_L_ip = AKisht(pHRTF_L_nm,false,ipSamplingGrid(:,1:2),'complex').';
        pHRTF_R_ip = AKisht(pHRTF_R_nm,false,ipSamplingGrid(:,1:2),'complex').';
        
        %Also interpolate TOAs if ppMethod = 'OBTA'
        if strcmp(ppMethod,'OBTA')
            
            %Transform TOAs to SH domain
            %Tikhonov regularization can be applied if weights in sampling grid are erased and tikhEps ~= 0
            toa_L_nm = AKsht(toa_L,false,samplingGrid,HRTFset.Nmax,'complex',fs,false,'real',tikhEps);
            toa_R_nm = AKsht(toa_R,false,samplingGrid,HRTFset.Nmax,'complex',fs,false,'real',tikhEps);

            %Apply SH interpolation to ipSamplingGrid
            toa_L_ip = AKisht(toa_L_nm, false, ipSamplingGrid(:,1:2), 'complex', false, false, 'real');
            toa_R_ip = AKisht(toa_R_nm, false, ipSamplingGrid(:,1:2), 'complex', false, false, 'real');
            
        end
        
        %Also interpolate Phase if ppMethod = 'MagPhase'
        if strcmp(ppMethod,'MagPhase')
            
            %Transform Phase to SH domain
            %Tikhonov regularization can be applied if weights in sampling grid are erased and tikhEps ~= 0
            pHRTF_L_ph_nm = AKsht(pHRTF_L_ph.',false,samplingGrid,HRTFset.Nmax,'complex',fs,false,'complex',tikhEps);
            pHRTF_R_ph_nm = AKsht(pHRTF_R_ph.',false,samplingGrid,HRTFset.Nmax,'complex',fs,false,'complex',tikhEps);

            %Apply SH interpolation to ipSamplingGrid
            pHRTF_L_ph_ip = AKisht(pHRTF_L_ph_nm,false,ipSamplingGrid(:,1:2),'complex').';
            pHRTF_R_ph_ip = AKisht(pHRTF_R_ph_nm,false,ipSamplingGrid(:,1:2),'complex').';
  
        end
        
    case 'NN'
        disp('HRTF Interpolation: NN')
        
        %Transform grids to cartesian coordinates
        [samplingGrid_cart(:,1), samplingGrid_cart(:,2), samplingGrid_cart(:,3)] = sph2cart(samplingGrid(:,1)/360*2*pi, pi/2-samplingGrid(:,2)/360*2*pi,1);
        [ipSamplingGrid_cart(:,1), ipSamplingGrid_cart(:,2), ipSamplingGrid_cart(:,3)] = sph2cart(ipSamplingGrid(:,1)/360*2*pi, pi/2-ipSamplingGrid(:,2)/360*2*pi,1);
        
        %NN interpolate preprocessed HRTFs to ipSamplingGrid
        pHRTF_L_ip = zeros(length(ipSamplingGrid),length(f));
        pHRTF_R_ip = zeros(length(ipSamplingGrid),length(f));
        
        for nPoint = 1:size(ipSamplingGrid,1)

            [idx_voronoi,w_voronoi] = findvoronoi(samplingGrid_cart,ipSamplingGrid_cart(nPoint,:));

            for n = 1:length(idx_voronoi) 
                pHRTF_L_ip(nPoint,:) = pHRTF_L_ip(nPoint,:) + w_voronoi(n)*pHRTF_L(idx_voronoi(n),:);  
                pHRTF_R_ip(nPoint,:) = pHRTF_R_ip(nPoint,:) + w_voronoi(n)*pHRTF_R(idx_voronoi(n),:);
            end
        end
        
        %Also interpolate TOAs if ppMethod = 'OBTA'
        if strcmp(ppMethod,'OBTA')
            
            toa_L_ip = zeros(1,length(ipSamplingGrid));
            toa_R_ip = zeros(1,length(ipSamplingGrid));
            
            for nPoint = 1:size(ipSamplingGrid,1)

                [idx_voronoi,w_voronoi] = findvoronoi(samplingGrid_cart,ipSamplingGrid_cart(nPoint,:));

                for n = 1:length(idx_voronoi) 
                    toa_L_ip(1,nPoint) = toa_L_ip(1,nPoint) + w_voronoi(n)*toa_L(1,idx_voronoi(n));  
                    toa_R_ip(1,nPoint) = toa_R_ip(1,nPoint) + w_voronoi(n)*toa_R(1,idx_voronoi(n));
                end
            end
        end
        
        %Also interpolate Phase if ppMethod = 'MagPhase'
        if strcmp(ppMethod,'MagPhase')  
            
            pHRTF_L_ph_ip = zeros(length(ipSamplingGrid),length(f));
            pHRTF_R_ph_ip = zeros(length(ipSamplingGrid),length(f));

            for nPoint = 1:size(ipSamplingGrid,1)

                [idx_voronoi,w_voronoi] = findvoronoi(samplingGrid_cart,ipSamplingGrid_cart(nPoint,:));

                for n = 1:length(idx_voronoi) 
                    pHRTF_L_ph_ip(nPoint,:) = pHRTF_L_ph_ip(nPoint,:) + w_voronoi(n)*pHRTF_L_ph(idx_voronoi(n),:);  
                    pHRTF_R_ph_ip(nPoint,:) = pHRTF_R_ph_ip(nPoint,:) + w_voronoi(n)*pHRTF_R_ph(idx_voronoi(n),:);
                end
            end 
        end

    case 'Bary'
        disp('HRTF Interpolation: Bary')
        
        %Transform grids to cartesian coordinates
        [samplingGrid_cart(:,1), samplingGrid_cart(:,2), samplingGrid_cart(:,3)] = sph2cart(samplingGrid(:,1)/360*2*pi, pi/2-samplingGrid(:,2)/360*2*pi,1);
        [ipSamplingGrid_cart(:,1), ipSamplingGrid_cart(:,2), ipSamplingGrid_cart(:,3)] = sph2cart(ipSamplingGrid(:,1)/360*2*pi, pi/2-ipSamplingGrid(:,2)/360*2*pi,1);
        
        %Get simplified convex hull (could also be done after delaunay-triangulation)
        convHullIdx = convhull(samplingGrid_cart(:,1), samplingGrid_cart(:,2), samplingGrid_cart(:,3), 'simplify', true);
        
        %Get HRTFs / directions assigned to convex hull triangles
        convHullHRTF_A = samplingGrid_cart(convHullIdx(:,1),:);
        convHullHRTF_B = samplingGrid_cart(convHullIdx(:,2),:);
        convHullHRTF_C = samplingGrid_cart(convHullIdx(:,3),:);
        
        %Interpolate with barycentric weights
        %Right triangle of convex hull picked by ray-triangle intersection test (intersectLinePolygon3d from geom3D toolbox should work too)
        %Function returns barycentric weights / coordinates u and v of the intersection point -> w = 1-u-v
        %P = w*A + u*B + v*C
        orig = [0 0 0];
        pHRTF_L_ip = zeros(length(ipSamplingGrid),length(f));
        pHRTF_R_ip = zeros(length(ipSamplingGrid),length(f));
        
        for nPoint = 1:size(ipSamplingGrid,1)

            [intersect, ~, u, v, ~] = TriangleRayIntersection(orig, ipSamplingGrid_cart(nPoint,:), convHullHRTF_A, convHullHRTF_B, convHullHRTF_C, 'border', 'inclusive');
            intersectIdx = find(intersect, 1, 'first'); %Evaluate first intersection point
            u = u(intersectIdx); v = v(intersectIdx); w = 1-u-v;
            bw = [w u v]; %Barycentric weights
            idx_bw = convHullIdx(intersectIdx,:); %Indizes of HRTFs of convex hull / triangles

            for n = 1:3 %Always 3 because of triangulization
                pHRTF_L_ip(nPoint,:) = pHRTF_L_ip(nPoint,:) + bw(n)*pHRTF_L(idx_bw(n),:);  
                pHRTF_R_ip(nPoint,:) = pHRTF_R_ip(nPoint,:) + bw(n)*pHRTF_R(idx_bw(n),:);
            end
        end
        
        %Also interpolate TOAs if ppMethod = 'OBTA'
        if strcmp(ppMethod,'OBTA')
            
            toa_L_ip = zeros(1,length(ipSamplingGrid));
            toa_R_ip = zeros(1,length(ipSamplingGrid));
            
            for nPoint = 1:size(ipSamplingGrid,1)

                [intersect, ~, u, v, ~] = TriangleRayIntersection(orig, ipSamplingGrid_cart(nPoint,:), convHullHRTF_A, convHullHRTF_B, convHullHRTF_C, 'border', 'inclusive');
                intersectIdx = find(intersect, 1, 'first'); %Evaluate first intersection point
                u = u(intersectIdx); v = v(intersectIdx); w = 1-u-v;
                bw = [w u v]; %Barycentric weights
                idx_bw = convHullIdx(intersectIdx,:); %Indizes of HRTFs of convex hull / triangles

                for n = 1:3 %Always 3 because of triangulization
                    toa_L_ip(1,nPoint) = toa_L_ip(1,nPoint) + bw(n)*toa_L(1,idx_bw(n));  
                    toa_R_ip(1,nPoint) = toa_R_ip(1,nPoint) + bw(n)*toa_R(1,idx_bw(n));
                end
            end
        end
        
        %Also interpolate Phase if ppMethod = 'MagPhase'
        if strcmp(ppMethod,'MagPhase')  
            
            pHRTF_L_ph_ip = zeros(length(ipSamplingGrid),length(f));
            pHRTF_R_ph_ip = zeros(length(ipSamplingGrid),length(f));

            for nPoint = 1:size(ipSamplingGrid,1)

                [intersect, ~, u, v, ~] = TriangleRayIntersection(orig, ipSamplingGrid_cart(nPoint,:), convHullHRTF_A, convHullHRTF_B, convHullHRTF_C, 'border', 'inclusive');
                intersectIdx = find(intersect, 1, 'first'); %Evaluate first intersection point
                u = u(intersectIdx); v = v(intersectIdx); w = 1-u-v;
                bw = [w u v]; %Barycentric weights
                idx_bw = convHullIdx(intersectIdx,:); %Indizes of HRTFs of convex hull / triangles

                for n = 1:3 %Always 3 because of triangulization
                    pHRTF_L_ph_ip(nPoint,:) = pHRTF_L_ph_ip(nPoint,:) + bw(n)*pHRTF_L_ph(idx_bw(n),:);  
                    pHRTF_R_ph_ip(nPoint,:) = pHRTF_R_ph_ip(nPoint,:) + bw(n)*pHRTF_R_ph(idx_bw(n),:);
                end
            end
        end

    case 'SARITA'
        disp('HRTF Interpolation: SARITA')     
        %SARITA works in time domain - Get pHRIRs from pHRTFs
        %Get mirror spectrum
        pHRIR_L = AKsingle2bothSidedSpectrum(pHRTF_L.');
        pHRIR_R = AKsingle2bothSidedSpectrum(pHRTF_R.');
        %Transform to time domain
        pHRIR_L = real(ifft(pHRIR_L));
        pHRIR_R = real(ifft(pHRIR_R));
    
        %Zero padding so that the output of SARITA is for sure containing all relevant samples
        pHRIR_L_tmp=[pHRIR_L.', zeros(size(pHRIR_L.'))]; 
        pHRIR_R_tmp=[pHRIR_R.', zeros(size(pHRIR_R.'))];
        
        %Perform interpolation using SARITA function
        [pHRIR_L_ip_tmp, ~] = Sarita_upsampling(pHRIR_L_tmp, samplingGrid(:,1:2)/(360/2/pi), ipSamplingGrid(:,1:2)/(360/2/pi), headRadius*1.2, 'frame_length', floor(size(pHRIR_L,1)/2), 'frame_overlap', 1);
        [pHRIR_R_ip_tmp, ~] = Sarita_upsampling(pHRIR_R_tmp, samplingGrid(:,1:2)/(360/2/pi), ipSamplingGrid(:,1:2)/(360/2/pi), headRadius*1.2, 'frame_length', floor(size(pHRIR_R,1)/2), 'frame_overlap', 1);
        
        %Remove zeros again
        pHRIR_L_ip=pHRIR_L_ip_tmp(:,1:size(pHRIR_L,1)).';  
        pHRIR_R_ip=pHRIR_R_ip_tmp(:,1:size(pHRIR_L,1)).';

        %Transform to frequency domain
        pHRTF_L_ip=fft(pHRIR_L_ip, [], 1);
        pHRTF_L_ip=pHRTF_L_ip(1:round(end/2)+1,:).';
        pHRTF_R_ip=fft(pHRIR_R_ip, [], 1);
        pHRTF_R_ip=pHRTF_R_ip(1:round(end/2)+1,:).';  

        clear pHRIR_L_tmp pHRIR_R_tmp pHRIR_L_ip_tmp pHRIR_R_ip_tmp
            
end

%Save intermediate result (interpolated pre-processed HRTFs) in output struct
interpHRTFset.p.pHRTF_L_ip = pHRTF_L_ip;
interpHRTFset.p.pHRTF_R_ip = pHRTF_R_ip;
if strcmp(ppMethod,'OBTA')
    interpHRTFset.p.toa_L_ip = toa_L_ip;
    interpHRTFset.p.toa_R_ip = toa_R_ip;
end
if strcmp(ppMethod,'MagPhase')
    interpHRTFset.p.pHRTF_L_ph_ip = pHRTF_L_ph_ip;
    interpHRTFset.p.pHRTF_R_ph_ip = pHRTF_R_ph_ip;
end

%% Perform postprocessing

switch ppMethod
    case 'SUpDEq'
        disp('Postprocessing: SUpDEq')
        %Apply SUpDEq - De-Equalize interpolated HRTFs

        %Get De-Equalization functions
        [deqTF_L,deqTF_R] = supdeq_getEqTF(eqDataset,ipSamplingGrid,'DEG',2,'ak',false);
        
        %De-Equalize
        HRTF_L_ip = pHRTF_L_ip.*deqTF_L;
        HRTF_R_ip = pHRTF_R_ip.*deqTF_R;
       
    case 'SUpDEq_Lim'
        disp('Postprocessing: SUpDEq_Lim')
        %Apply SUpDEq - De-Equalize interpolated HRTFs

        %Get De-Equalization functions of limited eqDataset
        [deqTF_L,deqTF_R] = supdeq_getEqTF(eqDataset,ipSamplingGrid,'DEG',2,'ak',false);
        
        %De-Equalize
        HRTF_L_ip = pHRTF_L_ip.*deqTF_L;
        HRTF_R_ip = pHRTF_R_ip.*deqTF_R;
        
    case 'SUpDEq_AP'
        disp('Postprocessing: SUpDEq_AP')
        %Apply SUpDEq_AP - De-Equalize interpolated HRTFs with rigid sphere allpass functions

        %Get De-Equalization functions
        [deqTF_L,deqTF_R] = supdeq_getEqTF(eqDataset,ipSamplingGrid,'DEG',2,'ak',true);
        
        %De-Equalize
        HRTF_L_ip = pHRTF_L_ip.*deqTF_L;
        HRTF_R_ip = pHRTF_R_ip.*deqTF_R;
        
    case 'SUpDEq_Lim_AP'
        disp('Postprocessing: SUpDEq_Lim_AP')
        %Apply SUpDEq_Lim_AP - De-Equalize interpolated HRTFs with limited rigid sphere allpass functions

        %Get De-Equalization functions of limited eqDataset
        [deqTF_L,deqTF_R] = supdeq_getEqTF(eqDataset,ipSamplingGrid,'DEG',2,'ak',true);
        
        %De-Equalize
        HRTF_L_ip = pHRTF_L_ip.*deqTF_L;
        HRTF_R_ip = pHRTF_R_ip.*deqTF_R;

    case 'PC'
        disp('Postprocessing: PC')
        %Apply Inverse Phase-Correction

        %Transform interpolation sampling grid to radiant
        sg = ipSamplingGrid * pi / 180;

        %Get inverse phase correction term for left/right ear according to Eq. 13 in [4]
        cosThetaL = cos(sg(:,2)')*cos(pi/2) + sin(sg(:,2)')*sin(pi/2) .* cos(sg(:,1)'-pi/2);
        %No minus before 1j to get inverse phase term. Could also be achieved by division instead of multiplication
        phaseCorL = exp(1j*headRadius * k .* cosThetaL);
        cosThetaR = cos(sg(:,2)')*cos(pi/2) + sin(sg(:,2)')*sin(pi/2) .* cos(sg(:,1)'+pi/2);
        %No minus before 1j to get inverse phase term. Could also be achieved by division instead of multiplication
        phaseCorR = exp(1j*headRadius * k .* cosThetaR);

        %Apply to interpolated HRTFs
        HRTF_L_ip = pHRTF_L_ip .* phaseCorL.'; %Just transpose, no conjugate complex
        HRTF_R_ip = pHRTF_R_ip .* phaseCorR.';

    case 'OBTA' 
        disp('Postprocessing: OBTA')
        %Re-Insert interpolated TOAs with fractional delays
        
        %Get interpolated HRIRs
        pHRIR_L_ip = AKsingle2bothSidedSpectrum(pHRTF_L_ip.');
        pHRIR_R_ip = AKsingle2bothSidedSpectrum(pHRTF_R_ip.');
        pHRIR_L_ip = real(ifft(pHRIR_L_ip));
        pHRIR_R_ip = real(ifft(pHRIR_R_ip));
        
        %Re-Insert TOAs with fractional delays
        HRIR_L_ip = AKfractionalDelayCyclic(pHRIR_L_ip, toa_L_ip);
        HRIR_R_ip = AKfractionalDelayCyclic(pHRIR_R_ip, toa_R_ip);
        
        %Transform back to frequency domain
        HRTF_L_ip = AKboth2singleSidedSpectrum(fft(HRIR_L_ip)).';
        HRTF_R_ip = AKboth2singleSidedSpectrum(fft(HRIR_R_ip)).';
        
    case 'MagPhase'
        disp('Postprocessing: MagPhase')
        
        %Combine magnitude and phase of interpolated HRTFs again
        HRTF_L_ip = pHRTF_L_ip.*exp(1i*pHRTF_L_ph_ip); 
        HRTF_R_ip = pHRTF_R_ip.*exp(1i*pHRTF_R_ph_ip); 
        
    case 'None' 
        disp('Postprocessing: None')

        %Simply copy pHRTF_L_ip to HRTF_L_ip - Easier for further processing
        HRTF_L_ip = pHRTF_L_ip;
        HRTF_R_ip = pHRTF_R_ip;
end

%% Correct magnitude if mc ~= nan

if ~isnan(mc)
    
    disp('MC postprocessing');
    
    %Save intermediate result (interpolated HRTFs without mc) in output struct
    interpHRTFset.p.HRTF_L_ip_noMC = HRTF_L_ip;
    interpHRTFset.p.HRTF_R_ip_noMC = HRTF_R_ip;
    
    %Get input HRTFs in auditory bands (ERB filter)
    fLowErb = 50;
    [cl,ferb] = AKerbErrorPersistent(HRIR_L(1:hrirLength,:),AKdirac(hrirLength),[fLowErb fs/2],fs);
    [cr] = AKerbErrorPersistent(HRIR_R(1:hrirLength,:),AKdirac(hrirLength),[fLowErb fs/2],fs);
    
    %Interpolate ERB filters to ipSamplingGrid with respective ipMethod
    switch ipMethod
        case 'SH'

            %Transform ERB filters to SH domain
            %Tikhonov regularization can be applied if weights in sampling grid are erased and tikhEps ~= 0
            cl_nm = AKsht(cl, false, samplingGrid, HRTFset.Nmax, 'complex', fs, false, 'real',tikhEps);
            cr_nm = AKsht(cr, false, samplingGrid, HRTFset.Nmax, 'complex', fs, false, 'real',tikhEps);

            %Apply SH interpolation to ipSamplingGrid
            cl_ip = AKisht(cl_nm, false, ipSamplingGrid(:,1:2), 'complex', false, false, 'real');
            cr_ip = AKisht(cr_nm, false, ipSamplingGrid(:,1:2), 'complex', false, false, 'real');

        case {'NN','SARITA'} %NN and SARITA handled the same, as SARITA uses NN interpolation
            
            %Transform grids to cartesian coordinates
            [samplingGrid_cart(:,1), samplingGrid_cart(:,2), samplingGrid_cart(:,3)] = sph2cart(samplingGrid(:,1)/360*2*pi, pi/2-samplingGrid(:,2)/360*2*pi,1);
            [ipSamplingGrid_cart(:,1), ipSamplingGrid_cart(:,2), ipSamplingGrid_cart(:,3)] = sph2cart(ipSamplingGrid(:,1)/360*2*pi, pi/2-ipSamplingGrid(:,2)/360*2*pi,1);

            %NN interpolate ERB filters to ipSamplingGrid
            cl_ip = zeros(length(ferb),length(ipSamplingGrid));
            cr_ip = zeros(length(ferb),length(ipSamplingGrid));

            for nPoint = 1:size(ipSamplingGrid,1)

                [idx_voronoi,w_voronoi] = findvoronoi(samplingGrid_cart,ipSamplingGrid_cart(nPoint,:));

                for n = 1:length(idx_voronoi) 
                    cl_ip(:,nPoint) = cl_ip(:,nPoint) + w_voronoi(n)*cl(:,idx_voronoi(n));  
                    cr_ip(:,nPoint) = cr_ip(:,nPoint) + w_voronoi(n)*cr(:,idx_voronoi(n));
                end
            end

        case 'Bary'
            
            %Transform grids to cartesian coordinates
            [samplingGrid_cart(:,1), samplingGrid_cart(:,2), samplingGrid_cart(:,3)] = sph2cart(samplingGrid(:,1)/360*2*pi, pi/2-samplingGrid(:,2)/360*2*pi,1);
            [ipSamplingGrid_cart(:,1), ipSamplingGrid_cart(:,2), ipSamplingGrid_cart(:,3)] = sph2cart(ipSamplingGrid(:,1)/360*2*pi, pi/2-ipSamplingGrid(:,2)/360*2*pi,1);

            %Get simplified convex hull (could also be done after delaunay-triangulation)
            convHullIdx = convhull(samplingGrid_cart(:,1), samplingGrid_cart(:,2), samplingGrid_cart(:,3), 'simplify', true);

            %Get HRTFs / directions assigned to convex hull triangles
            convHullHRTF_A = samplingGrid_cart(convHullIdx(:,1),:);
            convHullHRTF_B = samplingGrid_cart(convHullIdx(:,2),:);
            convHullHRTF_C = samplingGrid_cart(convHullIdx(:,3),:);

            %Interpolate with barycentric weights
            %Right triangle of convex hull picked by ray-triangle intersection test (intersectLinePolygon3d from geom3D toolbox should work too)
            %Function returns barycentric weights / coordinates u and v of the intersection point -> w = 1-u-v
            %P = w*A + u*B + v*C
            orig = [0 0 0];
            cl_ip = zeros(length(ferb),length(ipSamplingGrid));
            cr_ip = zeros(length(ferb),length(ipSamplingGrid));

            for nPoint = 1:size(ipSamplingGrid,1)

                [intersect, ~, u, v, ~] = TriangleRayIntersection(orig, ipSamplingGrid_cart(nPoint,:), convHullHRTF_A, convHullHRTF_B, convHullHRTF_C, 'border', 'inclusive');
                intersectIdx = find(intersect, 1, 'first'); %Evaluate first intersection point
                u = u(intersectIdx); v = v(intersectIdx); w = 1-u-v;
                bw = [w u v]; %Barycentric weights
                idx_bw = convHullIdx(intersectIdx,:); %Indizes of HRTFs of convex hull / triangles

                for n = 1:3 %Always 3 because of triangulization
                    cl_ip(:,nPoint) = cl_ip(:,nPoint) + bw(n)*cl(:,idx_bw(n));  
                    cr_ip(:,nPoint) = cr_ip(:,nPoint) + bw(n)*cr(:,idx_bw(n));
                end
            end
    end
    
    %Save interpolated ERBs in 3D-array for further processing
    c_ip(:,:,1) = cl_ip;
    c_ip(:,:,2) = cr_ip;
    
    %Get interpolated HRTFs in auditory bands (ERB filter)
    HRIR_L_ip = AKsingle2bothSidedSpectrum(HRTF_L_ip.');
    HRIR_R_ip = AKsingle2bothSidedSpectrum(HRTF_R_ip.');
    HRIR_ip(:,:,1) = real(ifft(HRIR_L_ip));
    HRIR_ip(:,:,2) = real(ifft(HRIR_R_ip));
    HRIR_ip = HRIR_ip(1:hrirLength,:,:);
    %Save in 3D-array for further processing (Spectrum c(hrir) interpolated)
    c_hrir_ip = AKerbErrorPersistent(HRIR_ip,AKdirac(size(HRIR_ip,1)),[fLowErb fs/2],fs);

    
    %ILD compensation <--- VAR Project code
    if ILDComp

        %Calculate ILDs post auditory smoothings from the conventional time aligned interpolation (CTAI)
        ILDs_hrir_ip = c_hrir_ip(:, :, 1) - c_hrir_ip(:, :, 2);

        %Calculate ILDs post auditory smoothings from the original sparse HRIRs
        ILDs_pre_ip = cl - cr;

        %Interpolation of the sparse ILDs to the ipSamplingGrid with respect to the chosen ipMethod
        % ILDs_ip = ...
        switch ipMethod
            case 'SH'

                %Transform ILDs to SH domain
                %Tikhonov regularization can be applied if weights in sampling grid are erased and tikhEps ~= 0
                ILDs_nm = AKsht(ILDs_pre_ip, false, samplingGrid, HRTFset.Nmax, 'complex', fs, false, 'real',tikhEps);

                %Apply SH interpolation to ipSamplingGrid
                ILDs_ip = AKisht(ILDs_nm, false, ipSamplingGrid(:,1:2), 'complex', false, false, 'real');

            case {'NN','SARITA'} %NN and SARITA handled the same, as SARITA uses NN interpolation

                %Transform grids to cartesian coordinates
                % --> already done previously

                %NN interpolate ERB filters to ipSamplingGrid
                ILDs_ip = zeros(length(ferb),length(ipSamplingGrid));

                for nPoint = 1:size(ipSamplingGrid,1)

                    [idx_voronoi,w_voronoi] = findvoronoi(samplingGrid_cart,ipSamplingGrid_cart(nPoint,:));

                    for n = 1:length(idx_voronoi)
                        ILDs_ip(:,nPoint) = ILDs_ip(:,nPoint) + w_voronoi(n)*ILDs_pre_ip(:,idx_voronoi(n));
                    end
                end

            case 'Bary'

                %Transform grids to cartesian coordinates
                % --> already done previously

                %Get simplified convex hull (could also be done after delaunay-triangulation)
                % --> already done previously

                %Get HRTFs / directions assigned to convex hull triangles
                % --> already done previously

                %Interpolate with barycentric weights
                %Right triangle of convex hull picked by ray-triangle intersection test (intersectLinePolygon3d from geom3D toolbox should work too)
                %Function returns barycentric weights / coordinates u and v of the intersection point -> w = 1-u-v
                %P = w*A + u*B + v*C
                orig = [0 0 0];
                ILDs_ip = zeros(length(ferb),length(ipSamplingGrid));

                for nPoint = 1:size(ipSamplingGrid,1)

                    [intersect, ~, u, v, ~] = TriangleRayIntersection(orig, ipSamplingGrid_cart(nPoint,:), convHullHRTF_A, convHullHRTF_B, convHullHRTF_C, 'border', 'inclusive');
                    intersectIdx = find(intersect, 1, 'first'); %Evaluate first intersection point
                    u = u(intersectIdx); v = v(intersectIdx); w = 1-u-v;
                    bw = [w u v]; %Barycentric weights
                    idx_bw = convHullIdx(intersectIdx,:); %Indizes of HRTFs of convex hull / triangles

                    for n = 1:3 %Always 3 because of triangulization
                        ILDs_ip(:,nPoint) = ILDs_ip(:,nPoint) + bw(n)*ILDs_pre_ip(:,idx_bw(n));
                    end
                end
        end

        %Division of the CTAI ILDs with the interpolated ILDs from the original sparse HRIRs
        ILDs = ILDs_hrir_ip - ILDs_ip;

        %Left or right ear ILDs selection

        %Generate filter...

        %Apply filter...



    end
    
    %Get correction filters for all HRTFs
    corrFilt = c_ip-c_hrir_ip;
    %Spline interpolate correction filters to 0-fs/2
    corrFilt_l = AKinterpolateSpectrum( squeeze(corrFilt(:,:,1)), ferb, NFFT, {'nearest' 'spline' 'nearest'}, fs);
    corrFilt_r = AKinterpolateSpectrum( squeeze(corrFilt(:,:,2)), ferb, NFFT, {'nearest' 'spline' 'nearest'}, fs);
    
    %Limit to f > fA (set to 0 below fA) if desired
    if limitMC
        
        %Get spatial alasing frequency
        fA = HRTFset.Nmax*c/2/pi/headRadius;
        
        disp(['MC postprocessing limited to f > fA, with fA = ',num2str(round(fA,2)),'Hz']);
        
        if strcmp(limFade,'fadeDown')
            
            disp('MC fade downwards');
            
            %Get third-octave below fA
            fAt = fA/2^(1/3);
            
            %Set to 0 (logarithmic) below fAt and apply linear fade from fAt to fA
            [~,fA_bin] = min(abs(HRTFset.f-fA));
            [~,fAt_bin] = min(abs(HRTFset.f-fAt));
            ramp = linspace(0,1,abs(fA_bin-fAt_bin)+1);
            corrFilt_l(1:fAt_bin-1,:) = 0;
            corrFilt_l(fAt_bin:fA_bin,:) = corrFilt_l(fAt_bin:fA_bin,:).*ramp.';
            corrFilt_r(1:fAt_bin-1,:) = 0;
            corrFilt_r(fAt_bin:fA_bin,:) = corrFilt_r(fAt_bin:fA_bin,:).*ramp.'; 
        end
        
        if strcmp(limFade,'fadeUp')
            
            disp('MC fade upwards');
            
            %Get third-octave above fA
            fAt = fA*2^(1/3);
            
            %Set to 0 (logarithmic) below fA and apply linear fade from fA to fAt
            [~,fA_bin] = min(abs(HRTFset.f-fA));
            [~,fAt_bin] = min(abs(HRTFset.f-fAt));
            ramp = linspace(0,1,abs(fAt_bin-fA_bin)+1);
            corrFilt_l(1:fA_bin-1,:) = 0;
            corrFilt_l(fA_bin:fAt_bin,:) = corrFilt_l(fA_bin:fAt_bin,:).*ramp.';
            corrFilt_r(1:fA_bin-1,:) = 0;
            corrFilt_r(fA_bin:fAt_bin,:) = corrFilt_r(fA_bin:fAt_bin,:).*ramp.'; 
        end
        
    end
    
    %Transform to linear values
    corrFilt_l = 10.^(corrFilt_l/20);
    corrFilt_r = 10.^(corrFilt_r/20);
    
    if ~isinf(mc)
        %Apply soft-limiting to correction filters according to mc if not inf
        disp(['MC maximum boost: ',num2str(mc),'dB'])
        disp(['MC soft-limiting knee: ',num2str(mcKnee),'dB'])

        corrFilt_l_lim = AKsoftLimit(corrFilt_l, mc, mcKnee,[0 fs/2], fs, true);
        corrFilt_r_lim = AKsoftLimit(corrFilt_r, mc, mcKnee,[0 fs/2], fs, true);
        
    else %No soft-limiting / Unlimited gain
        
        disp('MC maximum boost: inf')
        
        corrFilt_l_lim = corrFilt_l;
        corrFilt_r_lim = corrFilt_r;
        
    end
    
    %Design minimum phase filters and use for correction instead of zero
    %phase filters (if mcMinPhase = false)
    if mcMinPhase
        
        disp('MC phase: minimum');
        
        corrFilt_l_lim = AKsingle2bothSidedSpectrum(corrFilt_l_lim);
        corrFilt_r_lim = AKsingle2bothSidedSpectrum(corrFilt_r_lim);

        corrFilt_l_lim = real(ifft(corrFilt_l_lim));
        corrFilt_r_lim = real(ifft(corrFilt_r_lim));

        corrFilt_l_lim = AKphaseManipulation(corrFilt_l_lim,fs,'min',1,0);
        corrFilt_r_lim = AKphaseManipulation(corrFilt_r_lim,fs,'min',1,0);

        %Go back to frequency domain
        corrFilt_l_lim = AKboth2singleSidedSpectrum(fft(corrFilt_l_lim));
        corrFilt_r_lim = AKboth2singleSidedSpectrum(fft(corrFilt_r_lim));
        
    else
        
        disp('MC phase: zero');
        
    end
    
    %Write in 3D array
    corrFilt_lim(:,:,1) = corrFilt_l_lim;
    corrFilt_lim(:,:,2) = corrFilt_r_lim;
    
    %Save intermediate results (correction filter) in output struct
    interpHRTFset.p.corrFilt_lim = corrFilt_lim;    
        
    %Apply magnitude correction filters to HRTFs
    HRTF_L_ip = HRTF_L_ip.*corrFilt_lim(:,:,1).';
    HRTF_R_ip = HRTF_R_ip.*corrFilt_lim(:,:,2).';
    
end

%% Write result in struct

%Get final HRIRs
HRIR_L_ip = AKsingle2bothSidedSpectrum(HRTF_L_ip.');
HRIR_R_ip = AKsingle2bothSidedSpectrum(HRTF_R_ip.');
HRIR_L_ip = real(ifft(HRIR_L_ip));
HRIR_R_ip = real(ifft(HRIR_R_ip));
%Cut zeros depending on FFToversize
HRIR_L_ip = HRIR_L_ip(1:hrirLength,:);
HRIR_R_ip = HRIR_R_ip(1:hrirLength,:);

interpHRTFset.HRTF_L = HRTF_L_ip;
interpHRTFset.HRTF_R = HRTF_R_ip;
interpHRTFset.HRIR_L = HRIR_L_ip;
interpHRTFset.HRIR_R = HRIR_R_ip;
interpHRTFset.f = HRTFset.f;
interpHRTFset.fs = fs;
interpHRTFset.FFToversize = HRTFset.FFToversize;
interpHRTFset.samplingGrid = ipSamplingGrid;
interpHRTFset.headRadius = headRadius;
interpHRTFset.ppMethod = ppMethod;
interpHRTFset.ipMethod = ipMethod;
interpHRTFset.headRadius = headRadius;
if tikhEps ~= 0
    interpHRTFset.tikhEps = tikhEps;
end
if ~isnan(mc)
    interpHRTFset.mc = mc;
    interpHRTFset.mcKnee = mcKnee;
    interpHRTFset.mcMinPhase = mcMinPhase;
    if limitMC
        interpHRTFset.limitMC = 1;
        interpHRTFset.limitMC_fA = fA;
        interpHRTFset.limitMC_fAt = fAt;
        interpHRTFset.limitMC_fade = limFade;
    end
end
disp('Done...');

end

