%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function HRIRs_sfd_ext = supdeq_rangeExt(HRIRs_sfd, r1, regMethod, ampLimit, compTimeShift, compEnergy, rmin, r0)
%
% This function performs range extrapolation of HRTFs in SH-domain. This
% way, near-field HRTFs can be obtained from far-field HRTFs. In brief,
% radial propagation functions are applied to the SH-coefficients to 
% calculate the range-extrapolated SH-coefficients for the new distance r1.?
% As the radial functions (in case of an exterior problem the hankel
% quotient) increase dramatically towards lower kr, a regularization method
% like for example soft-limiting or the Tikhonov regularization needs to be
% applied. See references for more information about range extrapolation 
% and regularization approaches.
%
% Output:
% HRIRs_sfd_ext     - Struct with Spherical Harmonic coefficients 
%                     (SH-Coefficients) for the left (Hl_nm) and right (Hr_nm) 
%                     channel/ear of the range-extrapolated HRTFs, absolute 
%                     frequency scale f, transform order N, FFToversize,
%                     (new range-extrapolated) source distance r1,
%                     non-regularized radial functions (hankel quotient),
%                     regularized (and applied) hankel quotient, and info
%                     about applied regularization and range extrapolation
%
% Input:        
% HRIRs_sfd         - Struct with Spherical Harmonic coefficients 
%                     (SH-Coefficients) for the left (Hl_nm) and right (Hr_nm) 
%                     channel/ear, absolute frequency scale f, 
%                     transform order N, FFToversize, and source distance r0
% r1                - New range-extrapolated source distance r1 in m
% regMethod         - Regularization method (string parameter) for the 
%                     radial function (the hankel quotient). The following 
%                     methods are available:
%                     'none'     - no regularization
%                     'disc'     - simply discard contributions where the
%                                  radial filter magnitude exceeds ampLimit
%                     'disc_kr'  - discard contributions at higher orders
%                                  based on kr_min. Results are pretty
%                                  similar to normal 'disc'
%                     'hardLimit'- set a hard limit for all radial filter
%                                  components at ampLimit instead of 
%                                  dropping them
%                     'softLimit'- soft-limit the components 
%                     'Tikh'     - apply Tikhonov regularization
%                     See Bernschütz2011 and Rettberg2014 for an overview
%                     of most regularization methods. Implementation of most 
%                     of the regularization methods according to Rettberg2014
%                     Regularization 'disc_kr' implemented according to
%                     Duraiswami2004 and Pollow2012
%                     Default: 'softLimit'
% ampLimit          - Amplitude limit of radial functions in dB. Not needed 
%                     for all regularization methods, but since the default
%                     regularization is 'softLimiting', the amplitude limit
%                     can be freely chosen here.
%                     Default: 20 dB
% compTimeShift     - true/false - compensate time shift induced by
%                     range-extrapolation with an inverse phase term
%                     Default: true
% compEnergy        - true/false - compensate the in-/decrease in energy
%                     averaged over all orders and modes 
%                     (single value compensation)
%                     Default: true
% rmin              - Smallest radius in m, meaning the radius of the 
%                     smallest sphere enclosing the listener?s head 
%                     with all its significant scattering sources. This
%                     defines the truncation order N = [krmin] when
%                     regMethod 'disc_kr' according to Duraiswami2004 and
%                     Pollow2012 is used. Otherwise, rmin has no influence.
%                     Default: 0.1
% r0                - Actual source distance in m input HRIRs_sfd.
%                     The HRIRs_sfd structs have a field with the actual
%                     source distance and the function will use the value
%                     in this field as r0. However, if this field is
%                     missing, r0 can be given as an argument.
%                     Default: [] (use value in HRIRs_sfd.sourceDistance) 
%
% Dependencies: AKtools
%
% References:
% R. Duraiswami, D. N. Zotkin, and N. A. Gumerov, ?Interpolation and range 
% extrapolation of HRTFs,? in Proceedings of the IEEE International 
% Conference on Acoustics, Speech, and Signal Processing, 2004, pp. IV45-IV48.
%
% M. Pollow, K.-V. Nguyen, O. Warusfel, T. Carpentier, M. Müller-Trapet, 
% M. Vorländer, and M. Noisternig, ?Calculation of Head-Related Transfer 
% Functions for Arbitrary Field Points Using Spherical Harmonics 
% Decomposition,? Acta Acust. united Ac., vol. 98, no. 1, pp. 72?82, 2012. 
%
% T. Rettberg and S. Spors, ?Time-Domain Behaviour of Spherical Microphone 
% Arrays at High Orders,? in Proceedings of the 40th DAGA, 2014, pp. 598?599.
%
% B. Bernschutz, C. Pörschmann, S. Spors, and S. Weinzierl, ?Soft-Limiting 
% der modalen Amplitudenverstärkung bei sphärischen Mikrofonarrays im 
% Plane Wave Decomposition Verfahren,? 
% in Proceedings of the 37th DAGA, 2011, pp. 661?662.
%
% Boaz Rafaely: Fundamentals of spherical array processing. In. Springer 
% topics in signal processing. Benesty, J.; Kellermann, W. (Eds.), 
% Springer, Heidelberg et al. (2015).
%
% Benjamin Bernschütz: Microphone Arrays and Sound Field Decomposition 
% for Dynamic Binaural Recording. Ph.D. dissertation, Technical University
% Berlin (2016).
%   
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function HRIRs_sfd_ext = supdeq_rangeExt(HRIRs_sfd, r1, regMethod, ampLimit, compTimeShift, compEnergy, rmin, r0)

if nargin < 1 || isempty(HRIRs_sfd)
    error('Please specifiy HRIRs_sfd!');
end

if nargin < 2 || isempty(r1)
    error('Please specifiy r1!');
end

if nargin < 3 || isempty(regMethod)
    regMethod = 'softLimit';
end

if nargin < 4 || isempty(ampLimit)
    ampLimit = 20;
end

if nargin < 5 || isempty(compTimeShift)
    compTimeShift = 1;
end

if nargin < 6 || isempty(compEnergy)
    compEnergy = 1;
end

if nargin < 7 || isempty(rmin)
    rmin = 0.1;
end

if nargin < 7 || isempty(rmin)
    rmin = 0.1;
end

if nargin < 8 || isempty(r0)
    if isfield(HRIRs_sfd,'sourceDistance')
        r0 = HRIRs_sfd.sourceDistance;
    else
        error('Please specifiy r0!');
    end
end

%% Get kr and N

f   = HRIRs_sfd.f;
c   = 343; %Set speed of sound as default with 343 m/s
k   = 2*pi*f/c;
kr0 = k*r0;
kr1 = k*r1;
N   = HRIRs_sfd.N; 

%% Get hankel quotient for range extrapolation

%Get hankel functions
hn_r0 = zeros(N+1,length(f));
hn_r1 = zeros(N+1,length(f));
for kk = 0:N
    hn_r0(kk+1,:) = AKshRadial(kr0,'hankel',2,kk);
    hn_r1(kk+1,:) = AKshRadial(kr1,'hankel',2,kk);
end

%Calculate hankel quotient
hn_ext = hn_r1./hn_r0;

%Compensate time shift if intended
if compTimeShift
    t = ((r0-r1)/c);
    timeShift = exp((0-1i) * complex((2*pi*f*t)+0i));
    for kk = 0:N
        hn_ext(kk+1,:) = hn_ext(kk+1,:).*timeShift;
    end
end
hn_ext_noReg = hn_ext;

%% Apply regularization to hankel quotient

%Implementation according to Rettberg2014
%https://github.com/spatialaudio/sfa-numpy/blob/master/micarray/modal/radial.py
%Function "regularize"

if ~strcmp(regMethod,'disc_kr') %Regularization disc_kr implemented in a different way...
    a0_lin = 10^(ampLimit/20);
    %Get values where abs(hn_ext) > maxGain
    idx = abs(hn_ext) > a0_lin;

    %Get regularisation function hn_reg depending on method
    if strcmp(regMethod,'none')
        hn_reg = ones(size(hn_ext,1),size(hn_ext,2));

    elseif strcmp(regMethod,'disc')
        hn_reg = ones(size(hn_ext,1),size(hn_ext,2));
        hn_reg(idx) = 0;

    elseif strcmp(regMethod,'hardLimit')
        hn_reg = ones(size(hn_ext,1),size(hn_ext,2));
        hn_reg(idx) = a0_lin ./ abs(hn_ext(idx));

    elseif strcmp(regMethod,'softLimit')
        scaling = pi/2;
        hn_reg = a0_lin ./ abs(hn_ext);
        hn_reg = 2/pi * atan(scaling*hn_reg);

    elseif strcmp(regMethod,'Tikh')
        a0_lin = sqrt(a0_lin/2);
        alpha = (1 - sqrt(1 - 1/(a0_lin^2))) / (1 + sqrt(1 - 1/(a0_lin^2)));
        hn_reg = 1 ./ (1 + alpha^2 * abs(hn_ext).^2);  
    end

    %Apply hn_reg to hn_ext
    hn_ext = hn_ext .* hn_reg;
end


%% Get hankel quotient matrix with size of SH-coefficients matrix

%Set nan values at kr = 0 to 0
hn_ext(isnan(hn_ext)) = 0;

hn_ext_matrix = zeros(size(HRIRs_sfd.Hl_nm,1),size(HRIRs_sfd.Hl_nm,2));
modeCounter = 1;
for n = 1:N+1
    %Get "corner"-modes
    modeLeft    = n^2-modeCounter+1;
    modeRight   = n^2;
       
    %Define modes vector
    modes = modeLeft:modeRight;
    
    for mcr = 1:length(modes)
        hn_ext_matrix(modes(mcr),:) = hn_ext(n,:);
    end
    modeCounter = modeCounter + 2;
end

%% Apply hankel quotient matrix to SH-coefficients

Hl_nm_r0 = HRIRs_sfd.Hl_nm;
Hr_nm_r0 = HRIRs_sfd.Hr_nm;
Hl_nm_r1 = zeros(size(Hl_nm_r0,1),size(Hl_nm_r0,2));
Hr_nm_r1 = zeros(size(Hr_nm_r0,1),size(Hr_nm_r0,2));

if  strcmp(regMethod,'disc_kr') %With disc_kr regularization according to Duraiswami2004
    krmin = k*rmin;
    Nh = round(krmin)+1; %Duraiswami2004, Eq.8
    Nh(Nh>N) = N;
    
    for bin = 1:size(Hl_nm_r0,2)
        Nbin = Nh(bin);
        range = (Nbin+1)^2;
        Hl_nm_r1(1:range,bin) = Hl_nm_r0(1:range,bin) .* hn_ext_matrix(1:range,bin);
        Hr_nm_r1(1:range,bin) = Hr_nm_r0(1:range,bin) .* hn_ext_matrix(1:range,bin);
    end
    
else %With any other regularization already applied above to hn_ext
    
    for kk = 1:size(hn_ext_matrix,1)
        Hl_nm_r1(kk,:) = Hl_nm_r0(kk,:) .* hn_ext_matrix(kk,:);
        Hr_nm_r1(kk,:) = Hr_nm_r0(kk,:) .* hn_ext_matrix(kk,:);
    end 
end

%% Apply energy compensation if intended

if compEnergy
    %Get average energy over all orders and modes of SH-coefficients
    [~,eL_r0] = AKshEnergy(Hl_nm_r0);
    [~,eR_r0] = AKshEnergy(Hr_nm_r0);
    [~,eL_r1] = AKshEnergy(Hl_nm_r1);
    [~,eR_r1] = AKshEnergy(Hr_nm_r1);
    
    %Apply energy-compensation factor to SH-coefficients
    Hl_nm_r1 = Hl_nm_r1 * sqrt(eL_r0/eL_r1);
    Hr_nm_r1 = Hr_nm_r1 * sqrt(eR_r0/eR_r1); 
end

%% Write new output struct with additional info on range extrapolation

HRIRs_sfd_ext = HRIRs_sfd;
HRIRs_sfd_ext.Hl_nm             = Hl_nm_r1;
HRIRs_sfd_ext.Hr_nm             = Hr_nm_r1;
HRIRs_sfd_ext.sourceDistance    = r1;
HRIRs_sfd_ext.rangeExt.r0               = r0;
HRIRs_sfd_ext.rangeExt.r1               = r1;
HRIRs_sfd_ext.rangeExt.regMethod        = regMethod;
HRIRs_sfd_ext.rangeExt.ampLimit         = ampLimit;
HRIRs_sfd_ext.rangeExt.compTimeShift    = compTimeShift;
HRIRs_sfd_ext.rangeExt.compEnergy       = compEnergy;
HRIRs_sfd_ext.rangeExt.rmin             = rmin;
HRIRs_sfd_ext.rangeExt.hn_ext           = hn_ext;
HRIRs_sfd_ext.rangeExt.hn_ext_noReg     = hn_ext_noReg; 

end

