function [InternalNoise_L,InternalNoise_R,IN_MaskerLevelDependentPart_L,IN_MaskerLevelDependentPart_R, OHC_L_dBSPL, OHC_R_dBSPL, IHC_L_dBHL, IHC_R_dBHL, Gamma] = vicente2020_internalnoise(fc,f0,Aud_L,Aud_R,N,B,Nlim,eta,max_OHC)
%VICENTE2020_INTERNALNOISE calculates internal noise level
%   Usage: [IN_L,IN_R] = vicente2020_internalnoise(fc,f0,Aud_L,Aud_R,N,B,Nlim,eta,max_OHC)
%
%   Input parameters:
%     fc      : center frequency of the model (vector)
%     f0      : audiometric frequency vector 
%     Aud     : left or right hearing thresholds @ f0
%     N       : Noise level (scalar or vector) 
%     B       : Parameter of the internal noise defining the internal noise floor
%     Nlim    : Parameter of the internal noise defining at which the internal noise level starts to increase with the external noise level
%     eta     : defining the percentage of hearing threshold elevation due to OHC loss
%     max_OHC : maximum allowed for OHC loss
%
%   Output parameters:
%     InternalNoise : Left/Right internal noise level @ model center frequency
%     IN_MaskerLevelDependentPart : Left/Right external noise level dependent part of the internal noise level @ model center frequency
%     Gamma : Transformation to convert dB HL into dB SPL @ model center frequency
%
%   VICENTE2020_INTERNALNOISE simulates the internal noise based on the input audiograms 
%   by calculating the IHC and OHC loss, taking into account the
%   external noise level.
%
%   See also: vicente2020 exp_lavandier2022
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/vicente2020_internalnoise.php


%   #StatusDoc: Perfect
%   #StatusCode: Good
%   #Verification: Qualified
%   #Requirements: MATLAB
%   #Author: Matthieu Lavandier
%   #Author: Clara Hollomey (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% check variables
if isvector(Aud_L)
    if size(Aud_L,1) > 1 
        Aud_L = Aud_L'; 
    end
else
    error('Left audiograms is not a vector')
end
if isvector(Aud_R)
    if size(Aud_R,1) > 1 
        Aud_R = Aud_R'; 
    end
else
    error('Right audiograms is not a vector')
end
if isvector(N)
    if size(N,2) > 1 
        N = N'; 
    end
else
    error('Variable ''N'' must be a vector')
end

% Compute OHC loss and IHC loss 
Aud_L_OHC = min(max_OHC,Aud_L*eta);
Aud_R_OHC = min(max_OHC,Aud_R*eta);
Aud_L_IHC = Aud_L - Aud_L_OHC;
Aud_R_IHC = Aud_R - Aud_R_OHC;

%Interpolation of OHC Loss + converting from dB HL to dB SPL
[OHC_L_dBSPL,OHC_R_dBSPL,~,~,Gamma] = local_interextrapaud_hl2spl(fc,f0,Aud_L_OHC,Aud_R_OHC);
%Interpolation of IHC Loss
[~,~,IHC_L_dBHL,IHC_R_dBHL,~] = local_interextrapaud_hl2spl(fc,f0,Aud_L_IHC,Aud_R_IHC);

% Compute the external-noise-level dependent part of the Internal noise
IN_MaskerLevelDependentPart_L = 10*log10(10^(B/10)+10.^((N - Nlim + IHC_L_dBHL)/10));
IN_MaskerLevelDependentPart_R = 10*log10(10^(B/10)+10.^((N - Nlim + IHC_R_dBHL)/10));

%Compute the entire Internal Noise. 
InternalNoise_L = OHC_L_dBSPL + IN_MaskerLevelDependentPart_L; 
InternalNoise_R = OHC_R_dBSPL + IN_MaskerLevelDependentPart_R;

end

function [aud_L_dBSPL,aud_R_dBSPL,aud_L_dBHL,aud_R_dBHL,adapted_HL2SPL] = local_interextrapaud_hl2spl(fc,f_audTon,Aud_L,Aud_R)

%to be used with vicente2020.m
% fc = center frequency of the model
% f_audTon
% Aud_* = Left or Right Audiogram
%Check inputs
if length(f_audTon) ~= length(Aud_L)
    error('Frequency vector does not have the same length as the the left audiogram vector')
elseif length(f_audTon) ~= length(Aud_R)
    error('Frequency vector does not have the same length as the the right audiogram vector')
end
    
%Set the transformation to convert db HL into dB SPL
Gamma = [17.6 15.7 13 11.8 10.7 9.1 8.9 9.2 8.6 10.5 14.5 17.4 15.3 12.9 13.4 17.7]; % Transformation to convert dB HL into dB SPL
f_HL2SPL = [200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300]; % center frequency of the transformation

% Find the indices corresponding to the frequency range in which the audiogram is interpolated
infTon = find(fc >= f_audTon(1), 1 ); 
supTon = find(fc > f_audTon(end), 1 )-1;

% Do the interpolation / extrapolation for the audiograms
aud_L_dBHL = zeros(1,length(fc)); 
aud_R_dBHL = zeros(1,length(fc));

aud_L_dBHL(1:infTon-1) = Aud_L(1)*ones(1,infTon-1);
aud_R_dBHL(1:infTon-1) = Aud_R(1)*ones(1,infTon-1);
aud_L_dBHL(infTon:supTon) = interp1(log10(f_audTon), Aud_L(:), log10(fc(infTon:supTon))); % interpolation inside the frequency range defined above   
aud_R_dBHL(infTon:supTon) = interp1(log10(f_audTon), Aud_R(:), log10(fc(infTon:supTon)));
aud_L_dBHL(supTon+1:length(fc)) = Aud_L(end)*ones(1,length(fc)-supTon); % extrapolation outside the frequency range defined above
aud_R_dBHL(supTon+1:length(fc)) = Aud_R(end)*ones(1,length(fc)-supTon);

% Find the indices corresponding to the frequency range in which the transformation is interpolated
infTon = find(fc>=f_HL2SPL(1), 1 ); 
supTon = find(fc>f_HL2SPL(end), 1 )-1;

% Do the interpolation / extrapolation for the transform
adapted_HL2SPL = zeros(1,length(fc));
adapted_HL2SPL(1:infTon-1) = Gamma(1)*ones(1,infTon-1);
adapted_HL2SPL(infTon:supTon) = interp1(log10(f_HL2SPL), Gamma(:), log10(fc(infTon:supTon)));   
adapted_HL2SPL(supTon+1:length(fc)) = Gamma(end)*ones(1,length(fc)-supTon);

% Convert audiogram in dB HL into dB SPL
aud_L_dBSPL = adapted_HL2SPL + aud_L_dBHL; 
aud_R_dBSPL = adapted_HL2SPL + aud_R_dBHL;

end



