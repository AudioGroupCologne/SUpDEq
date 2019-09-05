%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function eqDataset_limited = supdeq_limitEqDataset(eqDataset, Nmax, radius)
%
% This function can be used to limit the equalization in frequency domain
% dependent on Nmax of the sparse HRTF set. Below the specific aliasing 
% frequency fA, the eqDataset will be set to 0, and thus the subsequent
% equalization will have no effect in this frequency range, resulting in 
% normal SH interpolation below fA.
% IMPORTANT (1): If limitEq is applied for equalization using supdeq_eq, 
% it must also be applied for de-equalization using supdeq_deq!
% IMPORTANT (2): As fA depends on the radius, and the transition between 
% normal SH interpolation and SUpDEq processing should be the same in the
% equalization and de-equalization, the radius (of the eqDataset) can be
% passed here. Thus, if a different radius is applied for de-equalization,
% the limiting of the deqDataset should be done based on the radius of the
% eqDataset.
%
% Output:
% eqDataset_limited     - Limited eqDataset      
%
% Input:        
% eqDataset             - eqDataset struct obtained with supdeq_getEqDataset
% Nmax                  - Maximum spatial order N of the sparse HRTF set
% radius                - Optionally, the radius can be passed if fA should
%                         be calculated based on another radius than 
%                         specified in the eqDataset struct. This can be 
%                         usefull, if the deqDataset has a different radius,
%                         but the limiting should affect the dataset in the 
%                         exact same frequency range (see important hints).
%                         If the radius is not passed, the values specified 
%                         in the eqDataset struct will be applied.
%
% Dependencies: -
%
% (C) 2019 by JMA, Johannes M. Arend
%             CP,  Christoph Pörschmann
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function eqDataset_limited = supdeq_limitEqDataset(eqDataset, Nmax, radius)

if nargin < 3 || isempty(radius)
    radius = eqDataset.radius;
end

if isfield(eqDataset,'limited')
    if eqDataset.limited
        error('Input eqDataset already limited!');
    end
end

%Just a variable to check if equalizing should be applied
noEq = false;

%Get aliasing frequency fA (N ~ kr)
fA = Nmax * eqDataset.c/ (2*pi*radius);
%Check if fA is >= fs/2
if fA >= eqDataset.f(end)
    noEq = true;
end
%Lower by 1/3 octave to ensure distance to aliasing frequency
fA_shift = fA/2^(1/3);
%Get frequency 1/8 octave above/below fA as start/end point of the
%amplitude fade
fAbove = fA_shift/2^(-1/8);
fBelow = fA_shift/2^(1/8);
%Get corresponding index of frequency bins
[~,fA_bin] = min(abs(eqDataset.f-fA));
[~,fA_shift_bin] = min(abs(eqDataset.f-fA_shift));
[~,fAbove_bin]=min(abs(eqDataset.f-fAbove));
[~,fBelow_bin]=min(abs(eqDataset.f-fBelow));
%Generate linear fade
amplitudeFade = 0:1/(fAbove_bin-fBelow_bin):1;

%Generate sigmoid function fade
linVec = linspace(40,-40,fA_shift_bin);
sigFade = 1 ./ (1 + exp(linVec));
sigFade = [sigFade,ones(1,length(eqDataset.f)-length(sigFade))];

%Get eqDataset_limited
eqDataset_limited = eqDataset;

if ~noEq %If SUpDEq should be partly applied
    %Linear fade
    %Set SH coefficients for N > 0 up to fBelow-1 to 0
    %eqDataset_limited.Hl_nm(2:end,1:fBelow_bin-1)=0;
    %eqDataset_limited.Hr_nm(2:end,1:fBelow_bin-1)=0;
    %Apply fade from fBelow_bin up to fAbove_bin
    %eqDataset_limited.Hl_nm(2:end,fBelow_bin:fAbove_bin)=eqDataset_limited.Hl_nm(2:end,fBelow_bin:fAbove_bin).*amplitudeFade;
    %eqDataset_limited.Hr_nm(2:end,fBelow_bin:fAbove_bin)=eqDataset_limited.Hr_nm(2:end,fBelow_bin:fAbove_bin).*amplitudeFade;
    
    %Sigmoid
    %eqDataset_limited.Hl_nm(2:end,:)=eqDataset_limited.Hl_nm(2:end,:).*sigFade;
    %eqDataset_limited.Hr_nm(2:end,:)=eqDataset_limited.Hr_nm(2:end,:).*sigFade;
    
    %Just to 0
    eqDataset_limited.Hl_nm(2:end,1:fA_shift_bin)=0;
    eqDataset_limited.Hr_nm(2:end,1:fA_shift_bin)=0;
    
else %If fA is > fs/2 and SUpDEq does not need to be applied
    %Note: Actually, in this case SH interpolation without any equalization is
    %perfectly fine anyway ;-).

    %Set all SH coefficients for N > 0 to 0
    eqDataset_limited.Hl_nm(2:end,1:end)=0;
    eqDataset_limited.Hr_nm(2:end,1:end)=0;
end

%Write some infos in new struct eqDataset_limited
eqDataset_limited.limited = true;
eqDataset_limited.limitInfo.Nmax = Nmax;
eqDataset_limited.limitInfo.appliedRadius = radius;
eqDataset_limited.limitInfo.fA = fA;
if ~noEq
    eqDataset_limited.limitInfo.fA_shift = fA_shift;
    %eqDataset_limited.limitInfo.fAbove = fAbove;
    %eqDataset_limited.limitInfo.fBelow = fBelow;
else
    eqDataset_limited.limitInfo.noEq = noEq;
end

fprintf('eqDataset limited depending on Nmax = %d.\n',Nmax);
