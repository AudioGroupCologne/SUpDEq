%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function eqDataset_limited = supdeq_limitEqDataset(eqDataset, Nmax, radius, safetyMargin, fade)
%
% This function can be used to limit the equalization in frequency domain
% dependent on Nmax of the sparse HRTF set. Below the specific aliasing 
% frequency fA, the eqDataset will be set to 0, and thus the subsequent
% equalization will have no effect in this frequency range, resulting in 
% normal SH interpolation below fA. N = 0 is not affected by the limiting!
%
% IMPORTANT (1): If limitEq is applied for equalization using supdeq_eq, 
% it must also be applied for de-equalization using supdeq_deq!
%
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
% safteyMargin          - Boolean to choose if safety margin (third-octave 
%                         below fA) is applied or not. 
%                         Default: true
% fade                  - Boolean to choose whether to apply a fade from fA
%                         to fA/(2^1/3) or not
%                         Default: false (to be compatabile with Arend et al. 2021)
%                         Only applies if safteyMargin = true
%
% Dependencies: -
%
% (C) 2019 by JMA, Johannes M. Arend
%             CP,  Christoph Pörschmann
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function eqDataset_limited = supdeq_limitEqDataset(eqDataset, Nmax, radius, safetyMargin, fade)

if nargin < 3 || isempty(radius)
    radius = eqDataset.radius;
end

if nargin < 4 || isempty(safetyMargin)
    safetyMargin = true;
end

if nargin < 5 || isempty(fade)
    fade = false;
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

if safetyMargin
    %Lower by 1/3 octave to ensure distance to aliasing frequency
    fA_shift = fA/2^(1/3);
else
    fA_shift = fA;
end
%Get corresponding index of frequency bins
[~,fA_shift_bin] = min(abs(eqDataset.f-fA_shift));

%Get eqDataset_limited
eqDataset_limited = eqDataset;

if ~noEq %If SUpDEq should be partly applied
    
    if safetyMargin && ~fade %Without fade
        %Set all SH coefficients for N > 0 to 0
        eqDataset_limited.Hl_nm(2:end,1:fA_shift_bin)=0;
        eqDataset_limited.Hr_nm(2:end,1:fA_shift_bin)=0;
    end
    
    if safetyMargin && fade %With fade from fA_shift to fA
        [~,fA_bin] = min(abs(eqDataset.f-fA));
        
        ramp = linspace(0,1,abs(fA_bin-fA_shift_bin)+1);
        
        eqDataset_limited.Hl_nm(2:end,1:fA_shift_bin-1)=0;
        eqDataset_limited.Hl_nm(2:end,fA_shift_bin:fA_bin)=eqDataset_limited.Hl_nm(2:end,fA_shift_bin:fA_bin).*ramp;
        eqDataset_limited.Hr_nm(2:end,1:fA_shift_bin-1)=0;
        eqDataset_limited.Hr_nm(2:end,fA_shift_bin:fA_bin)=eqDataset_limited.Hr_nm(2:end,fA_shift_bin:fA_bin).*ramp;
    end
    
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
    if safetyMargin && fade
        eqDataset_limited.limitInfo.fade = true;
    end
else
    eqDataset_limited.limitInfo.noEq = noEq;
end

fprintf('eqDataset limited depending on Nmax = %d\n',Nmax);
