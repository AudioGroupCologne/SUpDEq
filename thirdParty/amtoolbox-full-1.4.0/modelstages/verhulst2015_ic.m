function IC = verhulst2015_ic(CN,fs,Tex,Tin,dly,Aic,Sic)
% VERHULST2015_IC inferior-colliculus model
%
%   Usage: IC = verhulst2015_ic(CN,fs,Tex,Tin,dly,Aic,Sic)
%
%   Inferior Colliculus model based on the same-frequency excitatory-inhibitory
%   (SFEI) from Nelson and Carney (2004). This function implements
%   Eq. 13 from Verhulst et al. (2015).
%
%   License:
%   --------
%
%   This model is licensed under the UGent Academic License. Further usage details are provided 
%   in the UGent Academic License which can be found in the AMT directory "licences" and at 
%   <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.
%
%   References:
%     P. C. Nelson and L. Carney. A phenomenological model of peripheral and
%     central neural responses to amplitude-modulated tones. J. Acoust. Soc.
%     Am., 116(4), 2004.
%     
%     S. Verhulst, H. Bharadwaj, G. Mehraei, C. Shera, and
%     B. Shinn-Cunningham. Functional modeling of the human auditory
%     brainstem response to broadband stimulation. jasa, 138(3):1637--1659,
%     2015.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/verhulst2015_ic.php


%   #License: ugent
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal PYTHON C
%   #Author: Alejandro Osses (2020): primary implementation
%   #Author: Piotr Majdak (2021): adaptations for the AMT 1.0

% This file is licenced under the terms of the UGent Academic License, which details can be found in the AMT directory "licences" and at <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.
% For non-commercial academic research, you can use this file and/or modify it under the terms of that license. This file is distributed without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. 

%   #License: ugent
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal PYTHON C
%   #Author: Alejandro Osses (2020): primary implementation based on https://github.com/HearingTechnology/Verhulstetal2018Model
%   #Author: Piotr Majdak (2021): adaptations for the AMT 1.0


if nargin < 3
    Tex=0.5e-3;
end
if nargin < 4
    Tin=2e-3;
end
if nargin < 5
    dly = 2e-3; % default in Verhulst 2018
end
if nargin < 6
    Aic=1;
end
if nargin < 7
    Sic=1.5;
end

N_Tau = length(Tex);

for i = 1:N_Tau

    Tex_here = Tex(i);
    Tin_here = Tin(i);
    dly_IC   = dly(i);
    Aic_here = Aic(i);
    Sic_here = Sic(i);

    inhibition_delay_samples = round(dly_IC*fs);

    if iscell(CN)
        CN_here = CN{i};
    else
        CN_here = CN;
    end
    L = size(CN_here,1);

    delayed_inhibition = zeros(size(CN_here));
    delayed_inhibition(inhibition_delay_samples+1:end,:) = CN_here(1:L-inhibition_delay_samples,:);

	[bEx,aEx] = local_irrcoeff(Tex_here,fs);
	[bIn,aIn] = local_irrcoeff(Tin_here,fs);

    IC_excitatory = filter(bEx,aEx,CN_here);
    IC_inhibitory = filter(bIn,aIn,delayed_inhibition);

    IC{i} = Aic_here*(IC_excitatory-Sic_here*IC_inhibitory);
end

if N_Tau == 1
    IC = IC{i};
end

function [num,den] = local_irrcoeff(Tau,fs)
%   The filter characterised by the transfer function with num and den as the 
%   coefficients of the numerator and denominator, respectively, is obtained 
%   as normalised alpha functions. The design of this filter was done applying 
%   a bilinear transformation. This filter design is equivalent to the alpha 
%   functions described by Nelson and Carney (2004).

factor = 1/(2*fs*Tau+1)^2;
m = (2*fs*Tau-1)/(2*fs*Tau+1);
a0  = 1;
a1  = -2*m;
a2  = m^2;
b0 = 1;
b1 = 2;
b2 = 1;

num = [b0 b1 b2]*factor;
den = [a0 a1 a2];

