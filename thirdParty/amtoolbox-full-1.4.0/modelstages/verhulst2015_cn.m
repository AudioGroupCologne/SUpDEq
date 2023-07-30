function CN = verhulst2015_cn(summedAN,fs,Tex,Tin,dly,Acn,Scn)
%VERHULST2015_CN cochlear-nucleus model
%
%   Usage: CN = verhulst2015_cn(summedAN,fs,Tex,Tin,dly,Acn,Scn)
%
%   Cochlear nucleus model based on the same-frequency inhibitory-excitatory
%   (SFIE) design from Nelson and Carney (2004). This function implements
%   Eq. 12 from Verhulst et al. (2015).
%
%   Make sure that the variables Tex, Tin, dly, Acn, Scn have the same dimensions.
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/verhulst2015_cn.php


%   #License: ugent
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal PYTHON C
%   #Author: Alejandro Osses (2020): primary implementation based on https://github.com/HearingTechnology/Verhulstetal2018Model
%   #Author: Piotr Majdak (2021): adaptations for the AMT 1.0

% This file is licenced under the terms of the UGent Academic License, which details can be found in the AMT directory "licences" and at <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.
% For non-commercial academic research, you can use this file and/or modify it under the terms of that license. This file is distributed without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. 

if nargin < 3
    Tex=0.5e-3; % Tau excitation
end
if nargin < 4
    Tin=2e-3; % Tau inhibition
end
if nargin < 5
    dly = 1e-3;
end
if nargin < 6
    Acn=1.5;
end
if nargin < 7
    Scn=0.6;
end

N_Tau = length(Tex);
L = size(summedAN,1);

for i = 1:N_Tau
    
    Tex_here = Tex(i);
    Tin_here = Tin(i);
    dly_CN   = dly(i);
    Acn_here = Acn(i);
    Scn_here = Scn(i);
    
    inhibition_delay_samples = round(dly_CN*fs); 

    delayed_inhibition = zeros(size(summedAN)); 
    delayed_inhibition(inhibition_delay_samples+1:end,:) = summedAN(1:L-inhibition_delay_samples,:); 

    [bEx,aEx] = local_irrcoeff(Tex_here,fs);
    [bIn,aIn] = local_irrcoeff(Tin_here,fs);

    CN_excitatory = filter(bEx,aEx,summedAN);
    CN_inhibitory = filter(bIn,aIn,delayed_inhibition);
    CN{i} = Acn_here*(CN_excitatory-Scn_here*CN_inhibitory);
    
end

if N_Tau == 1
    CN = CN{1};
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


