function [Ldn,E,Cam,CF] = chen2011(inputF,inputLdB,varargin)
%CHEN2011  Fast excitation pattern estimation
%   Usage: [Ldn,E,Cam,CF] = chen2011(inputF,inputLdB)
%
%
%
%   Input parameters:
%     inputF           : vector of input frequency values
%     inputLdB         : vector of input level
%
%   Output parameters:
%     Ldn              : loudness
%     E                : excitation at the output of the passive filter
%     Cam              : ERB-number in Cams (Glasberg and Moore, 1990; Moore, 2003)
%     CF               : center frequency of auditory filters
%
%
%   CHEN2011 provides a fast way of loudness estimation.
%
%   Optional parameters:
%
%     'HLcf',HLcf                parameters to simulate hearing loss: audiogram frequencies
%
%     'HLohcdB0',HLohcdB0        parameters to simulate hearing loss: OHC loss at audiogram frequencies
%
%     'HLihcdB0',HLihcdB0        parameters to simulate hearing loss: IHC loss at audiogram frequencies
%
%     'cambin',cambin            spacing [ERB] between successive auditory filter CFs
%
%     'flow',flow                lowest center frequency of an auditory filter
%
%     'fhigh',fhigh              highest center frequency of an auditory filter
%
%   Moreover, the model supports 3 flags for outer ear correction ('FreeField', 'PDR10', 'Eardrum').
%
%
%
%   References:
%     Z. Chen, G. Hu, B. R. Glasberg, and C. Moore, Brian. A new model for
%     calculating auditory excitation patterns and loudness for cases of
%     cochlear hearing loss. J. Acoust. Soc. Am., 282(1), 2011.
%     
%     B. R. Glasberg and B. C. J. Moore. A Model of Loudness Applicable to
%     Time-Varying Sounds. J. Audio Eng. Soc, 50(5):331--342, 2002.
%     
%     B. C. J. Moore, B. R. Glasberg, and T. Baer. A Model for the Prediction
%     of Thresholds, Loudness, and Partial Loudness. J. Audio Eng. Soc,
%     45(4):224--240, 1997.
%     
%
%   See also: demo_chen2011 exp_chen2011
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/chen2011.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Author : Zhangli Chen: original code
%   #Author : Clara Hollomey (2020): integration in the AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


if length(inputF) ~= length(inputLdB)
    error('inputF and inputLdB should be dB/Hz and have same length');
end

definput.import = {'chen2011'}; % load defaults from arg_chen2011

[flags,kv]  = ltfatarghelper({},definput,varargin);

%calculation of auditory filter CF
Cam = (fc2erb(kv.flow):kv.cambin:fc2erb(kv.fhigh))';     % 40-17000Hz
CF = 1000*(10.^(Cam/21.4)-1)/4.37;   % Hz

if length(kv.HLcf) ~= length(kv.HLohcdB0) || length(kv.HLcf) ~= length(kv.HLihcdB0)
    error('HLcf, HLohcdB and HLihcdB should have same length');
elseif length(kv.HLcf) == 1
    HLohcdB = kv.HLohcdB0 * ones(length(CF),1);
    HLihcdB = kv.HLihcdB0 * ones(length(CF),1);
else
    HLohcdB = interp1(kv.HLcf,kv.HLohcdB0,CF,'linear','extrap');
    HLihcdB = interp1(kv.HLcf,kv.HLihcdB0,CF,'linear','extrap');
end
HLohcdB = max(HLohcdB,0);
HLihcdB = max(HLihcdB,0);


%%%%%%%%%%%%%%%%%%%%% step1: outer ear correction

if strcmp(flags.outerearcorrection,'FreeField')
    freefield_F = [0 20 25 31.5 40 50 63 80 100 ...
    125 160 200 250 315 400 500 630 ...
    750 800 1000 1250 1500 1600 2000 2500 ...
    3000 3150 4000 5000 6000 6300 8000 9000 ...
    10000 11200 12500 14000 15000 16000 20000];
    freefield_dB = [0 0 0 0 0 0 0 0 0 ...
    0.1 0.3 0.5 0.9 1.4 1.6 1.7 2.5 ...
    2.7 2.6 2.6 3.2 5.2 6.6 12 16.8 ...
    15.3 15.2 14.2 10.7 7.1 6.4 1.8 -0.9 ...
    -1.6 1.9 4.9 2 -2 2.5 2.5];

    inputLdB = inputLdB + interp1(freefield_F,freefield_dB,inputF,'linear','extrap');
    
elseif strcmp(flags.outerearcorrection,'PDR10')
    PDR10_F =  [100  200  300  500 700  1000 1500 2000 2500 3000 3500 4000 5000 6000 7000];
    PDR10_dB = [-19.0 -9.5 -4.5 0 0.5 0 1.5  3.0  3.0  3.0  4.0  3.0  1.5  -3.0  -28.0];
    
    inputLdB = inputLdB + interp1(PDR10_F,PDR10_dB,inputF,'linear','extrap');
    
elseif strcmp(flags.outerearcorrection,'Eardrum')
else
    error('No such correction');
end

%%%%%%%%%%%%%%%%%%%%% step2: middle ear correction
MidEar_F = [20 25 31.5 40 50 63 80 100 125 ...
            160 200 250 315 400 500 630 750 800 1000 ...
            1250 1500 1600 2000 2500 3000 3150 4000 5000 6000 ...
            6300 8000 9000 10000 11200 12500 14000 15000 16000 20000];
MidEar_dB = [-33.2 -28.2 -23.2 -19.4 -16.3 -13.3 -10.2 -8.0 -6.1 ...
            -4.7 -3.5 -2.8 -2.4 -1.9 -1.8 -2.1 -2.5 -2.3 -2.6 ...
            -3.7 -5.5 -6.7 -11.4 -14.5 -11.5 -11.0 -10.5 -10.8 -12.8 ...
            -13.6 -16.5 -15.8 -15.0 -16.9 -18.8 -20.7 -21.9 -22.3 -24.1];

inputLdB = inputLdB + interp1(MidEar_F,MidEar_dB,inputF,'linear','extrap');

inputL = 10.^(inputLdB/10);                 % Intensity in linear

%%%%%%%%%%%%%%%%%%%%% step3: passive filter
tl = CF./(0.1084*CF+2.3301);
tu = 15.0 * ones(length(CF),1);

E_pf = zeros(length(CF),1);
for i = 1:length(CF)
    g = inputF/CF(i)-1;
    indexl = find(g<0);
    gl = abs(g(indexl));
    E_pf(i) = sum((1+gl*tl(i)).*exp(-gl*tl(i)).*inputL(indexl));
    indexu = find(g>=0);
    gu = g(indexu);
    E_pf(i) = E_pf(i) + sum((1+gu*tu(i)).*exp(-gu*tu(i)).*inputL(indexu));
end
E_pf = max(E_pf,10^(-10));
EdB_pf = 10*log10(E_pf);

EdB_pf = max(EdB_pf,0);    

%%%%%%%%%%%%%%%%%%%% step4: gain decided by passive input
GdBmax = CF./(0.0191*CF+1.1) - HLohcdB;

GdB = GdBmax.*( 1 - 1./(1+exp(-0.05*(EdB_pf-(100-GdBmax)))) + 1./(1+exp(-0.05*(0-(100-GdBmax)))));   
index = find(EdB_pf>30);
GdB(index) = GdB(index) - 0.003*(EdB_pf(index)-30).^2;

GdB = min(max(GdB,-20),GdBmax);
G = 10.^(GdB/10); 

%%%%%%%%%%%%%%%%%%%% step5: active tip filter (af)
pl = CF./(0.0272*CF+5.4365);
pu = 27.9 * ones(length(CF),1);

E_af = zeros(length(CF),1);
for i = 1:length(CF)
    g = inputF/CF(i)-1;
    indexl = find(g<0);
    gl = abs(g(indexl));
    E_af(i) = G(i) * sum((1+gl*pl(i)).*exp(-gl*pl(i)).*inputL(indexl));
    indexu = find(g>=0);
    gu = g(indexu);
    E_af(i) = E_af(i) + G(i) * sum((1+gu*pu(i)).*exp(-gu*pu(i)).*inputL(indexu));
end

E = E_pf + E_af;
E= max(E,10^(-10));
EdB = 10*log10(E);

EdB = EdB- HLihcdB.*(1-0.5./(1+exp(-0.2*((EdB-52)-(HLihcdB+20))))); 
E = 10.^(EdB/10);

Ldn = sum(E) * kv.cambin * 1.525*1e-8;


