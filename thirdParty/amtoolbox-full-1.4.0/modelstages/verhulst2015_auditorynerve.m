function Pfiring = verhulst2015_auditorynerve(Vm,fs,kSR,kmax_or_cf)
%VERHULST2015_AUDITORYNERVE Auditory nerve models used by Verhulst et al. 2015
%
%   Auditory nerve models used in the model Verhulst et al. 2015 based on the 
%   three-store diffusion model of Westerman and Smith (1988).
%
%   License:
%   --------
%
%   This model is licensed under the UGent Academic License. Further usage details are provided 
%   in the UGent Academic License which can be found in the AMT directory "licences" and at 
%   <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.
%
%   References:
%     S. Verhulst, H. Bharadwaj, G. Mehraei, C. Shera, and
%     B. Shinn-Cunningham. Functional modeling of the human auditory
%     brainstem response to broadband stimulation. jasa, 138(3):1637--1659,
%     2015.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/verhulst2015_auditorynerve.php


%   #License: ugent
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal PYTHON C
%   #Author: Alejandro Osses (2020): primary implementation based on https://github.com/HearingTechnology/Verhulstetal2018Model
%   #Author: Piotr Majdak (2021): adaptations for the AMT 1.0

% This file is licenced under the terms of the UGent Academic License, which details can be found in the AMT directory "licences" and at <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.
% For non-commercial academic research, you can use this file and/or modify it under the terms of that license. This file is distributed without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. 


dt=1/fs;

%%% 1. Parameters
% 1.1 Fitted parameters:
if nargin < 3; kSR = 68.5;  end % Spontanous rate [spikes/s], param for HSR
if nargin < 4; kmax_or_cf = 3000; end % Peak exocytosis rate [spikes/s], param for HSR


% Memory allocation:
[N_samples,N_ch] = size(Vm);
Pfiring=zeros(N_samples,N_ch); % Memory allocation

% def anf_model(Vm):
cf = transpose(kmax_or_cf(:)); % makes sure 'cf' is a row array

Vth     = 50e-6; % [V] - Vihc threshold below which the AN remains at SR, for HSR
Vth_SR  = 2e-3 + Vth; % [V] - slightly more than 2e-3. Vihc at ~40 dB for LSR, MSR
Vsat_max= 1e-3; % [V] - Vihc yielding maximum PI

kSR = kSR/(1-0.75e-3*kSR)*ones(1,N_ch); % to compensate for the division at the end of the processing chain
A_SS = 150+(cf/100.0); % Only frequency dependence here, comes from Liberman1978 */

TauR = 2e-3;   % Rapid Time Constant eq.10 
TauS = 60e-3;  % Short Time Constant eq.10 
A_RS = kSR;     % Ratio of A_slow and A_fast

PTS=1+(6*kSR/(6+kSR));

AR  = (A_RS/(1+A_RS))*(PTS*A_SS-A_SS);
AST = (1./(1+A_RS)).*(PTS*A_SS-A_SS);
PI1 = kSR.*(PTS*A_SS-kSR)./(PTS*A_SS.*(1-kSR./A_SS)); % Verhulst2015, Eq. 9
PI2 = (PTS*A_SS-kSR)./(1-kSR./A_SS); % Verhulst2015, Eq. 10. From original Westerman
CG = 1;

gamma1 = CG./kSR;
gamma2 = CG./A_SS;
k1     = -1/TauR;
k2     = -1/TauS;

nume = (1-((PTS*A_SS)./kSR));
deno = gamma1.*((AR.*(k1-k2)./(CG*PI2))+(k2./(PI1.*gamma1))-(k2./(PI2.*gamma2)));
VI0=nume./deno;

% Same numerator, but slightly different denominator:
deno = gamma1.*((AST*(k2-k1)./(CG*PI2))+(k1./(PI1.*gamma1))-(k1./(PI2.*gamma2)));
VI1 = nume./deno;

VI=(VI0+VI1)/2;

alpha=(CG*TauR*TauS)./A_SS;
beta=(1/TauS+1/TauR)*alpha;
theta1=(alpha.*PI2)./VI;
theta2=VI./PI2;
theta3=1./A_SS-1./PI2;

PL=(((beta-theta2.*theta3)./theta1)-1).*PI2;
PG=1./(theta3-1./PL);
VL=theta1.*PL.*PG;
CI=kSR./PI1; % CI, or 'q(t)'- concentration of synaptic neurotransmitters in the immediate store
CL=CI.*(PI1+PL)./PL;

for i = 1:N_samples
	 PI=((PI2-PI1)./(Vsat_max)).*(Vm(i,:)-(Vth_SR./exp(kSR)))+PI1; % Eq. 8
	 
	 idx_PI_rest=find(Vm(i,:)<(Vth+(Vth_SR./exp(kSR))));
	 PI(idx_PI_rest)=PI1(idx_PI_rest); % Eq. 7

	 CI_prev = CI;
	 CI = CI + (dt./VI).*(-PI.*CI + PL.*(CL-CI));
	 CL = CL + (dt./VL).*(-PL.*(CL - CI_prev) + PG.*(CG - CL));
	 temp = 1./PG+1./PL+1./PI;
	 CI_n=find(CI<0);
	 CI(CI_n) = CG./(PI(CI_n).*temp(CI_n));

	 CL(CI_n) = CI(CI_n).*(PI(CI_n)+PL(CI_n))./PL(CI_n);
	 Pfiring(i,:) =PI.*CI; % Verhulst2015, Eq. 6
end


