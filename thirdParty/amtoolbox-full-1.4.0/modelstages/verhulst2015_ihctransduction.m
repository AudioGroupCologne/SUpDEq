function Vihc=verhulst2015_ihctransduction(Vin,fs,keyvals)
%VERHULST2015_IHCTRANSDUCTION ihc transduction in Verhulst et al. 2015
%
%   Inner hair cell transduction as used in Verhulst et al. models. 
%   Vihc represents the inner-hair-cell membrane potential in Volts.
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/verhulst2015_ihctransduction.php


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
    keyvals.ihc_scal_constant = [];
end

[N_samples,N_ch] = size(Vin); % 'size' in Python code

% Memory allocation:
Vihc = zeros([N_samples N_ch]);

if ~isempty(keyvals.ihc_scal_constant)
    scal_constant = keyvals.ihc_scal_constant;
else
	VBMmax =  41e-6; % m/s, former variable name: Mvel
	YBmax  = 200e-9; % m, or 200 nm
	scal_constant = YBmax/VBMmax; % see Verhulst2015, Table I, parameter 'G' = 0.0049
end
Vin = scal_constant*Vin;
    
VihcNF=zeros([N_samples N_ch]); % Memory allocation
% Some constants:
F_LPC=1000;
[b,a] = local_LowPass_coeff(F_LPC,fs);

%%% The following is the same as implemented in 'ol_NLogarithm':
A0=0.008; % %0.1 scalar in IHC nonlinear function
B=2000*6000; % #2000 par in IHC nonlinear function
C=0.33; % #1.74 par in IHC nonlinear function
D=200e-9;

idx_neg=find(Vin<0);  % negative indexes
idx_pos=find(Vin>=0); % positive indexes, null samples are 'unprocessed'

VihcNF(idx_pos)= A0*log(1+B*abs(Vin(idx_pos))); % Eq. 5(A)

% The '3' in the denominator was taken from the original code, but in
%   the paper it seems that that value should be 0.3 instead:
VihcNF(idx_neg)=-A0*(((abs(Vin(idx_neg)).^C)+D)./((3*abs(Vin(idx_neg)).^C)+D)).* ...
					log(1+B*abs(Vin(idx_neg))); % Eq. 5(B)
%%% End of 'ol_NLogarithm'

past_Vin     = zeros(1,N_ch);
past_output1 = zeros(1,N_ch);
past_output2 = zeros(1,N_ch);
for i=1:N_samples % check across which dimension I should go
	%%% Cascade of IIR filters vectorised (all sections at once)
	%       two identical first order filters
	y1       = (-a(2)*past_output1+b(1)*VihcNF(i,:)+ b(2)*past_Vin)/a(1); % intermediate output of the iir cascade
	Vihc(i,:)= (-a(2)*past_output2+b(1)*y1     + b(2)*past_output1)/a(1);
	
	% Update filters' past values
	past_Vin = VihcNF(i,:);
	past_output1=y1;
	past_output2= Vihc(i,:);
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b,a] = local_LowPass_coeff(Fc,fs)
%
% Returns the coefficients of a low-pass filter of first order with a cut-off
% frequency at Fc for a digital sampling frequency fs.

c = 2*fs;
C1LP = ( c - 2*pi*Fc ) / ( c + 2*pi*Fc );
C2LP = 2*pi*Fc / (2*pi*Fc + c);
   
b = [C2LP C2LP]; % filter coefficients (numerator)
a = [1   -C1LP]; % filter coefficients (denominator)

