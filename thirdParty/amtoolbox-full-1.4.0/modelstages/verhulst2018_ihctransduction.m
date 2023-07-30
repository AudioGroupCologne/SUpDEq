function Vihc=verhulst2018_ihctransduction(Vin,fs,version_year,keyvals)
%VERHULST2018_IHCTRANSDUCTION ihc transduction in Verhulst et al. 2018
%
%   Inner hair cell transduction as used in Verhulst et al. models. 
%   Vihc represents the inner-hair-cell membrane potential in Volts.
%
%   The following IHC transduction models are implemented here:
%
%   version_year = 2018 ('verhulst2018b')
%
%   version_year = 2015
% 
%   The 2018 model Vihc is obtained from the interplay between K+ ion channels
%   and mechanoelectrical transduction (MET) currents. This model is described
%   in verhulst2018a (Appendix A) and is based on the description/validation
%   described in altoe2018. The influence of calcium channels is not accounted
%   in this script but it is in verhulst2018_auditorynerve.m, where the effective
%   exocystosis rate is controlled by Ca-related parameters. In this script
%   an effect of the calcium channels could be added by setting a leak current
%   but in the current validation such Ikeak is not active (because Gleak is 0).
%   
%   License:
%   --------
%
%   This model is licensed under the UGent Academic License. Further usage details are provided 
%   in the UGent Academic License which can be found in the AMT directory "licences" and at 
%   <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.
%
%   References:
%     S. Verhulst, A. Altoè, and V. Vasilkov. Functional modeling of the
%     human auditory brainstem response to broadband stimulation.
%     hearingresearch, 360:55--75, 2018.
%     
%     A. Altoè, V. Pulkki, and S. Verhulst. The effects of the activation of
%     the inner-hair-cell basolateral k+ channels on auditory nerve
%     responses. hearingresearch, 364:68--80, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/verhulst2018_ihctransduction.php


%   #License: ugent
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal PYTHON C
%   #Author: Alejandro Osses (2020): primary implementation based on https://github.com/HearingTechnology/Verhulstetal2018Model
%   #Author: Piotr Majdak (2021): adaptations for the AMT 1.0

% This file is licenced under the terms of the UGent Academic License, which details can be found in the AMT directory "licences" and at <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.
% For non-commercial academic research, you can use this file and/or modify it under the terms of that license. This file is distributed without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. 

%   #License: ugent
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal PYTHON C
%   #Author: Alessandro Altoe (2018): original code of the model
%   #Author: Alejandro Osses (2020): primary integration in the AMT
%   #Author: Piotr Majdak (2021): adaptations for the AMT 1.0



if nargin < 4
    keyvals = [];
    keyvals.ihc_scal_constant = [];
end

if nargin < 3
    version_year = 2018;
end

[N_samples,N_ch] = size(Vin); % 'size' in Python code

% Memory allocation:
Vihc = zeros([N_samples N_ch]);

if ~isempty(keyvals.ihc_scal_constant)
    scal_constant = keyvals.ihc_scal_constant;
else
    switch version_year
        case 2018
            scal_constant = 0.118; % [s] ratio to bring BM velocity to IHC cilia displacement range
        case 2015
            VBMmax =  41e-6; % m/s, former variable name: Mvel
            YBmax  = 200e-9; % m, or 200 nm
            scal_constant = YBmax/VBMmax; % see Verhulst2015, Table I, parameter 'G' = 0.0049
    end
end
Vin = scal_constant*Vin;
    
switch version_year
    case 2018
        flag = 'verhulst2018b'; % Calcium channels not implemented here yet. It is like in Verhulst2018b
        %%% Loading parameters
        pars = il_IHC_param(flag); % run(filename)
        tau_kf   = pars.tkf1; % time constant, K, fast
        tau_ks   = pars.tks1; % time constant, K, slow
        tau_Met  = pars.tauMet;
        Gmet_max = pars.Gmet_max;
        x0 = pars.x0;
        s0 = pars.s0;
        Gleak  = pars.Gleak;
        Gk_max = pars.Gk_max;
        vk05 = pars.vk05;
        sk = pars.sk;
        EP = pars.EP;
        Ekf = pars.Ekf;
        Eks = pars.Eks;
        Cm = pars.Cm;
        Eca = pars.Eca;
        s1 = pars.s1;
        
        if Gleak ~= 0
            warning('The Verhulst2018 model has not been validated with a leaky channel');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        resting_potential = -0.05703; % 'Extra parameter'

        % Memory allocation (check if it is needed):
        Vm = zeros([1 N_ch])+resting_potential;
        nMet = zeros([1 N_ch]);
        %%% End of memory allocation

        dt=1/fs;
        alpha_Met = exp(-dt/tau_Met);
        alpha_kf  = exp(-dt/tau_kf);
        alpha_ks  = exp(-dt/tau_ks);
        factor1   = exp(x0/s1);
        factor0   = exp(x0/s0);
        mt0 = 1/(1+factor1*(1+factor0)); % similar to Altoe2018, Eq. A2
        nMet = nMet+mt0;

        factor1 = exp((x0-Vin)/s1);
        factor0 = exp((x0-Vin)/s0);
        nMet_inf= 1./(1+factor1.*(1+factor0)); % Altoe2018, Eq. A2, note change of sign: -(mu-x0) as written is the same as (x0-mu)
                                               % Also as Verhulst2018, Eq. A2
        nk_f = 1./(1+exp(-(Vm-vk05)/sk));
        nk_s = 1./(1+exp(-(Vm-vk05)/sk));

        %%%
        zero_sample = round(fs*50e-3);
        for i = 1:zero_sample
            Imet=(Gmet_max*nMet).*(Vm-EP);
            Ileak=Gleak*(Vm-Eca);
            nk_inf=1./(1+exp(-(Vm-vk05)/sk)); % Eq. A5
            nk_f=(1-alpha_kf)*nk_inf+alpha_kf*nk_f;
            nk_s=(1-alpha_ks)*nk_inf+alpha_ks*nk_s;
            
            Ikf=(Gk_max*nk_f).*(Vm-Ekf); % fast activating basolateral K+ current
            Iks=(Gk_max*nk_s).*(Vm-Eks); % slow activating basolateral K+ current
            dV=-(Ileak+Imet+Ikf+Iks)/Cm;
            Vm=Vm+dV*dt;
        end
        %%% Vm is now 'pre-charged'

        for i = 1:N_samples
            % Equations from Verhulst2018:
            nMet=(1-alpha_Met)*nMet_inf(i,:)+alpha_Met*nMet;
            
            Imet =(Gmet_max*nMet).*(Vm-EP); % Eq. A4
            Ileak=Gleak*(Vm-Eca); % In the Verhulst2018 implementation Ileak should be 0
            
            nk_inf=1./(1+exp(-(Vm-vk05)/sk)); % Eq. A5
            nk_f=(1-alpha_kf)*nk_inf+alpha_kf*nk_f;
            nk_s=(1-alpha_ks)*nk_inf+alpha_ks*nk_s;
            
            Ikf=(Gk_max*nk_f).*(Vm-Ekf); % Eq. A7 - fast activating basolateral K+ current
            Iks=(Gk_max*nk_s).*(Vm-Eks); % Eq. A7 - slow activating basolateral K+ current
            dV=-(Ileak+Imet+Ikf+Iks)/Cm; % Eq. A1 - Imet is an inward current, Ikf and Iks are outward currents
            Vm=Vm+dV*dt;

            Vihc(i,:)=Vm;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2015
        VihcNF=zeros([N_samples N_ch]); % Memory allocation
        % Some constants:
        F_LPC=1000;
        [b,a] = il_LowPass_coeff(F_LPC,fs);
        
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par = il_IHC_param(flag)
% function par = IHC_param(flag)
%
% Parameters according to Altoe2018, Table 1. Parameters that differ from
%   what it is defined in the table:
%             here    paper
%       s1 =  48 nm   35 nm
%     tks1 =  10 ms    8 ms
%     Gkfs =  var     19.8 nS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    flag = 'verhulst2018';
end

par.Cm=12.5e-12; % F (pF), IHC capacitance
par.EP=90e-3; % V (mV), Endocochlear potencial

%%% MET channels:
par.Gmet_max=30e-9; % S (nS), maximum conductance
par.x0=20e-9; % m (nm), offset - Verhulst2018: set to obtain about 13% of the MET channels to be open at rest
par.s0=16e-9; % m (nm), sensitivity
par.s1=par.s0*3; % m (nm), sensitivity, AO: deviates from Alessandro's paper, same as in Verhulst2018b
par.tauMet=50e-6; % s (us), activation time constant

%%% Fast and slow K+ channels:
par.Gk_max= 230e-9; % S (nS), maximum conductance
par.vk05  =-31e-3;  % V (mV), half action potential
par.sk    = 10.5e-3;% V (mV), sensitivity
par.Ekf   =-71e-3;  % V (mV), reversal potential, fast
par.Eks   =-78e-3;  % V (mV), reversal potential, slow
par.tkf1  = 0.3e-3; % s (ms), activation time constant, fast
par.tks1  = 10e-3;  % s (ms), activation time constant, slow, AO: deviates from Alessandro's paper, same as in Verhulst2018b

% Not indicated in the table:
par.Gleak=0e-9; % No leak
par.tkf2=0.1e-3;
par.tks2=2e-3;

%%% Ca2+ Channels:
par.GcaM=4.1e-9; % S (pS), maximum conductance
par.sca=7.5e-3; % V (mV), sensitivity
par.tau_Ca=0.2e-3; % s (ms), activation time constant
switch flag
    case 'verhulst2018b'
        par.xca=-30e-3; % V (mV), half action potential
    case 'altoe'
        par.xca=-25e-3; % V (mV), half action potential
end
par.Eca=45.0e-3; % V (mV), reversal potential

% Not indicated in the table:
switch flag
    case 'verhulst2018b'
        par.tauInS=0.5;
        par.CaInMax=0.4;
        par.xCaIn=-43e-3;
    case 'altoe'
        par.tauInS=1;
        par.CaInMax=0.0;
        par.xCaIn=-42e-3;
end
par.tauInF=50e-3;
par.sCaIn=6e-3;
par.Ileak=0.0e-9; % No leak

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b,a] = il_LowPass_coeff(Fc,fs)
% double [b,a] = il_LowPass_coeff(Fc,fs)
%
% Returns the coefficients of a low-pass filter of first order with a cut-off
% frequency at Fc for a digital sampling frequency fs.

c = 2*fs;
C1LP = ( c - 2*pi*Fc ) / ( c + 2*pi*Fc );
C2LP = 2*pi*Fc / (2*pi*Fc + c);
   
b = [C2LP C2LP]; % filter coefficients (numerator)
a = [1   -C1LP]; % filter coefficients (denominator)


