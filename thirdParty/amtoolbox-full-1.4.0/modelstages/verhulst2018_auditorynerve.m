function Pfiring = verhulst2018_auditorynerve(Vm,fs,kSR,kmax_or_cf,version_year)
%VERHULST2018_AUDITORYNERVE Auditory nerve models used by Verhulst et al. 2018 and 2015
%
%   Auditory nerve models used in the models by Verhulst et al. 2018 (if
%   version_year == 2018) or Verhulst et al. 2015 (if version_year == 2015).
%   Both versions are based on the three-store diffusion model of Westerman
%   and Smith (1988).
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
%     S. Verhulst, A. Altoè, and V. Vasilkov. Functional modeling of the
%     human auditory brainstem response to broadband stimulation.
%     hearingresearch, 360:55--75, 2018.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/verhulst2018_auditorynerve.php


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
%   #Author: Alejandro Osses (2020): primary implementation based on https://github.com/HearingTechnology/Verhulstetal2018Model
%   #Author: Piotr Majdak (2021): adaptations for the AMT 1.0


dt=1/fs;

%%% 1. Parameters
% 1.1 Fitted parameters:
if nargin < 3; kSR = 68.5;  end % Spontanous rate [spikes/s], param for HSR
if nargin < 4; kmax_or_cf = 3000; end % Peak exocytosis rate [spikes/s], param for HSR
if nargin < 5; version_year = 2018; end

% Memory allocation:
[N_samples,N_ch] = size(Vm);
Pfiring=zeros(N_samples,N_ch); % Memory allocation

switch version_year
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2018
        kmax = kmax_or_cf;
        alpha_l = 220; % reserve pool max. replenishment rate [spikes/s]

        L = 60;      % [Nr.] RP size: Max. vesicles in the second pool
        s = 1.5e-3;  % [mV]  Sensitivity of the Boltzmann function relating V_IHC and driven exocytosis rate

        % 1.2 Parameters from literature:
        alpha_q = 700;    % RRP maximum replenishment rate [spikes/s] (Pangrisc 2010,Chapocnikov 2014)
        M       = 14;     % Max. vesicles in the ready release pool or release sites (Meyer 2009). With the paramaters is 250 release/s, because of refractoriness the steady state spike rate goes around 200 spike/s
        tau_r   = 0.6e-3; % Time constant of relative refractoriness (Peterson and Heil 2014)
        tau_a   = 0.6e-3; % Time constant of absolute refractoriness (Peterson and Heil 2014)
        tau_Ca  = 0.2e-3; % Verhulst2018a, in Eq. (11). Time constant of the Ca^2+ channel (Johnson and Marcotti 2008)
        V_rest  = -0.05703; % resting_potential at equilibrium, from IHC model
        % peak_potential=-0.04; % peak resting potential at 100 dB 4 kHz (where nerve fibers saturate)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Memory allocation:
        eta= nan(N_samples,N_ch);
        k  = nan(N_samples,N_ch);
        qt = nan(N_samples,N_ch);
        lt = nan(N_samples,N_ch);
        rel_refract=nan(N_samples,N_ch);
        available=nan(N_samples,N_ch);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% Steady-state values:
        L_steady = L*(1         -kSR/alpha_l); % Initial Nr. of vescicles in RP 
        M_steady = M*(L_steady/L-kSR/alpha_q); % Initial Nr. of vescicles in RRP 
        %%% End of Steady-state values

        alpha_q_dt = alpha_q*dt; % [spikes per unit of time]
        alpha_l_dt = alpha_l*dt; % [spikes per unit of time]

        % Conversion from time constants
        alpha_a  = exp(-dt/tau_a);
        alpha_r  = exp(-dt/tau_r);
        alpha_Ca = exp(-dt/tau_Ca);

        % Calculations from the Equations in Verhulst2018a:
        V05_SR  = log((kmax-kSR)/kSR)*s+V_rest; % Eq. (13)
        eta_inf = sqrt(1./(1+exp(-(Vm-V05_SR)/s))); % Eq. (11) 

        %%% parameters for refractoriness

        rel_refract0 = 0; % how much the firing probability decreases due to relative refractoriness
        available0   = 1.0; % number of fibers not in a refractory state
        dly_a   = round(tau_a*fs); % length of buffer to store the firing history (in order to account for absolute refractoriness)
        buf_ref = zeros([dly_a N_ch]); % buffer to store the history of firing (1001, 60)
        idx_buf=1; 
        pp=kmax/M_steady; % multiplier relating the activation nonlinearity (between 0 and 1) with actual firing rate
        % parameters to relate Vm with firing rate

        eta0=sqrt(1/(1+exp(-(V_rest-V05_SR)/s)))+zeros([1 N_ch]); % driven exocytosis rate at rest
        %           take the square root, filter it with a first order filter and then square it. This is 
        %           equivalent to a second order activation of the ion channels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dl_description = 'replenishReserve';
        dr_description = 'replenishRRP';

        kt0 = pp*(eta0.^2); % vesicleReleaseRate, does not change
        kdt = kt0*dt;
        zero_time=round(50e-3*fs);
        %%% Filling in the buffer:
        qt_here = M_steady; % initial value
        lt_here = L_steady;
        for i=1:zero_time
            ydt=kdt.*qt_here; % 'ejected'
            dl = alpha_l_dt*(L-lt_here)/L; % dl - replenishReserve
            dr = alpha_q_dt*max(lt_here/L-qt_here/M,0); % dr - replenishRRP
            qt(i,:) = qt_here + dr - ydt;
            lt(i,:) = lt_here - dr + dl;
            
            firing0       = (available0-rel_refract0).*ydt;
            Pfiring_dly_a = buf_ref(idx_buf,:); % Pfiring_dly_a: 'recovered'
            rel_refract0 = alpha_r*rel_refract0+(1-alpha_a)*Pfiring_dly_a;
            available0   = available0-firing0+Pfiring_dly_a;
            buf_ref(idx_buf,:)=firing0;
            idx_buf=mod(idx_buf+1,dly_a); 
            if idx_buf == 0
                idx_buf = dly_a;
            end
        end

        for i=1:N_samples
            
            % suffix 'prev' stands for 'previous':
            if i ~= 1
                eta_prev = eta(i-1,:);
                qt_prev  = qt(i-1,:);
                lt_prev  = lt(i-1,:);
                available_prev = available(i-1,:);
                rel_refract_prev = rel_refract(i-1);
            else
                % Using initial values if i==1:
                eta_prev = eta0;
                qt_prev  = M_steady;
                lt_prev  = L_steady;
                available_prev = available0;
                rel_refract_prev = rel_refract0;
            end
            
            eta(i,:) = alpha_Ca*eta_prev+(1-alpha_Ca)*eta_inf(i,:); % Analytical solution of Eq. (11)
            k(i,:)   = pp*(eta(i,:).^2); % yt - effective exocytosis rate
            ydt = k(i,:).*qt_prev*dt; % Eq. (14.1) -- kdt: Probability that one vescicle is released, text after Eq. (14)
            dr  = alpha_q_dt*max(lt_prev/L-qt_prev/M,0); % Eq. (14.3), but derivated
            dl  = alpha_l_dt*(L-lt_prev)/L; % Eq. (14.4)
            
            qt(i,:) = qt_prev + dr - ydt; % Eq. (14.2), populating the ready-releasable pool (RRP) of neurotransmitter
            lt(i,:) = lt_prev - dr + dl;  % populating the reserve pool (RP) of neurotransmitter
            P_Rel   = available_prev-rel_refract_prev; % rel_refract: relative refractoriness
            Pfiring(i,:) = P_Rel.*ydt;
            Pfiring_dly_a(i,:)=buf_ref(idx_buf,:); % Contains 'Pfiring(i-1)'
            rel_refract(i,:)  = alpha_a*rel_refract_prev+(1-alpha_r)*Pfiring_dly_a(i,:);
            rel_refract(i,:)  = alpha_a*rel_refract_prev+(1-alpha_r)*Pfiring_dly_a(i,:);
            available(i,:)    = available_prev-(Pfiring(i,:)-Pfiring_dly_a(i,:));
            
            buf_ref(idx_buf,:)=Pfiring(i,:); % overwrites already used buffer with sample to be used later
            idx_buf = mod(idx_buf+1,dly_a);
            
            if idx_buf == 0
                idx_buf = dly_a;
            end
        end

        Pfiring = Pfiring/dt; % AO, Verhulst2018a, Eq. (16)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2015
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
        
        disp('')
end


