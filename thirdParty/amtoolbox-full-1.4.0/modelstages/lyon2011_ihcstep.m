function [ihc_out, state] = lyon2011_ihcstep(filters_out, coeffs, state)
%LYON2011_IHCSTEP update of inner-hair-cell (IHC) model
%
%   Usage: [ihc_out, state] = lyon2011_ihcstep(filters_out, coeffs, state);
%
%   LYON2011_IHCSTEP updates the inner-hair-cell (IHC) model, including the
%   detection nonlinearity and one or two capacitor state variables. This 
%   file is part of an implementation of Lyon's cochlear model:
%   "Cascade of Asymmetric Resonators with Fast-Acting Compression"
%
%
%   See also:   lyon2011_agcstep lyon2011_carstep
%               lyon2011_closeagcloop lyon2011_design
%               lyon2011_ihcstep lyon2011_init
%               lyon2011_spatialsmooth
%               demo_lyon2011
%
%   References:
%     R. F. Lyon. Cascades of two-pole–two-zero asymmetric resonators are
%     good models of peripheral auditory function. J. Acoust. Soc. Am.,
%     130(6), 2011.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lyon2011_ihcstep.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #License: Apache2
%   #Author: Richard F. Lyon (2013): original implementation (https://github.com/google/carfac)
%   #Author: Amin Saremi (2016): adaptations for the AMT
%   #Author: Clara Hollomey (2021): integration in the AMT 1.0
%   #Author: Richard Lyon (2022): bug fixes for AMT
%   #Author: Mihajlo Velimirovic (2022): implementation of the option ihc_potential

% This file is licensed unter the Apache License Version 2.0 which details can 
% be found in the AMT directory "licences" and at 
% <http://www.apache.org/licenses/LICENSE-2.0>. 
% You must not use this file except in compliance with the Apache License 
% Version 2.0. Unless required by applicable law or agreed to in writing, this 
% file is distributed on an "as is" basis, without warranties or conditions 
% of any kind, either express or implied.

%   #Author: Amin Saremi (2016) adaptations for the AMT (based on
%     <https://github.com/google/carfac>, Richard F. Lyon)
%   #Author: Richard Lyon (2023) further adaptations and new two-cap model.


% AC couple the filters_out, with 20 Hz corner
ac_diff = filters_out - state.ac_coupler;
state.ac_coupler = state.ac_coupler + coeffs.ac_coeff * ac_diff;

if coeffs.just_hwr
  ihc_out = min(2, max(0, ac_diff));  % limit it for stability
else
  %conductance = lyon2011_detect(ac_diff);
  % An IHC-like sigmoidal detection nonlinearity for the CARFAC.
  % Resulting conductance is in about [0...1.3405]
  a = 0.175;   % offset of low-end tail into neg x territory
  % this parameter is adjusted for the book, to make the 20% DC
  % response threshold at 0.1

  set = ac_diff > -a;
  z = ac_diff(set) + a;

  % Detection nonlinearity; zero is the final answer for many points:
  conductance = zeros(size(ac_diff));
  conductance(set) = z.^3 ./ (z.^3 + z.^2 + 0.1);

  if coeffs.one_cap;
    % Output comes from receptor current like in Hall and Allen's models.
    ihc_out = conductance .* state.cap_voltage;
    state.cap_voltage = state.cap_voltage - ...
      ihc_out .* coeffs.out_rate + ...
      (1 - state.cap_voltage) .* coeffs.in_rate;
    ihc_out = ihc_out * coeffs.output_gain;
    % Smooth it twice with LPF:
    state.lpf1_state = state.lpf1_state + coeffs.lpf_coeff * ...
      (ihc_out - state.lpf1_state);
    state.lpf2_state = state.lpf2_state + coeffs.lpf_coeff * ...
      (state.lpf1_state - state.lpf2_state);
    ihc_out = state.lpf2_state - coeffs.rest_output;
  else
    % Change to 2-cap version mediated by receptor potential at cap1:
    % Geisler book fig 8.4 suggests 40 to 800 Hz corner.
    receptor_current = conductance .* state.cap1_voltage;
    % "out" means charge depletion; "in" means restoration toward 1.
    state.cap1_voltage = state.cap1_voltage - ...
      receptor_current .* coeffs.out1_rate + ...
      (1 - state.cap1_voltage) .* coeffs.in1_rate;
    % Amount of depletion below 1 is receptor potential.
    receptor_potential = 1 - state.cap1_voltage;  % Already smooth.
    % Identity map from receptor potential to neurotransmitter conductance.
    ihc_out = receptor_potential .* state.cap2_voltage;  % Now a current.
    % cap2 represents transmitter store; another adaptive gain.
    % Deplete the transmitter store like in Meddis models:
    state.cap2_voltage = state.cap2_voltage - ...
      ihc_out .* coeffs.out2_rate + ...
      (1 - state.cap2_voltage) .* coeffs.in2_rate;
    % Awkwardly, gain needs to be here for the init states to be right.
    ihc_out = ihc_out * coeffs.output_gain;
    % smooth once more with LPF (receptor potential was already smooth):
    state.lpf1_state = state.lpf1_state + coeffs.lpf_coeff * ...
      (ihc_out - state.lpf1_state);
    ihc_out = state.lpf1_state - coeffs.rest_output;
  end
end

state.ihc_accum = state.ihc_accum + ihc_out;  % for where decimated output is useful



