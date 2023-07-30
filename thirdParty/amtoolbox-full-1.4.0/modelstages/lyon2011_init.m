function CF = lyon2011_init(CF)
%LYON2011_INIT initializes the CARFAC model
%
%   Usage:
%     CF = lyon2011_init(CF)
%
%   LYON2011_INIT allocates the state vector storage in the CF struct for one or more ears of CF
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lyon2011_init.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #License: Apache2
%   #Author: Richard F. Lyon (2013): original implementation (https://github.com/google/carfac)
%   #Author: Amin Saremi (2016): adaptations for the AMT
%   #Author: Clara Hollomey (2021): integration in the AMT 1.0
%   #Author: Richard Lyon (2022): bug fixes for AMT
%   #Author: Mihajlo Velimirovic (2022): implementation of the option ihc_potential

%   This file is licensed unter the Apache License Version 2.0 which details can 
%   be found in the AMT directory "licences" and at 
%   <http://www.apache.org/licenses/LICENSE-2.0>. 
%   You must not use this file except in compliance with the Apache License 
%   Version 2.0. Unless required by applicable law or agreed to in writing, this 
%   file is distributed on an "as is" basis, without warranties or conditions 
%   of any kind, either express or implied.



n_ears = CF.n_ears;

for ear = 1:n_ears
  % for now there's only one coeffs, not one per ear
  CF.ears(ear).CAR_state = carinitstate(CF.ears(ear).CAR_coeffs);
  CF.ears(ear).IHC_state = ihcinitstate(CF.ears(ear).IHC_coeffs);
  CF.ears(ear).AGC_state = agcinitstate(CF.ears(ear).AGC_coeffs);
end


function state = carinitstate(coeffs)
n_ch = coeffs.n_ch;
state = struct( ...
  'z1_memory', zeros(n_ch, 1), ...
  'z2_memory', zeros(n_ch, 1), ...
  'zA_memory', zeros(n_ch, 1), ...
  'zB_memory', coeffs.zr_coeffs, ...
  'dzB_memory', zeros(n_ch, 1), ...
  'zY_memory', zeros(n_ch, 1), ...
  'g_memory', coeffs.g0_coeffs, ...
  'dg_memory', zeros(n_ch, 1) ...
  );


function state = agcinitstate(coeffs)
n_ch = coeffs(1).n_ch;
n_AGC_stages = coeffs(1).n_AGC_stages;
state = struct([]);
for stage = 1:n_AGC_stages
  % Initialize state recursively...
  state(stage).AGC_memory = zeros(n_ch, 1);
  state(stage).input_accum = zeros(n_ch, 1);
  state(stage).decim_phase = 0;  % integer decimator phase
end


function state = ihcinitstate(coeffs)
n_ch = coeffs.n_ch;
if coeffs.just_hwr
  state = struct('ihc_accum', zeros(n_ch, 1), ...
                 'ac_coupler', zeros(n_ch, 1));
else
  if coeffs.one_cap
    state = struct( ...
      'ihc_accum', zeros(n_ch, 1), ...
      'cap_voltage', coeffs.rest_cap * ones(n_ch, 1), ...
      'lpf1_state', coeffs.rest_output * ones(n_ch, 1), ...
      'lpf2_state', coeffs.rest_output * ones(n_ch, 1), ...
      'ac_coupler', zeros(n_ch, 1) ...
      );
  else
    state = struct( ...
      'ihc_accum', zeros(n_ch, 1), ...
      'cap1_voltage', coeffs.rest_cap1 * ones(n_ch, 1), ...
      'cap2_voltage', coeffs.rest_cap2 * ones(n_ch, 1), ...
      'lpf1_state', coeffs.rest_output * ones(n_ch, 1), ...
      'ac_coupler', zeros(n_ch, 1) ...
      );
  end
end



