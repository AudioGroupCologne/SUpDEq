function CF = lyon2011_closeagcloop(CF)
%LYON2011_CLOSEAGCLOOP active gain control loop
%   Usage: [CF, decim_naps, naps, BM, ohc, agc] = lyon2011_closeagcloop(CF,input_waves, AGC_plot_fig_num, open_loop);
%
%
%   Input parameters:
%     CF              : The CF struct holds the filterbank design and
%                       state; if you want to break the input up into
%                       segments, you need to use the updated CF
%                       to keep the state between segments.
%
%   Output parameters:
%     CF              : processed CF struct
%
%
%
%   This file is part of an implementation of Lyon's cochlear model:
%   "Cascade of Asymmetric Resonators with Fast-Acting Compression"
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lyon2011_closeagcloop.php


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

% fastest decimated rate determines interp needed:
decim1 = CF.AGC_params.decimation(1);

for ear = 1:CF.n_ears
  undamping = 1 - CF.ears(ear).AGC_state(1).AGC_memory; % stage 1 result
  % degrade the OHC active undamping if the ear is less than healthy:
  undamping = undamping .* CF.ears(ear).CAR_coeffs.OHC_health;
  % Update the target stage gain for the new damping:
  new_g = lyon2011_stageg(CF.ears(ear).CAR_coeffs, undamping);
  % Return the stage gain g needed to get unity gain at DC

  r1 = CF.ears(ear).CAR_coeffs.r1_coeffs;  % at max damping
  a0 = CF.ears(ear).CAR_coeffs.a0_coeffs;
  c0 = CF.ears(ear).CAR_coeffs.c0_coeffs;
  h  = CF.ears(ear).CAR_coeffs.h_coeffs;
  zr = CF.ears(ear).CAR_coeffs.zr_coeffs;
  r  = r1 + zr .* undamping;
  g  = (1 - 2*r.*a0 + r.^2) ./ (1 - 2*r.*a0 + h.*r.*c0 + r.^2);

  % set the deltas needed to get to the new damping:
  CF.ears(ear).CAR_state.dzB_memory = ...
    (CF.ears(ear).CAR_coeffs.zr_coeffs .* undamping - ...
    CF.ears(ear).CAR_state.zB_memory) / decim1;
  CF.ears(ear).CAR_state.dg_memory = ...
    (new_g - CF.ears(ear).CAR_state.g_memory) / decim1;
end

