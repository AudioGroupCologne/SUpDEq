function [car_out, state] = lyon2011_carstep(x_in, CAR_coeffs, state)
%LYON2011_CARSTEP One sample-time update step for the filter part of the model
%
%   Usage: [car_out, state] = lyon2011_carstep(x_in, CAR_coeffs, state);
%
%   Input parameters:
%     x_in          : The CF struct holds the filterbank design and 
%                     state; if you want to break the input up into
%                     segments, you need to use the updated CF
%                     to keep the state between segments.
%     CAR_coeffs           : input_waves is a column vector if there's just one
%                            audio channel; more generally, it has a row per 
%                            time sample, a column per audio channel. The 
%                            input_waves are assumed to be sampled at the 
%                            same rate as the CARFAC is designed for. 
%                            A resampling may be needed before calling this.
%     state            : Plot automatic gain control figure. Default is 0.
%
%   Output parameters:
%     car_out           : The CF struct holds the filterbank design and 
%                         state; if you want to break the input up into
%                         segments, you need to use the updated CF
%                         to keep the state between segments.
%     state           : decim_naps is like naps but time-decimated by 
%                       the int CF.decimation.
%
%   LYON2011_CARSTEP updates the filter cascade state taking into account
%   outer hair cell feedback. This 
%   file is part of an implementation of Lyon's cochlear model:
%   "Cascade of Asymmetric Resonators with Fast-Acting Compression"
%
%   See also:  lyon2011_agcstep lyon2011_carstep
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lyon2011_carstep.php


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

%   #Author: Amin Saremi (2016) adaptations for the AMT (based on <https://github.com/google/carfac>, Richard F. Lyon)
%   #Author: Clara Hollomey (2021) adaptation for the AMT 1.0
%   #License: gpl3


% Most of the update is parallel; finally we ripple inputs at the end.
ohcFeedback = 1;

% do the DOHC stuff:
g = state.g_memory + state.dg_memory;  % interp g
zB = state.zB_memory + state.dzB_memory; % AGC interpolation state
% update the nonlinear function of "velocity", and zA (delay of z2):
zA = state.zA_memory;
v = state.z2_memory - zA;

%if ohcFeedback % widen v with feedback
%   nlf = 1 ./ (1 + ...
%  ((v. * widen) * CAR_coeffs.velocity_scale + CAR_coeffs.v_offset) .^ 2 );
% else
    nlf = 1 ./ (1 + ...
  (v * CAR_coeffs.velocity_scale + CAR_coeffs.v_offset) .^ 2 );
%endif

% zB * nfl is "undamping" delta r:
r = CAR_coeffs.r1_coeffs + zB .* nlf; 
zA = state.z2_memory;

% now reduce state by r and rotate with the fixed cos/sin coeffs:
z1 = r .* (CAR_coeffs.a0_coeffs .* state.z1_memory - ...
  CAR_coeffs.c0_coeffs .* state.z2_memory);
% z1 = z1 + inputs;
z2 = r .* (CAR_coeffs.c0_coeffs .* state.z1_memory + ...
  CAR_coeffs.a0_coeffs .* state.z2_memory);

zY = CAR_coeffs.h_coeffs .* z2;  % partial output

% Ripple input-output path, instead of parallel, to avoid delay...
% this is the only part that doesn't get computed "in parallel":
in_out = x_in;
for ch = 1:length(zY)
  % could do this here, or later in parallel:
  z1(ch) = z1(ch) + in_out;
  % ripple, saving final channel outputs in zY
  in_out = g(ch) * (in_out + zY(ch));
  zY(ch) = in_out;
end

% put new state back in place of old
% (z1 is a genuine temp; the others can update by reference in C)
state.z1_memory = z1;
state.z2_memory = z2;
state.zA_memory = zA;
state.zB_memory = zB;
state.zY_memory = zY;
state.g_memory = g;

car_out = zY;




