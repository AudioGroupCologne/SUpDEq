function CF = lyon2011_design(n_ears, fs, CF_CAR_params, CF_AGC_params, CF_IHC_params)
%LYON2011_DESIGN computes all the coefficients needed to run Lyon's CARFAC model
%
%   Usage:
%     CF = lyon2011_design(n_ears, fs, CF_CAR_params, CF_AGC_params, CF_IHC_params)
%
%
%   Input parameters:
%     n_ears        : number of input signals
%     fs            : sampling frequency [Hz]
%     CF_CAR_params : struct, pole-zero filter cascade parameters
%     CF_AGC_params : struct, automatic gain control parameters
%     CF_IHC_params : struct, inner hair cell parameters
%
%   Output parameters:
%     CF            : filter coefficients
%
%   References:
%     R. F. Lyon. Cascades of two-pole–two-zero asymmetric resonators are
%     good models of peripheral auditory function. J. Acoust. Soc. Am.,
%     130(6), 2011.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lyon2011_design.php


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


if nargin < 1
  n_ears = 1;  % if more than 1, make them identical channels;
  % then modify the design if necessary for different reasons
end

if nargin < 2
  fs = 22050;
end

if nargin < 3
  CF_CAR_params = struct( ...
    'velocity_scale', 0.1, ...  % for the velocity nonlinearity
    'v_offset', 0.04, ...  % offset gives a quadratic part
    'min_zeta', 0.10, ... % minimum damping factor in mid-freq channels
    'max_zeta', 0.35, ... % maximum damping factor in mid-freq channels
    'first_pole_theta', 0.85*pi, ...
    'zero_ratio', sqrt(2), ... % how far zero is above pole
    'high_f_damping_compression', 0.5, ... % 0 to 1 to compress zeta
    'ERB_per_step', 0.5, ... % assume G&M's ERB formula
    'min_pole_Hz', 30, ...
    'ERB_break_freq', 165.3, ...  % 165.3 is Greenwood map's break freq.
    'ERB_Q', 1000/(24.7*4.37));  % Glasberg and Moore's high-cf ratio
end

if nargin < 4
  CF_AGC_params = struct( ...
    'n_stages', 4, ...
    'time_constants', 0.002 * 4.^(0:3), ...
    'AGC_stage_gain', 2, ...  % gain from each stage to next slower stage
    'decimation', [8, 2, 2, 2], ...  % how often to update the AGC states
    'AGC1_scales', 1.0 * sqrt(2).^(0:3), ...   % in units of channels
    'AGC2_scales', 1.65 * sqrt(2).^(0:3), ... % spread more toward base
    'AGC_mix_coeff', 0.5);
end

if nargin < 5
  % HACK: these constants control the defaults
  one_cap = 1;         % bool; 1 for Allen model, as text states we use
  just_hwr = 0;        % bool; 0 for normal/fancy IHC; 1 for HWR

  CF_IHC_params = struct( ...
    'just_hwr', just_hwr, ...  % not just a simple HWR
    'one_cap', one_cap, ...    % bool; 0 for new two-cap hack
    'ac_corner_Hz', 20, ...    % AC couple at 20 Hz corner
    'tau_lpf', 0.000080, ...   % 80 microseconds smoothing twice
    'tau_out', 0.0005, ...     % depletion tau is pretty fast
    'tau_in', 0.010, ...       % recovery tau is slower
    'tau1_out', 0.000500, ...  % depletion tau is fast 500 us
    'tau1_in', 0.000200, ...   % recovery tau is very fast 200 us
    'tau2_out', 0.001, ...     % depletion tau is pretty fast 1 ms
    'tau2_in', 0.010)          % recovery tau is slower 10 ms
end


% first figure out how many filter stages (PZFC/CARFAC channels):
pole_Hz = CF_CAR_params.first_pole_theta * fs / (2*pi);
n_ch = 0;
while pole_Hz > CF_CAR_params.min_pole_Hz
  n_ch = n_ch + 1;
  % pole_Hz = pole_Hz - CF_CAR_params.ERB_per_step * ...
  %  lyon2011_erbhz(pole_Hz, CF_CAR_params.ERB_break_freq, CF_CAR_params.ERB_Q);
  pole_Hz = pole_Hz - CF_CAR_params.ERB_per_step * ...
    f2erb(pole_Hz, CF_CAR_params.ERB_break_freq, CF_CAR_params.ERB_Q);
end
% Now we have n_ch, the number of channels, so can make the array
% and compute all the frequencies again to put into it:
pole_freqs = zeros(n_ch, 1);
pole_Hz = CF_CAR_params.first_pole_theta * fs / (2*pi);
for ch = 1:n_ch
  pole_freqs(ch) = pole_Hz;
%   pole_Hz = pole_Hz - CF_CAR_params.ERB_per_step * ...
%     lyon2011_erbhz(pole_Hz, CF_CAR_params.ERB_break_freq, CF_CAR_params.ERB_Q);
  pole_Hz = pole_Hz - CF_CAR_params.ERB_per_step * ...
    f2erb(pole_Hz, CF_CAR_params.ERB_break_freq, CF_CAR_params.ERB_Q);
end
% Now we have n_ch, the number of channels, and pole_freqs array.

max_channels_per_octave = log(2) / log(pole_freqs(1)/pole_freqs(2));

% Convert to include an ear_array, each w coeffs and state...
CAR_coeffs = local_designfilters(CF_CAR_params, fs, pole_freqs);
AGC_coeffs = local_designagc(CF_AGC_params, fs, n_ch);
IHC_coeffs = local_designihc(CF_IHC_params, fs, n_ch);

% Copy same designed coeffs into each ear (can do differently in the
% future, e.g. for unmatched OHC_health).
for ear = 1:n_ears
  ears(ear).CAR_coeffs = CAR_coeffs;
  ears(ear).AGC_coeffs = AGC_coeffs;
  ears(ear).IHC_coeffs = IHC_coeffs;
end

CF = struct( ...
  'fs', fs, ...
  'max_channels_per_octave', max_channels_per_octave, ...
  'CAR_params', CF_CAR_params, ...
  'AGC_params', CF_AGC_params, ...
  'IHC_params', CF_IHC_params, ...
  'n_ch', n_ch, ...
  'pole_freqs', pole_freqs, ...
  'ears', ears, ...
  'n_ears', n_ears, ...
  'open_loop', 0, ...
  'linear_car', 0);


%% Design the filter coeffs:
function CAR_coeffs = local_designfilters(CAR_params, fs, pole_freqs)

n_ch = length(pole_freqs);

% the filter design coeffs:
% scalars first:
CAR_coeffs = struct( ...
  'n_ch', n_ch, ...
  'velocity_scale', CAR_params.velocity_scale, ...
  'v_offset', CAR_params.v_offset ...
  );

% don't really need these zero arrays, but it's a clue to what fields
% and types are needed in other language implementations:
CAR_coeffs.r1_coeffs = zeros(n_ch, 1);
CAR_coeffs.a0_coeffs = zeros(n_ch, 1);
CAR_coeffs.c0_coeffs = zeros(n_ch, 1);
CAR_coeffs.h_coeffs = zeros(n_ch, 1);
CAR_coeffs.g0_coeffs = zeros(n_ch, 1);

CAR_coeffs.OHC_health = ones(n_ch, 1);  % 0 to 1 to derate OHC activity.

% zero_ratio comes in via h.  In book's circuit D, zero_ratio is 1/sqrt(a),
% and that a is here 1 / (1+f) where h = f*c.
% solve for f:  1/zero_ratio^2 = 1 / (1+f)
% zero_ratio^2 = 1+f => f = zero_ratio^2 - 1
f = CAR_params.zero_ratio^2 - 1;  % nominally 1 for half-octave

% Make pole positions, s and c coeffs, h and g coeffs, etc.,
% which mostly depend on the pole angle theta:
theta = pole_freqs .* (2 * pi / fs);

c0 = sin(theta);
a0 = cos(theta);

% different possible interpretations for min-damping r:
% r = exp(-theta * CF_CAR_params.min_zeta).
% Compress theta to give somewhat higher Q at highest thetas:
ff = CAR_params.high_f_damping_compression;  % 0 to 1; typ. 0.5
x = theta/pi;
theta = pi * (x - ff * x.^3);  % when ff is 0, this is just theta,
%                       and when ff is 1 it goes to zero at theta = pi.
max_zeta = CAR_params.max_zeta;
CAR_coeffs.r1_coeffs = (1 - theta .* max_zeta);  % "r1" for the max-damping condition

min_zeta = CAR_params.min_zeta;
if min_zeta <= 0  % Use this to do a new design strategy
  local_low_level_q = pole_freqs ./ lyon2011_erbhz( ...
    pole_freqs, CAR_params.ERB_break_freq, CAR_params.ERB_Q);
  % Number of overlapping channels is about ERB_per_step^-1, so this:
  min_zetas = CAR_params.ERB_per_step^-0.5 ./ (2*local_low_level_q);
  min_zetas = min(min_zetas, 0.75*max_zeta);  % Keep some low CF action.
  % "r1" for the max-damping condition
  CAR_coeffs.r1_coeffs = exp(-theta .* max_zeta);
  r0_coeffs = exp(-theta .* min_zetas);  % min_damping condition.
  CAR_coeffs.zr_coeffs = r0_coeffs - CAR_coeffs.r1_coeffs;
else
  % Increase the min damping where channels are spaced out more, by pulling
  % toward lyon2011_erbhz/pole_freqs (close to 0.1 at high f)
  min_zetas = min_zeta + 0.25*(f2erb(pole_freqs, ...
    CAR_params.ERB_break_freq, CAR_params.ERB_Q) ./ pole_freqs - min_zeta);
  CAR_coeffs.r1_coeffs = (1 - theta .* max_zeta);  % "r1" for the max-damping condition
  CAR_coeffs.zr_coeffs = theta .* ...
    (max_zeta - min_zetas);  % how r relates to undamping
end

% undamped coupled-form coefficients:
CAR_coeffs.a0_coeffs = a0;
CAR_coeffs.c0_coeffs = c0;

% the zeros follow via the h_coeffs
h = c0 .* f;
CAR_coeffs.h_coeffs = h;

relative_undamping = CAR_coeffs.OHC_health;  % Typically just ones.
% this function needs to take CAR_coeffs even if we haven't finished
% constucting it by putting in the g0_coeffs:
CAR_coeffs.g0_coeffs = lyon2011_stageg(CAR_coeffs, relative_undamping);


%% the AGC design coeffs:
function AGC_coeffs = local_designagc(AGC_params, fs, n_ch)

n_AGC_stages = AGC_params.n_stages;

% AGC1 pass is smoothing from base toward apex;
% AGC2 pass is back, which is done first now (in double exp. version)
AGC1_scales = AGC_params.AGC1_scales;
AGC2_scales = AGC_params.AGC2_scales;

decim = 1;

total_DC_gain = 0;

%%
% Convert to vector of AGC coeffs
AGC_coeffs = struct([]);
for stage = 1:n_AGC_stages
  AGC_coeffs(stage).n_ch = n_ch;
  AGC_coeffs(stage).n_AGC_stages = n_AGC_stages;
  AGC_coeffs(stage).AGC_stage_gain = AGC_params.AGC_stage_gain;

  AGC_coeffs(stage).decimation = AGC_params.decimation(stage);
  tau = AGC_params.time_constants(stage);  % time constant in seconds
  decim = decim * AGC_params.decimation(stage);  % net decim to this stage
  % epsilon is how much new input to take at each update step:
  AGC_coeffs(stage).AGC_epsilon = 1 - exp(-decim / (tau * fs));

  % effective number of smoothings in a time constant:
  ntimes = tau * (fs / decim);  % typically 5 to 50

  % decide on target spread (variance) and delay (mean) of impulse
  % response as a distribution to be convolved ntimes:
  % TODO (dicklyon): specify spread and delay instead of scales???
  delay = (AGC2_scales(stage) - AGC1_scales(stage)) / ntimes;
  spread_sq = (AGC1_scales(stage)^2 + AGC2_scales(stage)^2) / ntimes;

  % get pole positions to better match intended spread and delay of
  % [[geometric distribution]] in each direction (see wikipedia)
  u = 1 + 1 / spread_sq;  % these are based on off-line algebra hacking.
  p = u - sqrt(u^2 - 1);  % pole that would give spread if used twice.
  dp = delay * (1 - 2*p +p^2)/2;
  polez1 = p - dp;
  polez2 = p + dp;
  AGC_coeffs(stage).AGC_polez1 = polez1;
  AGC_coeffs(stage).AGC_polez2 = polez2;

  % try a 3- or 5-tap FIR as an alternative to the double exponential:
  n_taps = 0;
  FIR_OK = 0;
  n_iterations = 1;
  while ~FIR_OK
    switch n_taps
      case 0
        % first attempt a 3-point FIR to apply once:
        n_taps = 3;
      case 3
        % second time through, go wider but stick to 1 iteration
        n_taps = 5;
      case 5
        % apply FIR multiple times instead of going wider:
        n_iterations = n_iterations + 1;
        if n_iterations > 16
          error('Too many n_iterations in lyon2011_designagc');
        end
      otherwise
        % to do other n_taps would need changes in lyon2011_spatialsmooth
        % and in Design_FIR_coeffs
        error('Bad n_taps in lyon2011_designagc');
    end
    [AGC_spatial_FIR, FIR_OK] = local_designFIRcoeffs( ...
      n_taps, spread_sq, delay, n_iterations);
  end
  % when FIR_OK, store the resulting FIR design in coeffs:
  AGC_coeffs(stage).AGC_spatial_iterations = n_iterations;
  AGC_coeffs(stage).AGC_spatial_FIR = AGC_spatial_FIR;
  AGC_coeffs(stage).AGC_spatial_n_taps = n_taps;

  % accumulate DC gains from all the stages, accounting for stage_gain:
  total_DC_gain = total_DC_gain + AGC_params.AGC_stage_gain^(stage-1);

  % TODO (dicklyon) -- is this the best binaural mixing plan?
  if stage == 1
    AGC_coeffs(stage).AGC_mix_coeffs = 0;
  else
    AGC_coeffs(stage).AGC_mix_coeffs = AGC_params.AGC_mix_coeff / ...
      (tau * (fs / decim));
  end
end

% adjust stage 1 detect_scale to be the reciprocal DC gain of the AGC filters:
AGC_coeffs(1).detect_scale = 1 / total_DC_gain;


%%
function [FIR, OK] = local_designFIRcoeffs(n_taps, delay_variance, ...
  mean_delay, n_iter)
% function [FIR, OK] = Design_FIR_coeffs(n_taps, delay_variance, ...
%   mean_delay, n_iter)
% The smoothing function is a space-domain smoothing, but it considered
% here by analogy to time-domain smoothing, which is why its potential
% off-centeredness is called a delay.  Since it's a smoothing filter, it is
% also analogous to a discrete probability distribution (a p.m.f.), with
% mean corresponding to delay and variance corresponding to squared spatial
% spread (in samples, or channels, and the square thereof, respecitively).
% Here we design a filter implementation's coefficient via the method of
% moment matching, so we get the intended delay and spread, and don't worry
% too much about the shape of the distribution, which will be some kind of
% blob not too far from Gaussian if we run several FIR iterations.

% reduce mean and variance of smoothing distribution by n_iterations:
mean_delay = mean_delay / n_iter;
delay_variance = delay_variance / n_iter;
switch n_taps
  case 3
    % based on solving to match mean and variance of [a, 1-a-b, b]:
    a = (delay_variance + mean_delay*mean_delay - mean_delay) / 2;
    b = (delay_variance + mean_delay*mean_delay + mean_delay) / 2;
    FIR = [a, 1 - a - b, b];
    OK = FIR(2) >= 0.25;
  case 5
    % based on solving to match [a/2, a/2, 1-a-b, b/2, b/2]:
    a = ((delay_variance + mean_delay*mean_delay)*2/5 - mean_delay*2/3) / 2;
    b = ((delay_variance + mean_delay*mean_delay)*2/5 + mean_delay*2/3) / 2;
    % first and last coeffs are implicitly duplicated to make 5-point FIR:
    FIR = [a/2, 1 - a - b, b/2];
    OK = FIR(2) >= 0.15;
  otherwise
    error('Bad n_taps in AGC_spatial_FIR');
end


%% the IHC design coeffs:
function IHC_coeffs = local_designihc(IHC_params, fs, n_ch)

if IHC_params.just_hwr
  IHC_coeffs = struct( ...
    'n_ch', n_ch, ...
    'just_hwr', 1);
else
  if IHC_params.one_cap
    gmax = lyon2011_detect(10);  % output conductance at a high level
    rmin = 1 / gmax;
    c = IHC_params.tau_out * gmax;
    ri = IHC_params.tau_in / c;
    % to get approx steady-state average, double rmin for 50% duty cycle
    saturation_current = 1 / (2/gmax + ri);
    % also consider the zero-signal equilibrium:
    g0 = lyon2011_detect(0);
    r0 = 1 / g0;
    rest_current = 1 / (ri + r0);
    cap_voltage = 1 - rest_current * ri;
    IHC_coeffs = struct( ...
      'n_ch', n_ch, ...
      'just_hwr', 0, ...
      'ac_coeff', 2 * pi * IHC_params.ac_corner_Hz / fs, ...
      'lpf_coeff', 1 - exp(-1/(IHC_params.tau_lpf * fs)), ...
      'out_rate', rmin / (IHC_params.tau_out * fs), ...
      'in_rate', 1 / (IHC_params.tau_in * fs), ...
      'one_cap', IHC_params.one_cap, ...
      'output_gain', 1 / (saturation_current - rest_current), ...
      'rest_output', rest_current / (saturation_current - rest_current), ...
      'rest_cap', cap_voltage);
    % one-channel state for testing/verification:
    IHC_state = struct( ...
      'cap_voltage', IHC_coeffs.rest_cap, ...
      'lpf1_state', 0, ...
      'lpf2_state', 0, ...
      'ihc_accum', 0);
  else
    g1max = lyon2011_detect(10);  % receptor conductance at high level
    r1min = 1 / g1max;
    c1 = IHC_params.tau1_out * g1max;  % capacitor for min depletion tau
    r1 = IHC_params.tau1_in / c1;  % resistance for recharge tau
    % to get approx steady-state average, double r1min for 50% duty cycle
    saturation_current1 = 1 / (2*r1min + r1);  % Approximately.
    % also consider the zero-signal equilibrium:
    g10 = lyon2011_detect(0);
    r10 = 1/g10;
    rest_current1 = 1 / (r1 + r10);
    cap1_voltage = 1 - rest_current1 * r1;  % quiescent/initial state

    % Second cap similar, but using receptor voltage as detected signal.
    max_vrecep = r1 / (r1min + r1);  % Voltage divider from 1.
    % Identity from receptor potential to neurotransmitter conductance:
    g2max = max_vrecep;  % receptor resistance at very high level
    r2min = 1 / g2max;
    c2 = IHC_params.tau2_out * g2max;  % capacitor for min depletion tau
    r2 = IHC_params.tau2_in / c2;  % resistance for recharge tau
    % to get approx steady-state average, double r2min for 50% duty cycle
    saturation_current2 = 1 / (2 * r2min + r2);
    % also consider the zero-signal equilibrium:
    rest_vrecep = r1 * rest_current1;
    g20 = rest_vrecep;
    r20 = 1 / g20;
    rest_current2 = 1 / (r2 + r20);
    cap2_voltage = 1 - rest_current2 * r2;  % quiescent/initial state

    IHC_coeffs = struct(...
      'n_ch', n_ch, ...
      'just_hwr', 0, ...
      'ac_coeff', 2 * pi * IHC_params.ac_corner_Hz / fs, ...
      'lpf_coeff', 1 - exp(-1/(IHC_params.tau_lpf * fs)), ...
      'out1_rate', r1min / (IHC_params.tau1_out * fs), ...
      'in1_rate', 1 / (IHC_params.tau1_in * fs), ...
      'out2_rate', r2min / (IHC_params.tau2_out * fs), ...
      'in2_rate', 1 / (IHC_params.tau2_in * fs), ...
      'one_cap', IHC_params.one_cap, ...
      'output_gain', 1 / (saturation_current2 - rest_current2), ...
      'rest_output', rest_current2 / (saturation_current2 - rest_current2), ...
      'rest_cap2', cap2_voltage, ...
      'rest_cap1', cap1_voltage);
    % one-channel state for testing/verification:
    IHC_state = struct( ...
      'cap1_voltage', IHC_coeffs.rest_cap1, ...
      'cap2_voltage', IHC_coeffs.rest_cap2, ...
      'lpf1_state', 0, ...
      'lpf2_state', 0, ...
      'ihc_accum', 0);
  end
end


% function g = lyon2011_stageg(CAR_coeffs, relative_undamping)
% % function g = lyon2011_stageg(CAR_coeffs, relative_undamping)
% % Return the stage gain g needed to get unity gain at DC
% 
% r1 = CAR_coeffs.r1_coeffs;  % at max damping
% a0 = CAR_coeffs.a0_coeffs;
% c0 = CAR_coeffs.c0_coeffs;
% h  = CAR_coeffs.h_coeffs;
% zr = CAR_coeffs.zr_coeffs;
% r  = r1 + zr .* relative_undamping;
% g  = (1 - 2*r.*a0 + r.^2) ./ (1 - 2*r.*a0 + h.*r.*c0 + r.^2);

