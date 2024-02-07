function [CF, decim_naps, naps, BM, ohc, agc, ihc_potential] = lyon2011(input_waves, CF)
%LYON2011 Cascade of asymmetric resonators with fast-acting compression (CARFAC) model
%   Usage: [CF, decim_naps, naps, BM, ohc, agc] = lyon2011(input_waves,CF);
%
%
%   Input parameters:
%     input_waves      : input_waves is a column vector if there's just one
%                        audio channel; more generally, it has a row per
%                        time sample, a column per audio channel. The
%                        input_waves are assumed to be sampled at the
%                        same rate as the CARFAC model is designed for.
%                        A resampling may be needed before calling this.
%     CF               : The CF struct holds the CARFAC design and state.
%
%   Output parameters:
%     CF               : The CF struct holds the CARFAC design and
%                        state; if you want to break the input up into
%                        segments, you need to use the updated CF
%                        to keep the state between segments.
%     decim_naps       : decim_naps is like naps but time-decimated by
%                        the int CF.seglen (defaults to about 20 ms long).
%     naps             : naps (neural activity patterns) has a row per
%                        time sample, a column per filterbank channel,
%                        and a layer per audio channel if n_ears > 1.
%     BM               : BM is basilar membrane motion (filter outputs
%                        before detection).
%     ohc              : optional extra output for diagnosing internals.
%     agc              : optional extra outputs for diagnosing internals.
%     ihc_potential    : optional extra output for the IHC potential equivalent.
%
%   LYON2011 runs the CARFAC model. That is, it filters a 1 or more channel
%   sound input to make one or more neural activity patterns (naps). This 
%   file is part of an implementation of Lyon's cochlear model:
%   "Cascade of Asymmetric Resonators with Fast-Acting Compression"
%
%   License
%   --------
%
%   This model is licensed under the Apache License Version 2.0. Further usage details 
%   are provided in the in the AMT directory "licences".
%
%   References:
%     R. F. Lyon. Cascades of two-pole–two-zero asymmetric resonators are
%     good models of peripheral auditory function. J. Acoust. Soc. Am.,
%     130(6), 2011.
%     
%
%
%   See also:   data_lyon2011 demo_lyon2011_impulseresponses demo_lyon2011
%               demo_lyon2011_compressivefunctions lyon2011_init
%               lyon2011_agcstep lyon2011_carstep
%               lyon2011_ihcstep lyon2011_crosscouple lyon2011_detect
%               lyon2011_stageg lyon2011_closeagcloop
%               lyon2011_spatialsmooth erbest f2erb
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/lyon2011.php


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


%  Probably a better scale for CARFAC; from 104 dB SPL at FS to 94.
input_waves = 0.316 * input_waves;

[n_samp, n_ears] = size(input_waves);
if n_ears ~= CF.n_ears
  error('bad number of input_waves channels passed to CARFAC_Run')
end
n_ch = CF.n_ch;

if ~isfield(CF, 'seglen')  % If no seglen specified...
  CF.seglen = round(CF.fs/50);  % anything should work; this is 20 ms.
end
seglen = CF.seglen;
n_segs = ceil(n_samp / seglen);

if nargout > 3
  BM = zeros(n_samp, n_ch, n_ears);
else
  BM = [];
end

if nargout > 4
  ohc = zeros(n_samp, n_ch, n_ears);
else
  ohc = [];
end

if nargout > 5
  agc = zeros(n_samp, n_ch, n_ears);
else
  agc = [];
end

if nargout > 6
  for ear = 1:n_ears
    if ~isfield(CF.ears(ear).IHC_state, 'cap1_voltage')
      error('To get the IHC potential, please use the two-cap version.')
    end
  end

  ihc_potential = zeros(n_samp, n_ch, n_ears);
else
  ihc_potential = [];
end

if nargout > 1
  % make decimated detect output:
  decim_naps = zeros(n_segs, CF.n_ch, CF.n_ears);
else
  decim_naps = [];
end

if nargout > 2
  % make decimated detect output:
  naps = zeros(n_samp, CF.n_ch, CF.n_ears);
else
  naps = [];
end

for seg_num = 1:n_segs
  if seg_num == n_segs
    % The last segement may be short of seglen, but do it anyway:
    k_range = (seglen*(seg_num - 1) + 1):n_samp;
  else
    k_range = seglen*(seg_num - 1) + (1:seglen);
  end
  % Process a segment to get a slice of decim_naps, and plot AGC state:
  if ~isempty(ihc_potential)
    [seg_naps, CF, seg_BM, seg_ohc, seg_agc, seg_ihc_potential] = CARFAC_Run_Segment(...
      CF, input_waves(k_range, :));
  elseif ~isempty(BM)
    [seg_naps, CF, seg_BM, seg_ohc, seg_agc] = CARFAC_Run_Segment(...
      CF, input_waves(k_range, :));
  else
    [seg_naps, CF] = CARFAC_Run_Segment(CF, input_waves(k_range, :));
  end

  if ~isempty(BM)
    for ear = 1:n_ears
      % Accumulate segment BM to make full BM
      BM(k_range, :, ear) = seg_BM(:, :, ear);
    end
  end

  if ~isempty(naps)
    for ear = 1:n_ears
      % Accumulate segment naps to make full naps
      naps(k_range, :, ear) = seg_naps(:, :, ear);
    end
  end

  if ~isempty(ohc)
    for ear = 1:n_ears
      % Accumulate segment naps to make full naps
      ohc(k_range, :, ear) = seg_ohc(:, :, ear);
    end
  end

  if ~isempty(agc)
    for ear = 1:n_ears
      % Accumulate segment naps to make full naps
      agc(k_range, :, ear) = seg_agc(:, :, ear);
    end
  end

  if ~isempty(ihc_potential)
    for ear = 1:n_ears
      % Accumulate segment IHC capacitor voltage to make full IHC vector
      ihc_potential(k_range, :, ear) = seg_ihc_potential(:, :, ear);
    end
  end

  if ~isempty(decim_naps)
    for ear = 1:n_ears
      decim_naps(seg_num, :, ear) = CF.ears(ear).IHC_state.ihc_accum / ...
        seglen;
      CF.ears(ear).IHC_state.ihc_accum = zeros(n_ch,1);
    end
  end

end  % segment loop
return


function [naps, CF, BM, seg_ohc, seg_agc, seg_ihc_potential] = CARFAC_Run_Segment(...
  CF, input_waves)
% function [naps, CF, BM, seg_ohc, seg_agc, seg_ihc_potential] = CARFAC_Run_Segment(...
%   CF, input_waves)
%
% This function runs the CARFAC; that is, filters a 1 or more channel
% sound input segment to make one or more neural activity patterns (naps);
% it can be called multiple times for successive segments of any length,
% as long as the returned CF with modified state is passed back in each
% time.
%
% input_waves is a column vector if there's just one audio channel;
% more generally, it has a row per time sample, a column per audio channel.
%
% naps has a row per time sample, a column per filterbank channel, and
% a layer per audio channel if more than 1.
% BM is basilar membrane motion (filter outputs before detection).
%
% the input_waves are assumed to be sampled at the same rate as the
% CARFAC is designed for; a resampling may be needed before calling this.
%
% The function works as an outer iteration on time, updating all the
% filters and AGC states concurrently, so that the different channels can
% interact easily.  The inner loops are over filterbank channels, and
% this level should be kept efficient.
%
% seg_ohc seg_agc are optional extra outputs useful for seeing what the
% ohc nonlinearity and agc are doing; both in terms of extra damping.
%
% seg_ihc_potential is also optional, and is used for getting the internal IHC
% state, comparable to the IHC potential.

if nargout > 2
  do_BM = 1;
else
  do_BM = 0;
end

[n_samp, n_ears] = size(input_waves);

if n_ears ~= CF.n_ears
  error('bad number of input_waves channels passed to CARFAC_Run_Segment')
end

n_ch = CF.n_ch;
naps = zeros(n_samp, n_ch, n_ears);  % allocate space for result
if do_BM
  BM = zeros(n_samp, n_ch, n_ears);
  seg_ohc = zeros(n_samp, n_ch, n_ears);
  seg_agc = zeros(n_samp, n_ch, n_ears);
end

if nargout > 5
  seg_ihc_potential = zeros(n_samp, n_ch, n_ears);
else
  seg_ihc_potential = [];
end

%detects = zeros(n_ch, n_ears);
for k = 1:n_samp
  % at each time step, possibly handle multiple channels
  for ear = 1:n_ears

    % This would be cleaner if we could just get and use a reference to
    % CF.ears(ear), but Matlab doesn't work that way...

    [car_out, CF.ears(ear).CAR_state] = lyon2011_carstep( ...
      input_waves(k, ear), CF.ears(ear).CAR_coeffs, ...
      CF.ears(ear).CAR_state);

    % update IHC state & output on every time step, too
    [ihc_out, CF.ears(ear).IHC_state] = lyon2011_ihcstep( ...
      car_out, CF.ears(ear).IHC_coeffs, CF.ears(ear).IHC_state);

    % run the AGC update step, decimating internally,
    [CF.ears(ear).AGC_state, updated] = lyon2011_agcstep( ...
      ihc_out, CF.ears(ear).AGC_coeffs, CF.ears(ear).AGC_state);

    % save some output data:
    naps(k, :, ear) = ihc_out;  % output to neural activity pattern
    if do_BM
      BM(k, :, ear) = car_out;
      state = CF.ears(ear).CAR_state;
      seg_ohc(k, :, ear) = state.zA_memory;
      seg_agc(k, :, ear) = state.zB_memory;
    end

    if ~isempty(seg_ihc_potential)
      seg_ihc_potential(k, :, ear) = 1 - CF.ears(ear).IHC_state.cap1_voltage;
    end
  end

  % connect the feedback from AGC_state to CAR_state when it updates;
  % all ears together here due to potential mixing across them:
  if updated
    if n_ears > 1
      % do multi-aural cross-coupling:
      CF.ears = lyon2011_croscouple(CF.ears);
    end
    if CF.open_loop
      for ear = 1:CF.n_ears
        CF.ears(ear).CAR_state.dzB_memory = 0;  % To stop intepolating.
        CF.ears(ear).CAR_state.dg_memory = 0;  % To stop intepolating.
      end
    else
      CF = lyon2011_closeagcloop(CF);
    end
  end
end

