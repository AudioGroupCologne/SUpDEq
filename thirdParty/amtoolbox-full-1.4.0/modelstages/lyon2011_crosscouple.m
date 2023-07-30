function ears = lyon2011_crosscouple(ears);
%LYON2011_CROSSCOUPLE adjust the intensity of the ear signals
%   Usage: [CF, decim_naps, naps, BM, ohc, agc] = lyon2011(CF, input_waves, AGC_plot_fig_num, open_loop);
%
%
%   Input parameters:
%     detects             : The CF struct holds the filterbank design and 
%                           state; if you want to break the input up into
%                           segments, you need to use the updated CF
%                           to keep the state between segments.
%     coeffs              : input_waves is a column vector if there's just one
%                           audio channel; more generally, it has a row per 
%                           time sample, a column per audio channel. The 
%                           input_waves are assumed to be sampled at the 
%                           same rate as the CARFAC is designed for. 
%                           A resampling may be needed before calling this.
%     state               : Plot automatic gain control figure. Default is 0.
%
%   Output parameters:
%     state               : The CF struct holds the filterbank design and 
%                           state; if you want to break the input up into
%                           segments, you need to use the updated CF
%                           to keep the state between segments.
%     update              : decim_naps is like naps but time-decimated by 
%                           the int CF.decimation.
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lyon2011_crosscouple.php


%   #Author: Amin Saremi (2016) adaptations for the AMT (based on <https://github.com/google/carfac>, Richard F. Lyon)
%   #Author: Clara Hollomey (2021) adaptation for the AMT 1.0
%   #License: gpl3



n_ears = length(ears);
if n_ears > 1
  n_stages = ears(1).AGC_coeffs(1).n_AGC_stages;
  % now cross-ear mix the stages that updated (leading stages at phase 0):
  for stage = 1:n_stages
    if ears(1).AGC_state(stage).decim_phase > 0
      break  % all recently updated stages are finished
    else
      mix_coeff = ears(1).AGC_coeffs(stage).AGC_mix_coeffs;
      if mix_coeff > 0  % Typically stage 1 has 0 so no work on that one.
        this_stage_sum = 0;
        % sum up over the ears and get their mean:
        for ear = 1:n_ears
          stage_state = ears(ear).AGC_state(stage).AGC_memory;
          this_stage_sum = this_stage_sum + stage_state;
        end
        this_stage_mean = this_stage_sum / n_ears;
        % now move them all toward the mean:
        for ear = 1:n_ears
          stage_state = ears(ear).AGC_state(stage).AGC_memory;
          ears(ear).AGC_state(stage).AGC_memory = ...
            stage_state +  mix_coeff * (this_stage_mean - stage_state);
        end
      end
    end
  end
end

