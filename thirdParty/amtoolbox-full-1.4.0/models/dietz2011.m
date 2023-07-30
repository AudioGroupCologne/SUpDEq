function [fine,fc,ild,env] = dietz2011(insig,fs,varargin)
%DIETZ2011  Sound lateral direction
%   Usage: [fine,fc,ild,env] = dietz2011(insig,fs);
%
%   Input parameters:
%       insig       : binaural signal for which values should be calculated
%       fs          : sampling rate (Hz)
%
%   Output parameters:
%       fine        : Information about the fine structure (see below)
%       fc          : center frequencies of gammatone filterbank
%       ild         : interaural level difference in dB
%       env         : Information about the envelope (see below)
%
%   DIETZ2011(insig,fs) calculates interaural phase, time and level
%   differences of fine- structure and envelope of the signal, as well as
%   the interaural coherence, which can be used as a weighting function.
%   The _lp values are not returned if the 'nolowpass' flag is set.
%
%   The steps of the binaural model to calculate the result are the
%   following (see also Dietz et al., 2011):
%
%   1) Middle ear filtering (500-2000 Hz 1st order bandpass)
%
%   2) Auditory bandpass filtering on the basilar membrane using a
%      4th-order all-pole gammatone filterbank, employing 23 filter
%      bands between 200 and 5000 Hz, with a 1 ERB spacing. The filter
%      width was set to correspond to 1 ERB.
%
%   3) Cochlear compression was simulated by power-law compression with
%      an exponent of 0.4.
%
%   4) The transduction process in the inner hair cells was modelled
%      using half-wave rectification followed by filtering with a 770-Hz
%      5th order lowpass.
%
%   5) Modulationfilterbank with three different filters applied to every
%      frequency channel. One 2nd order gammatone filter for the fine structure
%      centered at the center frequency of the frequency channel. One 2nd order
%      gammatone filter for the envelope of the signal centered at 135 Hz.
%      And a 2nd order lowpass filter with a cutoff frequency of 30 Hz to
%      extract the ILD of the signal.
%
%   6) Calculation of binaural parameters such as IPD, ITD, IC for fine
%      structure and envelope filter signals and ILD for the ILD filter.
%
%
%   The output structures fine and env have the following fields:
%
%     .itf        transfer function
%
%     .ipd        phase difference in rad
%
%     .itd        interaural time difference based on instantaneous frequency
%
%     .itd_C      interaural time difference based on center frequency
%
%     .f_inst_1   instantaneous frequencies of left ear signal
%
%     .f_inst_2   instantaneous frequencies of right ear canal signal
%
%     .f_inst     instantaneous frequencies (average of f_inst1 and 2)
%
%     .ic         interaural coherence
%
%     .rms        rms value of frequency channels for weighting
%
%     .ild_lp     based on low passed-filtered insig, level difference in dB
%
%     .ipd_lp     based on lowpass-filtered itf, phase difference in rad
%
%     .itd_lp     based on lowpass-filtered itf, interaural time difference
%
%     .itd_C_lp   based on lowpass-filtered itf, interaural time difference
%
%     .f_inst_lp  lowpass instantaneous frequencies
%
%
%
%   DIETZ2011 accepts the following optional parameters:
%
%     'flow',flow    Set the lowest frequency in the filterbank to
%                    flow. Default value is 200 Hz.
%
%     'fhigh',fhigh  Set the highest frequency in the filterbank to
%                    fhigh. Default value is 5000 Hz.
%
%     'basef',basef  Ensure that the frequency basef is a center frequency
%                    in the filterbank. The default value  is 1000.
%
%     'filters_per_ERB',filters_per_erb
%                    Filters per erb. The default value is 1.
%
%     'middle_ear_thr',r
%                    Bandpass frequencies for middle ear transfer. The
%                    default value is [500 2000].
%
%     'middle_ear_order',n 
%                    Order of middle ear filter. Only even numbers are
%                    possible. The default value is 2.
%
%     'compression_power',cpwr
%                    Applied compression of the signal on the cochlea with
%                    ^compression_power. The default value is 0.4.
%
%     'alpha',alpha  Internal noise strength. Convention 65dB = 0.0354.
%                    The default value is 0.
%
%     'int_randn'    Internal noise by adding random noise with rms = alpha.
%                    This is the default.
%
%     'int_mini'     Internal noise by setting all values < alpha to alpha.
%
%     'filter_order',fo
%                    Filter order for the two gammatone filter used for the fine
%                    structure and envelope of the modulation filter bank. The
%                    default value is 2.
%
%     'filter_attenuation_db',fadb
%                    Filter attenuation for the two gammatone filter used for the fine
%                    structure and envelope of the modulation filter bank. The
%                    default value is 10.
%                     
%     'fine_filter_finesse',fff
%                    Filter finesse (determines the bandwidth with fc/finesse)
%                    for the fine structure gammatone filter. The default value
%                    is 3.
%
%     'mod_center_frequency_hz',mcf_hz
%                    Center frequency of the gammatone envelope filter. The
%                    default value is 135.
%
%     'mod_filter_finesse',mff
%                    Filter finesse (determines the bandwidth with fc/finesse)
%                    for the envelope gammatone filter. The default value is 8.
% 
%     'level_filter_cutoff_hz',lfc_hz
%                    Cutoff frequency off the low pass filter used for ILD
%                    calculation. The default value is 30.
%
%     'level_filter_order',lforder
%                    Order of low pass filter for the ILD calculation. The
%                    default value is 2.
%
%     'tau_cycles',tau_cycles
%                    Temporal resolution of binaural processor in terms of
%                    cycles per frequency channel. The default value is 5.
%
%     'signal_level_dB_SPL',signal_level
%                    Sound pressure level of left channel. Used for data
%                    display and analysis. Default value is 70.
%
%     'lowpass'      Calculate the interaural parameters of the lowpassed
%                    signal/ITF (_lp return values). This is the default.
%
%     'nolowpass'    Don't calculate the lowpass based interaural parameters.
%                    The _lp values are not returned.
%
%     'debug'        Display what is happening.
%
%
%   See also: dietz2011_interauralfunctions, dietz2011_filterbank,
%     ihcenvelope, auditoryfilterbank itd2angle itd2angle_lookuptable
%     sig_competingtalkers dietz2011_filterbank dietz2011_interauralfunctions
%     wierstorf2013_estimateazimuth dietz2011_unwrapitd
%     exp_steidle2019 exp_dietz2011 breebaart2001 wierstorf2013 hohmann2002
%
%   References:
%     M. Dietz, S. D. Ewert, and V. Hohmann. Auditory model based direction
%     estimation of concurrent speakers from binaural signals. Speech
%     Communication, 53(5):592--605, 2011. [1]http ]
%     
%     References
%     
%     1. http://www.sciencedirect.com/science/article/pii/S016763931000097X
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/dietz2011.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: M-Signal
%   #Author: Tobias Peters (2002)
%   #Author: Mathias Dietz (2006-2009)
%   #Author: Martin Klein-Hennig (2011)
%   #Author: Martin Klein-Hennig 
%   #Author: Hagen Wierstorf (2013): for AMT 

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% This model got a major update, and the figures of exp_dietz2011() look now
% slightly different than before. If you want to restore the original version of
% the model please checkout version 29a048b with git or install AMToolbox 0.9.5

  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(insig) || min(size(insig))~=2
    error('%s: insig has to be a numeric two channel signal!',upper(mfilename));
end

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs has to be a positive scalar!',upper(mfilename));
end

% import default arguments from other functions
definput.import={'auditoryfilterbank','ihcenvelope','dietz2011_filterbank','dietz2011_interauralfunctions'};
% changing default parameters
definput.importdefaults = { ...
    'flow',200, ...     % gammatone lowest frequency / Hz
    'fhigh',5000, ...   % gammatone highest frequency / Hz
    'basef',1000, ...   % auditory filter should be centered at basef / Hz
    'ihc_breebaart2001' ... % use haircell parameters as in Breebarts model
};
% Preprocessing parameters
definput.keyvals.middle_ear_thr = [500 2000]; % Bandpass freqencies for middle ear transfer
definput.keyvals.middle_ear_order = 2;        % Only even numbers possible
definput.keyvals.alpha = 0;                   % Internal noise strength
                                              % 65dB = 0.0354
% randn: add random noise with rms = alpha
% mini: set all values < alpha to alpha
definput.flags.int_noise_case = {'int_randn','int_mini'};

% debugging messages
definput.flags.disp = {'no_debug','debug'};

[flags,kv]  = ltfatarghelper({},definput,varargin);



%% Model processing starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amt_disp('DIETZ2011:',flags.disp);

%% ---- middle ear band pass filtering ------
% Compare Puria et al. (1997)
amt_disp('  band pass filtering of input according to middle ear transfer charact.',flags.disp);
[b,a] = butter(kv.middle_ear_order,kv.middle_ear_thr(2)/(fs/2),'low');
inoutsig = filter(b,a,insig);
[b,a] = butter(kv.middle_ear_order,kv.middle_ear_thr(1)/(fs/2),'high');
inoutsig = filter(b,a,inoutsig);

%% ---- inner ear ------
% gammatone filterbank
amt_disp('  splitting signal into frequency channels',flags.disp);
[inoutsig,fc] = auditoryfilterbank(inoutsig,fs,'argimport',flags,kv);
% cochlea compression
inoutsig = sign(inoutsig).*abs(inoutsig).^kv.compression_power;
% rectification and lowpass filtering of filtered signals
amt_disp('  haircell processing of frequency bands',flags.disp);
inoutsig = ihcenvelope(inoutsig,fs,'argimport',flags,kv);

%% ---- internal noise ------
% additive white noise
if flags.do_int_randn
  amt_disp('  adding internal random noise',flags.disp);
  inoutsig = inoutsig + kv.alpha*randn(size(inoutsig));
end
% replace values<kv.alpha with kv.alpha
if flags.do_int_mini
  amt_disp('  adding internal noise via minimum',flags.disp);
  inoutsig = max(inoutsig,kv.alpha);
end

%% ---- modulation filterbank ------
% filter signals with three different filters to get fine structure, envelope
% and ILD low pass
amt_disp('  apply second filterbank',flags.disp)
[inoutsig_fine,fc_fine,inoutsig_env,fc_env,inoutsig_ild] = ...
  dietz2011_filterbank(inoutsig,fs,fc,'argimport',flags,kv);

%% ---- binaural processor ------
% calculate interaural parameters for fine structure and envelope and calculate
% ILD
% -- fine structure
amt_disp('  calculating interaural functions from haircell fine structure',flags.disp);
fine = dietz2011_interauralfunctions(inoutsig_fine,fs,fc_fine,'argimport',flags,kv);
% --envelope
amt_disp('  calculating interaural functions from haircell modulation',flags.disp);
env = dietz2011_interauralfunctions(inoutsig_env,fs, ...
  kv.mod_center_frequency_hz+0*fc_env,'argimport',flags,kv);
% -- ILD
% interaural level difference, eq. 5 in Dietz (2011)
% max(sig,1e-4) avoids division by zero
amt_disp('  determining ILD',flags.disp);
ild = 20/kv.compression_power*log10(max(inoutsig_ild(:,:,2),1e-4)./max(inoutsig_ild(:,:,1),1e-4));

end


