function stimulus = sig_vanderheijden1999(c_itd,c_sphase,c_noise_mode,c_tone_level,fs,dur,cos_rise_time,bw,fc,spl)
%sig_vanderheijden1999   Tone masked by a constant-ITD inner-band noise
%   Usage: stimulus = sig_vanderheijden1999(c_itd,c_sphase,c_noise_mode,c_tone_level,fs,dur,cos_rise_time,bw,fc,spl)
%
%   See also: eurich2022 exp_eurich2022
%
%   References:
%     B. Eurich, J. Encke, S. D. Ewert, and M. Dietz. Lower interaural
%     coherence in off-signal bands impairs binaural detection. The Journal
%     of the Acoustical Society of America, 151(6):3927--3936, 06 2022.
%     [1]arXiv | [2]http ]
%     
%     References
%     
%     1. http://arxiv.org/abs/https://pubs.aip.org/asa/jasa/article-pdf/151/6/3927/16528275/3927\_1\_online.pdf
%     2. https://doi.org/10.1121/10.0011673
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_vanderheijden1999.php


%   #Author: Bernhard Eurich (2022): original implementation

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

spar.itd = c_itd;
spar.noise_mode = c_noise_mode;
spar.fs = fs;
spar.dur = dur;
spar.cos_rise_time = cos_rise_time;
spar.bw = bw;
spar.spl = spl;
spar.fc = fc;
spar.sphase = c_sphase;
spar.tone_level = c_tone_level;
spar.nlcut = 50; % lower freq limit of masker
spar.nhcut = 950; % upper
% this function generates the stimulus for the experiment from
% van der Heijden & Trahiotis 1999 (as simulated in Fig. 3 in Eurich et al. 2022)


% tone
tone_mono = gen_tone(spar.fc,spar.dur, spar.fs,0);
tone      = [tone_mono spar.sphase*tone_mono];
tone_spl  = set_dbspl(tone,spar.tone_level);

[window] = cosine_fade_window([1:spar.dur*spar.fs]', spar.cos_rise_time, spar.fs);


% generate noises
[mDDnoise, mSDnoise]     = gen_DDNoise(spar.dur*fs, spar.nlcut, spar.nhcut, spar.itd, spar.fs);

% noise switch
if spar.noise_mode == 1
    noise = mSDnoise;
elseif spar.noise_mode == 2
    noise = mDDnoise;
end

[noise_spl]   = set_dbspl(noise,spar.spl + 10*log10(spar.bw));





% add tone to noise
stimulus = tone_spl + noise_spl;


end


function [DDnoise, SDnoise] = gen_DDNoise(len, flcut, fhcut, itd, fs)

% generate two independent noises

% bw =  (fhcut - flcut);

% ====  create shifted noises in f domain ====

% nyquist freq bin, depending on even or odd number of samples
fftpts = fs; % fft points for 1s length;
spacing = 1; % Hz per bin
nyBin = floor(fs/2) + 1;

% vector of freq per bin up to Nyquist freq
freqVec = [0:nyBin-1]' * spacing;

% lowest and highest bins
% lower flanking noise
lbin = max( round(flcut/spacing) + 1, 2);
hbin = min( round(fhcut/spacing) + 1, nyBin );

a1 = zeros(fftpts,1);
b1 = a1;
a2 = zeros(fftpts,1);
b2 = a2;

% two independent noises with each real and imag part
a1(lbin:hbin) = randn(hbin-lbin+1,1);
a2(lbin:hbin) = randn(hbin-lbin+1,1);
b1(lbin:hbin) = randn(hbin-lbin+1,1);
b2(lbin:hbin) = randn(hbin-lbin+1,1);


% complex noise spectra
fspec1 = a1+ i*b1;
fspec2 = a2+ i*b2;

fspec1_shift = fspec1;
fspec2_shift = fspec2;


% introduce opposing ITDs
itd_shift = 2*pi*freqVec(lbin:hbin)*itd;
fspec1_shift(lbin:hbin) = fspec1(lbin:hbin) .* exp(1i*itd_shift);
fspec2_shift(lbin:hbin) = fspec2(lbin:hbin) .* exp(1i*-itd_shift);

SDN1_l = (2*real(fftpts*ifft(fspec1))) / std(2*real( fftpts * ifft(fspec1) ));
SDN1_r = (2*real(fftpts*ifft(fspec1_shift)))/std(2*real(fftpts*ifft(fspec1_shift)));

SDN2_l = (2*real(fftpts*ifft(fspec2))) / std(2*real( fftpts * ifft(fspec2) ));
SDN2_r = (2*real(fftpts*ifft(fspec2_shift)))/std(2*real(fftpts*ifft(fspec2_shift)));

spec_DDN_l = fspec1 + fspec2;
spec_DDN_r = fspec1_shift + fspec2_shift;

DDN_l = (2*real(fftpts*ifft(spec_DDN_l)))/std(2*real(fftpts*ifft(spec_DDN_l)));
DDN_r = (2*real(fftpts*ifft(spec_DDN_r)))/std(2*real(fftpts*ifft(spec_DDN_r)));
    
SDnoise     = [SDN1_l(1:len) SDN1_r(1:len)];
DDnoise     = [DDN_l(1:len) DDN_r(1:len)];

end

function [signal_out] = set_dbspl(signal_in,dbspl_val)
%set_dbspl set SPL for each signal channel
% SIGNAL_OUT = set_dbspl(SIGNAL_IN, DBSPL_VAL)
% SIGNAL_IN     preassure waveform, multi channels in colums
% DBSPL_VAL     level in dB SPL (ref value is 20e-6 pa), per chanel 
%
% see also get_dbspl, add_dbgain, audio_signal_info


p0 = 20e-6; %ref value
    
val = sqrt(mean(signal_in.^2));   

factor = (p0 * 10.^(dbspl_val / 20)) ./ val;

signal_out = signal_in .* factor;

end

function [window] = cosine_fade_window(signal, rise_time, fs)
%cosine_fade_window returns a weighting vector for windowing (fade-in + fade-out) a signal
%WINDOW = cosine_fade_window(SIGNAL, RISE_TIME, FS)
% SIGNAL    preassure waveform, multi channels in colums
% RISE_TIME time in seconds
% FS        sampling frequency
% WINDOW    vector with the length of SIGNAL
%   
%           ................
%        .                    .
%     .                          .
%  .                                .
% |rise_time|              |rise_time|
% |          signal_time             | 
%
% EXAMPLE:
%       fs = 48e3;
%       sig = generate_tone(100,.5,fs);
%       window = cosine_fade_window(sig,.1,fs);
%       ramped_sig = sig .* window;

n_ramp = round(rise_time * fs);
n_signal = size(signal,1);

window = ones(1, n_signal);
flank = 0.5 * (1 + cos(pi / n_ramp * (-n_ramp:-1)));
window(1:n_ramp) = flank;
window(end-n_ramp+1:end) = fliplr(flank);
window = window';

end


function [sine,time] = gen_tone(frequency, duration, fs, start_phase)
% gen_tone returns a cosine wave
% [SINE, TIME ] = gen_tone(FREQUENCY, DURATION, FS, START_PHASE)
% FREQUENCY     frequency in Hz
% DURATION      duration in seconds
% FS            sampling frequenc
% START_PHASE   default is 0
% SINE          sine wave
% TIME          time vector for the cosine-wave
%
% see also getAMS gen_sam audio_signal_info
    if nargin < 4
        start_phase = 0;
    end

    nsamp = round(duration * fs);
    time = get_time(nsamp, fs);

    sine = cos(2 * pi * frequency * time + start_phase);

end

function [time] = get_time(signal, fs)
%get_time returns time-vector for a given 
% signal or number of samples
%[TIME] = get_time(SIGNAL, FS)
% SIGNAL    signal or number of samples
% FS        sampling frequency
% TIME      time vector (start with 0)
%
% see also nsamples

dt = 1. / fs;

if length(signal) > 1
    nsamp = length(signal);
else
    nsamp = signal;
end
    
max_time = nsamp * dt;

time = 0:dt:(max_time - dt);
    
time = time';
end

