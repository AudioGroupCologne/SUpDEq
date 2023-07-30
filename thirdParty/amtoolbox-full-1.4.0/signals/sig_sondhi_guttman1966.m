function [stimulus] = sig_sondhi_guttman1966(c_innerbw,c_flanking_phase,c_noise_mode,c_tone_level,c_itd,c_sphase,fs,dur,cos_rise_time,bw,fc,spl)

spar.innerbw = c_innerbw;
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
spar.innerbw = c_innerbw;
spar.flanking_phase = c_flanking_phase;
spar.nlcut = 0; % lower freq limit of masker
spar.nhcut = 2500; % upper
% this function generates the stimulus for the experiment from
% Sondhi & Guttman 1966  (as simulated in Fig. 4 in Eurich et al. 2022)
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_sondhi_guttman1966.php


%%%%%

% tone
tone_mono = gen_tone(spar.fc,spar.dur, spar.fs,0);
tone      = [tone_mono spar.sphase*tone_mono]; 
tone_spl  = set_dbspl(tone,spar.tone_level);

len = spar.dur * spar.fs;

[window] = cosine_fade_window([1:spar.dur*spar.fs]', spar.cos_rise_time, spar.fs);


spar.nlcut = floor(spar.fc-spar.bw/2); % noise
spar.nhcut = ceil(spar.fc+spar.bw/2);
bw_flank =  (spar.nhcut - spar.nlcut- spar.innerbw)/2;


% freq spacing
fftpts = spar.fs; % length 1ms for now, truncated later
spacing = spar.fs/fftpts; % 1 Hz spacing

% nyquist freq bin, depending on even or odd number of samples
nyBin = floor(fftpts/2) + 1;

% vector of freq per bin up to Nyquist freq
freqVec = [0:nyBin-1]' * spacing;

% lowest and highest bins
% lower flanking noise
lbin1 = max( round(spar.nlcut/spacing) + 1, 2);
hbin1 = min( round((spar.nlcut+bw_flank)/spacing) + 1, nyBin );

% higher flanking noise

correction = spar.innerbw==0; % so that in case of no inner band we don't double the middle bin
lbin2 = max( round((spar.nhcut-bw_flank)/spacing) + 1+correction, 2); % + 1 so that the edge bin isn't contained twice
hbin2 = min( round(spar.nhcut/spacing) + 1, nyBin );

% inner noise
lbin3 = max( round((spar.nlcut+bw_flank)/spacing) + 1, 2);
hbin3 = min( round(spar.nhcut-bw_flank/spacing) + 1, nyBin );

% create Gaussian white noises
a1 = zeros(fftpts,1);
b1 = a1;
a2 = zeros(fftpts,1);
b2 = a2;
a3 = a2;
b3 = b2;
a4 = a2;
b4 = b2;

% noises for each real and imag part: lower/upper flank, twice SDN
a1(lbin1:hbin1) = randn(hbin1-lbin1+1,1);
a2(lbin2:hbin2) = randn(hbin2-lbin2+1,1);
a3(lbin3:hbin3) = randn(hbin3-lbin3+1,1);
%     a4(lbin3:hbin3) = randn(hbin3-lbin3+1,1);

b1(lbin1:hbin1) = randn(hbin1-lbin1+1,1);
b2(lbin2:hbin2) = randn(hbin2-lbin2+1,1);
b3(lbin3:hbin3) = randn(hbin3-lbin3+1,1);


% complex noise spectra
fspec1 = a1+ 1i*b1;
fspec2 = a2+ 1i*b2;
fspec3 = a3+ 1i*b3;


fspec1_shift = fspec1;
fspec2_shift = fspec2;
fspec3_shift = fspec3;

% phase shifts flanking bands
if (hbin1-lbin1 >= 1) && (hbin2-lbin2 >= 1)
    fspec1_shift(lbin1:hbin1) = fspec1(lbin1:hbin1) .* exp(1i*spar.flanking_phase);
    fspec2_shift(lbin2:hbin2) = fspec2(lbin2:hbin2) .* exp(1i*spar.flanking_phase);
end


% phase shifts inner band: if flanking phase = pi, inner band phase = 0 and vice versa
if (hbin3-lbin3 >= 1)
    fspec3_shift(lbin3:hbin3) = fspec3(lbin3:hbin3) .* exp(1i*abs(spar.flanking_phase-pi));
end


flank1_l = (2*real(fftpts*ifft(fspec1)))/std(2*real(fftpts*ifft(fspec1)));
flank1_r =  (2*real(fftpts*ifft(fspec1_shift)))/std(2*real(fftpts*ifft(fspec1_shift)));

flank2_l = (2*real(fftpts*ifft(fspec2)))/std(2*real(fftpts*ifft(fspec2)));
flank2_r = (2*real(fftpts*ifft(fspec2_shift)))/std(2*real(fftpts*ifft(fspec2_shift)));

inner_l = (2*real(fftpts*ifft(fspec3))) / std(2*real( fftpts * ifft(fspec3) ));
inner_r = (2*real(fftpts*ifft(fspec3_shift)))/std(2*real(fftpts*ifft(fspec3_shift)));


flank_band_1 = [flank1_l(1:len) flank1_r(1:len)];
flank_band_2 = [flank2_l(1:len) flank2_r(1:len)];
inner_band   = [inner_l(1:len) inner_r(1:len)];


% calibration
[inner_band_spl]   = set_dbspl(inner_band,spar.spl + 10*log10(spar.innerbw));
[flank_band_1_spl]   = set_dbspl(flank_band_1,spar.spl + 10*log10(bw_flank));
[flank_band_2_spl]   = set_dbspl(flank_band_2,spar.spl + 10*log10(bw_flank));


noise_flanked_spl = inner_band_spl +flank_band_1_spl +flank_band_2_spl;


% add tone to noise
stimulus = tone_spl + noise_flanked_spl;

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

