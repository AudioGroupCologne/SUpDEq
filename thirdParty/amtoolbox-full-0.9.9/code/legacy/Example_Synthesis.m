% This Example demonstrates how to create and how to use the combined
% analysis-synthesis Filterbank system.
%
% FIXME: This example seems to be a demo for hohmann2002. Check it.
% copyright: Universitaet Oldenburg
% author   : tp
% date     : Jan, Mar 2002, Nov 2003, Nov 2006
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/legacy/Example_Synthesis.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%%% First, create a filterbank analyzer as in Example_Filterbank.m %%%

flow =    70;
fhigh =  6700;
base_frequency_hz         =  1000;
sampling_rate_hz          = 16276;
filters_per_ERB           =     1.0;
desired_delay_in_seconds  =     0.004;
filter_order              =     4;
bandwidth_factor          =     1.0;

amt_disp(['Building analysis filterbank']);
analyzer = hohmann2002(sampling_rate_hz, flow, ...
                            base_frequency_hz, fhigh,...
			    filters_per_ERB, filter_order, bandwidth_factor);


%%% Now create a synthesizer that can resynthesize the analyzer's output %%%

amt_disp(['Building synthesizer for an analysis-synthesis delay of ', ...
      num2str(desired_delay_in_seconds), ' seconds']);
synthesizer = hohmann2002_synth(analyzer, desired_delay_in_seconds);

%%% Extract the synthesizer's parameters %%%
amt_disp(['The synthesizers parameters:';...
      '----------------------------']);
delay = synthesizer.delay;
mixer = synthesizer.mixer;

bands = length(mixer.gains);

fprintf(1,'%3s|%7s | %22s | %5s\n\n', ...
        '# ', 'delay ', 'phase factor    ', 'gain / dB');

for band = 1:bands
  fprintf(1,'%3d|%7d | %9f + %9fi | %5.2f\n', ...
	  band,                            ...
	  delay.delays_samples(band),      ...
	  real(delay.phase_factors(band)), ...
	  imag(delay.phase_factors(band)), ...
	  20*log10(mixer.gains(band)));
end

%%%  plot the resynthesized impulse and the frequency response of the  %%%
%%%  analysis-synthesis system                                         %%%

impulse = [1, zeros(1,8191)];                                          
[analyzed_impulse, analyzer] = hohmann2002_process(analyzer, impulse);
[resynthesized_impulse, synthesizer] = ...
    hohmann2002_process(synthesizer, analyzed_impulse);

figure(1);
plot([0:8191]/sampling_rate_hz*1e3, resynthesized_impulse);
axis([40/sampling_rate_hz*1e3, 120/sampling_rate_hz*1e3, -1, 1]);
title('impulse response of the analysis-synthesis system');
xlabel('time / ms');
ylabel('system output');

frequency = [0:8191] * sampling_rate_hz / 8192;
figure(2)
plot(frequency, 20 * log10(abs(fft(resynthesized_impulse'))));
axis([0, sampling_rate_hz/2, -40, 5]);
title('frequency response of the analysis-synthesis-system');
xlabel('frequency / Hz');
ylabel('system response level / dB');

amt_disp('Figure 1 shows the impulse response of the analysis-synthesis system');
amt_disp('in the time domain; figure 2 shows its frequency response.');

%OLDFORMAT

