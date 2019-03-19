% DEMO_GAMMATONE Demo for gammatone.m 
%
%   DEMO_GAMMATONE demonstrates the usage of various Gammatone filter
%   implementations in the AMT. 
%   
%   Figure 1: Gammatone IRs obtained with the `gammatone` implementation with the option classic (Patterson et al., 1987).
%
%   This figure shows IRs of Gammatone filters in 24 erb-spaced channels 
%   derived from gammatone with parameters classic and real (Patterson et al., 1987).
%   Left panel shows the causalphase option. Right panel shows the
%   peakphase option, which aligns the IRs to the maxima of the corresponding
%   envelope maxima.
%
%   Figure 2: Gammatone IRs obtained with the `gammatone` implementation with the option complex (Hohmann, 2002).
%
%   This figure shows IRs of Gammatone filters in 24 erb-spaced channels 
%   derived from gammatone with parameters allpole and complex (Hohmann,
%   2002). Left panel shows the causalphase option. Right panel shows the
%   peakphase option, which aligns the IRs to the maxima of the corresponding
%   envelope maxima.
%
%   Figure 3: Gammatone IRs obtained with the `gfb` implementation corresponding to complex-valued allpole filters (Hohmann, 2002).
%
%   This figure shows IRs of Gammatone filters in 24 erb-spaced channels 
%   derived from gfb_analyzer_new and other gfb function. This
%   implementation corresponds to complex-valued allpole Gammatone filters
%   (Hohmann, 2002). Left panel shows the IRs as output by
%   gfb_analyzer_new, the right panel shows those IRs manually delayed
%   such that their maxima are aligned across frequency channels.
%
%
%   See also: exp_gammatone gammatone exp_hohmann2002 demo_hohmann2002 hohmann2002
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/demos/demo_gammatone.php

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


% AUTHOR: Christian Klemenschitz, 2014

%% Figure gammatone, classic, Real

    % Parameters; 
    flow = 100;                     % Lowest center frequency in Hz;
    fhigh = 4000;                   % Highest center frequency in Hz;
    fs = 28000;                     % Sampling rate in Hz;
    fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
    nchannels = length(fc);         % Number of channels;
    N = 8192;                       % Number of samples;
    insig = [1, zeros(1,8191)];     % Impulse signal;
    treal = (1:N)/fs*1000;          % Time axis;
    
    %---- classic real ----
    % Derive filter coefficients and filter impulse responses; 
    [b,a] = gammatone(fc,fs,'classic','causalphase','real');
    outsig1 = 2*real(ufilterbankz(b,a,insig));
    outsig1 = permute(outsig1,[3 2 1]);
    
    % Derive filter coefficients and filter impulse responses with option 'peakphase';
    [b,a] = gammatone(fc,fs,'classic','peakphase','real');
    outsig2 = 2*real(ufilterbankz(b,a,insig));
    outsig2 = permute(outsig2,[3 2 1]);
    
    % Figure 1;
    type1   = 'classic / causalphase / real';    % Type of implemantion for headline;
    type2   = 'classic / peakphase / real';      % Type of implemantion for headline;
    
    % Plot
    % Figure 1 - classic real;
    figure('units','normalized','outerposition',[0 0.525 1 0.475]);
    set(gcf,'Name','gammatone, classic, real');
  %  subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig1(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type1)])
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    print(gcf,'-dpng','Gammatone')   
       
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig2(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type2)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    
%% Figure gammatone, classic, complex: This parameter configuration does not work and should be avoided
%   .. figure::
%
%      Classic gammatone implementation with complex-valued filter coefficients derived from gammatone.m.
%
%   This figure shows in the first plot an array of 24 erb-spaced channels of
%   classic gammatone implementation (Patterson et al., 1987) derived from
%   gammatone.m with complex-valued filter coefficients and in the second
%   plot this implementation with option 'peakphase', which makes the phase
%   of each filter be zero when the envelope of the impulse response of the
%   filter peaks (both do not scale correctly).
%
%     % Parameters; 
%     flow = 100;                     % Lowest center frequency in Hz;
%     fhigh = 4000;                   % Highest center frequency in Hz;
%     fs = 28000;                     % Sampling rate in Hz;
%     fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
%     nchannels = length(fc);         % Number of channels;
%     N = 8192;                       % Number of samples;
%     insig = [1, zeros(1,8191)];     % Impulse signal;
%     treal = (1:N)/fs*1000;          % Time axis;
%     
%     %---- classic complex ----
%     amt_disp('Classic, complex-valued, causal-phased gammtone implementation:')
%     % Derive filter coefficients and filter impulse responses;
%     [b,a] = gammatone(fc,fs,'classic','causalphase','complex');
%     outsig1 = 2*real(ufilterbankz(b,a,insig));
%     outsig1 = permute(outsig1,[3 2 1]);
% 
%     amt_disp('Classic, complex-valued, peak-phased gammtone implementation:')
%     % Derive filter coefficients and filter impulse responses with option 'peakphase';
%     [b,a] = gammatone(fc,fs,'classic','peakphase','complex');
%     outsig2 = 2*real(ufilterbankz(b,a,insig));
%     outsig2 = permute(outsig2,[3 2 1]);
%     
%     % Figure 2;
%     type1   = 'classic / causalphased / complex'; % Type of implemantion for headline;
%     type2   = 'classic / peakphased / complex';   % Type of implemantion for headline;
%     
%     % Plot
%     % Figure 2 - classic complex;
%     figure('units','normalized','outerposition',[0 0.05 1 0.475])
%     set(gcf,'Name','gammatone, classic, complex');
%     subplot(1,2,1)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,40*real(outsig1(ii,:)) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type1)])
%     xlabel 'Time (ms)'
%     ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
%     ylabel('# Frequency Channel (ERB): Frequency (Hz)');
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, 0, 26]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
%         ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
%         ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
%         ['#24: ' num2str(round(fc(24))) ' Hz']})
%     box on
%     hold off
%         
%     subplot(1,2,2)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,40*real(outsig2(ii,:)) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type2)])
%     xlabel 'Time (ms)'
%     ylabel('# Frequency Channel (ERB): Frequency (Hz)');
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, 0, 26]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
%         ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
%         ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
%         ['#24: ' num2str(round(fc(24))) ' Hz']})
%     box on
%     hold off
    
%% Figure gammatone, allpole, complex

    % Parameters; 
    flow = 100;                     % Lowest center frequency in Hz;
    fhigh = 4000;                   % Highest center frequency in Hz;
    fs = 28000;                     % Sampling rate in Hz;
    fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
    nchannels = length(fc);         % Number of channels;
    N = 8192;                       % Number of samples;
    insig = [1, zeros(1,8191)];     % Impulse signal;
    treal = (1:N)/fs*1000;          % Time axis;
    
    %---- allpole complex ----
    % Derive filter coefficients and filter impulse responses;
    [b,a] = gammatone(fc,fs,'allpole','causalphase','complex');
    outsig1 = 2*real(ufilterbankz(b,a,insig));
    outsig1 = permute(outsig1,[3 2 1]);

    % Derive filter coefficients and filter impulse responses with option 'peakphase';
    [b,a] = gammatone(fc,fs,'allpole','peakphase','complex');
    outsig2 = 2*real(ufilterbankz(b,a,insig));
    outsig2 = permute(outsig2,[3 2 1]);
    
    % Figure 3;
    type1   = 'allpole / causalphase / complex'; % Type of implemantion for headline;
    type2   = 'allpole / peakphase / complex';   % Type of implemantion for headline;
    
    % Plot
    % Figure 3 - allpole complex;
    figure('units','normalized','outerposition',[0 0.05 1 0.475]);
    set(gcf,'Name','gammatone, allpole, complex');
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig1(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type1)])
    xlabel 'Time (ms)'
    ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
        
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig2(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Gammatone impulse responses of type: ', num2str(type2)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    
%% Figure gfb implementation

    % Hohmann implementation;
    % Parameters;
    fs = 28000;                 % Sampling rate in Hz;
    flow = 100;                 % Lowest center frequency in Hz;
    basef = 888.44;             % Base center frequency in Hz;
    fhigh  = 4000;              % Highest center frequency in Hz;
    filters_per_ERBaud = 1;     % Filterband density on ERB scale;     
    
    % Construct new analyzer object;
    analyzer = hohmann2002(fs,flow, basef, fhigh,filters_per_ERBaud);
    % Impulse signal;
    impulse = [1, zeros(1,8191)];
    % Filter signal;
    [impulse_response, analyzer] = hohmann2002_process(analyzer, impulse);
       
    % Find peak at envelope maximum for lowest channel and add one sample;
    delay_samples = find(abs(impulse_response(1,:)) == max(abs(impulse_response(1,:)))) + 1;
    
    % 
    delay = hohmann2002_delay(analyzer, delay_samples);
     
    [outsig, delay] = hohmann2002_process(delay, impulse_response);
    
    % Figure 4;
    type1   = 'plotted as they are';   % Type of implemantion for headline;
    type2   = 'manually delayed to align their max';   % Type of implemantion for headline;
    
    % Plot
    % Figure 4 - gfb implementation, corresponds to allpole complex;
    figure ('units','normalized','outerposition',[0 0.05 1 0.475])
    set(gcf,'Name','gfb implementation');
    subplot(1,2,1)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(impulse_response(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['Impulse responses (IRs) from the gfb implementation: ', num2str(type1)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
    
    subplot(1,2,2)
    hold on
    dy = 1;
    for ii = 1:nchannels
       plot(treal,40*real(outsig(ii,:)) + dy)
       dy = dy + 1;
    end;
    clear dy;
    title (['IRs from the gfb implementation: ', num2str(type2)])
    xlabel 'Time (ms)'
    ylabel('# Frequency Channel (ERB): Frequency (Hz)');
    xt = 0:5:25;
    yt = [1 4 8 12 16 20 24];
    axis([0, 25, 0, 26]);
    set(gca,'XTick',xt, 'YTick', yt)
    set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
        ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
    box on
    hold off
   
%% Figure, gammatone classic, real, but with 6dbperoctave
%   .. figure::
%
%      Classic gammatone implementation with real-valued filter coefficients derived from gammatone.m with otpion '6dBperoctave'.
%
%   Figure 5 shows in the first plot an array of 24 erb-spaced channels of
%   classic gammatone implementation (Patterson et al., 1987) with real-valued
%   filter coefficients and option '6dBperoctave', which scales the amplitude
%   +/- 6 dB per octave with 0 dB at 4000 Hz and in the second plot this
%   implementation with option 'exppeakphase', which makes the phase of each
%   filter be zero when the envelope of the impulse response of the filter peaks.

%     % Parameters; 
%     flow = 100;                     % Lowest center frequency in Hz;
%     fhigh = 4000;                   % Highest center frequency in Hz;
%     fs = 28000;                     % Sampling rate in Hz;
%     fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
%     nchannels = length(fc);         % Number of channels;
%     N = 8192;                       % Number of samples;
%     insig = [1, zeros(1,8191)];     % Impulse signal;
%     treal = (1:N)/fs*1000;          % Time axis;
%     
%     
%     %---- classic real 6dBperoctave ----    
%     % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
%     [b,a] = gammatone(fc,fs,'classic','real','6dBperoctave');
%     outsig1 = 2*real(ufilterbankz(b,a,insig));
%     outsig1 = permute(outsig1,[3 2 1]);
%     
%     %---- classic real peakphase(new) 6dBperoctave ----
%     % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
%     [b,a] = gammatone(fc,fs,'classic','real','exppeakphase','6dBperoctave');
%     outsig2 = 2*real(ufilterbankz(b,a,insig));
%     outsig2 = permute(outsig2,[3 2 1]);
%     
%     % Figure 5;
%     type1   = 'classic / casualphased / real / +/-6dB';   % Type of implemantion for headline;
%     type2   = 'classic / peakphased(new) / real / +/-6dB';  % Type of implemantion for headline;
%     
%     % Plot
%     % Figure 5 - classic / causalphased / real / +/-6dB;
%     figure('units','normalized','outerposition',[0 0.525 1 0.475])
%     subplot(1,2,1)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,8*real(outsig1(ii,:)) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type1)])
%     xlabel 'Time (ms)'
%     ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
%     ylabel('# Frequency Channel (ERB): Frequency (Hz)');
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, -4, 26]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
%         ['#08: ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
%         ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
%         ['#24: ' num2str(round(fc(24))) ' Hz']})
%     box on
%     hold off
%         
%     subplot(1,2,2)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,8*real(outsig2(ii,:)) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type2)])
%     xlabel 'Time (ms)'
%     ylabel('# Frequency Channel (ERB): Frequency (Hz)');
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, -4, 26]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
%         ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
%         ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
%         ['#24: ' num2str(round(fc(24))) ' Hz']})
%     box on
%     hold off
% 

%% Figure gammatone allpole, complex, but with 6dbperoctave
%   .. figure::
%
%      Allpole gammatone implementation with real-valued filter coefficients derived from gammatone.m with otpion '6dBperoctave'.
%
%   Figure 6 shows in the first plot an array of 24 erb-spaced channels of
%   allpole gammatone implementation (Hohmann, 2002) with complex-valued
%   filter coefficients and option '6dBperoctave', which scales the amplitude
%   +/- 6 dB per octave with 0 dB at 4000 Hz and in the second plot this
%   implementation with option 'exppeakphase', which makes the phase of each
%   filter be zero when the envelope of the impulse response of the filter peaks.

% 
%     % Parameters; 
%     flow = 100;                     % Lowest center frequency in Hz;
%     fhigh = 4000;                   % Highest center frequency in Hz;
%     fs = 28000;                     % Sampling rate in Hz;
%     fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
%     nchannels = length(fc);         % Number of channels;
%     N = 8192;                       % Number of samples;
%     insig = [1, zeros(1,8191)];     % Impulse signal;
%     treal = (1:N)/fs*1000;          % Time axis;
%    
%     %---- allpole complex 6dBperoctave ----    
%     % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
%     [b,a] = gammatone(fc,fs,'allpole','complex','6dBperoctave');
%     outsig1 = 2*real(ufilterbankz(b,a,insig));
%     outsig1 = permute(outsig1,[3 2 1]);
%     
%     %---- allpole peakphase(new) complex 6dBperoctave ----
%     % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
%     [b,a] = gammatone(fc,fs,'allpole','complex','exppeakphase','6dBperoctave');
%     outsig2 = 2*real(ufilterbankz(b,a,insig));
%     outsig2 = permute(outsig2,[3 2 1]);
%     
%     % Figure 6;
%     type1   = 'allpole / causalphased  / complex / +/-6dB';   % Type of implemantion for headline;
%     type2   = 'allpole / peakphased(new) / complex / +/-6dB';  % Type of implemantion for headline;
%     
%     % Plot
%     % Figure 6 - classic / causalphased / real / +/-6dB;
%     figure('units','normalized','outerposition',[0 0.05 1 0.475])
%     subplot(1,2,1)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,8*real(outsig1(ii,:)) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type1)])
%     xlabel 'Time (ms)'
%     ylabel ([ num2str(nchannels),' channels / ', num2str(flow), ' - ',num2str(fhigh),' Hz'])
%     ylabel('# Frequency Channel (ERB): Frequency (Hz)');
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, -2, 26]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
%         ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
%         ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
%         ['#24: ' num2str(round(fc(24))) ' Hz']})
%     box on
%     hold off
%         
%     subplot(1,2,2)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,8*real(outsig2(ii,:)) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type2)])
%     xlabel 'Time (ms)'
%     ylabel('# Frequency Channel (ERB): Frequency (Hz)');
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, -2, 26]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
%         ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
%         ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
%         ['#24: ' num2str(round(fc(24))) ' Hz']})
%     box on
%     hold off
%     
%% Figure demo of manual delay
%   .. figure::
%
%      Classic with real-valued and allpole with complex-valued filter coefficients of gammatone implementation with otpion '6dBperoctave', peakphased and delayed.
%
%   Figure 7 shows in the first plot an array of 24 erb-spaced channels of
%   classic gammatone implementation (Patterson et al., 1987) with real-valued
%   filter coefficients and option '6dBperoctave', which scales the amplitude
%   +/- 6 dB per octave with 0 dB at 4000 Hz and option 'exppeakphase', which
%   makes the phase of each filter be zero when the envelope of the impulse
%   response of the filter peaks. Further the filters are delayed so the peaks
%   of each filter are arranged above each other. The second plot shows the
%   same but with allpole implementation (Hohmann, 2002) with complex-valued
%   filter coefficients derived from gammatone.m.   

% 
%     % Parameters; 
%     flow = 100;                     % Lowest center frequency in Hz;
%     fhigh = 4000;                   % Highest center frequency in Hz;
%     fs = 28000;                     % Sampling rate in Hz;
%     fc = erbspacebw(flow,fhigh);    % 24 erb-spaced channels;
%     nchannels = length(fc);         % Number of channels;
%     N = 8192;                       % Number of samples;
%     insig = [1, zeros(1,8191)];     % Impulse signal;
%     treal = (1:N)/fs*1000;          % Time axis;
%     
%     %---- classic real 6dBperoctave peakphased and delayed ----    
%     % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
%     [b,a] = gammatone(fc,fs,'classic','real','exppeakphase','6dBperoctave');
%     outsig1 = 2*real(ufilterbankz(b,a,insig));
%     outsig1 = permute(outsig1,[3 2 1]);
%     
%     % Find peak at maximum
%     envmax1 = zeros(1,nchannels);
%     for ii = 1:nchannels
%         % maximum per channel
%         envmax1(ii) = find(abs(outsig1(ii,:)) == max(abs(outsig1(ii,:))));
%     end;
%         
%     % Delay signal
%     for ii = 1:nchannels
%         % Time to delay
%         delay =  zeros(1,max(envmax1)+1 - envmax1(ii));
%         % Add delay
%         outsig1(ii,:) = [delay outsig1(ii,1:N - length(delay))]; 
%     end;
%     
%     % Type of implemantion for headline;
%     type1   = 'classic / exppeakphased / real / 6dbperoctave / manually delayed'; 
% 
%     % Plot 
%     figure ('units','normalized','outerposition',[0 0.05 1 0.475])
%     subplot(1,2,1)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,8*real(outsig1(ii,:)) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type1)])
%     xlabel 'Time (ms)'
%     ylabel('# Frequency Channel (ERB): Frequency (Hz)');
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, -4, 26]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
%         ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
%         ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
%         ['#24: ' num2str(round(fc(24))) ' Hz']})
%     box on
%     hold off
%     
%     %---- allpole complex 6dBperoctave peakphased and delayed ----
%     % Derive filter coefficients and filter impulse responses with option '6dBperoctave';
%     [b,a] = gammatone(fc,fs,'allpole','complex','exppeakphase','6dBperoctave');
%     outsig2 = 2*real(ufilterbankz(b,a,insig));
%     outsig2 = permute(outsig2,[3 2 1]);
%     
%     % Find peak at envelope maximum
%     envmax2 = zeros(1,nchannels);
%     for ii = 1:nchannels
%         % Envelope maximum per channel
%         envmax2(ii) = find(abs(outsig2(ii,:)) == max(abs(outsig2(ii,:))));
%     end;
%         
%     % Delay signal
%     for ii = 1:nchannels
%         % Time to delay
%         delay =  zeros(1,max(envmax2)+1 - envmax2(ii));
%         % Add delay
%         outsig2(ii,:) = [delay outsig2(ii,1:N - length(delay))]; 
%     end;
%     
%     % Type of implemantion for headline;
%     type2   = 'allpole / peakphased(new) / complex / +/-6dB / delayed';   
%     
%     % Plot;    
%     subplot(1,2,2)
%     hold on
%     dy = 1;
%     for ii = 1:nchannels
%        plot(treal,8*real(outsig2(ii,:)) + dy)
%        dy = dy + 1;
%     end;
%     clear dy;
%     title (['Gammatone impulse responses of type: ', num2str(type2)])        
%     xlabel 'Time (ms)'
%     ylabel('# Frequency Channel (ERB): Frequency (Hz)');
%     xt = 0:5:25;
%     yt = [1 4 8 12 16 20 24];
%     axis([0, 25, -4, 26]);
%     set(gca,'XTick',xt, 'YTick', yt)
%     set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
%         ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))),' Hz'], ... 
%         ['#16: ' num2str(round(fc(16))),' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
%         ['#24: ' num2str(round(fc(24))) ' Hz']})
%     box on
%     hold off
%     
