function exp_gammatone(varargin)
%EXP_GAMMATONE  Creates various figures related to the Gammatone filters
%   Usage: exp_gammatone(flags)
%
%   EXP_GAMMATONE(flags) reproduces some figures from  
%   Patterson et al. (1987), Lyon (1997), and Hohmann (2002)
%
%   The following flags can be specified;
%
%     'fig1_patterson1987'
%
%                   Reproduce Fig.1 from Patterson et al. (1987):
%                   The amplitude characteristics of three roex(p) filters
%                   centered at 427 Hz, 1000 Hz and 2089 Hz. The lower and
%                   upper filters are centered 6 ERBs below and above the
%                   1000 Hz filter respectively. In each case, the range of
%                   the abscissa extends from an octave below to an octave
%                   above the centre frequency of the filter, on a linear
%                   frequency scale. The range of the ordinate is 40dB.
%
%
%     'fig2_patterson1987'
%
%                   Reproduce Fig.2 from Patterson et al. (1987):
%                   An array of 24 impulse responses for roex(p) filters
%                   whose centre frequencies ranges from 100 to 4000 Hz.
%                   The linear-phase assumption leads to symmetric impulse
%                   responses which have been aligned at their temporal
%                   mid-points.
%
%
%     'fig3_patterson1987'
%
%                   Reproduce Fig.3 similiar to Patterson et al. (1987):
%                   A cochleagram of four cycles of the [ae] in "past"
%                   produced by by a gammatone filterbank without phase 
%                   compensation. The triangular objects are the upper
%                   three formants of the vowel. The duration of each
%                   period is ~8 ms. The ordinate is filter centre
%                   frequency of 189 equally spaced channels from 100 to
%                   4000 Hz. Note the strong rightward skew induced  by the
%                   phase lags of the low-frequency filters in the lower
%                   half of the figure.
%                       
%
%     'fig4_patterson1987'
%
%                   Reproduce Fig.4 similiar to Patterson et al. (1987):
%                   A cochleagram of four cycles of the [ae] in "past"
%                   produced by a gammatone filterbank with phase compensation.
%                   The coordinates are the same as for Figure 3. Note that
%                   the strong rightward skew produced by the phase lags of
%                   the lower frequency filters has been removed now.
%
%
%     'fig5_patterson1987'
%
%                   Reproduce Fig.5 from Patterson et al. (1987):
%                   An array of gamma impulse responses for a 24-channel
%                   auditory filterbank.(lower portion), and the equvalent
%                   array of gamma envelopes (upper portion). The range of
%                   abcissa is 25 ms; the filter centre frequencies range 
%                   from 100 - 4,000 Hz.
%
%
%     'fig6_patterson1987'
%
%                   Reproduce Fig.6 from Patterson et al. (1987):
%                   A comparison of the gammatone(4,b) and roex(p) filters
%                   at three centre frequencies, 0.43, 1.00 and 2.09 kHz.
%                   In this case, the gammatone filter has been matched to
%                   the roex filter by equating the ERB, thus minimising
%                   the difference in the area under the curves. The range
%                   of the ordinate is 60 dB, the abscissa ranges from an
%                   octave below to an octave above the centre frequency in
%                   each case.
%
%
%     'fig7_patterson1987'
%
%                   Reproduce Fig.7 from Patterson et al. (1987):
%                   The comparison of the gammatone(4,b) and roex(p)
%                   filters at three centre frequencies (0.43, 1.00 and
%                   2.09 kHz). In this case, the bandwidth of the gammatone
%                   filter has been increased by 10% to minimise the
%                   decibel difference between it and the roex filter. The
%                   range of the ordinate is 60 dB, and the abscissa ranges
%                   from an octave below to an octave above the centre
%                   frequency of the filter in each case.
%
%
%     'fig8_patterson1987'
%
%                   Reproduce Fig.8 from Patterson et al. (1987):
%                   A comparison of the gammatone(2,b) and the roex(p,w,t)
%                   filters at three centre frequencies (0.43, 1.00 and
%                   2.09 kHz). The parameters for the roex filter are
%                   taken from Patterson et al (1982). The gammatone filter
%                   has been fitted to the roex by equating their ERBs.
%                   The range of the ordinate is 50dB, and the abscissa ranges
%                   from an octave below to an octave above the centren frequency,
%                   in each case.
%
%
%     'fig9_patterson1987'
%
%                   Reproduce Fig.9 from Patterson et al. (1987):
%                   A comparison of the gammatone(2,b) and the roex(p,w,t)
%                   filters at three centre frequencies (0.43, 1.00 and
%                   2.09 kHz). In this case the roex parameters, w and t,
%                   have been adjusted to improve the fit to the
%                   gammatone(2,b) filter to show that the discrepancy can
%                   easily be minimised. The range of the ordinate is 50 dB,
%                   the abscissa shows a range from an octave below to an octave
%                   above the centre frequency of the filter in each case.
%
%                   
%     'fig10a_patterson1987'
%                   
%                   Reproduce Fig.10a from Patterson et al. (1987):
%                   The impulse responses for a gammatone auditory filterbank
%                   without phase compensation. The filterbank contains 37 
%                   channels ranging from 100 to 5,000 Hz. The range of abcissa
%                   is 25 ms.
%
%
%     'fig10b_patterson1987'
%                   
%                   Reproduce Fig.10b from Patterson et al. (1987):
%                   The impulse responses for a gammatone auditory
%                   filterbank with envelope phase-compensation; that is,
%                   the peaks of the impulse-response envelopes have been
%                   aligned vertically. The filterbank contains 37 channels
%                   ranging from 100 to 5,000 Hz. The range of the abscissa
%                   is 25 ms.
%
%
%     'fig10c_patterson1987'
%                   
%                   Reproduce Fig.10c from Patterson et al. (1987):
%                   The impulse responses for a gammatone auditory filterbank
%                   with envelope and fine structure phase-compensation; that is,
%                   the envelope peaks have been aligned and then a fine structure
%                   peak has been aligned with the envelope peak. The filterbank
%                   contains 37 channels ranging from 100 to 5,000 Hz.
%                   The range of abcissa is 25 ms.
%
%
%     'fig11_patterson1987'
%
%                   Reproduce similiar to Fig.11 from Patterson et al. (1987):
%                   The impulse responses of a gammatone filterbank from a
%                   pulsetrain.
%
%
%     'fig1_lyon1997'
%
%                   Reproduce parts of Fig.1 from Lyon (1997) (no DAPGF):
%                   Comparison of GTF and APGF transfer function for two different
%                   values of the real part of the pole position. Notice that
%                   the ordering of gain near the peak for the classic gammatone
%                   filter is not maintained in the tail. The variation in the
%                   GTF tail is not accounted for by the usual phase-independent
%                   symetric GTF approximation.
%
%
%     'fig2_lyon1997'
%
%                   Reproduce parts of Fig.2 from Lyon (1997) (no DAPGF):
%                   Impulse responses of the APGF and sine-phased GTF from 
%                   Figure 1. Note that the GTF's zero crossings are equally
%                   spaced in time, while those of the APGF are stretched out
%                   in early cycles.
%
%
%     'fig1_hohmann2002'
%
%                   Reproduce Fig.l from Hohmann (2002):
%                   Impulse response of the example Gammatone filter
%                   (center frequency fc = 1000 Hz; 3-db bandwidth fb = 100 Hz;
%                   sampling frequency fs = 10kHz). Solid and dashed lines
%                   show the real and imaginary part of the filter output,
%                   respectively. The absolute value of the filter output
%                   (dashdotted line) clearly represents the envelope.
%
%
%     'fig2_hohmann2002'
%
%                   Reproduce Fig.2 from Hohmann (2002):
%                   Frequency response of the example Gammatone filter
%                   (upper two panels) and of the real-to-imaginary
%                   response(lower two panels). Pi/2 was added to the phase
%                   of the latter (see text). The frequency axis goes up to
%                   half the sampling rate (z=pi).
%       
%
%     'fig3_hohmann2002'
%
%                   Reproduce Fig.3 from Hohmann (2002):
%                   Magnitude frequency response of the Gammatone
%                   filterbank. In this example, the filter channel density
%                   is 1 on the ERB scale and the filter bandwidth is 1
%                   ERBaud. The sampling frequency was 16276Hz and the
%                   lower and upper boundary for the center frequencies
%                   were 70Hz and 6.7kHz, respectively.
%
%
%     'fig4_hohmann2002'
%
%                   Treatment of an impulse response with envelope maximum
%                   to the left of the desired group delay (at sample 65).
%                   The original complex impulse response (real part and
%                   envelope plotted in the upper panel) is multiplied with
%                   a complex factor and delayed so that the envelope maximum 
%                   and the maximum of the real part coincide with the desired
%                   group delay (lower panel).
%
%
%   Examples:
%   ---------
%
%   To display Fig. 1 from Patterson et al. (1987) use :
%
%     exp_gammatone('fig1_patterson1987');
%
%   To display Fig. 2 Patterson et al. (1987) use :
%
%     exp_gammatone('fig2_patterson1987');
%
%   To display Fig. 3 Patterson et al. (1987) use :
%
%     exp_gammatone('fig3_patterson1987');
%
%   To display Fig. 4 Patterson et al. (1987) use :
%
%     exp_gammatone('fig4_patterson1987');
%
%   To display Fig. 5 Patterson et al. (1987) use :
% 
%     exp_gammatone('fig5_patterson1987');
%
%   To display Fig. 6 Patterson et al. (1987) use :
% 
%     exp_gammatone('fig6_patterson1987');
%
%   To display Fig. 7 Patterson et al. (1987) use :
% 
%     exp_gammatone('fig7_patterson1987');
%
%   To display Fig. 8 Patterson et al. (1987) use :
% 
%     exp_gammatone('fig8_patterson1987');
%
%   To display Fig. 9 Patterson et al. (1987) use :
% 
%     exp_gammatone('fig9_patterson1987');
%
%   To display Fig. 10a Patterson et al. (1987) at use :
%
%     exp_gammatone('fig10a_patterson1987');
%
%   To display Fig. 10b Patterson et al. (1987) at use :
%
%     exp_gammatone('fig10a_patterson1987');
%
%   To display Fig. 10c Patterson et al. (1987) at use :
%
%     exp_gammatone('fig10c_patterson1987');
%
%   To display Fig. 11 Patterson et al. (1987) use :
%
%     exp_gammatone('fig11_patterson1987');
%
%   To display Fig. 1 Lyon (1997) use :
%
%     exp_gammatone('fig1_lyon1997');
%
%   To display Fig. 2 Lyon (1997) use :
%
%     exp_gammatone('fig2_lyon1997');
%
%   To display Fig. 1 Hohmann (2002) use :
%
%     exp_gammatone('fig1_hohmann2002');
%
%   To display Fig. 2 Hohmann (2002) use :
%
%     exp_gammatone('fig2_hohmann2002');
%
%   To display Fig. 3 Hohmann (2002) use :
%
%     exp_gammatone('fig3_hohmann2002');
%
%   To display Fig. 4 Hohmann (2002) use :
%
%     exp_gammatone('fig4_hohmann2002');
%
%   References:
%     V. Hohmann. Frequency analysis and synthesis using a gammatone
%     filterbank. Acta Acustica united with Acoustica, 88(3):433-442, 2002.
%     
%     R. Lyon. All pole models of auditory filtering. In E. R. Lewis, G. R.
%     Long, R. F. Lyon, P. M. Narins, C. R. Steele, and E. Hecht-Poinar,
%     editors, Diversity in auditory mechanics, pages 205-211. World
%     Scientific Publishing, Singapore, 1997.
%     
%     B. Moore and B. Glasberg. Suggested formulae for calculating
%     auditory-filter bandwidths and excitation patterns. J. Acoust. Soc.
%     Am., 74:750-753, 1983.
%     
%     R. Patterson, I. Nimmo-Smith, J. Holdsworth, and P. Rice. An efficient
%     auditory filterbank based on the gammatone function. APU report, 2341,
%     1987.
%     
%
%   See also: gammatone exp_hohmann2002 demo_gammatone demo_hohmann2002 hohmann2002 hohmann2002_process
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/experiments/exp_gammatone.php

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

% AUTHOR: CK, 2014

%% ------ Check input options --------------------------------------------
    definput.import={'amtredofile'};
    definput.keyvals.FontSize = 12;
    definput.keyvals.MarkerSize = 6;
    definput.flags.type = {'missingflag', 'fig1_patterson1987', 'fig2_patterson1987', ...
    'fig3_patterson1987', 'fig4_patterson1987', 'fig5_patterson1987', 'fig6_patterson1987', ...
    'fig7_patterson1987', 'fig8_patterson1987', 'fig9_patterson1987', ...
    'fig10a_patterson1987', 'fig10b_patterson1987', 'fig10c_patterson1987', ...
    'fig11_patterson1987', 'fig1_lyon1997', 'fig2_lyon1997',...
    'fig1_hohmann2002', 'fig2_hohmann2002', 'fig3_hohmann2002', 'fig4_hohmann2002'};

    % Parse input options
    [flags,~]  = ltfatarghelper({'FontSize','MarkerSize'},definput,varargin);

    if flags.do_missingflag
        flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
                   sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
        error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
    end;

%% Pattersons Paper - Roex(p) Filter
% Figure 1; 

   if flags.do_fig1_patterson1987
       
      % Global parameters:
      fs = 25000;           % Sampling frequency in Hz; 
 
      % Roex(p) at 1000 Hz;
      % Parameters:
      fc = 1000;                                % Center frequency in Hz; 
      df = 1/fc;                                % Distance on frequency scale in Hz;
      g = (0:df:1-df);                          % Frequency scale;
      erb = audfiltbw(fc);                      % Equivalent rectangular bandwidth;
      % Equation from paper Moore, 1983 found in text after Equation 6.
      p = 4*fc/erb;                             % Filter shape parameter;

      % Equation 6 from paper Moore, 1983. (slopes are mirrored);
      weight_l = (1+p.*g).*exp(-p.*g);          % for lower slope;
      weight_u = (1+p.*g).*exp(-p.*g);          % for upper slope;
      rp_y = [fliplr(weight_l) 1 weight_u];     % Filteramplitude;  
      rp_x = ([-fliplr(g) 0 g]+1);              % Frequency-axis;

      % Plot;
      figure('units','normalized','outerposition',[0 0.05 1 0.95])
      subplot(1,3,2)
      plot(rp_x, 10*log10(rp_y), 'b')
      hold on
      title('Roex(p) filters')
      xlabel(['Roex(p) filter at ', num2str(fc), ' Hz center frequency'])
      ylabel('Attenuation (dB) ')
      xt = 0:0.2:2;
      yt = -40:10:0;
      axis([0 2 -40 0])
      set(gca,'XTick',xt)
      set(gca,'YTick',yt)
      box on
      hold off


      % Roex(p) at 427 Hz;
      % Parameters:
      fc = 427;                                 % Center frequency in Hz; 
      df = 1/fc;                                % Distance on frequency scale in Hz;
      g = (0:df:1-df);                          % Frequency scale;
      erb = audfiltbw(fc);                      % Equivalent rectangular bandwidth;
      % Equation from paper Moore, 1983 found in text after Equation 6.
      p = 4*fc/erb;                             % Filter shape;

      % Equation 6 from paper Moore, 1983 (slopes are mirrored);
      weight_l = (1+p.*g).*exp(-p.*g);          % for lower slope;
      weight_u = (1+p.*g).*exp(-p.*g);          % for upper slope;
      rp_y = [fliplr(weight_l) 1 weight_u];     % Filteramplitude;
      rp_x = ([-fliplr(g) 0 g]+1);              % Frequency-axis;

      % Plot;
      subplot(1,3,1)
      plot(rp_x, 10*log10(rp_y), 'b')
      hold on
       xlabel(['Roex(p) filter at ', num2str(fc), ' Hz center frequency'])
      ylabel('Attenuation (dB) ')
      xt = 0:0.2:2;
      yt = -40:10:0;
      axis([0 2 -40 0])
      set(gca,'XTick',xt)
      set(gca,'YTick',yt)
      box on
      hold off


      % Roex(p) at 2089 Hz;
      % Parameters:
      fc = 2089;                                % Center frequency in Hz; 
      df = 1/fc;                                % Distance on frequency scale in Hz;
      g = (0:df:1-df);                          % Frequency scale;
      erb = audfiltbw(fc);                      % Equivalent rectangular bandwidth;
      % Equation from paper Moore, 1983 found in text after Equation 6.
      p = 4*fc/erb;                             % Filter shape;


      % Equation 6 from paper Moore, 1983. (slopes are mirrored)
      weight_l = (1+p.*g).*exp(-p.*g);          % for lower slope;
      weight_u = (1+p.*g).*exp(-p.*g);          % for upper slope;
      rp_y = [fliplr(weight_l) 1 weight_u];     % Filteramplitude;
      rp_x = ([-fliplr(g) 0 g]+1);              % Frequency-axis;

      % Plot;
      subplot(1,3,3)
      plot(rp_x, 10*log10(rp_y), 'b')
      hold on
      xlabel(['Roex(p) filter at ', num2str(fc), ' Hz center frequency'])
      ylabel('Attenuation (dB) ')
      xt = 0:0.2:2;
      yt = -40:10:0;
      axis([0 2 -40 0])
      set(gca,'XTick',xt)
      set(gca,'YTick',yt)
      box on
      hold off
   end;
       
%% Pattersons Paper - Roex(p) Filter
% Figure 2 

  if flags.do_fig2_patterson1987
      
      % Parameters:  
      fs = 20000;                       % Sampling frequency in Hz;
      treal = 1:250;                    % Samples x-axis
      fn = fs/2;                        % Nyquist frequency;
      f = 0:fn;                         % Frequency vector
      flow = 100;                       % ERB lowest center frequency Hz;
      fhigh = 4000;                     % ERB highest center frequency Hz;
      fc = erbspacebw(flow,fhigh);      % 24 erb spaced channels;
      nchannels = length(fc);           % Number of channels
       
      g = zeros(nchannels,fn+1);
      erb = zeros(1,nchannels);
      p = zeros(1,nchannels);
      weight = zeros(nchannels,fn+1);
      yfft = zeros(nchannels,fs);
      yifft = zeros(nchannels,fs);
      yshift = zeros(nchannels,fs);
      for ii = 1:nchannels
          g(ii,:) = abs(f(:)-fc(ii))./fc(ii);   % Normalized distance from filter center
          % Equation from paper Moore, 1983 found in text after Equation 6
          erb(ii) = audfiltbw(fc(ii));          % Equivalent rectangular bandwidth;
          p(ii) = 4*fc(ii)/erb(ii);             % Filter shape
          % Equation 6 from paper Moore, 1983 (slopes are mirrored);
          weight(ii,:) = (1+p(ii).*g(ii,:)).*exp(-p(ii).*g(ii,:));      % Filter slope;
          yfft(ii,:) = [weight(ii,:) fliplr(weight(ii,(2:end-1)))];     % Filter read for IFFT
          yifft(ii,:) = (ifft(yfft(ii,:)));                             % IFFT
          yshift(ii,:)= circshift(yifft(ii,:),[0 size(yifft,2)/2]);     % Circular shift 
      end;

      % Plot
      figure ('units','normalized','outerposition',[0 0.05 1 0.95])
      dy = 1;
      for ii = 1:nchannels
          plot(treal,25*yshift(ii,9875:10124)+dy);
          hold on
          dy = dy+1;
      end;
      title('Roex(p) filters')
      xlabel('Sample')
      ylabel('# Frequency Channel (ERB): Frequency (Hz)');
      yt = [1 4 8 12 16 20 24];
      axis([0, 250, 0, 26]);
      set(gca,'YTick', yt)
      set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))) ' Hz'], ... 
        ['#16: ' num2str(round(fc(16))) ' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
      box on
      hold off
  end;
  
%% Pattersons Paper - Classic Gammatone Filter
% Figure 3 
% gammatone 'classic' (causalphase, real)

   if flags.do_fig3_patterson1987
        % Input parameters
        x=amt_load('gammatone','a_from_past_uk.mat'); % Load mat-file with vowel a in variable insig and sampling frequency in variable fs;
        fs=x.fs;
        ts = 1/fs;                  % Time between sampling points in s;
        insig = x.data(:,1).';        % Extract one channel from signal;
        insig = insig(610:963);     % Extract ~8 ms;
        insig = insig-mean(insig);  % Removes DC;
        insig = [insig insig insig insig insig insig insig insig]; % 8 cycles;
        nchannels = 189;            % Number of channels in filterbank;
        N=length(insig);            % Length of signal;
        treal = (1:N)*ts*1000;      % Time axis in ms;
        flow = 150;                 % Lowest center frequency in Hz;
        fhigh = 4000;               % Highest center frequency in Hz;
        fc = erbspace(flow,fhigh,nchannels);   % 189 erb-spaced channels; 
        
        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,'classic');
        % Filters impulse signal with filter coefficients from above.
        outsig = real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Plot;
        figure ('units','normalized','outerposition',[0 0.05 1 0.95])
        hold on
        dy = 1;
        for ii = 1:nchannels
           plot(treal,24*outsig(ii,:)+dy)
           dy = dy+1;
        end;
        title ([num2str(nchannels), ' channel array of gammatone impulse responses of vowel a from past'])
        xlabel 'Time (ms)'
        ylabel('# Frequency Channel (ERB): Frequency (Hz)');
        set(gca, 'Xlim', [20 56], 'Ylim', [1 189])
        yt = [1 20 40 60 80 100 120 140 160 180 189];
        set(gca, 'YTick', yt )
        set(gca,'YTickLabel', {['#001:   '  num2str(round(fc(1))) ' Hz'] ,['#020:   '  num2str(round(fc(20))) ' Hz'] , ...
        ['#040:   '  num2str(round(fc(40))) ' Hz'] , ['#060:   ' num2str(round(fc(60))) ' Hz'], ...
        ['#080:   ' num2str(round(fc(80))) ' Hz'], ['#100: ' num2str(round(fc(100))) ' Hz'], ... 
        ['#120: ' num2str(round(fc(120))) ' Hz'], ['#140: ' num2str(round(fc(140))) ' Hz'], ...
        ['#160: ' num2str(round(fc(160))) ' Hz'], ['#180: ' num2str(round(fc(180))) ' Hz'],['#189: ' num2str(round(fc(189))) ' Hz']})
        box on
        hold off
   end;
%% Pattersons Paper - Classic Gammatone Filter
% Figure 4 
% gammatone 'classic' (causalphase, real)

   if flags.do_fig4_patterson1987
        % Input parameters
        x=amt_load('gammatone','a_from_past_uk.mat'); % Load mat-file with vowel a in variable insig and sampling frequency in variable fs;
        fs=x.fs;
        ts = 1/fs;                  % Time between sampling points in s;
        insig = x.data(:,1).';        % Extract one channel from signal;
        insig = insig(610:963);     % Extract ~8 ms;
        insig = insig-mean(insig);  % Removes DC;
        insig = [insig insig insig insig insig insig insig insig]; % 8 cycles;
        nchannels = 189;            % Number of channels in filterbank;
        N=length(insig);            % Length of signal;
        treal = (1:N)*ts*1000;      % Time axis in ms;
        flow = 150;                 % Lowest center frequency in Hz;
        fhigh = 4000;               % Highest center frequency in Hz;
        fc = erbspace(flow,fhigh,nchannels);   % 189 erb-spaced channels; 
        
        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,'classic','exppeakphase');
        % Filters impulse signal with filter coefficients from above.
        outsig = real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Impulse response of filter;
        impulse = [1 zeros(1,352)];
        outimp = real(ufilterbankz(b,a,impulse));
        outimp = permute(outimp,[3 2 1]);
        
        % Find impulse response envelopes maximums;
        envmax = zeros(1,nchannels);
        for ii = 1:nchannels
            envmax(ii) = find( abs(outimp(ii,:)) == max(abs(outimp(ii,:))));
        end;
        
        % Delay output so the impulse response envelopes peak above each other;
        for ii = 1:nchannels
            % Time to delay
            delay =  zeros(1,(max(envmax(:))+1) - envmax(ii));
            % Add delay
            outsig(ii,:) = [delay outsig(ii,1:N - length(delay))]; 
        end;
        
        % Plot;
        figure ('units','normalized','outerposition',[0 0.05 1 0.95])
        hold on
        dy = 1;
        for ii = 1:nchannels
           plot(treal,24*outsig(ii,:)+dy)
           dy = dy+1;
        end;
        title ([num2str(nchannels), ' channel array of gammatone impulse responses of vowel a from past'])
        xlabel 'Time (ms)'
        ylabel('# Frequency Channel (ERB): Frequency (Hz)');
        set(gca, 'Xlim', [20 56], 'Ylim', [1 189])
        yt = [1 20 40 60 80 100 120 140 160 180 189];
        set(gca, 'YTick', yt )
        set(gca,'YTickLabel', {['#001:   '  num2str(round(fc(1))) ' Hz'] ,['#020:   '  num2str(round(fc(20))) ' Hz'] , ...
        ['#040:   '  num2str(round(fc(40))) ' Hz'] , ['#060:   ' num2str(round(fc(60))) ' Hz'], ...
        ['#080:   ' num2str(round(fc(80))) ' Hz'], ['#100: ' num2str(round(fc(100))) ' Hz'], ... 
        ['#120: ' num2str(round(fc(120))) ' Hz'], ['#140: ' num2str(round(fc(140))) ' Hz'], ...
        ['#160: ' num2str(round(fc(160))) ' Hz'], ['#180: ' num2str(round(fc(180))) ' Hz'],['#189: ' num2str(round(fc(189))) ' Hz']})
        box on
        hold off
   end;   
%% Pattersons Paper - Classic Gammatone Filter
% Figure 5 
% gammatone 'classic' (causalphase, real)

    if flags.do_fig5_patterson1987
        % Input parameters
        fs = 25000;                     % Sampling frequency in Hz;
        ts = 1/fs;                      % Time between sampling points in s
        N=4096;                         % Length of signal;
        treal = (1:N)*ts*1000;          % Time axis in ms;
        flow = 100;                     % ERB lowest center frequency Hz;
        fhigh = 4000;                   % ERB highest center frequency Hz;
        fc = erbspacebw(flow,fhigh);    % Array of centre frequencies in Hz;
        nchannels = length(fc);         % Number of channels
        insig = zeros(1,N);             % Impulse as input signal;
        insig(1) = 1;
        
        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,'classic');
        % Filters impulse signal with filter coefficients from above.
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Lower figure;
        figure('units','normalized','outerposition',[0.5 0.1 0.5 0.9])
        subplot(2,1,2)
        hold on
        dy = 1;
        for ii = 1:size(outsig,1)
           plot(treal,14*outsig(ii,:) + dy)
           dy = dy + 1;
        end;
        clear dy;
        title ([ num2str(nchannels),' channel array of gamma impulse responses (classic)'])
        xlabel 'Time (ms)'
        ylabel('# Frequency Channel (ERB): Frequency (Hz)');
        xt = 0:5:27;
        yt = [1 4 8 12 16 20 24];
        axis([0, 25, 0, 26]);
        set(gca,'XTick',xt, 'YTick', yt)
        set(gca,'YTickLabel', {['#01:  '  num2str(round(fc(1))) ' Hz'] , ['#04:  ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:  ' num2str(round(fc(8))) ' Hz'], ['#12:  ' num2str(round(fc(12))) ' Hz'], ... 
        ['#16: ' num2str(round(fc(16))) ' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz']})
        box on
        hold off
        
        % Upper figure;
        subplot(2,1,1)
        hold on
        dy = 1;
        for ii = 1:size(outsig,1)
           env = abs(hilbert(outsig(ii,:),N));
           plot(treal,(14*env)+dy)
           dy = dy+1;
        end;
        title ([num2str(nchannels),' channel array of gamma impulse responses envelopes (classic)'])
        xlabel 'Time (ms)'
        ylabel('# Frequency Channel (ERB): Frequency (Hz)');
        xt = 0:5:25;
        yt = [1 4 8 12 16 20 24];
        axis([0, 25, 0, 26]);
        set(gca,'XTick',xt, 'YTick', yt)
        set(gca,'YTickLabel', {['#01:  '  num2str(fc(1),'%6.2f') ' Hz'] , ['#04:  ' num2str(fc(4),'%6.2f') ' Hz'], ...
        ['#08:  ' num2str(fc(8),'%6.2f') ' Hz'], ['#12:  ' num2str(fc(12),'%6.2f') ' Hz'], ... 
        ['#16: ' num2str(fc(16),'%6.2f') ' Hz'], ['#20: ' num2str(fc(20),'%6.2f') ' Hz'], ...
        ['#24: ' num2str(fc(24),'%6.2f') ' Hz']})
        box on
        hold off
    end;
    
%% Pattersons Paper - Roex(pwt) Filter
% Figure 6 

   if flags.do_fig6_patterson1987
       
      % Global parameters:
      fs = 25000;               % Sampling frequency in Hz;
      df = 1/1000;              % Frequency resolution roex(pwt);
      g = 0:df:1-df;            % Frequency axis roex(pwt);
      N = 4178;                 % Signallength for gammatone filter;
      f = (0:N-1)*fs/N;         % Frequency axis gammatone;

      
      % Roex(p) at 1000 Hz;
      % Parameters:
      fc = 1000;                % Center frequency in Hz; 
      erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
      
      % Equation from paper Moore, 1983 found in text after Equation 6.
      p = 4*fc/erb;             % Filter shape;
      w = 0.002;                % Realtive weight for second exponatial;
      t = 0.4*p;                % Filter shape tail;     

      % Equation A8 from Patterson 1982 (slopes are mirrored).
      weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for lower slope;
      weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for upper slope;
      rpwt_y = [fliplr(weight_l) 1 weight_u];                           % Filteramplitude;
      rpwt_x = ([-fliplr(g) 0 g]+1);                                    % Frequency axis;
          
      % Gammatone filter at 1000 Hz center frequency;
      % Parameters:
      fc = 1000;                % Center frequency in Hz;
      ind = find(f/fc > 2,1);   % Index where frequency axis is 2;
      insig = zeros(1,N);       % Input impulse signal;
      insig(1) = 1;             % with peak at 1;
      n = 4;                    % Filter order;
      betamul = 1.02;           % Filter bandwidth;
                   
      % Derive gammatone filter coefficients; 
      [b,a] = gammatone(fc,fs,n,betamul,'classic');
      % Filter impulse signal;
      outsig = 2*real(ufilterbankz(b,a,insig));
     
      % FFT;
      outsigfft = fft(outsig);
      % Normalize to 0;
      outsigfft = (outsigfft)/max(outsigfft);

      % Plot;
      figure('units','normalized','outerposition',[0 0.05 1 0.95])
      subplot(1,3,2)
      plot(rpwt_x, 10*log10(rpwt_y), 'b')
      hold on
      plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
      title('Roex(pwt) filters in blue vs. Gammatone classic in red')
      xlabel(['Center frequency at ', num2str(fc), ' Hz'])
      ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
      xt = 0:0.2:2;
      yt = -60:10:0;
      axis([0 2 -60 0])
      set(gca,'XTick',xt)
      set(gca,'YTick',yt)
      box on
      hold off

      
      % Roex(p) at 427 Hz
      % Parameters:
      fc = 427;                 % Center frequency in Hz; 
      erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
      
      % Equation from paper Moore, 1983 found in text after Equation 6.
      p = 4*fc/erb;             % Filter shape;
      w = 0.002;                % Realtive weight for second exponatial;
      t = 0.4*p;                % Filter shape tail;

      % Equation A8 from Patterson 1982 (slopes are mirrored).
      weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for lower slope;
      weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for upper slope;
      rpwt_y = [fliplr(weight_l) 1 weight_u];                           % Filteramplitude;
      rpwt_x = ([-fliplr(g) 0 g]+1);                                    % Frequency axis;
      
      % Gammatone filter at 427 Hz center frequency;
      % Parameters:
      fc = 427;                 % Center frequency in Hz;
      ind = find(f/fc > 2,1);   % Index where Frequency axis is 2;
      insig = zeros(1,N);       % Input impulse signal;
      insig(1) = 1;             % with peak at 1;
      n = 4;                    % Filter order;
      betamul = 1.02;           % Filter bandwidth;
                   
      % Derive gammatone filter coefficients; 
      [b,a] = gammatone(fc,fs,n,betamul,'classic');
      % Filter impulse signal;
      outsig = 2*real(ufilterbankz(b,a,insig));
     
      % FFT;
      outsigfft = fft(outsig);
      % Normalize to 0;
      outsigfft = (outsigfft)/max(outsigfft);

      % Plot;
      subplot(1,3,1)
      plot(rpwt_x, 10*log10(rpwt_y), 'b')
      hold on
      plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
      xlabel(['Center frequency at ', num2str(fc), ' Hz'])
      ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
      xt = 0:0.2:2;
      yt = -60:10:0;
      axis([0 2 -60 0])
      set(gca,'XTick',xt)
      set(gca,'YTick',yt)
      box on
      hold off


      % Roex(p) at 2089 Hz;
      % Parameters:
      fc = 2089;                % Center frequency in Hz; 
      erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
      
      % Equation from paper Moore, 1983 found in text after Equation 6.
      p = 4*fc/erb;             % Filter shape;
      w = 0.002;                % Realtive weight for second exponatial;
      t = 0.4*p;                % Filter shape tail;

      % Equation A8 from Patterson 1982 (slopes are mirrored).
      weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for lower slope;
      weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for upper slope;
      rpwt_y = [fliplr(weight_l) 1 weight_u];                           % Filteramplitude;
      rpwt_x = ([-fliplr(g) 0 g]+1);                                    % Frequency axis;
      
      % Gammatone filter at 2089 Hz center frequency;
      % Parameters:
      fc = 2089;                % Center frequency in Hz;
      ind = find(f/fc > 2,1);   % Index where Frequency axis is 2;
      insig = zeros(1,N);       % Input impulse signal;
      insig(1) = 1;             % with peak at 1;
      n = 4;                    % Filter order;
      betamul = 1.02;           % Filter bandwidth;
                   
      % Derive gammatone filter coefficients; 
      [b,a] = gammatone(fc,fs,n,betamul,'classic');
      % Filter impulse signal;
      outsig = 2*real(ufilterbankz(b,a,insig));
     
      % FFT;
      outsigfft = fft(outsig);
      % Normalize to 0;
      outsigfft = (outsigfft)/max(outsigfft);

      % Plot;
      subplot(1,3,3)
      plot(rpwt_x, 10*log10(rpwt_y), 'b')
      hold on
      plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
      xlabel(['Center frequency at ', num2str(fc), ' Hz'])
      ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
      xt = 0:0.2:2;
      yt = -60:10:0;
      axis([0 2 -60 0])
      set(gca,'XTick',xt)
      set(gca,'YTick',yt)
      box on
      hold off
   end; 
   
%% Pattersons Paper - Roex(pwt) Filter
% Figure 7 

   if flags.do_fig7_patterson1987
       
      % Global parameters:
      fs = 25000;               % Sampling frequency in Hz;
      df = 1/1000;              % Frequency resolution roex(pwt);
      g = 0:df:1-df;            % Frequency axis roex(pwt);
      N = 4178;                 % Signallength for gammatone filter;
      f = (0:N-1)*fs/N;         % Frequency axis gammatone;

   
      % Roex(p) at 1000 Hz;
      % Parameters:
      fc = 1000;                % Center frequency in Hz; 
      erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
      
      % Equation from paper Moore, 1983 found in text after Equation 6.
      p = 4*fc/erb;             % Filter shape;
      w = 0.002;                % Realtive weight for second exponatial;
      t = 0.4*p;                % Filter shape tail;     

      % Equation A8 from Patterson 1982 (slopes are mirrored).
      weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for lower slope;
      weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for upper slope;
      rpwt_y = [fliplr(weight_l) 1 weight_u];                           % Filteramplitude;
      rpwt_x = ([-fliplr(g) 0 g]+1);                                    % Frequency axis;
          
      % Gammatone filter at 1000 Hz center frequency;
      % Parameters:
      fc = 1000;                % Center frequency in Hz;
      ind = find(f/fc > 2,1);   % Index where frequency axis is 2;
      insig = zeros(1,N);       % Input impulse signal;
      insig(1) = 1;             % with peak at 1;
      n = 4;                    % Filter order;
      betamul = 1.02*1.1;       % Filter bandwidth;      
      
      % Derive gammatone filter coefficients; 
      [b,a] = gammatone(fc,fs,n,betamul,'classic');
      % Filter impulse signal;
      outsig = 2*real(ufilterbankz(b,a,insig));
     
      % FFT;
      outsigfft = fft(outsig);
      % Normalize to 0;
      outsigfft = (outsigfft)/max(outsigfft);

      % Plot;
      figure('units','normalized','outerposition',[0 0.05 1 0.95])
      subplot(1,3,2)
      plot(rpwt_x, 10*log10(rpwt_y), 'b')
      hold on
      plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
      title('Roex(pwt) filters in blue vs. Gammatone classic in red')
      xlabel(['Center frequency at ', num2str(fc), ' Hz'])
      ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
      xt = 0:0.2:2;
      yt = -60:10:0;
      axis([0 2 -60 0])
      set(gca,'XTick',xt)
      set(gca,'YTick',yt)
      box on
      hold off

      
      % Roex(p) at 427 Hz;
      % Parameters:
      fc = 427;                 % Center frequency in Hz; 
      erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
      
      % Equation from paper Moore, 1983 found in text after Equation 6.
      p = 4*fc/erb;             % Filter shape;
      w = 0.002;                % Realtive weight for second exponatial;
      t = 0.4*p;                % Filter shape tail;

      % Equation A8 from Patterson 1982 (slopes are mirrored).
      weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for lower slope;
      weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for upper slope;
      rpwt_y = [fliplr(weight_l) 1 weight_u];                           % Filteramplitude;
      rpwt_x = ([-fliplr(g) 0 g]+1);                                    % Frequency axis;
      
      % Gammatone filter at 427 Hz center frequency;
      % Parameters:
      fc = 427;                 % Center frequency in Hz;
      ind = find(f/fc > 2,1);   % Index where frequency axis is 2;
      insig = zeros(1,N);       % Input impulse signal;
      insig(1) = 1;             % with peak at 1;
      n = 4;                    % Filter order;
      betamul = 1.02*1.1;       % Filter bandwidth;
                   
      
      % Derive gammatone filter coefficients; 
      [b,a] = gammatone(fc,fs,n,betamul,'classic');
      outsig = 2*real(ufilterbankz(b,a,insig));
     
      % FFT;
      outsigfft = fft(outsig);
      % Normalize to 0;
      outsigfft = (outsigfft)/max(outsigfft);

      % Plot;
      subplot(1,3,1)
      plot(rpwt_x, 10*log10(rpwt_y), 'b')
      hold on
      plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
      xlabel(['Center frequency at ', num2str(fc), ' Hz'])
      ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
      xt = 0:0.2:2;
      yt = -60:10:0;
      axis([0 2 -60 0])
      set(gca,'XTick',xt)
      set(gca,'YTick',yt)
      box on
      hold off


      % Roex(p) at 2089 Hz;
      % Parameters:
      fc = 2089;                % Center frequency in Hz; 
      erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
      
      % Equation from paper Moore, 1983 found in text after Equation 6.
      p = 4*fc/erb;             % Filter shape;
      w = 0.002;                % Realtive weight for second exponatial;
      t = 0.4*p;                % Filter shape tail;

      % Equation A8 from Patterson 1982 (slopes are mirrored).
      weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for lower slope;
      weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);    % for upper slope;
      rpwt_y = [fliplr(weight_l) 1 weight_u];                           % Filteramplitude;
      rpwt_x = ([-fliplr(g) 0 g]+1);                                    % Frequency axis;
      
      % Gammatone filter at 2089 Hz center frequency;
      % Parameters:
      fc = 2089;                % Center frequency in Hz;
      ind = find(f/fc > 2,1);   % Index where frequency axis is 2;
      insig = zeros(1,N);       % Input impulse signal;
      insig(1) = 1;             % with peak at 1;
      n = 4;                    % Filter order;
      betamul = 1.02*1.1;       % Filter bandwidth;
                    
      % Derive gammatone filter coefficients; 
      [b,a] = gammatone(fc,fs,n,betamul,'classic');
      outsig = 2*real(ufilterbankz(b,a,insig));
     
      % FFT;
      outsigfft = fft(outsig);
      % Normalize to 0;
      outsigfft = (outsigfft)/max(outsigfft);

      % Plot;
      subplot(1,3,3)
      plot(rpwt_x, 10*log10(rpwt_y), 'b')
      hold on
      plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
      xlabel(['Center frequency at ', num2str(fc), ' Hz'])
      ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
      xt = 0:0.2:2;
      yt = -60:10:0;
      axis([0 2 -60 0])
      set(gca,'XTick',xt)
      set(gca,'YTick',yt)
      box on
      hold off
   end;
 
%% Pattersons Paper - Roex(pwt)
% Figure 8

    if flags.do_fig8_patterson1987
        
        % Global parameters:
        fs = 25000;               % Sampling frequency in Hz;
        df = 1/1000;              % Frequency resolution roex(pwt);
        g = 0:df:1-df;            % Frequency axis roex(pwt);
                  
        % Roex(pwt) filters at 427 Hz;
        % Parameters:
        fc = 427;                 % Center frequency at 427 Hz;
        erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
        
        % Equation from paper Moore, 1983 found in text after Equation 6.
        p = 4*fc/erb;             % Filter shape;
        w = 0.0025;               % Realtive weight for second exponatial;
        t = 0.2*p;                % Filter shape tail;                      
        
        % Equation A8 from Patterson 1982 (slopes are mirrored).
        weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for lower slope;
        weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for upper slope;
        rpwt_y = [fliplr(weight_l) 1 weight_u];                         % Filteramplitude;
        rpwt_x = ([-fliplr(g) 0 g]+1);                                  % Frequency-axis;
                
        % Gammatone filter at 427 Hz center frequency;
        % Parameters:
        fc = 427;                 % Center frequency at 427 Hz;
        N = 4178;                 % Signallength for gammatone filter;
        f = (0:N-1)*fs/N;         % Frequency axis;
        ind = find(f/fc > 2,1);   % Index where frequency axis is 2;
        insig = zeros(1,N);       % Input impulse signal;
        insig(1) = 1;             % with peak at 1;
        n = 2;                    % Filter order;
        betamul = 0.6;            % ERB multiplicator;
        
        % Derive gammatone filter coefficients; 
        [b,a] = gammatone(fc,fs,n,betamul,'classic');
        % Filter impulse signal;
        outsig = 2*real(ufilterbankz(b,a,insig));
        
        % FFT;
        outsigfft = fft(outsig);
        %Normalize; 
        outsigfft = (outsigfft)/max((outsigfft));
                         
        % Plot;
        figure('units','normalized','outerposition',[0 0.05 1 0.95])
        subplot(1,3,1)
        plot(rpwt_x,10*log10(rpwt_y), 'b')
        hold on
        plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
        xlabel(['Center frequency at ', num2str(fc), ' Hz'])
        ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
        xt = 0:0.2:2;
        yt = -50:10:0;
        axis([0 2 -50 0])
        set(gca,'XTick',xt)
        set(gca,'YTick',yt)
        box on
        hold off
        
          
        % Filters at 1000 Hz;
        % Parameters:
        fc = 1000;                % Center frequency at 427 Hz;
        erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
        
        % Equation from paper Moore, 1983 found in text after Equation 6.
        p = 4*fc/erb;             % Filter shape;
        w = 0.0025;               % Realtive weight for second exponatial;
        t = 0.2*p;                % Filter shape tail;                      
        
        % Equation A8 from Patterson 1982 (slopes are mirrored).
        weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for lower slope;
        weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for upper slope;
        rpwt_y = [fliplr(weight_l) 1 weight_u];                         % Filteramplitude;
        rpwt_x = ([-fliplr(g) 0 g]+1);                                  % Frequency-axis;
        
        % Gammatone filter at 1000 Hz center frequency;
        % Parameters:
        fc = 1000;                % Center frequency at 427 Hz;
        N = 4178;                 % Signallength for gammatone filter;
        f = (0:N-1)*fs/N;         % Frequency axis;
        ind = find(f/fc > 2,1);   % Index where frequency axis is 2;
        insig = zeros(1,N);       % Input impulse signal;
        insig(1) = 1;             % with peak at 1;
        n = 2;                    % Filter order;
        betamul = 0.6;            % ERB multiplicator;
        
        % Derives gammatone filter coefficients;
        [b,a] = gammatone(fc,fs,n,betamul,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
         
        % FFT;
        outsigfft = fft(outsig);
        outsigfft = (outsigfft)/max((outsigfft));

        % Plot;
        subplot(1,3,2)
        plot(rpwt_x,10*log10(rpwt_y), 'b')
        hold on
        plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
        title('Roex(pwt) filters in blue vs. Gammatone classic in red')
        xlabel(['Center frequency at ', num2str(fc), ' Hz'])
        ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
        xt = 0:0.2:2;
        yt = -50:10:0;
        axis([0 2 -50 0])
        set(gca,'XTick',xt)
        set(gca,'YTick',yt)
        box on
        hold off
          
          
        % Filters at 2089 Hz;
        % Parameters:
        fc = 2089;                % Center frequency at 427 Hz;
        erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
        
        % Equation from paper Moore, 1983 found in text after Equation 6.
        p = 4*fc/erb;             % Filter shape;
        w = 0.0025;               % Realtive weight for second exponatial;
        t = 0.2*p;                % Filter shape tail;                      
        
        % Equation A8 from Patterson 1982 (slopes are mirrored).
        weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for lower slope;
        weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for upper slope;
        rpwt_y = [fliplr(weight_l) 1 weight_u];                         % Filteramplitude;
        rpwt_x = ([-fliplr(g) 0 g]+1);                                  % Frequency-axis;
        
        % Gammatone filter at 2089 Hz center frequency;
        % Parameters:
        fc = 2089;                % Center frequency at 427 Hz;
        N = 4178;                 % Signallength for gammatone filter;
        f = (0:N-1)*fs/N;         % Frequency axis;
        ind = find(f/fc > 2,1);   % Index where frequency axis is 2;
        insig = zeros(1,N);       % Input impulse signal;
        insig(1) = 1;             % with peak at 1;
        n = 2;                    % Filter order;
        betamul = 0.6;            % ERB multiplicator;
        
        % Derives gammatone filter coefficients;
        [b,a] = gammatone(fc,fs,n,betamul,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
        
        % FFT;
        outsigfft = fft(outsig);
        outsigfft = (outsigfft)/max((outsigfft));
        
        % Plot;
        subplot(1,3,3)
        plot(rpwt_x,10*log10(rpwt_y), 'b')
        hold on
        plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
        xlabel(['Center frequency at ', num2str(fc), ' Hz'])
        ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
        xt = 0:0.2:2;
        yt = -50:10:0;
        axis([0 2 -50 0])
        set(gca,'XTick',xt)
        set(gca,'YTick',yt)
        box on
        hold off
    end;
    
%% Pattersons Paper - Roex(pwt)
% Figure 9

    if flags.do_fig9_patterson1987
        
        % Global parameters:
        fs = 25000;               % Sampling frequency in Hz;
        df = 1/1000;              % Frequency resolution roex(pwt);
        g = 0:df:1-df;            % Frequency axis roex(pwt);
                  
        % Roex(pwt) filters at 427 Hz;
        % Parameters:
        fc = 427;                 % Center frequency at 427 Hz;
        erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
        
        % Equation from paper Moore, 1983 found in text after Equation 6.
        p = 4*fc/erb;             % Filter shape;
        w = 0.005;                % Realtive weight for second exponatial;
        t = 0.24*p;               % Filter shape tail;                      
        
        % Equation A8 from Patterson 1982 (slopes are mirrored).
        weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for lower slope;
        weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for upper slope;
        rpwt_y = [fliplr(weight_l) 1 weight_u];                         % Filteramplitude;
        rpwt_x = ([-fliplr(g) 0 g]+1);                                  % Frequency-axis;
                
        % Gammatone filter at 427 Hz center frequency;
        % Parameters:
        fc = 427;                 % Center frequency at 427 Hz;
        N = 4178;                 % Signallength for gammatone filter;
        f = (0:N-1)*fs/N;         % Frequency axis;
        ind = find(f/fc > 2,1);   % Index where frequency axis is 2;
        insig = zeros(1,N);       % Input impulse signal;
        insig(1) = 1;             % with peak at 1;
        n = 2;                    % Filter order;
        betamul = 0.6;            % ERB multiplicator;
        
        % Derive gammatone filter coefficients; 
        [b,a] = gammatone(fc,fs,n,betamul,'classic');
        % Filter impulse signal;
        outsig = 2*real(ufilterbankz(b,a,insig));
        
        % FFT;
        outsigfft = fft(outsig);
        %Normalize; 
        outsigfft = (outsigfft)/max((outsigfft));
                         
        % Plot;
        figure('units','normalized','outerposition',[0 0.05 1 0.95])
        subplot(1,3,1)
        plot(rpwt_x,10*log10(rpwt_y), 'b')
        hold on
        plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
        xlabel(['Center frequency at ', num2str(fc), ' Hz'])
        ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
        xt = 0:0.2:2;
        yt = -50:10:0;
        axis([0 2 -50 0])
        set(gca,'XTick',xt)
        set(gca,'YTick',yt)
        box on
        hold off
        
          
        % Filters at 1000 Hz;
        % Parameters:
        fc = 1000;                % Center frequency at 427 Hz;
        erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
        
        % Equation from paper Moore, 1983 found in text after Equation 6.
        p = 4*fc/erb;             % Filter shape;
        w = 0.005;                % Realtive weight for second exponatial;
        t = 0.24*p;               % Filter shape tail;                      
        
        % Equation A8 from Patterson 1982 (slopes are mirrored).
        weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for lower slope;
        weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for upper slope;
        rpwt_y = [fliplr(weight_l) 1 weight_u];                         % Filteramplitude;
        rpwt_x = ([-fliplr(g) 0 g]+1);                                  % Frequency-axis;
        
        % Gammatone filter at 1000 Hz center frequency;
        % Parameters:
        fc = 1000;                % Center frequency at 427 Hz;
        N = 4178;                 % Signallength for gammatone filter;
        f = (0:N-1)*fs/N;         % Frequency axis;
        ind = find(f/fc > 2,1);   % Index where frequency axis is 2;
        insig = zeros(1,N);       % Input impulse signal;
        insig(1) = 1;             % with peak at 1;
        n = 2;                    % Filter order;
        betamul = 0.6;            % ERB multiplicator;
        
        % Derives gammatone filter coefficients;
        [b,a] = gammatone(fc,fs,n,betamul,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
         
        % FFT;
        outsigfft = fft(outsig);
        outsigfft = (outsigfft)/max((outsigfft));

        % Plot;
        subplot(1,3,2)
        plot(rpwt_x,10*log10(rpwt_y), 'b')
        hold on
        plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
        title('Roex(pwt) filters in blue vs. Gammatone classic in red')
        xlabel(['Center frequency at ', num2str(fc), ' Hz'])
        ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
        xt = 0:0.2:2;
        yt = -50:10:0;
        axis([0 2 -50 0])
        set(gca,'XTick',xt)
        set(gca,'YTick',yt)
        box on
        hold off
          
          
        % Filters at 2089 Hz;
        % Parameters:
        fc = 2089;                % Center frequency at 427 Hz;
        erb = audfiltbw(fc);      % Equivalent rectangular bandwidth;
        
        % Equation from paper Moore, 1983 found in text after Equation 6.
        p = 4*fc/erb;             % Filter shape;
        w = 0.005;                % Realtive weight for second exponatial;
        t = 0.24*p;               % Filter shape tail;                      
        
        % Equation A8 from Patterson 1982 (slopes are mirrored).
        weight_l = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for lower slope;
        weight_u = (1-w).*(1+p.*g).*exp(-p.*g)+w*(1+t.*g).*exp(-t.*g);  % for upper slope;
        rpwt_y = [fliplr(weight_l) 1 weight_u];                         % Filteramplitude;
        rpwt_x = ([-fliplr(g) 0 g]+1);                                  % Frequency-axis;
        
        % Gammatone filter at 2089 Hz center frequency;
        % Parameters:
        fc = 2089;                % Center frequency at 427 Hz;
        N = 4178;                 % Signallength for gammatone filter;
        f = (0:N-1)*fs/N;         % Frequency axis;
        ind = find(f/fc > 2,1);   % Index where frequency axis is 2;
        insig = zeros(1,N);       % Input impulse signal;
        insig(1) = 1;             % with peak at 1;
        n = 2;                    % Filter order;
        betamul = 0.6;            % ERB multiplicator;
        
        % Derives gammatone filter coefficients;
        [b,a] = gammatone(fc,fs,n,betamul,'classic');
        outsig = 2*real(ufilterbankz(b,a,insig));
        
        % FFT;
        outsigfft = fft(outsig);
        outsigfft = (outsigfft)/max((outsigfft));
        
        % Plot;
        subplot(1,3,3)
        plot(rpwt_x,10*log10(rpwt_y), 'b')
        hold on
        plot(f(1:ind)/fc,20*log10(abs(outsigfft(1:ind))),'r')
        xlabel(['Center frequency at ', num2str(fc), ' Hz'])
        ylabel('Attenuation (dB) roex: 10*log10 vs. gammatone: 20*log10')
        xt = 0:0.2:2;
        yt = -50:10:0;
        axis([0 2 -50 0])
        set(gca,'XTick',xt)
        set(gca,'YTick',yt)
        box on
        hold off
    end;
    
%% Pattersons Paper 
% Figure 10a
% gammatone 'classic' (causalphase, real)

    if flags.do_fig10a_patterson1987
        % Input parameters
        fs = 25000;                 % Sampling frequency in Hz;
        ts = 1/fs;                  % Time between sampling points in s;
        N=4096;                     % Length of signal;
        treal = (1:N)*ts*1000;      % Time axis in ms;
        nchannels = 37;             % Number of channels in filterbank;
        flow = 100;                 % ERB lowest center frequency in Hz;
        fhigh = 5000;               % ERB highest center frequency in Hz;
        fc = erbspace(flow,fhigh,nchannels);   % 37 erb-spaced channels 
        insig = zeros(1,N);         % Impulse as input signal 
        insig(1) = 1;
        
        % Derives filter coefficients for 37 erb spaced channels
        [b,a] = gammatone(fc,fs,'classic');
        % Filters impulse signal with filter coefficients from above
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Plot
        figure
        hold on
        dy = 1;
        for ii = 1:nchannels
           plot(treal,14*outsig(ii,:) + dy)
           dy = dy+1;
        end;
        clear dy;
        title 'Gammatone impulse response filterbank without phase-compensation'
        xlabel 'Time (ms)'
        ylabel('# Frequency Channel (ERB): Frequency (Hz)');
        xt = 0:5:25;
        yt = [1 4 8 12 16 20 24 28 32 37];
        axis([0, 25, 0, 40]);
        set(gca,'XTick',xt, 'YTick', yt)
        set(gca,'YTickLabel', {['#01:   '  num2str(round(fc(1))) ' Hz'] , ['#04:   ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:   ' num2str(round(fc(8))) ' Hz'], ['#12:   ' num2str(round(fc(12))) ' Hz'], ... 
        ['#16:   ' num2str(round(fc(16))) ' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz'], ['#28: ' num2str(round(fc(28))) ' Hz'], ...
        ['#32: ' num2str(round(fc(32))) ' Hz'], ['#37: ' num2str(round(fc(37))) ' Hz'], })
        box on
        hold off
    end;

%% Pattersons Paper
% Figure 10b
% gammatone 'classic','peakphase','complex'

    if flags.do_fig10b_patterson1987
        % Input parameters
        fs = 25000;                 % Sampling frequency;
        ts = 1/fs;                  % Time between sampling points in s
        N=4096;                     % Length of signal
        treal = (1:N)*ts*1000;      % Time axis in ms
        nchannels = 37;             % Number of channels in filterbank;
        flow = 100;                 % ERB lowest center frequency in Hz;
        fhigh = 5000;               % ERB highest center frequency
        fc = erbspace(flow,fhigh,nchannels);   % 37 erb-spaced channels 
        insig = zeros(1,N);         % Impulse as input signal 
        insig(1) = 1;

        % Derives filter coefficients for 37 erb-spaced channels
        [b,a] = gammatone(fc,fs,'classic');
        % Filters impulse signal with filter coefficients from above
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Find peak at envelope maximum
        envmax = zeros(1,nchannels);
        for ii = 1:nchannels
            % Envelope maximum per channel
            envmax(ii) = find(abs(outsig(ii,:)) == max(abs(outsig(ii,:))));
        end;
        
        % Adding Delay
        desiredpeak  = max(envmax)+1;
        for ii = 1:nchannels
            % Time to delay
            delay =  zeros(1,desiredpeak - envmax(ii));
            % Add delay
            outsig(ii,:) = [delay outsig(ii,1:size(outsig,2) - length(delay))];
        end;
        
        %Plot
        figure
        hold on
        dy = 1;
        for ii = 1:nchannels
           plot(treal,14*real(outsig(ii,:))+dy)
           dy = dy+1;
        end;
        clear dy;
        title 'Gammatone filterbank with envelope phase-compensation'
        xlabel 'Time (ms)'
        ylabel('# Frequency Channel (ERB): Frequency (Hz)');
        xt = 0:5:25;
        yt = [1 4 8 12 16 20 24 28 32 37];
        axis([0, 25, 0, 39]);
        set(gca,'XTick',xt, 'YTick', yt)
        set(gca,'YTickLabel', {['#01:   '  num2str(round(fc(1))) ' Hz'] , ['#04:   ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:   ' num2str(round(fc(8))) ' Hz'], ['#12:   ' num2str(round(fc(12))) ' Hz'], ... 
        ['#16:   ' num2str(round(fc(16))) ' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz'], ['#28: ' num2str(round(fc(28))) ' Hz'], ...
        ['#32: ' num2str(round(fc(32))) ' Hz'], ['#37: ' num2str(round(fc(37))) ' Hz'], })
        box on
        hold off
    end;
    
%% Pattersons Paper
% Figure 10c
% gammatone 'classic','peakphase','complex'

    if flags.do_fig10c_patterson1987
        % Input parameters
        fs = 25000;                 % Sampling frequency;
        ts = 1/fs;                  % Time between sampling points in s
        N=4096;                     % Length of signal
        treal = (1:N)*ts*1000;      % Time axis in ms
        nchannels = 37;             % Number of channels in filterbank;
        flow = 100;                 % ERB lowest center frequency in Hz;
        fhigh = 5000;               % ERB highest center frequency
        fc = erbspace(flow,fhigh,nchannels);   % 37 erb-spaced channels 
        insig = zeros(1,N);         % Impulse as input signal 
        insig(1) = 1;

        % Derives filter coefficients for 37 erb-spaced channels
        [b,a] = gammatone(fc,fs,'classic', 'exppeakphase');
        % Filters impulse signal with filter coefficients from above
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Find peak at envelope maximum
        envmax = zeros(1,nchannels);
        for ii = 1:nchannels
            % Envelope maximum per channel
            envmax(ii) = find(abs(outsig(ii,:)) == max(abs(outsig(ii,:))));
        end;
        
        % Adding Delay
        desiredpeak  = max(envmax)+1;
        for ii = 1:nchannels
            % Time to delay
            delay =  zeros(1,desiredpeak - envmax(ii));
            % Add delay
            outsig(ii,:) = [delay outsig(ii,1:size(outsig,2) - length(delay))];
        end;
        
        %Plot
        figure
        hold on
        dy = 1;
        for ii = 1:nchannels
           plot(treal,14*real(outsig(ii,:))+dy)
           dy = dy+1;
        end;
        clear dy;
        title 'Gammatone filterbank with envelope phase-compensation and delayed channel'
        xlabel 'Time (ms)'
        ylabel('# Frequency Channel (ERB): Frequency (Hz)');
        xt = 0:5:25;
        yt = [1 4 8 12 16 20 24 28 32 37];
        axis([0, 25, 0, 39]);
        set(gca,'XTick',xt, 'YTick', yt)
        set(gca,'YTickLabel', {['#01:   '  num2str(round(fc(1))) ' Hz'] , ['#04:   ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:   ' num2str(round(fc(8))) ' Hz'], ['#12:   ' num2str(round(fc(12))) ' Hz'], ... 
        ['#16:   ' num2str(round(fc(16))) ' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz'], ['#28: ' num2str(round(fc(28))) ' Hz'], ...
        ['#32: ' num2str(round(fc(32))) ' Hz'], ['#37: ' num2str(round(fc(37))) ' Hz'], })
        box on
        hold off
    end;
    
%% Pattersons Paper 
% Figure 11 
% gammatone 'classic' (causalphase, real) pulsetrain
% Warning: Not sure if b=3 or b=1/3 // b=2 or b=1/2  
  
    if flags.do_fig11_patterson1987
        
    % Input parameters
        fs = 10000;                 % Sampling frequency in Hz;
        ts = 1/fs;                  % Time between sampling points in s;
        N=4096;                     % Length of signal;
        treal = (1:N)*ts*1000;      % Time axis in ms;
        nchannels = 37;             % Number of channels in filterbank;
        flow = 100;                 % ERB lowest center frequency in Hz;
        fhigh = 5000;               % ERB highest center frequency in Hz;
        fc = erbspace(flow,fhigh,nchannels);   % 37 erb-spaced channels 
        insig = zeros(1,N);         % Impulse as input signal 
        insig(160) = 1;             % at 16 ms,
        insig(240) = 1;             % at 24 ms,
        insig(320) = 1;             % at 32 ms;

        % Derives filter coefficients for erb spaced channels
        [b,a] = gammatone(fc,fs,'classic');
        % Filters impulse signal with filter coefficients from above
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
              
        % Plot
        figure('units','normalized','outerposition',[0 0.05 1 0.95])
        hta1 = subplot(3,1,1);
        posAxes1 = get(hta1, 'Position');
        posAxes1(2) = posAxes1(2) + (1 - posAxes1(2))*5/8;
        posAxes1(4) = posAxes1(4)/4;
        set(hta1, 'Position',posAxes1); 
        plot(treal,insig);
        title 'Pulsetrain'
        ylabel 'Amplitude'
        set(gca, 'XLim',[14.5 39],'YLim',[0 1])
        dy = 1;
        hta2 = subplot(3,1,2);
        posAxes2 = get(hta2, 'Position');
        posAxes2(2) = posAxes1(2)- posAxes2(2);
        posAxes2(4) = posAxes1(4)*6;
        set(hta2, 'Position',posAxes2);
        hold on
        for ii = 1:size(outsig,1)
           plot(treal,4*outsig(ii,:)+dy)
           dy = dy+1;
        end;
        title 'Gammatone impulse response filterbank but without envelope phase-compensation'
        ylabel('# Channel (ERB) : Frequency (Hz)');
        xt = 0:5:25;
        yt = [1 4 8 12 16 20 24 28 32 37];
        axis([14.5, 39, 0, 38]);
        set(gca,'XTick',xt, 'YTick', yt)
        set(gca,'YTickLabel', {['#01:   '  num2str(round(fc(1))) ' Hz'] , ['#04:   ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:   ' num2str(round(fc(8))) ' Hz'], ['#12:   ' num2str(round(fc(12))) ' Hz'], ... 
        ['#16:   ' num2str(round(fc(16))) ' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz'], ['#28: ' num2str(round(fc(28))) ' Hz'], ...
        ['#32: ' num2str(round(fc(32))) ' Hz'], ['#37: ' num2str(round(fc(37))) ' Hz'], })
        box on
        hta3 = subplot(3,1,3);
        posAxes3 = get(hta3, 'Position');
        posAxes3(2) =  posAxes3(2) - (1 - posAxes1(2))*1/8; % *3/8
        posAxes3(4) = posAxes1(4)*6;
        set(hta3, 'Position',posAxes3);
 %       title 'test'
        xlabel 'Time (ms)'
         ylabel('# Channel (ERB) : Frequency (Hz)');
        xt = 15:5:40;
        yt = [1 4 8 12 16 20 24 28 32 37];
        axis([14.5, 39, 0, 38]);
        set(gca,'XTick',xt, 'YTick', yt)
        set(gca,'YTickLabel', {['#01:   '  num2str(round(fc(1))) ' Hz'] , ['#04:   ' num2str(round(fc(4))) ' Hz'], ...
        ['#08:   ' num2str(round(fc(8))) ' Hz'], ['#12:   ' num2str(round(fc(12))) ' Hz'], ... 
        ['#16:   ' num2str(round(fc(16))) ' Hz'], ['#20: ' num2str(round(fc(20))) ' Hz'], ...
        ['#24: ' num2str(round(fc(24))) ' Hz'], ['#28: ' num2str(round(fc(28))) ' Hz'], ...
        ['#32: ' num2str(round(fc(32))) ' Hz'], ['#37: ' num2str(round(fc(37))) ' Hz'], })
        hold off
    end
    
%% Lyons Paper - Allpole Gammatone Filter
% Figure 1
% gammatone (allpole,causalphase,real')
% gammatone 'classic' (causalphase, real)
% Warning: Not sure if b should be: b=3 or b=1/3 // b=2 or b=1/2;  
    if flags.do_fig1_lyon1997
                
        % Input parameters
        fs = 24000;                     % Sampling frequency in Hz;
        flow = 100;                     % ERB lowest center frequency in Hz;
        fhigh = 4000;                   % ERB highest center frequency in Hz;
        fc = erbspacebw(flow,fhigh);    % 24 ERB spaced channels;
        channel = 20;                   % Channel
        n = 6;                          % Filter order;
        betamul = 3;                    % Betamul
        
        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,n,betamul); % default = APGF
        % Returns frequency response vector 'h' and the corresponding
        % angular frequency vector 'w' for channel 13.
        [h,w] = freqz(b(channel,:), a(channel,:));
        % Scale first value to zero.
        h = h(:)/h(1);
        % For scaling of following transfer function (peak on peak).
        nor = real(max(h));
             
        % Plot
        figure
        semilogx(w/(fc(channel)/fs*2*pi), 10*log10(abs(h)),'r-');
        title 'Transfer functions of classic and allpole gammatone filters for b = \omega_r/3 and b = \omega_r/2'
        xlabel '\omega/\omega_r'
        ylabel 'Magnitude gain (db)'
        set(gca,'XLim',[0.03 3],'YLim',[-40, 35]) 
        box on
        hold on
        
        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,n,betamul,'classic'); % GTF
        % Returns frequency response vector 'h' and the corresponding
        % angular frequency vector 'w' for channel 13.
        [h,w] = freqz(b(20,:), a(20,:));
        % Scale tranfer function to previous transfer fuction.
        h=abs(h) * nor;
        
        % Plot 
        semilogx(w/(fc(channel)/fs*2*pi), 10*log10(abs(h)),'b--');
        
        
        % Change betamul;
        betamul = 2;
        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,n,betamul);
        % Returns frequency response vector 'h' and the corresponding
        % angular frequency vector 'w' for channel 13.
        [h,w] = freqz(b(channel,:), a(channel,:));
        % Scale first value to zero;
        h = h(:)/h(1);
        % For scaling of following transfer function (peak on peak).
        nor = real(max(h));
       
        % Plot
        semilogx(w/(fc(channel)/fs*(2*pi)), 10*log10(abs(h)),'r-');

        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,n,betamul,'classic');
        % Returns frequency response vector 'h' and the corresponding
        % angular frequency vector 'w' for channel 13.
        [h,w] = freqz(b(channel,:), a(channel,:));
        % Scale tranfer function to previous transfer fuction.
        h=abs(h) * nor;
        
        % Plot
        semilogx(w/(fc(channel)/fs*2*pi), 10*log10(abs(h)),'b--');
        legend('APGF, b=3', 'GTF, b=3', 'APGF, b=2', 'GTF, b=2', 'Location', 'NorthWest')  
        box on
        hold off
        warning('Not sure if b should be: b=3 or b=1/3 // b=2 or b=1/2');  
    end;

%% Lyons Paper 
% Figure 2
% gammatone 'classic' (causalphase, real)
% gammatone (allpole,causalphase) 'complex'
% Warning: Not sure if b should be: b=3 or b=1/3 // b=2 or b=1/2;  

    if flags.do_fig2_lyon1997
        % Input parameters
        fs = 24000;             % Sampling frequency in Hz;
        flow = 100;             % ERB lowest center frequency in Hz;
        fhigh = 4000;           % ERB highest center frequency in Hz;
        fc = erbspacebw(flow,fhigh);   % 24 ERB spaced channels;
        N = 512;                % Signal length;
        insig = zeros(1,N);     % Input signal;
        insig(1) = 1;           % Impulse at sample 1024;
        channel = 12;           % Channel
        n = 6;                  % Filter order 6;
        betamul = 1/3;            % Betamul;
        treal = (0:N-1)*((2*pi*fc(channel))/fs/pi);
        
        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,n,betamul,'classic');
        % Filters impulse signal with filter coefficients from above.
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
   
        % Plot
        figure
        subplot(2,1,1)
        plot(treal, outsig(channel,:),'b--')
        title 'Impulse responses for b = \omega_r / 3 (low damping, low BW, high gain)'
        xlabel 't * \omega_r / pi'
        ylabel 'Amplitude'
        set(gca,'XLim',[0 16], 'YLim',[-0.25 0.25])
        box on
        hold on
        %'XLim',[1020 1084],'
        % Derives filter coefficients for erb spaced channels.
        %[b,a] = gammatone(erb,fs,n,betamul);
        [b,a] = gammatone(fc,fs,n,betamul);
        % Filters impulse signal with filter coefficients from above.
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Plot
        plot(treal,outsig(channel,:),'r-')
        legend('GTF','APGF')
        hold off

        % Change filter order 
        n = 6;
        % Change betamul
        betamul = 1/2;
        
        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,n,betamul,'classic');
        % Filters impulse signal with filter coefficients from above.
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Plot
        subplot(2,1,2)
        plot(treal, outsig(channel,:),'b--')
        title 'Impulse responses for b = \omega_r / 2 (high damping, high BW, low gain)'
        xlabel 't * \omega_r / pi'
        ylabel 'Amplitude'
        set(gca, 'XLim',[0 16],'YLim',[-0.25 0.25])
        box on
        hold on

        % Derives filter coefficients for erb spaced channels.
        %[b,a] = gammatone(erb,fs,n,betamul);
        [b,a] = gammatone(fc,fs,n,betamul);
        % Filters impulse signal with filter coefficients from above.
        outsig = 2*real(ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Plot
        plot(treal,outsig(channel,:),'r-')
        legend('GTF', 'APGF')
        box on
        hold off
        warning('Not sure if b should be: b=3 or b=1/3 // b=2 or b=1/2');  
     end;
    
%% Hohmanns paper
% Figure 1
% gammatone (allpole,causalphase) 'complex'

    if flags.do_fig1_hohmann2002
        % Input parameters
        fs = 10000;             % Sampling frequency in Hz;
        fc = 1000;              % Center frequency at 1000 Hz;
        N = 4096;               % Signal length;
%         n = 4;                  % Filter order
%         betamul = 0.75;         % Bandwidth multiplicator
        insig = zeros(1,N);     % Input signal with
        insig(1) = 1;           % impulse at 1;

        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,'complex');
        % Filters impulse signal with filter coefficients from above.
        outsig = (ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Envelope of impulse response of gammatone filter at 1000 Hz
        outsigenv = hilbert(real(outsig),fs);

        % Plot
        figure  
        plot(abs(outsigenv),'b-.')
        xlabel('Sample')
        ylabel('Amplitude')
        set(gca, 'XLim',[0 200],'YLim',[-0.04 0.04])
        box on
        hold on
        plot(real(outsig),'r-')
        plot(imag(outsig),'g--')
        hold off
    end;
    
%% Hohmanns paper
% Figure 2
% gammatone (allpole,causalphase) 'complex'

    if flags.do_fig2_hohmann2002
        % Input parameters
        fs = 10000;             % Sampling frequency in Hz;
        fc = 1000;              % Center frequency at 1000 Hz;
        fn = fs/2;              % Nyquist frequency in Hz   
        Tn = 1/fn;              % Time between samples for Nyquist frequency
        N = 4096;               % Signal length;
        df = fs/N;              % Distance frequency of fft 
        xfn = (0:df:fn-df)*Tn;  % Frequency vector for fft
        insig = zeros(1,N);     % Input signal with
        insig(1) = 1;           % impulse at one;

        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,'complex');
        % Filters impulse signal with filter coefficients from above.
        outsig = (ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);
        
        % Upper two panels
        % FFT of real part of output signal
        outfftr = fft(2*real(outsig));
        % Returns unwrapped phase response vectors evaluated at n equally-spaced points
        % around the unit circle from 0 to 2? radians/sample. 
        [phi] = phasez(b,a,2048);
        phi = phi+2*pi;
        phi = phi.';
        
        % Lower two panels
        % The real-to-imaginary transfer function is calculated by
        % dividing the FFT-spectra of imaginary and real part of the
        % impulse response.
        outffti = fft(2*imag(outsig));
        outsigrtoi = outffti./outfftr;
        % Extract phase angle of signal.
        theta = angle(outsigrtoi);
        % Add Pi/2 to phase.
        phirtoi = theta+ pi/2;
        
        % Plot
        figure
        subplot(4,1,1)
        plot(xfn, 20*log10(abs(outfftr(1:N/2))))
        set(gca,'Xlim',[0 1], 'YLim',[-70 0])
        ylabel('Magnitude/db')
        box on
        subplot(4,1,2)
        plot(xfn,phi)
        ylabel('Phase/rad')
        set(gca,'Xlim',[0 1],'Ylim',[-5 5])
        box on
        subplot(4,1,3)
        plot(xfn, 20*log10(abs(outsigrtoi(1:N/2))))
        set(gca,'Xlim',[0 1],'Ylim',[-10 10])
        ylabel('Magnitude/db')
        box on
        subplot(4,1,4)
        plot(xfn, phirtoi(1:N/2))
        set(gca,'Xlim',[0 1],'Ylim',[-2 2])
        xlabel('Frequency / \pi')
        ylabel('Phase + \pi / 2 / rad')
        box on
        grid
    end;
    
%% Hohmanns paper
% Figure 3
% gammatone (allpole,causalphase) 'complex'

     if flags.do_fig3_hohmann2002
        fs = 16276;                     % Sampling frequency in Hz;
        flow = 70;                      % ERB lowest center frequency in Hz;
        fhigh = 6700;                   % ERB highest center frequency in Hz;
        fc = erbspacebw(flow,fhigh);    % 24 ERB spaced channels;
        fn = fs/2;                      % Nyquist frequency in Hz
        N = 2048;                       % Signal length

        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,'complex');
        
        % Frequency response per channel
        h = zeros(length(fc),N/4);
        w = zeros(length(fc),N/4);
        for ii=1:length(fc)
           [h(ii,:),w(ii,:)] = freqz(b(ii,:),a(ii,:));
        end;
        % Normalize peak to zero
        h(:) = abs(h(:));
        
        % Plot
        figure
        hold on
        for ii = 1:length(fc)
            plot(w(ii,:)/(2*pi)*fs, 20*log10(h(ii,:)))
        end;
        xlabel('Frequency / Hz')
        ylabel('Level / dB')
        set(gca,'Xlim',[0 fn],'Ylim',[-40 0])
        box on
        hold off
     end;
     
%% Hohmanns paper
% Figure 4
% gammatone (allpole,causalphase) 'complex'
% Result is delayed and peakphased.
% Case 2a: envelope maximum before desired peak;

    if flags.do_fig4_hohmann2002
        fs = 16276;                     % Sampling frequency in Hz;
        flow = 70;                      % ERB lowest center frequency in Hz;
        fhigh = 6700;                   % ERB highest center frequency in Hz;
        fc = erbspacebw(flow,fhigh);    % 24 ERB spaced channels;
        insig = zeros(1,fs);            % Input signal with
        insig(1) = 1;                   % impulse at one
       
        % The envelope maximum envmax is earlier in time than the desired
        % peak at channels 13 to 30.
        channel = 20;
        desiredpeak = 65;

        % Derives filter coefficients for erb spaced channels.
        [b,a] = gammatone(fc,fs,'complex');
        % Filters impulse signal with filter coefficients from above. 
        outsig = (ufilterbankz(b,a,insig));
        outsig = permute(outsig,[3 2 1]);

        % Envelope maximum
        outenv = hilbert(real(outsig(channel,:)),fs);
        envmax = find(abs(outenv(1,:)) == max(abs(outenv(1,:))));
        % Signal maximum
        sigmax = find(real(outsig(channel,:)) == max(real(outsig(channel,:))));

        % Phase delay between envelope maximum and signal maximum
        % Equation 18 Phi_k
        phi = fc(channel)*(-2*pi)*(envmax - sigmax)/fs;
        
        % Derives signal with envelope maximum at signal maximum 
        % Equation 19
        outsignew = real(outsig(channel,:)) * cos(phi) - imag(outsig(channel,:)) * sin(phi);
          
          % For testing purpose 
%         outenvnew = hilbert(real(outsignew(1,:)),fs);
%         envmaxnew = find(abs(outenvnew(1,:)) == max(abs(outenvnew(1,:))));
%         sigmaxnew = find(real(outsignew(1,:)) == max(real(outsignew(1,:))));
%         
%         warning(['FIXME: The maximum of the real part of the impulse response ' ... 
%                  'misses the maximum of the envelope in channel '  num2str(channel)...
%                  ' by ' num2str(envmaxnew - sigmaxnew) ' sample(s).']);
        
        % Derives time delay        
        % Equation 20
        delay = outsignew(end-desiredpeak+envmax+1:end);
        
        % Delay signal
        % Equation 21
        outsigdelay= [delay outsignew(1:length(outsignew) - delay)];
        outsigdelay = outsigdelay(1:length(outsignew) - delay);
        
        % Envelope of delayed signal
        outsigdelayenv = hilbert(real(outsigdelay),fs);
        
        % Plot
        figure
        subplot(2,1,1)
        plot(real(outsig(channel,:)),'b-')
        hold on
        set(gca, 'XLim',[0 200],'YLim',[-0.04 0.04])
        line([65 65], [-0.05 0.05])
        xlabel('Sample')
        ylabel('Amplitude')
        box on
        plot(abs(outenv),'r-')
        hold off
        subplot(2,1,2)
        plot (real(outsigdelay),'b-')
        hold on
        set(gca, 'XLim',[0 200],'YLim',[-0.04 0.04])
        xlabel('Sample')
        ylabel('Amplitude')
        line([65 65], [-0.05 0.05])
        plot(abs(outsigdelayenv),'r-')
        box on
        hold off
    end;
     
end
