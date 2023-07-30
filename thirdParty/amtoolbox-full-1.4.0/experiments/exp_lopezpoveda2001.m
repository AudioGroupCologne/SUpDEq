function output=exp_lopezpoveda2001(varargin)
%EXP_LOPEZPOVEDA2001   Figures from Lopez-Poveda and Meddis (2001)
%   Usage: output = exp_lopezpoveda2001(flag)
%
%   EXP_LOPEZPOVEDA2001(flags,... ) reproduces experiments from the Lopez
%   & Poveda (2001) paper.
%
%   The following flags can be specified;
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'no_plot'  Don't plot, only return data.
%
%     'fig2'    Reproduce Fig. 2 from Lopez & Poveda (2001)
%
%               Fig. 2a represents the outer ear filter - pressure gain (dB) over
%               frequency with data points from Pralong and Carlile (1996).
%
%               Fig. 2b represents the middle ear filter - stapes velocity at
%               0dB over frequency in one plot Fig. 2b shows for default 
%               fs = 22050 Hz: 
%               
%                Data points directly derived from Goode et al. 1994 
%
%                FIR filter with data points from Goode et al. 1994
%
%                Data points of figure 2b from Lopez-Poveda and Meddis
%                 2001 (read from fig 2b, actually also derived from Goode et
%                 al. 1994) The output are the data points of the respective
%                 figure.
%
%               The dimensions of the output are: frequency values x data
%               points  x figure no.
%
%     'fig2a'   Reproduce just Fig. 2a.
%
%     'fig2b'   Reproduce just Fig. 2b.
%
%     'fig3bc'  Reproduce Fig. 3b and c from Lopez & Poveda
%               (2001). Isointensity response of the linear, nonlinear and
%               summed response of the lopezpoveda2001 filter for an input level of
%               30dB (fig 3b) and 85dB (fig 3c) SPL The output is the output
%               of the lopezpoveda2001 filter for the different input levels. The
%               dimensions of the output are: input frequency x [frequency
%               values, linear output, nonlinear output, summed lopezpoveda2001 output]
%               x input level.
%
%     'fig3b'   Reproduce just Fig. 3b.
%
%     'fig3c'   Reproduce just fig. 3c.
%
%     'fig4'    Reproduce Fig. 4 from Lopez & Poveda (2001) - pulsation threshold 
%               data (just the model results, not the experimental data)
%               The output is the model result for the different parameter sets.
%               The dimensions of the output are: signal level x masker level x signal frequency
%               masker level consists of 4 columns:
%
%               1) Signal level in dB SPL (x axis in the plots)
%               
%               2) Results for parameter set of YO, table I
%
%               3) Results for average parameter set, table II
% 
%               4) Results for regression lines, table III
%
%   See also: lopezpoveda2001, data_lopezpoveda2001, data_pralong1996, data_goode1994
%
%   Examples:
%   ---------
%
%   To display Figure 2 use :
%
%     exp_lopezpoveda2001('fig2');
%
%   To display Figure 3b and 3c use :
%
%     exp_lopezpoveda2001('fig3bc');
%
%   To display Figure 4 use :
%
%     exp_lopezpoveda2001('fig4');
%
%   References:
%     R. Goode, M. Killion, K. Nakamura, and S. Nishihara. New knowledge
%     about the function of the human middle ear: development of an improved
%     analog model. The American journal of otology, 15(2):145--154, 1994.
%     
%     E. Lopez-Poveda and R. Meddis. A human nonlinear cochlear filterbank.
%     J. Acoust. Soc. Am., 110:3107--3118, 2001.
%     
%     D. Pralong and S. Carlile. The role of individualized headphone
%     calibration for the generation of high fidelity virtual auditory space.
%     J. Acoust. Soc. Am., 100:3785--3793, 1996.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_lopezpoveda2001.php


%  #Author: Katharina Egger

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% ------ Check input options --------------------------------------------

  definput.flags.type = {'missingflag','fig2','fig2a','fig2b','fig3bc','fig3b','fig3c','fig4'};
  definput.flags.plot = {'plot','no_plot'};
  definput.keyvals.predrnl = {};
  definput.keyvals.postdrnl = {};

  % Parse input options
  [flags,kv]  = ltfatarghelper({},definput,varargin);
        
if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% parameter set of YO, table I
% The data is specified in this way, because the data for figure 4 is
% not specified as the two parameter fit, but instead specified
% directly. The parameters below accomplish this by removing all
% frequency dependence in the lopezpoveda2001. The parameters can then be
% specified exactly for a single frequency, but only for one frequency
% at a time.

f250 = {...
  'flow',250,...
  'fhigh',250,...
  'lin_fc', [log10(235) 0],...
  'lin_bw', [log10(115) 0],...
  'lin_gain', [log10(1400) 0],...
  'lin_lp_cutoff', [log10(235) 0],...      
  'nlin_fc_before', [log10(250) 0],...
  'nlin_bw_before', [log10(84) 0],...
  'nlin_lp_cutoff', [log10(250) 0],...      
  'nlin_a', [log10(2124) 0],...
  'nlin_b', [log10(0.45) 0] };

f500 = {...
  'flow',500,...
  'fhigh',500,...
  'lin_fc', [log10(460) 0],...
  'lin_bw', [log10(150) 0],...
  'lin_gain', [log10(800) 0],...
  'lin_lp_cutoff', [log10(460) 0],...      
  'nlin_fc_before', [log10(500) 0],...
  'nlin_bw_before', [log10(103) 0],...
  'nlin_lp_cutoff', [log10(500) 0],...      
  'nlin_a', [log10(4609) 0],...
  'nlin_b', [log10(0.28) 0] };

f1000 = {...
  'flow',1000,...
  'fhigh',1000,...
  'lin_fc', [log10(945) 0],...
  'lin_bw', [log10(240) 0],...
  'lin_gain', [log10(520) 0],...
  'lin_lp_cutoff', [log10(945) 0],...      
  'nlin_fc_before', [log10(1000) 0],...
  'nlin_bw_before', [log10(175) 0],...
  'nlin_lp_cutoff', [log10(1000) 0],...      
  'nlin_a', [log10(4598) 0],...
  'nlin_b', [log10(0.130) 0]};

f2000 = {...
  'flow',2000,...
  'fhigh',2000,...
  'lin_fc', [log10(1895) 0],...
  'lin_bw', [log10(390) 0],...
  'lin_gain', [log10(400) 0],...
  'lin_lp_cutoff', [log10(1895) 0],...      
  'nlin_fc_before', [log10(2000) 0],...
  'nlin_bw_before', [log10(300) 0],...
  'nlin_lp_cutoff', [log10(2000) 0],...      
  'nlin_a', [log10(9244) 0],...
  'nlin_b', [log10(0.078) 0]};

f4000 = {...
  'flow',4000,...
  'fhigh',4000,...
  'lin_fc', [log10(3900) 0],...
  'lin_bw', [log10(620) 0],...
  'lin_gain', [log10(270) 0],...
  'lin_lp_cutoff', [log10(3900) 0],...      
  'nlin_fc_before', [log10(4000) 0],...
  'nlin_bw_before', [log10(560) 0],...
  'nlin_lp_cutoff', [log10(4000) 0],...      
  'nlin_a', [log10(30274) 0],...
  'nlin_b', [log10(0.06) 0]};

f8000 = {...
  'flow',8000,...
  'fhigh',8000,...
  'lin_fc', [log10(7450) 0],...
  'lin_bw', [log10(1550) 0],...
  'lin_gain', [log10(250) 0],...
  'lin_lp_cutoff', [log10(7450) 0],...      
  'nlin_fc_before', [log10(8000) 0],...
  'nlin_bw_before', [log10(1100) 0],...
  'nlin_lp_cutoff', [log10(8000) 0],...      
  'nlin_a', [log10(76354) 0],...
  'nlin_b', [log10(0.035) 0]};

expparsYO = [f250; f500; f1000; f2000; f4000; f8000];

%% average parameter set, table II
f250avg = {...
  'flow',250,...
  'fhigh',250,...
  'lin_fc', [log10(244) 0],...
  'lin_bw', [log10(100) 0],...
  'lin_gain', [log10(1150) 0],...
  'lin_lp_cutoff', [log10(244) 0],...      
  'lin_ngt', 3,...
  'nlin_fc_before', [log10(250) 0],...
  'nlin_bw_before', [log10(84) 0],...
  'nlin_lp_cutoff', [log10(250) 0],...      
  'nlin_a', [log10(2194) 0],...
  'nlin_b', [log10(0.45) 0]};

f500avg = {...
  'flow',500,...
  'fhigh',500,...
  'lin_fc', [log10(480) 0],...
  'lin_bw', [log10(130) 0],...
  'lin_gain', [log10(850) 0],...
  'lin_lp_cutoff', [log10(480) 0],...      
  'lin_ngt', 3,...
  'nlin_fc_before', [log10(500) 0],...
  'nlin_bw_before', [log10(103) 0],...
  'nlin_lp_cutoff', [log10(500) 0],...      
  'nlin_a', [log10(5184) 0],...
  'nlin_b', [log10(0.28) 0]};

f1000avg = {...
  'flow',1000,...
  'fhigh',1000,...
  'lin_fc', [log10(965) 0],...
  'lin_bw', [log10(240) 0],...
  'lin_gain', [log10(520) 0],...
  'lin_lp_cutoff', [log10(965) 0],...      
  'lin_ngt', 3,...
  'nlin_fc_before', [log10(1000) 0],...
  'nlin_bw_before', [log10(175) 0],...
  'nlin_lp_cutoff', [log10(1000) 0],...      
  'nlin_a', [log10(7558) 0],...
  'nlin_b', [log10(0.130) 0]};

f2000avg = {...
  'flow',2000,...
  'fhigh',2000,...
  'lin_fc', [log10(1925) 0],...
  'lin_bw', [log10(400) 0],...
  'lin_gain', [log10(410) 0],...
  'lin_lp_cutoff', [log10(1925) 0],...      
  'lin_ngt', 3,...
  'nlin_fc_before', [log10(2000) 0],...
  'nlin_bw_before', [log10(300) 0],...
  'nlin_lp_cutoff', [log10(2000) 0],...      
  'nlin_a', [log10(9627) 0],...
  'nlin_b', [log10(0.078) 0]};

f4000avg = {...
  'flow',4000,...
  'fhigh',4000,...
  'lin_fc', [log10(3900) 0],...
  'lin_bw', [log10(660) 0],...
  'lin_gain', [log10(320) 0],...
  'lin_lp_cutoff', [log10(3900) 0],...      
  'lin_ngt', 3,...
  'nlin_fc_before', [log10(4000) 0],...
  'nlin_bw_before', [log10(560) 0],...
  'nlin_lp_cutoff', [log10(4000) 0],...      
  'nlin_a', [log10(22288) 0],...
  'nlin_b', [log10(0.045) 0]};

f8000avg = {...
  'flow',8000,...
  'fhigh',8000,...
  'lin_fc', [log10(7750) 0],...
  'lin_bw', [log10(1450) 0],...
  'lin_gain', [log10(220) 0],...
  'lin_lp_cutoff', [log10(7750) 0],...      
  'lin_ngt', 3,...
  'nlin_fc_before', [log10(8000) 0],...
  'nlin_bw_before', [log10(1100) 0],...
  'nlin_lp_cutoff', [log10(8000) 0],...      
  'nlin_a', [log10(43584) 0],...
  'nlin_b', [log10(0.030) 0]};

expparsavg = [f250avg; f500avg; f1000avg; f2000avg; f4000avg; f8000avg];

%% Lopez-Poveda and Meddis 2001, Figure 2

  %% Lopez-Poveda and Meddis 2001, Figure 2, a)

  if flags.do_fig2a || flags.do_fig2
    fs=22050;
    hpdata=data_pralong1996;
    bout=headphonefilter(fs);

    % Manually calculate the frequency response.
    fout = 20*log10(abs(fftreal(bout)));

    % Half the filter length.
    n2=length(fout);
    output(:,:,1) = hpdata;
  end


  %% Lopez-Poveda and Meddis 2001, Figure 2, b)

  if flags.do_fig2b || flags.do_fig2
    fs = 22050;
    
    % data points directly derived from Goode et al. 1994
    gde = middleearfilter;           

    % control data points directly read from Lopez-Poveda and Meddis 2001
    stapes_data = data_lopezpoveda2001('fig2b', 'no_plot', 'lopezpoveda');    
    
    % Calculate the filters.
    bmid = middleearfilter(fs);      % Goode et al. 1994
    
    % Manually calculate the frequency response for an input of 0dB SPL
    fmid = abs(fftreal(bmid*20e-6));
    % Half the filter length.
    n2=length(fmid);
    % x-values for plotting.
    xplot=linspace(0,fs/2,n2);
    outB = gde;
    if exist('output','var')
      output(end+1:end+(length(outB)-length(output)),:,1) = 0;
      output(:,:,2) = outB;
    else
      output(:,:,1) = outB;
    end

  end
    
%% plots    
if flags.do_plot
    if flags.do_fig2       
        figure
        set(gcf,'Position',[50,50,500,760])
        subplot(2,1,1)
        hold on;
        % Plot the measured data
        x=hpdata(:,1);
        freqresp=20*log10(hpdata(:,2));
        plot(x,freqresp,'ro','MarkerFaceColor', 'r');
    
        % Plot the filter
        x_filter=linspace(0,fs/2,n2);
        plot(x_filter,fout);
        axis([100 10000 -30 20])
        set(gca,'XScale','log','YTick',[-30,-20,-10,0,10,20])
        set(gca,'Position',[0.15,0.55,0.8,0.4]);
        leg1=legend('Pralong and Carlile (1996) + extrapolated points', ...
            'FIR filter');
        xlabel('Frequency (Hz)')
        ylabel('Pressure gain (dB)')
        title('Lopez-Poveda and Meddis 2001, Figure 2')
        set(leg1,'Position',[0.1887, 0.5991, 0.666, 0.0553]);

        subplot(2,1,2)
        p = loglog (stapes_data(:,1),stapes_data(:,2),':ok', 'MarkerFaceColor', 'k');
        hold on
        g = loglog (gde(:,1),gde(:,2),':or', 'MarkerFaceColor', 'r');
        firG = loglog(xplot,fmid,'r');
        axis([100 10000 1E-10 1E-07])
        set(gca,'Position',[0.15,0.07,0.8,0.4]);
        xlabel('Frequency (Hz)')
        ylabel('Stapes velocity (m/s) at 0dB SPL')
        leg2=legend([g,firG,p],'directly derived from Goode et al. 1994', ...
                    'FIR filter with data points from Goode et al. 1994', ...
                    'Control points as in Lopez-Poveda and Meddis 2001');        
        set(leg2,'Position',[0.1637, 0.085, 0.716, 0.07]);
    
    elseif flags.do_fig2a
        hold on;
        % Plot the measured data
        x=hpdata(:,1);
        freqresp=20*log10(hpdata(:,2));
        plot(x,freqresp,'ro','MarkerFaceColor', 'r');
    
        % Plot the filter
        x_filter=linspace(0,fs/2,n2);
        plot(x_filter,fout);
        axis([100 10000 -30 20])
        set(gca,'XScale','log','YTick',[-30,-20,-10,0,10,20])
        legend('Pralong and Carlile (1996) + extrapolated points', ...
            'FIR filter');
        title('Lopez-Poveda and Meddis 2001, Figure 2a) - Pressure gain (dB) as a function of frequency')
        xlabel('Frequency (Hz)')
        ylabel('Pressure gain (dB)')
        hold off
        
    elseif flags.do_fig2b
        p = loglog (stapes_data(:,1),stapes_data(:,2),':ok', 'MarkerFaceColor', 'k');
        hold on
        g = loglog (gde(:,1),gde(:,2),':or', 'MarkerFaceColor', 'r');
        firG = loglog(xplot,fmid,'r');
        axis([100 10000 1E-10 1E-07])
        title('Lopez-Poveda and Meddis 2001, Figure 2b) - Stapes velocity as a function of frequency')
        xlabel('Frequency (Hz)')
        ylabel('Stapes velocity (m/s) at 0dB SPL')
        legend([g,firG,p],'directly derived from Goode et al. 1994', 'FIR filter with data points from Goode et al. 1994', ...
        'Control points as in Lopez-Poveda and Meddis 2001')
        hold off
    end
end


%% Lopez-Poveda and Meddis 2001, Figure 3

if flags.do_fig3b || flags.do_fig3bc
     
  %% Lopez-Poveda and Meddis 2001, Figure 3, b)
  % input signal: 50ms pure tones, sampled at 100kHz
  fs = 100000;
  T = 0.05;       
  t = (0:1/fs:T - 1/fs)';
  fsig = 250:25:1750;    
  
  result3b = zeros(1,length(fsig));
  lin3b = zeros(1,length(fsig));
  nlin3b = zeros(1,length(fsig));

  level = 20e-6 * 10^(30/20);
  
  for ii = 1:length(fsig)
    
    insig = sin(2*pi*fsig(ii).*t)*(2^0.5) * level;  
    
    hp_fir = headphonefilter(fs);
    insig = filter(hp_fir,1,insig);
    
    [y_lin, ~] = lopezpoveda2001(insig, fs, kv.predrnl{:}, f1000{:},'linonly', kv.postdrnl{:});    
    [y_nlin, ~] = lopezpoveda2001(insig, fs, kv.predrnl{:},f1000{:},'nlinonly', kv.postdrnl{:});
    
    outsig = y_lin + y_nlin;
    
    result3b(1,ii) = rms(outsig(floor(length(insig)/2):end));
    lin3b(1,ii) = rms(y_lin(floor(length(insig)/2):end));
    nlin3b(1,ii) = rms(y_nlin(floor(length(insig)/2):end));
  end
  
  output(:,:,1) = [fsig', lin3b', nlin3b', result3b'];
end;

if flags.do_fig3c || flags.do_fig3bc
  %% Lopez-Poveda and Meddis 2001, Figure 3, c)
  % input signal: 50ms pure tones, sampled at 100kHz
  fs = 100000;
  T = 0.05;       
  t = (0:1/fs:T - 1/fs)';
  fsig = 250:25:1750;   
  
  result3c = zeros(1,length(fsig));
  lin3c = zeros(1,length(fsig));
  nlin3c = zeros(1,length(fsig));
  
  level = 20e-6 * 10^(85/20);
  
  for ii = 1:length(fsig)
    
    insig = sin(2*pi*fsig(ii).*t)*(2^0.5)* level;
    hp_fir = headphonefilter(fs);
    insig = filter(hp_fir,1,insig);
    
    [y_lin, ~] = lopezpoveda2001(insig, fs, kv.predrnl{:}, f1000{:},'linonly', kv.postdrnl{:});    
    [y_nlin, ~] = lopezpoveda2001(insig, fs, kv.predrnl{:}, f1000{:},'nlinonly', kv.postdrnl{:});

    outsig = y_lin + y_nlin;
        
    result3c(1,ii) = rms(outsig(floor(length(insig)/2):end));
    lin3c(1,ii) = rms(y_lin(floor(length(insig)/2):end));
    nlin3c(1,ii) = rms(y_nlin(floor(length(insig)/2):end));
    
  end
  
  outB = [fsig', lin3c', nlin3c', result3c'];
  if exist('output','var')
      output(:,:,2) = outB;
  else
      output(:,:,1) = outB;
  end
end;

%% plots    
if flags.do_plot
    if flags.do_fig3bc       
        figure
        set(gcf,'Position',[50,50,400,760])
        subplot(2,1,1)
        plot(fsig,result3b)
        hold on
        plot(fsig,lin3b,'-.g')
        plot(fsig,nlin3b,':r')
        set(gca,'YScale','log')
        % grid on
        set(gca,'XLim',[250 1750],'Layer','top')
        set(gca,'YLim',[1e-07 1e-03],'Layer','top')
        set(gca,'Position',[0.285,0.5838,0.62,0.3412]);
        title('30 dB SPL')
        xlabel('Frequency (Hz)')
        ylabel('lopezpoveda2001 filter output (m/s)')

        subplot(2,1,2)
        plot(fsig,result3c)
        hold on
        plot(fsig,lin3c,'-.g')
        plot(fsig,nlin3c,':r')
        set(gca,'YScale','log')
        % grid on
        set(gca,'XLim',[250 1750],'Layer','top')
        set(gca,'YLim',[1e-05 1e-01],'Layer','top')
        set(gca,'Position',[0.285,0.11,0.62,0.3412]);
        title('85 dB SPL')
        xlabel('Frequency (Hz)')
        ylabel('lopezpoveda2001 filter output (m/s)')
        leg=legend('lopezpoveda2001 output', 'Linear path output', 'Nonlinear path output');
        set(leg,'Position',[0.0133, 0.4759, 0.4525, 0.0798]);
        
    elseif flags.do_fig3b
        plot(fsig,result3b)
        hold on
        plot(fsig,lin3b,'-.g')
        plot(fsig,nlin3b,':r')
        set(gca,'YScale','log')
        set(gca,'XLim',[250 1750],'Layer','top')
        set(gca,'YLim',[1e-07 1e-03],'Layer','top')
        title('30 dB SPL')
        xlabel('Frequency (Hz)')
        ylabel('lopezpoveda2001 filter output (m/s)')
        leg=legend('lopezpoveda2001 output', 'Linear path output', 'Nonlinear path output');
        set(gcf,'Position',[150,150,400,400])
        set(leg,'Position',[0.2333, 0.1467, 0.4525, 0.1517]);
        hold off
        
    elseif flags.do_fig3c
        plot(fsig,result3c)
        hold on
        plot(fsig,lin3c,'-.g')
        plot(fsig,nlin3c,':r')
        set(gca,'YScale','log')
        set(gca,'XLim',[250 1750],'Layer','top')
        set(gca,'YLim',[1e-05 1e-01],'Layer','top')
        title('85 dB SPL')
        xlabel('Frequency (Hz)')
        ylabel('lopezpoveda2001 filter output (m/s)')
        leg=legend('lopezpoveda2001 output', 'Linear path output', 'Nonlinear path output');
        set(gcf,'Position',[150,150,400,400])
        set(leg,'Position',[0.2333, 0.1467, 0.4525, 0.1517]);
        hold off
    end
end

%% Lopez-Poveda and Meddis 2001, Figure 4
if flags.do_fig4

  %% input signal
  fs=64000;
  fsig = [250 500 1000 2000 4000 8000];
  T = 0.084;       
  t = (0:1/fs:T - 1/fs)';
  %basef = fsig;
  hp_fir = headphonefilter(fs);
  Tramp = 0.002;          % duration of ramps up and down
  ramp = (0:1/fs:Tramp - 1/fs)';
  sig=zeros(length(t),length(fsig));
  mask=zeros(length(t),length(fsig));
  
  %% experiment
  LSDB = 30:0.5:85;           % Signal level
  n=1;
  for jj = 30:0.5:85
    levelS(n) = 20e-6 * 10^(jj/20);
    n = n+1;
  end
  
  LMDB = 30:0.5:100;          % Masker level
  n=1;
  for jj = 30:0.5:100
    levelM(n) = 20e-6 * 10^(jj/20);
    n = n+1;
  end
  
  
  OMavg = zeros(length(levelM),length(fsig));
  OM = zeros(length(levelM),length(fsig));
  OSavg = zeros(length(levelS),length(fsig));
  OS = zeros(length(levelS),length(fsig));
  ratio = zeros(length(levelM),length(levelS),length(fsig));
  ratioavg = zeros(length(levelM),length(levelS),length(fsig));
  indx = zeros(length(levelS),length(fsig));
  indxavg = zeros(length(levelS),length(fsig));
  
  output = zeros(length(LSDB),4,length(fsig));

  for ii = 1:length(fsig)
    
    % first calculate the model's response to the masker for every
    % possible masker level
    for kk = 1:length(levelM)
      % rampsignal sine window equivalent to cosine ramps are used
      mask(:,ii) = rampsignal(sin(2*pi*fsig(ii)*0.6.*t),length(ramp),'sine').*(2^0.5);
      insig = mask(:,ii) * levelM(kk);
      outsig = filter(hp_fir,1,insig);

      % average parameter set, table II
      outsigavg = lopezpoveda2001(outsig, fs, kv.predrnl{:}, expparsavg{ii,:}, kv.postdrnl{:});         

      OMavg(kk,ii) = max(outsigavg(floor(length(insig)/2):end));    
      
      % parameter set of YO, table I
      outsig = lopezpoveda2001(outsig, fs, kv.predrnl{:}, expparsYO{ii,:}, kv.postdrnl{:});             
      OM(kk,ii) = max(outsig(floor(length(insig)/2):end));            
    end
    
    % then calculate model's response to the signal and find for every
    % signal level the masker level such that the ratio signal/masker
    % is equal to a value of one
    for mm = 1:length(levelS)
      sig(:,ii) = rampsignal(sin(2*pi*fsig(ii).*t),length(ramp),'sine').*(2^0.5);
      insig = sig(:,ii) * levelS(mm);
      outsig = filter(hp_fir,1,insig);

      % average parameter set, table II
      outsigavg = lopezpoveda2001(outsig, fs, kv.predrnl{:}, expparsavg{ii,:}, kv.postdrnl{:});         

      OSavg(mm,ii) = max(outsigavg(floor(length(insig)/2):end)); 

      % parameter set of YO, table I
      outsig = lopezpoveda2001(outsig, fs, kv.predrnl{:}, expparsYO{ii,:}, kv.postdrnl{:});

      OS(mm,ii) = max(outsig(floor(length(insig)/2):end)); 
      
      % ratio for parameter set of YO, table I
      ratio(:,mm,ii) = OS(mm,ii) ./ OM(:,ii);                 

      [~, indx(mm,ii)] = min(abs(1-ratio(:,mm,ii)),[],1);     

      % ratio for average parameter set, table II           
      ratioavg(:,mm,ii) = OSavg(mm,ii) ./ OMavg(:,ii);        
      [~, indxavg(mm,ii)] = min(abs(1-ratioavg(:,mm,ii)),[],1);          
    end
  
  output(:,1,ii) = LSDB;
  output(:,2,ii) = LMDB(indx(:,ii));
  output(:,3,ii) = LMDB(indxavg(:,ii));   
  end
  
  % same procedure for lopezpoveda2001 parameters calculated with regression lines, table III
  OMrl = zeros(length(levelM),length(fsig));
  OSrl = zeros(length(levelS),length(fsig));
  ratiorl = zeros(length(levelM),length(levelS),length(fsig));
  indxrl = zeros(length(levelS),length(fsig));
  
  for ii = 1:length(fsig)
    for kk = 1:length(levelM)
      insig = mask(:,ii) * levelM(kk);
      outsig = filter(hp_fir,1,insig);
      outsigrl = lopezpoveda2001(outsig, fs, kv.predrnl{:}, 'flow', fsig(ii), 'fhigh', fsig(ii), 'lin_ngt', 3, kv.postdrnl{:});
      OMrl(kk,ii) = max(outsigrl(floor(length(insig)/2):end));                       
    end
    
    for mm = 1:length(levelS)
      insig = sig(:,ii) * levelS(mm);
      outsig = filter(hp_fir,1,insig);
      outsigrl = lopezpoveda2001(outsig, fs, kv.predrnl{:}, 'flow', fsig(ii), 'fhigh', fsig(ii), 'lin_ngt', 3, kv.postdrnl{:});
      OSrl(mm,ii) = max(outsigrl(floor(length(insig)/2):end)); 
      
      ratiorl(:,mm,ii) = OSrl(mm,ii) ./ OMrl(:,ii);           
      [~, indxrl(mm,ii)] = min(abs(1-ratiorl(:,mm,ii)),[],1);                  
    end
    
  output(:,4,ii) = LMDB(indxrl(:,ii));
    
  end
  
  
  %% plots    
  if flags.do_plot
    subplot(2,3,1)
    plot(LSDB,LMDB(indx(:,1)),'-', 'LineWidth', 2)
    hold on
    plot(LSDB,LMDB(indxavg(:,1)),'-')
    plot(LSDB,LMDB(indxrl(:,1)),'--')
    plot(LSDB,LSDB,':')
    grid on
    set(gca,'XLim',[20 90],'Layer','top','YTick',[55,65,75,85,95])
    axis([20,90,55,100]);
    title('250 Hz')
    xlabel('Signal level (dB SPL)')
    ylabel('Masker level (dB SPL)')
    
    subplot(2,3,2)
    plot(LSDB,LMDB(indx(:,2)),'-', 'LineWidth', 2)
    hold on
    plot(LSDB,LMDB(indxavg(:,2)),'-')
    plot(LSDB,LMDB(indxrl(:,2)),'--')
    plot(LSDB,LSDB,':')
    grid on
    set(gca,'XLim',[20 90],'Layer','top','YTick',[55,65,75,85,95])
    axis([20,90,55,100]);
    title('500 Hz')
    xlabel('Signal level (dB SPL)')
    ylabel('Masker level (dB SPL)')
    
    subplot(2,3,3)
    plot(LSDB,LMDB(indx(:,3)),'-', 'LineWidth', 2)
    hold on
    plot(LSDB,LMDB(indxavg(:,3)),'-')
    plot(LSDB,LMDB(indxrl(:,3)),'--')
    plot(LSDB,LSDB,':')
    grid on
    set(gca,'XLim',[20 90],'Layer','top','YTick',[55,65,75,85,95])
    axis([20,90,55,100]);
    title('1000 Hz')
    xlabel('Signal level (dB SPL)')
    ylabel('Masker level (dB SPL)')
    
    subplot(2,3,4)
    plot(LSDB,LMDB(indx(:,4)),'-', 'LineWidth', 2)
    hold on
    plot(LSDB,LMDB(indxavg(:,4)),'-')
    plot(LSDB,LMDB(indxrl(:,4)),'--')
    plot(LSDB,LSDB,':')
    grid on
    set(gca,'XLim',[20 90],'Layer','top','YTick',[55,65,75,85,95])
    axis([20,90,55,100]);
    title('2000 Hz')
    xlabel('Signal level (dB SPL)')
    ylabel('Masker level (dB SPL)')
    
    subplot(2,3,5)
    plot(LSDB,LMDB(indx(:,5)),'-', 'LineWidth', 2)
    hold on
    plot(LSDB,LMDB(indxavg(:,5)),'-')
    plot(LSDB,LMDB(indxrl(:,5)),'--')
    plot(LSDB,LSDB,':')
    grid on
    set(gca,'XLim',[20 90],'Layer','top','YTick',[55,65,75,85,95])
    axis([20,90,55,100]);
    title('4000 Hz')
    xlabel('Signal level (dB SPL)')
    ylabel('Masker level (dB SPL)')
   
    subplot(2,3,6)
    plot(LSDB,LMDB(indx(:,6)),'-', 'LineWidth', 2)
    hold on
    plot(LSDB,LMDB(indxavg(:,6)),'-')
    plot(LSDB,LMDB(indxrl(:,6)),'--')
    plot(LSDB,LSDB,':')
    grid on
    set(gca,'XLim',[20 90],'Layer','top','YTick',[55,65,75,85,95])
    axis([20,90,55,100]);
    title('8000 Hz')
    xlabel('Signal level (dB SPL)')
    ylabel('Masker level (dB SPL)')
    legend('parameter set of YO, table I','average parameter set, table II','regression lines, table III','linear behavior') 
  end      
  
end;


