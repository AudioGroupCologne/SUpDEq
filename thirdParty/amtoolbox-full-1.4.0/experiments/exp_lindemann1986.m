function output = exp_lindemann1986(varargin)
%EXP_LINDEMANN1986 Figures from Lindemann (1986)
%   Usage: output = exp_lindemann1986(flag)
%
%   EXP_LINDEMANN1986(flag) reproduces the results for the figure given
%   by flag from the Lindemann (1986) paper. It will also plot the
%   results.  The format of its output depends on the chosen figure. If not
%   otherwise stated, the cross-correlation of a pure tone sinusoids with
%   f=500 Hz and different ITDs, ILDs or a combination of both is
%   calculated. Because of the stationary character of the input signals,
%   T_{int} = inf is used to produce only one time step in the crosscorr
%   output.
%   
%   The following flags can be specified;
%
%     'plot'     plot the output of the experiment. This is the default.
%
%     'no_plot'   Don't plot, only return data.
%
%     'auto'     Re-calculate the file if it does not exist. Return 1 if the
%                file exist, otherwise 0. This is the default
%
%     'redo'     Always recalculate the results.
%
%     'cached'   Always use the cached version. Default.
%
%     'fig6'  Reproduce Fig.6 from Lindemann (1986).  The cross-correlation is
%             calculated for different ITDs and different inhibition factors
%             c_s=0,0.3,1. Afterwards for every c_s the correlation is
%             plotted for every used ITD dependend on the correlation-time
%             delay. The output is cross-correlation result of the figure.
%             The output has dimensions: number of c_s conditions x
%             nitds x delay line length
%
%     'fig7'  Reproduce Fig.7 from Lindemann (1986).  The cross-correlation
%             is calculated for different ITDs and different inhibition
%             factors c_s=0,0.2,0.4,0.6,0.8,1.0. Afterwards for every
%             c_s the displacement of the centroid of the auditory image
%             is calculated and plotted depending on the ITD. The output is
%             the displacement of the centroid for different c_s values.
%             The output has dimensions: c_s x nitds
%
%     'fig8'  Reproduce Fig.8 from Lindemann (1986).  The cross-correlation
%             is calculated for different ILDs and different inhibition factors
%             c_s = 0.3, 1. Afterwards for every c_s the ILD is plotted
%             depending on the correlation time. The output is the
%             cross-correlation result of the figure. The dimensions of the
%             output are: number of c_s conditions x nilds x delay line length.
%
%     'fig10'  Reproduce Fig.10 from Lindemann (1986). The
%              cross-correlation is calculated for different ILDs with an
%              inhibition factor of c_s = 0.3 and a monaural detector
%              factor w_f = 0.035. Afterwards the ILD is plotted depending
%              on the correlation time.  The output is the cross-correlation
%              result of the figure.  The output has dimensions:
%              number of c_s conditions x nilds x delay line length
%
%     'fig11'  Reproduce Fig.11 from Lindemann (1986).  The centroid
%              position is calculated for different ILDs, an inhibition
%              factor c_s = 0.3 and a monaural detector factor w_f =
%              0.035. Afterwards for every c_s the displacement of the
%              centroid is plotted dependend on the ILD. The output is the
%              cross-correlation result of the to figure. The dimensions of
%              the output are: number of c_s conditions x nilds x delay
%              line length
%
%     'fig12'  Reproduce Fig.12 from Lindemann (1986). The centroids are
%              calculated for combinations of ITDs and ILDs.  After the
%              calculation the values for the centroids of the stimuli are
%              searched to find the nearest value to 0. The corresponding ILD
%              value is stored in output and plotted depending on the
%              ITD. The output is the simulated results for a trading
%              experiment. The ILD value for getting a centroid near the
%              center for an combined ITD, ILD stimulus with a given ITD
%              value.  The output has dimensions: number of ITDs x 1
%
%     'fig13'  Reproduce Fig.13 from Lindemann (1986). The centroids are
%              calculated for ILD only and ITD/ILD combination
%              stimuli. After the calculation the values for the centroids
%              of the ILD only stimuli are searched to find the nearest
%              value to the centroid of a given combined stimulus. The
%              resulting ILD value is stored for different combinaition
%              values and plotted dependend on ITD. The output is the ILD
%              value for getting the same lateralization with an ILD only
%              stimulus compared to a stimulus with both ITD and ILD.
%              The output has dimensions: number of ILDs x number of ITDs
%
%     'fig14a'  Reproduce fig.14 (a) from Lindemann (1986). The
%               cross-correlations for a combination of ILDs and a ITD of
%               ~1ms are calculated. This is done for different ILDs and
%               different inhibition factors c_s = 0,0.2,0.4,0.6,1.
%               Afterwards for every c_s the centroid of the auditory image
%               is calculated and plotted dependend on the ILD. The output is
%               the displacement of the centroid for different c_s values
%               and ILDs. The output has dimensions: c_s x nilds
%
%     'fig14b'  Reproduce Fig.14 (b) from Lindemann (1986). The
%               cross-correlations for a combination of ILDs and a ITD of ~1ms
%               are calculated. This is done for different small ILDs with a
%               standard deviation of 0-5. Afterwards for every standard
%               deviation the mean centroid displacement is calculated and
%               plotted dependend on the ITD. The output is the displacement
%               of the centroid as a function of the ITD averaged over
%               different small ILD with a standard deviation of 0,1,2,3,4,5.
%               The output has dimensions: nilds x nitds
%               THIS FIGURE DOES NOT SEEM TO MATCH THE PAPER!
%
%     'fig15'  Reproduce Fig.15 from Lindemann (1986). The cross-correlation
%              for an ITD of -0.5ms is calculated. This is done for
%              different ILDs, an inhibition factor c_s = 0.3 and a
%              monaural detector factor w_f = 0.035. Afterwards for every
%              c_s the ILD is plotted depending on the correlation
%              time. The output is the cross-correlation result. The output
%              has dimensions: number of c_s conditions x nilds x delay line
%              length
%
%     'fig16'  Reproduces Fig.16 from Lindemann (1986). The
%              cross-correlations for combinations of ITDs and ILDs are
%              calculated. Afterwards the combinations of ILD and ITD are
%              looked for the ones that have two peaks in its
%              cross-correlation which have the nearly same height. The
%              corresponding ILD value is than stored in output and plotted
%              dependend on the ITD. The output is the ILD value for which
%              the cross-correlation for a combined ITD, ILD stimulus has
%              two peaks with the same (nearest) height, depending on the
%              ITD. The output has dimensions: number of ITDs x 1
%
%     'fig17'  Reproduce Fig.17 from Lindemann (1986). The
%              cross-correlation for ITD/ILD combination and ITD only
%              stimuli is calculated. Afterwards the values for the
%              centroids and maxima of the ITD only stimuli are searched to
%              find the nearest value to the centroid and maxima of a given
%              combined stimulus. The resulting ITD value is stored for
%              different combinaition values. The output is the ITD
%              value for getting the same lateralization with an ITD only
%              stimulus compared to a stimulus with both ITD and ILD using
%              the centroid of the cross-correlation.  The output has
%              dimensions: number of ITDs x number of ILDs. 
%
%     'fig18'  Reproduce Fig.18 from Lindemann (1986). The cross-correlation
%              of pink noise with different interaural coherence values is
%              calculated.  Afterwards for every interaural coherence value
%              the correlation is plotted dependend on the correlation-time
%              delay. The output is the cross-correlation result of the
%              figure. The output has dimensions: number of interaural
%              coherences x delay line length
%
%   If no flag is given, the function will print the list of valid flags.
%
%   Examples:
%   ---------
%
%   To display Figure 6 use :
%
%     exp_lindemann1986('fig6');
%
%   To display Figure 7 use :
%
%     exp_lindemann1986('fig7');
%
%   To display Figure 8 use :
%
%     exp_lindemann1986('fig8');
%
%   To display Figure 10 use :
%
%     exp_lindemann1986('fig10');
%
%   To display Figure 11 use :
%
%     exp_lindemann1986('fig11');
%
%   To display Figure 12 use :
%
%     exp_lindemann1986('fig12');
%
%   To display Figure 13 use :
%
%     exp_lindemann1986('fig13');
%
%   To display Figure 14a use :
%
%     exp_lindemann1986('fig14a');
%
%   To display Figure 15 use :
%
%     exp_lindemann1986('fig15');
%
%   To display Figure 16 use :
%
%     exp_lindemann1986('fig16');
%
%   To display Figure 17 use :
%
%     exp_lindemann1986('fig17');
%
%   To display Figure 18 use :
%
%     exp_lindemann1986('fig18');
%
%   References:
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. I. Simulation of lateralization for
%     stationary signals. J. Acoust. Soc. Am., 80:1608--1622, 1986.
%     
%
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_lindemann1986.php


%  #Author: Hagen Wierstorf

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.import={'amt_cache'};
definput.flags.type={'missingflag','fig6','fig7','fig8','fig10','fig11',...
                    'fig12','fig13','fig14a','fig14b','fig15',...
                    'fig16','fig17','fig18'};

definput.flags.plot={'plot','no_plot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% ------ FIG 6 -----------------------------------------------------------
if flags.do_fig6
    
  % Sampling rate
  fs = 44100;
  % Frequency of the sinusoid
  f = 500;
  T = 1/f;
  fc = round(freqtoerb(f));   % corresponding frequency channel
  
  % Model parameter
  T_int = inf;
  w_f = 0;
  M_f = 6; % not used, if w_f==0
  c_s = [0,0.3,1];
  
  % NOTE: the longer the signal, the more time we need for computation. On the
  % other side N_1 needs to be long enough to eliminate any onset effects.
  % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
  % results for this demo.
  N_1 = ceil(25*T*fs);
  siglen = ceil(30*T*fs);
  
  % Calculate crosscorrelations for 21 ITD points between 0~ms and 1~ms
  nitds = 21; % number of used ITDs
  ndl = round(fs/1000)+1;   % length of the delay line (see lindemann1986_bincorr.m)
  itd = linspace(0,1,nitds);
  
  output=amt_cache('get','fig6',flags.cachemode);
  
  if isempty(output)
    output = zeros(length(c_s),nitds,ndl);
    for ii = 1:nitds; 
        % Generate ITD shifted sinusoid
        sig = sig_itdsin(f,itd(ii),fs);
        % Use only the beginning of the signal to generate only one time instance of
        % the cross-correlation and apply onset window
        sig = sig(1:siglen,:);
        sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
        % Calculate cross-correlation for different inhibition factor c_s
        for jj = 1:length(c_s)
            % Calculate cross-correlation (and squeeze due to T_int==inf)
            tmp = squeeze(lindemann1986(sig,fs,c_s(jj),w_f,M_f,T_int,N_1));
            % Store the needed frequency channel. NOTE: the cross-correlation
            % calculation starts with channel 5, so we have to subtract 4.
            output(jj,ii,:) =  tmp(:,fc-4);
        end
    end

    amt_cache('set','fig6',output);
  end;
    
  if flags.do_plot
    % ------ Plotting ------
    % Generate time axis
    tau = linspace(-1,1,ndl);
    % Plot figure for every c_s condition
    for jj = 1:length(c_s)
      figure;
      mesh(tau,itd(end:-1:1),squeeze(output(jj,:,:)));
      view(0,57);
      xlabel('correlation-time tau (ms)');
      ylabel('interaural time difference (ms)');
      set(gca,'YTick',0:0.2:1);
      set(gca,'YTickLabel',{'1','0.8','0.6','0.4','0.2','0'});
      tstr = sprintf('c_s = %.1f\nw_f = 0\nf = %i Hz\n',c_s(jj),f);
      title(tstr);
    end
  end;

end;


%% ------ FIG 7 -----------------------------------------------------------
if flags.do_fig7
    
    % Sampling rate
    fs = 44100;
    % Frequency of the sinusoid
    f = 500;
    T = 1/f;
    fc = round(freqtoerb(f));   % corresponding frequency channel

    % Model parameter
    T_int = inf;
    w_f = 0;
    M_f = 6; % not used, if w_f==0
    c_s = 0:0.2:1;

    % NOTE: the longer the signal, the more time we need for computation. On the
    % other side N_1 needs to be long enough to eliminate any onset effects.
    % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
    % results for this demo.
    N_1 = ceil(25*T*fs);
    siglen = ceil(30*T*fs);

    % Calculate crosscorrelations for 21 ITD points between 0~ms and 1~ms
    nitds = 21; % number of used ITDs
    itd = linspace(0,1,nitds);
    
    output=amt_cache('get','fig7',flags.cachemode);
    if isempty(output)
      
      output = zeros(length(c_s),nitds);
      for ii = 1:nitds 
        % Generate ITD shifted sinusoid
        sig = sig_itdsin(f,itd(ii),fs);
        % Use only the beginning of the signal to generate only one time instance of
        % the cross-correlation and apply an onset window
        sig = sig(1:siglen,:);
        sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
        % Calculate cross-correlation for different inhibition factor c_s 
        for jj = 1:length(c_s)
          % Calculate cross-correlation (and squeeze due to T_int==inf)
          tmp = squeeze(lindemann1986(sig,fs,c_s(jj),w_f,M_f,T_int,N_1));
          % Store the needed frequency channel. NOTE: the cross-correlation
          % calculation starts with channel 5, so we have to subtract 4.
          cc = tmp(:,fc-4);
          % Calculate the position of the centroid
          output(jj,ii) = lindemann1986_centroid(cc);
        end
      end
      amt_cache('set','fig7',output);
    end;
    
    if flags.do_plot
      % ------ Plotting ------
      figure;
      for jj = 1:length(c_s)
        plot(itd,output(jj,:));
        hold on;
      end
      xlabel('interaural time difference (ms)');
      ylabel('displacement of the centroid d');
      tstr = sprintf('w_f = 0\nf = %i Hz\n',f);
      title(tstr);
    end;
  
end;


%% ------ FIG 8 -----------------------------------------------------------
if flags.do_fig8
  
  % Sampling rate
  fs = 44100;
  % Frequency of the sinusoid
  f = 500;
  T = 1/f;
  fc = round(freqtoerb(f));   % corresponding frequency channel
  
  % Model parameter
  T_int = inf;
  w_f = 0;
  M_f = 6; % not used, if w_f==0
  c_s = [0.3,1];

  % NOTE: the longer the signal, the more time we need for computation. On the
  % other side N_1 needs to be long enough to eliminate any onset effects.
  % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
  % results for this demo.
  N_1 = ceil(25*T*fs);
  siglen = ceil(30*T*fs);
  
  % Calculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
  nilds = 26; % number of used ILDs
  ndl = 2*round(fs/2000)+1;   % length of the delay line (see bincorr.m)
  ild = linspace(0,25,nilds);
  
  output=amt_cache('get','fig8',flags.cachemode);
  if isempty(output)
    
    output = zeros(2,nilds,ndl);
    for ii = 1:nilds 
      % Generate sinusoid with given ILD
      sig = sig_ildsin(f,ild(ii),fs);
      % Use only the beginning of the signal to generate only one time instance of
      % the cross-correlation and apply onset window
      sig = sig(1:siglen,:);
      sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
      % Calculate cross-correlation for different inhibition factor c_s 
      for jj = 1:length(c_s)
        % Calculate cross-correlation (and squeeze due to T_int==inf)
          tmp = squeeze(lindemann1986(sig,fs,c_s(jj),w_f,M_f,T_int,N_1));
          % Store the needed frequency channel. NOTE: the cross-correlation
          % calculation starts with channel 5, so we have to subtract 4.
          output(jj,ii,:) = tmp(:,fc-4);
      end
    end
    
    amt_cache('set','fig8',output);
  end;
  
  if flags.do_plot
    % ------ Plotting ------
    % Generate time axis
    tau = linspace(-1,1,ndl);
    % Plot figure for every c_s condition
    for jj = 1:length(c_s) 
      figure;
      mesh(tau,ild(end:-1:1),squeeze(output(jj,:,:)));
      view(0,57);
      xlabel('correlation-time delay (ms)');
      ylabel('interaural level difference (dB)');
      set(gca,'YTick',0:5:25);
      set(gca,'YTickLabel',{'25','20','15','10','5','0'});
      tstr = sprintf('c_{inh} = %.1f\nw_f = 0\nf = %i Hz\n',c_s(jj),f);
      title(tstr);
    end
  end;
  
end;


%% ------ FIG 10 ----------------------------------------------------------
if flags.do_fig10

  % Sampling rate
  fs = 44100;
  % Frequency of the sinusoid
  f = 500;
  T = 1/f;
  fc = round(freqtoerb(f));   % corresponding frequency channel
  
  % Model parameter
  T_int = inf;
  w_f = 0.035;
  M_f = 6; % not used, if w_f==0
  c_s = 0.3;
  
  % NOTE: the longer the signal, the more time we need for computation. On the
  % other side N_1 needs to be long enough to eliminate any onset effects.
  % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
  % results for this demo.
  N_1 = ceil(25*T*fs);
  siglen = ceil(30*T*fs);
  
  % Calculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
  nilds = 26; % number of used ILDs
  ndl = 2*round(fs/2000)+1;   % length of the delay line (see bincorr.m)
  ild = linspace(0,25,nilds);
  
  output=amt_cache('get','fig10',flags.cachemode);
  if isempty(output)

    output = zeros(length(c_s),nilds,ndl);
    for ii = 1:nilds 
      % Generate sinusoid with given ILD
      sig = sig_ildsin(f,ild(ii),fs);
      % Use only the beginning of the signal to generate only one time instance of
      % the cross-correlation and apply onset window
      sig = sig(1:siglen,:);
      sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
      % Calculate cross-correlation for different inhibition factor c_s 
      for jj = 1:length(c_s)
        % Calculate cross-correlation (and squeeze due to T_int==inf)
        tmp = squeeze(lindemann1986(sig,fs,c_s(jj),w_f,M_f,T_int,N_1));
        % Store the needed frequency channel. NOTE: the cross-correlation
        % calculation starts with channel 5, so we have to subtract 4.
        output(jj,ii,:) = tmp(:,fc-4);
      end
    end
    
    amt_cache('set','fig10',output);
  end;
  
  
  % ------ Plotting ------
  if flags.do_plot
    % Generate time axis
    tau = linspace(-1,1,ndl);
    % Plot figure for every c_s condition
    for jj = 1:length(c_s)
      figure;
      mesh(tau,ild(end:-1:1),squeeze(output(jj,:,:)));
      view(0,57);
      xlabel('correlation-time delay (ms)');
      ylabel('interaural level difference (dB)');
      set(gca,'YTick',0:5:25);
      set(gca,'YTickLabel',{'25','20','15','10','5','0'});
      tstr = sprintf('c_{inh} = %.1f\nw_f = 0.035\nf = %i Hz\n',c_s(jj),f);
      title(tstr);
    end
  end;
  
end;


%% ------ FIG 11 ----------------------------------------------------------
if flags.do_fig11

    % Sampling rate
    fs = 44100;
    % Frequency of the sinusoid
    f = 500;
    T = 1/f;
    fc = round(freqtoerb(f));   % corresponding frequency channel

    % Model parameter
    T_int = inf;
    w_f = [0,0,0.035];
    M_f = 6; % not used, if w_f==0
    c_s = [0.3,1,0.3];

    % NOTE: the longer the signal, the more time we need for computation. On the
    % other side N_1 needs to be long enough to eliminate any onset effects.
    % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
    % results for this demo.
    N_1 = ceil(25*T*fs);
    siglen = ceil(30*T*fs);

    % Caelculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
    nilds = 26; % number of used ILDs
    ild = linspace(0,25,nilds);
          
    output=amt_cache('get','fig11',flags.cachemode);
    if isempty(output)
      
      output = zeros(length(c_s),nilds);
      for ii = 1:nilds 
        % Generate sinusoid with given ILD
        sig = sig_ildsin(f,ild(ii),fs);
        % Use only the beginning of the signal to generate only one time instance of
        % the cross-correlation and apply onset window
        sig = sig(1:siglen,:);
        sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
        % Calculate cross-correlation for different inhibition factor c_s 
        for jj = 1:length(c_s)
          % Calculate cross-correlation (and squeeze due to T_int==inf)
          tmp = squeeze(lindemann1986(sig,fs,c_s(jj),w_f(jj),M_f,T_int,N_1));
          % Store the needed frequency channel
          cc = tmp(:,fc-4);
          % Calculate the position of the centroid
          output(jj,ii) = lindemann1986_centroid(cc);
        end
      end
      
      amt_cache('set','fig11',output);
    end;


    if flags.do_plot
      % ------ Plotting ------
      figure;
      % Plot data from experiments
      data = data_lindemann1986('fig11_yost');
      plot(data(:,1),data(:,2)*0.058,'+');
      hold on;
      data = data_lindemann1986('fig11_sayers');
      plot(data(:,1),data(:,2)*0.058,'*');
      % Plot line for every condition
      for jj = 1:length(c_s) 
        plot(ild,output(jj,:));
      end
      axis([0 25 0 0.8]);
      legend('500 Hz, Yost','600 Hz, Sayers','500 Hz, model');
      xlabel('interaural level difference (dB)');
      ylabel('displacement of the centroid d');
    end;

end;


%% ------ FIG 12 ----------------------------------------------------------
if flags.do_fig12

    % Sampling rate
    fs = 44100;
    % Frequency of the sinusoid
    f = 500;
    T = 1/f;
    fc = round(freqtoerb(f));   % corresponding frequency channel

    % Model parameter
    T_int = inf;
    w_f = 0.035;
    M_f = 6; % not used, if w_f==0
    c_s = 0.3;

    % NOTE: the longer the signal, the more time we need for computation. On the
    % other side N_1 needs to be long enough to eliminate any onset effects.
    % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
    % results for this demo.
    N_1 = ceil(25*T*fs);
    siglen = ceil(30*T*fs);

    % Calculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
    nilds = 21; % number of used ILDs for the ILD only stimuli
    nitds = 21; % number of used ITDs for the combined stimuli
    ild = linspace(0,10,nilds);
    itd = linspace(-1,0,nitds);

    output=amt_cache('get','fig12',flags.cachemode);
    if isempty(output)

    
      % Calculate the centroids for ILD+ITD stimuli
      cen = zeros(nitds,nilds);
      for ii = 1:nitds
        for jj = 1:nilds
          % Generate sinusoid with given ILD
          sig = sig_itdildsin(f,itd(ii),ild(jj),fs);
          % Use only the beginning of the signal to generate only one time 
          % instance of the cross-correlation and apply a linear onset window
          sig = sig(1:siglen,:);
          sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
          % Calculate cross-correlation (and squeeze due to T_int==inf)
          tmp = squeeze(lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1));
          % Store the needed frequency channel
          cc = tmp(:,fc-4);
          % Calculate the position of the centroid
          cen(ii,jj) = lindemann1986_centroid(cc);
        end
      end
      
      
      % ------ Fiting ------
      % For every ITD find the ILD with gives a centroid near 0
      output = zeros(nitds,1);
      for ii = 1:nitds
        % Find centroid nearest to 0
        [val,idx] = min(abs(cen(ii,:)));
        output(ii) = ild(idx);
      end
      
      amt_cache('set','fig12',output);
    end;

    if flags.do_plot
      % ------ Plotting ------
      figure;
      data = data_lindemann1986('fig12_400');
      plot(data(:,1),data(:,2),'+');
      hold on;
      data = data_lindemann1986('fig12_600');
      plot(data(:,1),data(:,2),'*');
      plot(itd,output);
      axis([-1 0 0 10]);
      legend('400 Hz','600 Hz','500 Hz, model');
      xlabel('interaural time difference (ms)');
      ylabel('interaural level difference (dB)');
    end;

end;


%% ------ FIG 13 ----------------------------------------------------------
if flags.do_fig13

    % Sampling rate
    fs = 44100;
    % Frequency of the sinusoid
    f = 500;
    T = 1/f;
    fc = round(freqtoerb(f));   % corresponding frequency channel

    % Model parameter
    T_int = inf;
    N_1 = 1;
    w_f = 0.035;
    M_f = 6; % not used, if w_f==0
    c_s = 0.3;

    % NOTE: the longer the signal, the more time we need for computation. On the
    % other side N_1 needs to be long enough to eliminate any onset effects.
    % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
    % results for this demo.
    N_1 = ceil(25*T*fs);
    siglen = ceil(30*T*fs);

    % Calculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
    nilds_p = 21; % number of used ILDs for the ILD only stimuli
    nitds_t = 21; % number of used ITDs for the combined stimuli
    nilds_t = 6;  % number of used ILDs for the combined stimuli
    ndl = 2*round(fs/2000)+1;     % length of the delay line (see bincorr.m)
    ild_p = linspace(-10,40,nilds_p);
    itd_t = linspace(-1,1,nitds_t);
    ild_t = [-3,0,3,9,15,25];

    output=amt_cache('get','fig13',flags.cachemode);
    if isempty(output)
      
      % Calculate the centroids for the ILD only stimuli
      cen_p = zeros(1,nilds_p);
      for ii = 1:nilds_p
        % Generate sinusoid with given ILD
        sig = sig_ildsin(f,ild_p(ii),fs);
        % Use only the beginning of the signal to generate only one time instance of
        % the cross-correlation and apply a linear onset window
        sig = sig(1:siglen,:);
        sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
        % Calculate cross-correlation (and squeeze due to T_int==inf)
        tmp = squeeze(lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1));
        % Store the needed frequency channel
        cc = tmp(:,fc-4);
        % Calculate the position of the centroid
        cen_p(ii) = lindemann1986_centroid(cc);
      end
      
      % Calculate the centroids for the combined stimuli
      cen_t = zeros(nilds_t,nitds_t);
      for ii = 1:nitds_t
        for jj = 1:nilds_t
          sig = sig_itdildsin(f,itd_t(ii),ild_t(jj),fs);
          sig = sig(1:siglen,:);
          sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
          tmp = squeeze(lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1));
          cc = tmp(:,fc-4);
          cen_t(jj,ii) = lindemann1986_centroid(cc);
        end
      end

      % ------ Fitting ------
      % For the results for the combined centroids find the nearest centroid for the
      % ILD only stimuli
      output = zeros(nilds_t,nitds_t);
      for ii = 1:nitds_t
        for jj = 1:nilds_t
          idx = findnearest(cen_p,cen_t(jj,ii));
          output(jj,ii) = ild_p(idx);
        end
      end
      
      amt_cache('set','fig13',output);
    end;
    
    if flags.do_plot
      % ------ Plotting ------
      % First plot the only data from experiments
      figure;
      d = data_lindemann1986('fig13');
      plot(d(:,1),d(:,2),'x-r', ...   % -3dB
           d(:,1),d(:,3),'x-b', ...   %  3dB
           d(:,1),d(:,4),'x-g', ...   %  9dB
           d(:,1),d(:,5),'x-b', ...   % 15dB
           d(:,1),d(:,6),'x-r')       % 25dB
      legend('25dB','15dB','9dB','3dB','-3dB');
      axis([-1 1 -10 40]);
      set(gca,'XTick',-1:0.4:1);
      xlabel('interaural time difference (ms)');
      ylabel('interaural level difference (dB)');
      % Then plot model results
      figure;
      % Plot line for every condition
      for jj = 1:nilds_t
        plot(itd_t,output(jj,:));
        hold on;
      end
      axis([-1 1 -10 40]);
      set(gca,'XTick',-1:0.4:1);
      xlabel('interaural time difference (ms)');
      ylabel('interaural level difference (dB)');
    end;

end;


%% ------ FIG 14a ---------------------------------------------------------
if flags.do_fig14a

    % Sampling rate
    fs = 44100;
    % Frequency of the sinusoid
    f = 500;
    T = 1/f;
    fc = round(freqtoerb(f));   % corresponding frequency channel

    % --- FIG 14 (a) ---
    % Model parameter
    T_int = inf;
    w_f = 0;
    M_f = 6; % not used, if w_f==0
    c_s = [0.0,0.2,0.4,0.6,1.0];

    % NOTE: the longer the signal, the more time we need for computation. On the
    % other side N_1 needs to be long enough to eliminate any onset effects.
    % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
    % results for this demo.
    N_1 = ceil(25*T*fs);
    siglen = ceil(30*T*fs);

    % Calculate crosscorrelations for 21 ITD points between 0~ms and 1~ms
    nilds = 21; % number of used ITDs
    ild = linspace(-3,3,nilds);
    itd = 2000/f;
    
    output=amt_cache('get','fig14a',flags.cachemode);
    if isempty(output)
     
      output = zeros(length(c_s),nilds);
      for ii = 1:nilds 
        % Generate ITD shifted sinusoid
        sig = sig_itdildsin(f,itd,ild(ii),fs);
        % Use only the beginning of the signal to generate only one time instance of
        % the cross-correlation
        sig = sig(1:siglen,:);
        % Apply a linear onset window with length N_1/2 to minimize onset effects
        % (see lindemann1986 p. 1614)
        sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
        % Calculate cross-correlation for different inhibition factor c_s 
        for jj = 1:length(c_s)
          % Calculate cross-correlation (and squeeze due to T_int==inf)
          tmp = squeeze(lindemann1986(sig,fs,c_s(jj),w_f,M_f,T_int,N_1));
          % Store the needed frequency channel. NOTE: the cross-correlation
          % calculation starts with channel 5, so we have to subtract 4.
            cc = tmp(:,fc-4);
            % Calculate the position of the centroid
            output(jj,ii) = lindemann1986_centroid(cc);
        end
      end
      
      amt_cache('set','fig14a',output);
    end;
    
    if flags.do_plot
      % ------ Plotting ------
      figure;
      for jj = 1:length(c_s)
        plot(ild,output(jj,:));
        hold on;
      end
      xlabel('interaural level difference (dB)');
      ylabel('displacement of the centroid d');
      tstr = sprintf('w_f = 0\nf = %i Hz\n',f);
      title(tstr);
    end;

end;


%% ------ FIG 14b ---------------------------------------------------------
if flags.do_fig14b

    % Sampling rate
    fs = 44100;
    % Frequency of the sinusoid
    f = 500;
    T = 1/f;
    fc = round(freqtoerb(f));   % corresponding frequency channel

    % Model parameter
    T_int = inf;
    w_f = 0;
    M_f = 6; % not used, if w_f==0
    c_s = 1.0;

    % NOTE: the longer the signal, the more time we need for computation. On the
    % other side N_1 needs to be long enough to eliminate any onset effects.
    % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
    % results for this demo.
    N_1 = ceil(25*T*fs);
    siglen = ceil(30*T*fs);

    % Calculate crosscorrelations for 21 ITD points between 0~ms and 1~ms
    nitds = 21; % number of used ITDS
    nilds = 201; % number of used ILDs
    itd = linspace(0,1,nitds);
    ild_std = [1,2,3,4,5];
    
    output=amt_cache('get','fig14b',flags.cachemode);
    if isempty(output)

      amt_disp('NOTE: this test function will need a lot of time!');

      output = zeros(length(ild_std)+1,nitds);
      centmp = zeros(nilds,nitds);
      for ii = 1:nitds
        % Show progress
        amt_disp([num2str(ii) ' of ' num2str(nitds)]);
        % First generate the result for std(ILD) == 0
        sig = sig_itdsin(f,itd(ii),fs);
        sig = sig(1:siglen,:);
        sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
        tmp = squeeze(lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1));
        cc = tmp(:,fc-4);
        cen(1,ii) = lindemann1986_centroid(cc);
        % Generate results for std(ILD) ~= 0
        for nn = 1:length(ild_std)
          % Generate normal distributed ILDs with mean ~ 0 and std = 1
          tmp = randn(nilds,1);
          tmp = tmp/std(tmp);
          % Generate desired ILD distribution
          ild = tmp * ild_std(nn);
          % For all distributed ILD values calculate the centroid
          for jj = 1:nilds
            sig = sig_itdildsin(f,itd(ii),ild(jj),fs);
            sig = sig(1:siglen,:);
            sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
            tmp = squeeze(lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1));
            cc = tmp(:,fc-4);
            centmp(jj,ii) = lindemann1986_centroid(cc);
          end
          % Calculate the mean centroid above the ILD distribution
          output(nn+1,ii) = mean(centmp(:,ii));
        end
      end
      
      amt_cache('set','fig14b',output);
    end;

    if flags.do_plot
      figure;
      for nn = 1:length(ild_std)+1
        plot(itd,output(nn,:)); hold on;
      end
      axis([0 1 0 1]);
      xlabel('interaural level difference (dB)');
      ylabel('displacement of the centroid d');
      tstr = sprintf('w_f = 0\nc_{inh} = 1\nf = %i Hz\n',f);
      title(tstr);
    end;

end;


%% ------ FIG 15 ----------------------------------------------------------
if flags.do_fig15

    % Sampling rate
    fs = 44100;
    % Frequency of the sinusoid
    f = 500;
    T = 1/f;
    fc = round(freqtoerb(f));   % corresponding frequency channel

    % Model parameter
    T_int = inf;
    w_f = 0.035;
    M_f = 6; % not used, if w_f==0
    c_s = 0.3;

    % NOTE: the longer the signal, the more time we need for computation. On the
    % other side N_1 needs to be long enough to eliminate any onset effects.
    % Lindemann uses N_1 = 17640 (200*T). Here I uses only N_1 = 2205 (25*T) 
    % which gives the same results for this demo.
    N_1 = ceil(25*T*fs);
    siglen = ceil(30*T*fs);

    % Calculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
    nilds = 26; % number of used ILDs
    ndl = 2*round(fs/2000)+1;   % length of the delay line (see bincorr.m)
    ild = linspace(0,25,nilds);
    itd = -0.5;
    
    output=amt_cache('get','fig15',flags.cachemode);
    if isempty(output)
        
      output = zeros(length(c_s),nilds,ndl);
      for ii = 1:nilds 
        % Generate sinusoid with given ILD
        sig = sig_itdildsin(f,itd,ild(ii),fs);
        % Use only the beginning of the signal to generate only one time instance of
        % the cross-correlation and apply onset window
        sig = sig(1:siglen,:);
        sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
        % Calculate cross-correlation for different inhibition factor c_s 
        for jj = 1:length(c_s)
          % Calculate cross-correlation (and squeeze due to T_int==inf)
            tmp = squeeze(lindemann1986(sig,fs,c_s(jj),w_f,M_f,T_int,N_1));
            % Store the needed frequency channel. NOTE: the cross-correlation
            % calculation starts with channel 5, so we have to subtract 4.
            output(jj,ii,:) = tmp(:,fc-4);
        end
      end

      amt_cache('set','fig15',output);
    end;

    if flags.do_plot
      % ------ Plottinqg --------------------------------------------------------
      % Generate time axis
      tau = linspace(-1,1,ndl);
      % Plot figure for every c_s condition
      for jj = 1:length(c_s)
        figure;
        mesh(tau,ild(end:-1:1),squeeze(output(jj,:,:)));
        view(0,57);
        xlabel('correlation-time delay (ms)');
        ylabel('interaural level difference (dB)');
        set(gca,'YTick',0:5:25);
        set(gca,'YTickLabel',{'25','20','15','10','5','0'});
        tstr = sprintf('c_{inh} = %.1f\nw_f = 0.035\nf = %i Hz\n',c_s(jj),f);
        title(tstr);
      end
    end;

end;


%% ------ FIG 16 ---------------------------------------------------------
if flags.do_fig16

  
  s = [mfilename('fullpath'),'_fig16.mat'];

    % Sampling rate
    fs = 44100;
    % Frequency of the sinusoid
    f = 500;
    T = 1/f;
    fc = round(freqtoerb(f));   % corresponding frequency channel

    % Model parameter
    T_int = inf;
    w_f = 0.035;
    M_f = 6; % not used, if w_f==0
    c_s = 0.3;

    % NOTE: the longer the signal, the more time we need for computation. On the
    % other side N_1 needs to be long enough to eliminate any onset effects.
    % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
    % results for this demo.
    N_1 = ceil(25*T*fs);
    siglen = ceil(30*T*fs);

    % Calculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
    nilds = 21; % number of used ILDs for the ILD only stimuli
    nitds = 21; % number of used ITDs for the combined stimuli
    ndl = 2*round(fs/2000)+1;     % length of the delay line (see bincorr.m)
    ild = linspace(0,10,nilds);
    itd = linspace(-1,0,nitds);

    output=amt_cache('get','fig16',flags.cachemode);
    if isempty(output)
      
      % Calculate the centroids for ILD+ITD stimuli
      cc = zeros(nitds,nilds,ndl);
      for ii = 1:nitds
        for jj = 1:nilds
          % Generate sinusoid with given ILD
          sig = sig_itdildsin(f,itd(ii),ild(jj),fs);
          % Use only the beginning of the signal to generate only one time 
          % instance of the cross-correlation and apply a linear onset window
          sig = sig(1:siglen,:);
          sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
          % Calculate cross-correlation (and squeeze due to T_int==inf)
          tmp = squeeze(lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1));
          % Store the needed frequency channel
          cc(ii,jj,:) = tmp(:,fc-4);
        end
      end

      
      % ------ Fiting ------
      % For every ITD find the ILD with gives a centroid near 0
      output = zeros(nitds,1);
      for ii = 1:nitds
        % Find two peaks of same height
        idx = findpeaks(squeeze(cc(ii,:,:)));
        output(ii) = ild(idx);
      end

      amt_cache('set','fig16',output);
    end;
    
      
    if flags.do_plot
      % ------ Plotting ------
      figure;
      data = data_lindemann1986('fig16');
      plot(data(:,1),data(:,2),'+');
      hold on;
      plot(itd,output);
      axis([-1 0 0 15]);
      xlabel('interaural time difference (ms)');
      ylabel('interaural level difference (dB)');
    end;

end;


%% ------ FIG 17 ----------------------------------------------------------
if flags.do_fig17

  
    % Sampling rate
    fs = 44100;
    % Frequency of the sinusoid
    f = 500;
    T = 1/f;
    fc = round(freqtoerb(f));   % corresponding frequency channel

    % Model parameter
    T_int = inf;
    w_f = 0.035;
    M_f = 6; % not used, if w_f==0
    c_s = 0.3;

    % NOTE: the longer the signal, the more time we need for computation. On the
    % other side N_1 needs to be long enough to eliminate any onset effects.
    % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
    % results for this demo.
    N_1 = ceil(25*T*fs);
    siglen = ceil(30*T*fs);

    % Calculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
    nitds_p = 21; % number of used ITDs for the ITD only stimuli
    nitds_t = 4; % number of used ITDs for the combined stimuli
    nilds_t = 21;  % number of used ILDs for the combined stimuli
    itd_p = linspace(-0.36,0.72,nitds_p);
    ild_t = linspace(-9,9,nilds_t);
    itd_t = [0,0.09,0.18,0.27];

    output=amt_cache('get','fig17',flags.cachemode);
    if isempty(output)
      
      % Calculate the centroids for the ITD only stimuli
      cen_p = zeros(1,nitds_p);
      max_p = zeros(1,nitds_p);
      for ii = 1:nitds_p
        % Generate sinusoid with given ILD
        sig = sig_itdsin(f,itd_p(ii),fs);
        % Use only the beginning of the signal to generate only one time instance of
        % the cross-correlation and apply a linear onset window
        sig = sig(1:siglen,:);
        sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
        % Calculate cross-correlation (and squeeze due to T_int==inf)
        tmp = squeeze(lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1));
        % Store the needed frequency channel
        cc = tmp(:,fc-4);
        % Find the maximum position
        max_p(ii) = findmax(cc);
        % Calculate the position of the centroid
        cen_p(ii) = lindemann1986_centroid(cc);
      end
      
      % Calculate the centroids for the combined stimuli
      cen_t = zeros(nitds_t,nilds_t);
      max_t = zeros(nitds_t,nilds_t);
      for ii = 1:nilds_t
        for jj = 1:nitds_t
          sig = sig_itdildsin(f,itd_t(jj),ild_t(ii),fs);
          sig = sig(1:siglen,:);
          sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
          tmp = squeeze(lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1));
          cc = tmp(:,fc-4);
          max_t(jj,ii) = findmax(cc);
          cen_t(jj,ii) = lindemann1986_centroid(cc);
        end
      end
      
      % ------ Fiting ------
      % For the results for the combined centroids find the nearest centroid for the
      % ITD only stimuli
      output.rmax = zeros(nitds_t,nilds_t);
      output.rcen = output.rmax;
      for ii = 1:nilds_t
        for jj = 1:nitds_t
          idx = findnearest(cen_p,cen_t(jj,ii));
          output.rcen(jj,ii) = itd_p(idx);
          idx = findnearest(max_p,max_t(jj,ii));
          output.rmax(jj,ii) = itd_p(idx);
        end
      end
      
      amt_cache('set','fig17',output);
    end;
    
      
    if flags.do_plot
      % ------ Plotting ------
      % First plot the only data from experiments
      figure; % fig 17 (a)
      d = data_lindemann1986('fig17');
      plot(d(:,1),d(:,5),'x-r', ...   %  0 ms
           d(:,1),d(:,4),'x-b', ...   %  0.09 ms
           d(:,1),d(:,3),'x-g', ...   %  0.18 ms
           d(:,1),d(:,2),'x-b')       %  0.27 ms
      legend('0.27 ms','0.18 ms','0.09 ms','0 ms');
      axis([-9 9 -0.36 0.72]);
      set(gca,'XTick',-9:3:9);
      xlabel('interaural level difference (dB)');
      ylabel('interaural time difference (ms)');
      % Then plot model results
      figure; % fig 17 (b)
              % Plot line for every condition
      for jj = 1:nitds_t
        plot(ild_t,output.rcen(jj,:));
        hold on;
      end
      axis([-9 9 -0.36 0.72]);
      set(gca,'XTick',-9:3:9);
      xlabel('interaural level difference (dB)');
      ylabel('interaural time difference (ms)');
      %
      figure; % fig 17 (c)
      for jj = 1:nitds_t
        plot(ild_t,output.rmax(jj,:));
        hold on;
      end
      axis([-9 9 -0.36 0.72]);
      set(gca,'XTick',-9:3:9);
      xlabel('interaural level difference (dB)');
      ylabel('interaural time difference (ms)');
    end;

end;


%% ------ FIG 18 ----------------------------------------------------------
if flags.do_fig18

  
    % Sampling rate
    fs = 44100;
    % Frequency of the sinusoid
    f = 500;
    T = 1/f;
    fc = round(freqtoerb(f));   % corresponding frequency channel

    % Model parameter
    T_int = inf;
    w_f = 0.035;
    M_f = 6; % not used, if w_f==0
    c_s = 0.11;

    % NOTE: the longer the signal, the more time we need for computation. On the
    % other side N_1 needs to be long enough to eliminate any onset effects.
    % Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
    % results for this demo.
    N_1 = ceil(25*T*fs);

    % Calculate crosscorrelations for 21 ITD points between 0~ms and 1~ms
    niacs = 11; % number of used interaural coherence values
    ndl = 2*round(fs/2000)+1;   % length of the delay line (see bincorr.m)
    iac = linspace(0,1,niacs);
    
    output=amt_cache('get','fig18',flags.cachemode);
    if isempty(output)

      output = zeros(niacs,ndl);
      for ii = 1:niacs; 
        % Generate ITD shifted sinusoid
        sig = sig_bincorrnoise(fs,iac(ii),'pink');
        % Aplly onset window
        sig = rampsignal(sig,[round(N_1/2)-1 0],'tria');
        % Calculate cross-correlation (and squeeze due to T_int==inf)
        tmp = squeeze(lindemann1986(sig,fs,c_s,w_f,M_f,T_int,N_1));
        % Store the needed frequency channel. NOTE: the cross-correlation
        % calculation starts with channel 5, so we have to subtract 4.
        output(ii,:) =  tmp(:,fc-4);
      end
      amt_cache('set','fig18',output);
    end;

    if flags.do_plot
      % ------ Plotting ------
      % Generate time axis
      tau = linspace(-1,1,ndl);
      % Plot figure for every c_s condition
      figure;
      mesh(tau,iac,output);
      view(0,57);
      xlabel('correlation-time tau (ms)');
      ylabel('degree of interaural coherence');
    end;

end;


%% ------ Subfunctions ----------------------------------------------------
function idx = findnearest(list,val)
[dif,idx] = min(abs(list-val));

function idx = findpeaks(cc)
% FINDPEAKS Finds two peaks of the same height in a given cross-correlation
h = zeros(size(cc,1),1);
for ii = 1:size(cc,1)
    % Find a peak in the two halfes of the cross-correlation and subtract them
    h(ii) = max(cc(ii,1:30)) - max(cc(ii,31:end));
end
% Find the index with the smallest distance between the two peaks
[val,idx] = min(abs(h));

function m = findmax(list)
[m,idx] = max(list);
d = linspace(-1,1,length(list));
m = d(idx);


