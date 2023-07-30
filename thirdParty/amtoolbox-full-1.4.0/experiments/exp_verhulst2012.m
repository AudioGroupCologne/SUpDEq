function output=exp_verhulst2012(varargin)
%EXP_VERHULST2012 Figures from Verhulst et al. (2012)
%
%   Usage: output = exp_verhulst2012(flag)
%
%   This script reproduces figures 2a and 2c from Verhulst et al. (2012). "Nonlinear 
%   time-domain cochlear model for transient stimulation and human otoacoustic emission."
%   The Journal of the Acoustical Society of America 132(6). pages 3842-3848.
%
%   Requirements and installation: 
%   ------------------------------
%
%   1) Python 3 is required with numpy and scipi packages. 
%      On Linux, use sudo apt-get install python-scipy python-numpy
%      On Windows, use 
% 
%   2) Compiled files with a C-compiler, e.g. gcc. In amtbase/bin/verhulst2012 start make (Linux) or make.bat (Windows)
%
%   3) On linux, when problems with GFORTRAN lib appear, try sudo ln -sf /usr/lib64/libgfortran.so.3.0.0 /mymatlabroot/sys/os/glnxa64/libgfortran.so.3 (mymatlabroot is usually /usr/local/MATLAB/version
%               
%   Examples:
%   ---------
%
%   To display Figure 2a from the Verhulst et al. (2012) use:
%
%     exp_verhulst2012('fig2a');
%
%   To display Figure 2c from the Verhulst et al. (2012) use:
%
%     exp_verhulst2012('fig2c');
%
%   References:
%     S. Verhulst, T. Dau, and C. A. Shera. Nonlinear time-domain cochlear
%     model for transient stimulation and human otoacoustic emission. J.
%     Acoust. Soc. Am., 132(6):3842 -- 3848, 2012.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_verhulst2012.php


%   #Author: Alessandro Altoe
%   #Author: Piotr Majdak (2021)
%   #Author: Alejandro Osses (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.import={'amt_cache'};
definput.flags.type={'missingflag','fig2a','fig2c'};

definput.flags.plot={'plot','no_plot'};

[flags,~]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%% ------ FIG 2a -----------------------------------------------------------
%BM displacement simulated for the 1-kHz cochlear CF location for clicks with intensities between 0 and 90 dB peSPL. The displacements were normalized 
%by the pressure at the stapes of the cochlea such that compression is observed as a reduction of the IR amplitude 
if flags.do_fig2a

    fs=48000; % Hz -- same sampling frequency for both figures

    dur = 50e-3; % 50 ms, duration of the test signals
    dt = 1/fs;
    t=0:dt:dur-dt;

    spl=10:10:90;
    Nr_signals=length(spl); % number of test signals  
  
    %%% Preparing the test signals (clicks):
    insig = zeros(Nr_signals,length(t)); % memory allocation
    insig(:,8:9)=2; % impulse (the value 2 is to have 2 peak-to-peak value)
    
    p0 = 2e-5;
    for i = 1:length(spl)
        insig(i,:) = p0*10^(spl(i)/20)*insig(i,:);
    end  

    %%% 2. Runs the model:
    fc = 1000; % Hz, one probe frequency will be requested
    [V,Y,OAE,CF]=verhulst2012(insig,fs,fc,spl);    
    output.y = Y;
    output.OAE = OAE;
    output.v = V;
    output.cf = CF;
    output.description = '''y, v'' are the BM displacement and velocity, respectively. ''OAE'' is the simulated middle-ear sound pressure';
    
    %%% 4. Plots the results:
    if flags.do_plot
        for i = 1:Nr_signals
          outsigs(i).y = Y(:,i);
          outsigs(i).OAE = OAE(:,i);
          outsigs(i).v = V(:,i);
        end
        leg_str=repmat('dB',Nr_signals,1);
        figure;
        for i=1:Nr_signals
            plot(t*1e3,outsigs(i).y./max(abs(outsigs(i).y)));
            hold all;
        end
        legend((horzcat(int2str(spl'),leg_str)));
        grid off;
        xlabel('Time (ms)');
        ylabel('Normalised y_{BM}');
        axis 'tight';
        title('Displacement in response to a click');
        xlim([0 25]);
    end
end

%% ------ FIG 2c -----------------------------------------------------------
%cochlear excitation patterns calculated as the rms level of displacement per cochlear section,
%for stimulation with a pure tone of 1 kHz with stimulus intensities between 10 and 90 dB SPL
%if flags.do_fig2c

    fs=48000; % Hz -- same sampling frequency for both figures

    dur = 50e-3; % 50 ms, duration of the test signals
    dt = 1/fs;
    t=0:dt:dur-dt;

    spl=10:10:90;
    Nr_signals=length(spl); % number of test signals
    
    %%% Preparing the test signals:
    % 1. Calibration of the test signals (pre-ramp):
   f0 = 1000; % Hz, centre frequency of the sinudoids
   insig=ones(Nr_signals,1)*(sin(2*pi*f0*t)); % base input signal
   for i = 1:Nr_signals
       dBFS = 94; % i.e., amplitude 1 is equal to 94 dB SPL
       insig(i,:) = scaletodbspl(insig(i,:),spl(i),dBFS);
   end
    % 2. Creating a linear ramp (fade in):
   dur_ramp_ms = 10; % ms. Ramp duration or 'onset duration'
   dur_ramp = round((dur_ramp_ms*1e-3)*fs); % duration ramp in samples
   rp = ones(1,size(insig,2)); 
   rp(1:dur_ramp+1) = 0:1/dur_ramp:1;  % fade in 10ms
    
    % 3. Applying the linear ramp:
   insig = insig.*repmat(rp,Nr_signals,1); % Applying the ramp up
    %%%
    
    %%% 2. Runs the model:
   fc_flag='all';
   norm_Rms=ones(Nr_signals,1);
   irr = ones(Nr_signals,1)'; % Zweig irregularities are turned on
   [~,Y,~,cf]=verhulst2012(insig,fs,fc_flag,spl,'normalize',norm_Rms,'irr',irr);
    %%% 3. Processes the model outputs:
   disp_ref = 0.01; % reference displacement in m
   for ii=1:Nr_signals
       Y_rms(ii,:) = squeeze(20*log10(rms(Y(:,:,ii))/disp_ref));
   end
   output.cf   = cf;
   output.Yrms = Y_rms;
   output.Yrms_description = 'Average (RMS) cochlear displacement [dB re 0.1 m]';
      
    %%% 4. Plots the results:
   if flags.do_plot  
       leg_str=repmat('dB',Nr_signals,1);
       figure;
       semilogx(cf(2:length(cf)),(Y_rms(:,2:length(cf))));
       hold all;
            
       legend((horzcat(int2str(spl'),leg_str)));
       grid on;
       set(gca, 'xdir','reverse');
       set(gca,'Xtick',[250 500 1000 2000 4000 8000 16000]);
       set(gca,'XTickLabel',{'0,25','0.5','1','2','4','8','16'});
       xlabel('Centre frequency (KHz)');
       ylabel('Y_{rms} (dB re .01)');
       axis([250 8000 -250 -100]);
       title('BM displacement in response to 1 KHz sinsusoid');
   end
end


