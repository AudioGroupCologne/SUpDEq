function output=exp_verhulst2012(varargin)
%EXP_VERHULST2012 Compute figurs from the Verhulst paper
%
%   Usage: output = exp_verhulst2012(flag)
%
%   This script reproduces figures 2a and 2c from Verhulst et al.  "Nonlinear time-domain cochlear
%   model for transient stimulation and human otoacoustic emission."  The
%   Journal of the Acoustical Society of America 132.6 (2012): 3842-3848.
%
%   Requirements and installation: 
%   ------------------------------
%
%   1) Python >2.6 is required with numpy and scipi packages. On Linux, use sudo apt-get install python-scipy python-numpy
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
%     Acoust. Soc. Am., 132(6):3842 - 3848, 2012.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/experiments/exp_verhulst2012.php

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

%   AUTHOR: Alessandro Altoe, Piotr Majdak

definput.import={'amtredofile'};
definput.flags.type={'missingflag','fig2a','fig2c'};

definput.flags.plot={'plot','noplot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;



%% ------ FIG 2c -----------------------------------------------------------
%cochlear excitation patterns calculated as the rms level of displacement per cochlear section,
%for stimulation with a pure tone of 1 kHz with stimulus intensities between 10 and 90 dB SPL
if flags.do_fig2c

  fs=48000;
  f0=1e3;
  t=0:1/fs:0.05;
  levels_n0=9;
  sig=ones(levels_n0,1)*(sin(2*pi*f0*t));
  onset_duration=1e-2*fs;
  sig(:,1:round(onset_duration)+1)=sig(:,1:round(onset_duration)+1).*(ones(levels_n0,1)*(0:1.0/round(onset_duration):1)); %fade in 10ms
  spl=10:80/(levels_n0-1):90;
  fc='all';
  norm_Rms=ones(levels_n0,1);
  irr=ones(levels_n0,1)';
  [v,y,e,cf]=verhulst2012(sig,fs,fc,spl,norm_Rms,1,irr);
  Y_rms=20.*log10(rms(y(:,:,:)./0.01));
  output=Y_rms;
  
  if flags.do_plot  
    leg_str=repmat('dB',levels_n0,1);
    figure;
    for i=1:levels_n0
    semilogx(cf(2:length(cf)),(Y_rms(:,2:length(cf),i)));
    hold all;
    end
    legend((horzcat(int2str(spl'),leg_str)));
    grid on;
    set(gca, 'xdir','reverse');
    set(gca,'Xtick',[250 500 1000 2000 4000 8000 16000]);
    set(gca,'XTickLabel',{'0,25','0.5','1','2','4','8','16'});
    xlabel('Center frequency (KHz)');
    ylabel('Y_{rms} (dB re .01)');
    axis([250 8000 -250 -100]);
    title('BM displacement in response to 1 KHz sinsusoid');
  end
end

%% ------ FIG 2a -----------------------------------------------------------
%BM displacement simulated for the 1-kHz cochlear CF location for clicks with intensities between 0 and 90 dB peSPL. The displacements were normalized 
%by the pressure at the stapes of the cochlea such that compression is observed as a reduction of the IR amplitude 
if flags.do_fig2a
  fs=48000;
  f0=1e3;
  t=0:1/fs:0.05;
  levels_n0=9;
  fc=[1e3];
  norm_Rms=ones(levels_n0,1);
  norm_Rms=norm_Rms.*0;
  sig=ones(levels_n0,1)*(sin(2*pi*f0*t));
  onset_duration=1e-2*fs;
  sig(:,1:round(onset_duration)+1)=sig(:,1:round(onset_duration)+1).*(ones(levels_n0,1)*(0:1.0/round(onset_duration):1)); %fade in 10ms  
  sig=sig.*0;  
  sig(:,8:9)=2; %impulse (the value 2 is to have 2 peak-to-peak value)
  spl=10:80/(levels_n0-1):90;
  [v,y,e,cf]=verhulst2012(sig,fs,fc,spl,norm_Rms);
  output.v=v;
  output.y=y;
  output.e=e;
  output.cf=cf;
  
  if flags.do_plot
    leg_str=repmat('dB',levels_n0,1);
    figure;
    for i=1:levels_n0
      plot(t*1e3,y(:,1,i)./max(abs(y(:,i))));
      hold all;
    end
    legend((horzcat(int2str(spl'),leg_str)));
    grid off;
    xlabel('Time (ms)');
    ylabel('Normalized y_{BM}');
    axis 'tight';
    title('Displacement in response to a click');
    xlim([0 25]);
  end
end
