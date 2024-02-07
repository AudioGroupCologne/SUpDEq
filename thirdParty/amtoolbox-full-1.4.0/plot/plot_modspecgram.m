function plot_modspecgram(f,fs,varargin)
%PLOT_MODSPECGRAM  Modulation spectrogram
%   Usage:  plot_modspecgram(f,fs);
%           plot_modspecgram(f,fs,...);
%
%   PLOT_MODSPECGRAM(f,fs) plot the modulation spectogram of the signal f sampled
%   at a sampling frequency of fs Hz.
%
%   C=PLOT_MODSPECGRAM(f,fs, ... ) returns the image to be displayed as a matrix. Use
%   this in conjunction with imwrite etc.
%
%   The function takes the following additional arguments
%
%     'win',g      Use the window g. See the help on gabwin for some
%                  possiblities. Default is to use a Gaussian window
%                  controlled by the 'thr' or 'wlen' parameters listed below.
%   
%     'tfr',v      Set the ratio of frequency resolution to time resolution.
%                  A value v=1 is the default. Setting v>1 will give better
%                  frequency resolution at the expense of a worse time
%                  resolution. A value of 0<v<1 will do the opposite.
%   
%     'wlen',s     Window length. Specifies the length of the window
%                  measured in samples. See help of PGAUSS on the exact
%                  details of the window length.
%   
%     'image'      Use imagesc to display the spectrogram. This is the
%                  default.
%   
%     'clim',clim  Use a colormap ranging from clim(1) to clim(2). These
%                  values are passed to imagesc. See the help on imagesc.
%   
%     'dynrange',r  Use a colormap in the interval [chigh-r,chigh], where
%                   chigh is the highest value in the plot.
%   
%     'fmax',fmax  Display fmax as the highest frequency.
%   
%     'mfmax',mfmax  
%                  Display mfmax as the highest modulation frequency.
%   
%     'xres',xres  Approximate number of pixels along x-axis / time.
%   
%     'yres',yres  Approximate number of pixels along y-axis / frequency
%   
%     'contour'    Do a contour plot to display the spectrogram.
%            
%     'surf'       Do a surf plot to display the spectrogram.
%   
%     'mesh'       Do a mesh plot to display the spectrogram.
%   
%     'colorbar'   Display the colorbar. This is the default.
%   
%     'no_colorbar'  Do not display the colorbar.
%   
%     'interp'     Interpolate the image to get the desired
%                  x-resolution. Turn this off by using 'no_interp'
%
%   The parameters 'dynrange' and 'mfmax' may be speficied first on the
%   argument line, in that order.
%
%   Examples:
%   ---------
%
%   The first example shows a Modulation spectrogram of modulated wide band
%   noise. The modulation frequency is 50 Hz and the signal is sampled at
%   44.1 kHz:
%
%     fm = 50;    % Modulation frequency
%     l = 2;      % Length of the signal in seconds
%     fs = 44100; % Sampling frequency
%     t = 0:1/fs:l;
%     n = length(t);
%     noise = 1-2*randn(1,n);
%     modnoise = noise.*(1+cos(2*pi*t*fm));
%     plot_modspecgram(modnoise,fs,90)
%     title('Sinusoidally Modulated Noise')
%
%   The second example shows a modulation spectrogram of a speech signal
%   sampled at 16 kHz:
%
%     plot_modspecgram(greasy,16000,60,500)
%     title('Greasy')
%
%   The third example shows a modulation spectrogram of a modulated sinusoid
%   with a carrier frequency of 5 kHz. FIXME: What is the modulation
%   frequency:
%
%     fm = 50;    % Modulation frequency
%     l = 2;      % Length of the signal in seconds
%     fs = 44100; % Sampling frequency
%     fc = 5000;  % Carrier frequency
%     t = 0:1/fs:l;
%     s = sin(2*pi*t*fc);
%     smod = s.*(1+0.5*cos(2*pi*t*fm));
%     plot_modspecgram(s,fs,50,2*fm,'fmax',2*fc)
%
%   See also:  plot_audspecgram
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_modspecgram.php


%   #Author: Peter L. Sondergaard

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
  
if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;


% Define initial value for flags and key/value pairs.
definput.flags.wlen={'no_wlen','wlen'};
definput.flags.plottype={'image','contour','mesh','pcolor'};

definput.flags.clim={'no_clim','clim'};
definput.flags.log={'db','lin'};
definput.flags.colorbar={'colorbar','no_colorbar'};
definput.flags.interp={'interp','no_interp'};

definput.keyvals.tfr=1;
definput.keyvals.win=[];
definput.keyvals.wlen=0;
definput.keyvals.thr=0;
definput.keyvals.clim=[0,1];
definput.keyvals.climsym=1;
definput.keyvals.fmax=[];
definput.keyvals.mfmax=[];
definput.keyvals.dynrange=[];
definput.keyvals.xres=800;
definput.keyvals.yres=600;

[flags,kv]=ltfatarghelper({'dynrange','mfmax'},definput,varargin,mfilename);

  
l = (length(f)-1)/fs;   % Length of the signal (in seconds)

% Downsample
if ~isempty(kv.fmax)
  resamp=kv.fmax*2/fs;

  f=fftresample(f,round(length(f)*resamp));
  fs=kv.fmax*2;
end;

M=kv.yres*2;

if ~isempty(kv.mfmax)
  % Choose "a" such that the subband sampling rate mathes the desired
  % mfmax.
  a=max(round((fs/2)/kv.mfmax),1);
else
  % Choose "a" such that the number of pixels is matched.
  a=max(round(length(f)/2/kv.xres),1);
end;

% Set an explicit window length, if this was specified.
if flags.do_wlen
  kv.tfr=kv.wlen^2/L;
end;

if isempty(kv.win)
  g={'gauss',kv.tfr};  
end;

% Discrete Gabor transform
fdgt = dgtreal(f,g,a,M);
L=size(fdgt,2)*a/fs;    % Length of the transform, in seconds

nwin = size(fdgt,2);    % Number of time windows
nfbin = size(fdgt,1);   % Number of frequency bins

% Do fft along second dimension (time). Normalize so different values of
% "a" does not change the magnitude.
gram=fftreal(abs(fdgt),[],2);

% Go to dB
gram = 20*log10(abs(gram));

if flags.do_interp
  gram=dctresample(gram,kv.xres,2);
end;

subbandfs=fs/a;
mfmax=subbandfs/2;
mfbins=size(gram,2);

% 'dynrange' parameter is handled by thresholding the coefficients.
if ~isempty(kv.dynrange)
  maxclim=max(gram(:));
  gram(gram<maxclim-kv.dynrange)=maxclim-kv.dynrange;
end;

%% Plotting the modulation spectrum

% Frequency axis
yr = linspace(0,fs/2,nfbin);
xr = linspace(0,mfmax,mfbins);

switch(flags.plottype)
  case 'image'
    if flags.do_clim
      imagesc(xr,yr,gram,kv.clim);
    else
      imagesc(xr,yr,gram);
    end;
  case 'contour'
    contour(xr,yr,gram);
  case 'surf'
    surf(xr,yr,gram);
  case 'pcolor'
    pcolor(xr,yr,gram);
end;

if flags.do_colorbar
  colorbar;
end;

axis('xy');

title('Temporal modulation spectrogram')
xlabel('Modulation Frequency (Hz)')
ylabel('Frequency (Hz)')



