function plot_stmodspecgram(f,fs,varargin)
%PLOT_STMODSPECGRAM  Spectro-Temporal Modulation spectrogram
%   Usage:  plot_stmodspecgram(f,fs);
%           plot_stmodspecgram(f,fs,...);
%
%   PLOT_STMODSPECGRAM(f,fs) plots the modulation spectogram of the signal f sampled
%   at a sampling frequency of fs Hz.
%
%   C=PLOT_STMODSPECGRAM(f,fs, ... ) returns the image to be displayed as a matrix. Use
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
%     'smfmax',smfmax
%                  Display smfmax as the highest spectral modulation
%                  frequency.  
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
%     'nocolorbar'  Do not display the colorbar.
%   
%     'interp'     Interpolate the image to get the desired
%                  x-resolution. Turn this off by using 'no_interp'
%
%   The parameters 'dynrange', 'mfmax' and 'smfmax' may be speficied
%   first on the argument line, in that order.
%
%   See also:  plot_modspecgram
%
%   References:
%     T. Elliott and F. Theunissen. The modulation transfer function for
%     speech intelligibility. PLoS Computational Biology, 5(3), 2009.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_stmodspecgram.php


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
definput.flags.plottype={'image','contour','mesh','pcolor'};

definput.flags.clim={'no_clim','clim'};
definput.flags.log={'db','lin'};
definput.flags.colorbar={'colorbar','no_colorbar'};
definput.flags.interp={'interp','no_interp'};

definput.keyvals.win=[];
definput.keyvals.wlen=fs*.01; % Default window, 10 ms.
definput.keyvals.thr=0;
definput.keyvals.clim=[0,1];
definput.keyvals.climsym=1;
definput.keyvals.fmax=[];
definput.keyvals.mfmax=[];
definput.keyvals.smfmax=[];
definput.keyvals.dynrange=[];
definput.keyvals.xres=800;
definput.keyvals.yres=600;

[flags,kv]=ltfatarghelper({'dynrange','mfmax','smfmax'},definput,varargin,mfilename);

  
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

if isempty(kv.win)
  % Authors use a Gaussian with a width in time of 10 ms and 16 Hz in
  % frequency.
  g={'gauss','width',kv.wlen};
end;

% Discrete Gabor transform, use the complex version, so we avoid having to deal with
% boundary conditions in frequency
c = dgt(f,g,a,M);
L=size(c,2)*a/fs;    % Length of the transform, in seconds

nwin = size(c,2);    % Number of time windows
nfbin = size(c,1);   % Number of frequency bins

% Get the log-spectrogram
s=log(abs(c)+eps);

% Convert to spectro-temporal modulation domain
% Use regular fft along time and fftreal along frequency.
st = fft(fftreal(s),[],2);

% XXX Verify that different hopsizes does not change the absolute
% values. Perhaps a proper scaling is needed.

% Go to dB. XXX Use 20 or 10? How do we interpret this?
st = 20*log10(abs(st)+eps);

if flags.do_interp
  st=fftresample(st,kv.xres,2);
end;

% HACK: The DC-term (which is just a single pixel) can have a value which is
% much larger than all other pixel, so a large part of the dynamic range
% would be spent on just this one pixel (often 10-20 dB). To fix this, set
% the value of this pixel to the avarage of its neighbours.
st(1,1)=(st(2,1)+st(1,2))/2;

% 'dynrange' parameter is handled by thresholding the coefficients.
if ~isempty(kv.dynrange)
  maxclim=max(st(:));
  st(st<maxclim-kv.dynrange)=maxclim-kv.dynrange;
end;

% Center the plot such that "0 temporal modulations" is in the middle
st=fftshift(st,2);

% Determine ticks for the x-axis.
subbandfs=fs/a;
mfmax=subbandfs/2;
mfbins=size(st,2);
xr = linspace(-mfmax,mfmax,mfbins);

% Determine ticks for the y-axis
spectral_fs=M/fs;
smfmax=spectral_fs/2;
smfbins=size(st,1);
% *1000, because we use samples/kHz.
yr = linspace(0,smfmax*1000,smfbins);

switch(flags.plottype)
  case 'image'
    if flags.do_clim
      imagesc(xr,yr,st,kv.clim);
    else
      imagesc(xr,yr,st);
    end;
  case 'contour'
    contour(xr,yr,st);
  case 'surf'
    surf(xr,yr,st);
  case 'pcolor'
    pcolor(xr,yr,st);
end;

if flags.do_colorbar
  colorbar;
end;

axis('xy');

title('Spectro-Temporal modulation spectrogram')
xlabel('Modulation Frequency (Hz)')
ylabel('Spectral Modulation Frequency (Cycles/kHz)')


