function varargout=plot_relanoiborra2019(insig,fs,fc,varargin)
%PLOT_RELANOIBORRA2019  Auditory spectrogram.
%   Usage: plot_relanoiborra2019(insig,fs,op1,op2, ... );
%          C=plot_relanoiborra2019(insig,fs, ... );
%
%   AUDSPECGRAM(insig,fs) plots an auditory spectrogram of the signal insig,
%   which has been sampled at a sampling rate of fs Hz. The output is
%   low-pass modulation filtered before presentation.
%
%   The frequency axis is diplayed on a erb-scale, but labelled in
%   Hz. Using the mouse to get plot coordinates will reveal the real
%   value in erb's. Use ERBTOFREQ to convert to Hz.
%
%   C=AUDSPECGRAM(insig,fs, ... ) returns the image to be displayed as a
%   matrix. Use this in conjunction with IMWRITE etc. Do NOT use this as a
%   method to compute an auditory representation. Use some of the model
%   preprocessing functions for this.
%
%   Be carefull with long signals, as the routine may lock up the
%   interpreter.
%
%   Additional arguments can be supplied like this:
%   AUDSPECGRAM(insig,fs,'dynrange',30). The arguments must be character
%   strings possibly followed by an argument:
%
%    'adapt'   - Model adaptation. This is the default. This options also
%                sets the output to be displayed on a linear scale.
%
%    'no_adt' - Do not model adaptation. This option also sets a Db scale to
%                display the output.
%
%    'rectify' - Use half-wave rectification followed by low-pass filtering to
%                extract the envelope of the auditory filters. This is the
%                default.
%
%    'hilbert' - Use the Hilbert transform to extract the envelope of the
%                auditory filters.
%
%    'classic' - Display a classic spectrogram. This option is equal to
%               'hilbert', 'no_adt', 'no_mf'
%
%    'mlp',f   - Modulation low-pass filter to frequency f. Default is to
%                low-pass filter to 50 Hz.
%
%    'no_mf'    - No modulation filtering of any kind.
%
%    'image'   - Use 'imagesc' to display the spectrogram. This is the default.
%
%    'clim',[clow,chigh] - Use a colormap ranging from clow to chigh. These
%               values are passed to IMAGESC. See the help on IMAGESC.
%
%    'dynrange',r - Limit the displayed dynamic range to r. This option
%                is especially usefull when displaying on a log/Db scale.
%
%    'fullrange' - Use the full dynamic range.
%
%    'ytick'   - A vector containing the frequency in Hz of the yticks.
%
%
%    'thr',r   - Keep only the largest fraction r of the coefficients, and
%                set the rest to zero.
%
%    'frange',[flow,fhigh] - Choose a frequency scale ranging from flow to
%                fhigh, values are eneterd in Hz. Default is to display from
%                0 to 8000 Hz.
%
%    'xres',xres - Approximate number of pixels along x-axis / time.
%
%    'yres',yres - Approximate number of pixels along y-axis / frequency If
%                 only one of 'xres' and 'yres' is specified, the default
%                 aspect ratio will be used.
%
%    'displayratio',r - Set the default aspect ratio.
%
%    'contour' - Do a contour plot to display the spectrogram.
%          
%    'surf'    - Do a surf plot to display the spectrogram.
%
%    'mesh'    - Do a mesh plot to display the spectrogram.
%
% 
%
%   See also:  erbtofreq, dau96
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_relanoiborra2019.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: M-Stats M-Signal M-Control M-Communication Systems
%   #Author: Helia Relano Iborra (March 2019): v4.0 provided to the AMT team
%   #Author: Clara Hollomey (2021): adapted to the AMT
%   #Author: Piotr Majdak (2021): adapted to the AMT 1.0

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
  
if nargin<3
  error('Too few input arguments.');
end;

if ~isnumeric(insig) || ~isvector(insig)
  error('%s: Input must be a vector.',upper(mfilename));
end;

global CASP_CONF;

% Get the default values.
defaults=CASP_CONF.plotdefaults;

% Switches to turn specific behaviours on and off.
doadapt=1;   % Model adaptation
dorectify=1; % Use half-wave rectification.
dolog=1;     % Use a Db scale.
domask=0;    % Remove small coefficients.
doclim=0;    % Limit the colorscale, see IMAGESC for help.
dorange=0;   % Limit the colorscale.
doxres=0;
doyres=0;
domlp=0;     % Modulation low-pass
dombp=0;     % Modulation band-pass

% Parse default and optional arguments
arguments={defaults{:},varargin{:}};
startii=1;

ii=startii-1;
while ii<length(arguments)   
  
  ii=ii+1;
  
  optarg=arguments{ii};
  
  if ~ischar(optarg)
    error('Argument %i must be a character string.',ii+1);
  end;
  
  switch lower(optarg)
   case 'rectify'
    dorectify=1;
   case 'hilbert'
    dorectify=0;
   case 'adapt'
    doadapt=1;
    dolog=0;
    dorange=0;
   case 'no_adt'
    doadapt=0;
    dolog=1;
    dorange=1;
    range=100;
   case 'classic'
    doadapt=0;
    dorectify=0;
    domlp=0;
    dolog=1;
    dorange=1;
    range=100;
   case 'contour'
    plottype='contour';
   case 'surf'
    plottype='surf';
   case 'mesh'
    plottype='mesh'; 
   case 'lin'
    dolog=0;
   case 'thr'
    domask=1;
   case 'image'
    plottype='image';
   case 'thr'
    domask=1;
    mask_val=arguments{ii+1};
    ii=ii+1;
   case 'clim'
    doclim=1;
    clim=arguments{ii+1};
    ii=ii+1;
   case 'dynrange'
    dorange=1;
    range=arguments{ii+1};
    ii=ii+1;
   case 'fullrange'
    dorange=0;
   case 'frange'
    frange=arguments{ii+1};
    ii=ii+1;
   case 'xres'
    xres=arguments{ii+1};
    ii=ii+1;
    doxres=1;
   case 'yres'
    yres=arguments{ii+1};
    ii=ii+1;
    doyres=1;
   case 'displayratio'
    displayratio=arguments{ii+1};
    ii=ii+1;
   case 'ytick'
    ytick=arguments{ii+1};
    ii=ii+1;
   case 'mlp'
    mlp=arguments{ii+1};
    ii=ii+1;
    domlp=1;              
   case 'no_mf'
    domlp=0;
    dombp=0;
   otherwise
    error([optarg,' : Unknown optional argument 1']);
  end;  
end;

if ~doxres
  xres=floor(yres/displayratio);
end;

if ~doyres
  yres=ceil(xres*displayratio);
end;

siglen=length(insig);

fhigh=frange(2);
flow =frange(1);

audlimits=freqtoaud('erb',frange);

% fhigh can at most be the Nyquest frequency
% fhigh=min(fhigh,fs/2);
% 
% % Downsample this signal if it is sampled at a much higher rate than
% % 2*fhigh. This reduces memory consumption etc. 1.5 and 1.2 are choosen as a
% % safeguard to not loose information.
% if fs>2*1.5*fhigh
%   
%   fsnew=round(fhigh*2*1.2);
%   siglennew=round(siglen/fs*fsnew);
%   
%   % Do the resampling using an FFT based method, as this is more flexible
%   % than the 'resample' method included in Matlab
%   %
%   % This is a reimplementation of the fftresample function from LTFAT.
%   
%   ffs=fft(insig);
%   if mod(siglennew,2)==0
% 
%     % even. Use average of endpoints.
%     ffs=[ffs(1:siglennew/2,:);
%        (ffs(siglennew/2+1,:)+ffs(siglen-siglennew/2+1,:))/2;
%        ffs(siglen-siglennew/2+2:siglen,:)];
%     
%   else
%     
%     % No problem, just cut.
%     ffs=[ffs(1:(siglennew+1)/2,:);
%            ffs(siglen-(siglennew-1)/2+1:siglen,:)];
%     
%   end;     
%   
%   insig=ifft(ffs)/siglen*siglennew;
% 
%   % Switch to new values
%   siglen=siglennew;
%   fs=fsnew;
%   
% end;

% Determine the hopsize
% Using a hopsize different from 1 is currently not possible because all
% the subsequent filters fail because of a to low subband sampling rate.
%hopsize=max(1,floor(siglen/xres));

hopsize=1;
outsig = insig;
% % find the center frequencies used in the filterbank
% fc = erbspace(flow,fhigh,yres);
% 
% % Calculate filter coefficients for the gammatone filter bank.
% [gt_b, gt_a]=gammatone(fc, fs);
% 
% % Apply the Gammatone filterbank
% outsig = filterbank(gt_b,gt_a,insig,hopsize);
% 
% % The subband are now (possibly) sampled at a lower frequency than the
% % original signal.
% fssubband=round(fs/hopsize);
% 
% if dorectify && (fssubband>2000)
%   outsig = envextract(2*real(outsig),fssubband);
% else
%   outsig = abs(outsig);
% end;
% 
% if doadapt
%   % non-linear adaptation loops
%   outsig = adaptloop(outsig, fssubband,10);
% end;
%   
% if domlp
%   % Calculate filter coefficients for the 50 Hz modulation lowpass filter.
%   mlp_a = exp(-mlp/fssubband);
%   mlp_b = 1 - mlp_a;
%   mlp_a = [1, -mlp_a];
%   
%   % Apply the low-pass modulation filter.
%   outsig = filter(mlp_b,mlp_a,outsig);
% end;
%   
% if domask
%   % keep only the largest coefficients.
%   outsig=largestr(outsig,mask_val);
% end
  
% Apply transformation to coefficients.
% if dolog
%   % This is a safety measure to avoid log of negative numbers if the
%   % users chooses an incorrect combination of options (e.g. 'adapt' and 'log').
%   outsig(:)=max(outsig(:),realmin);
% 
%   outsig=20*log10(outsig);
% end;

% % 'dynrange' parameter is handled by threshholding the coefficients.
% if dorange
%   maxclim=max(outsig(:));
%   outsig(outsig<maxclim-range)=maxclim-range;
% end;

% Set the range for plotting
xr=(0:hopsize:siglen-1)/fs;
yr=linspace(audlimits(1),audlimits(2),length(fc));

% Determine the labels and position for the y-label.
ytickpos=freqtoerb(ytick);

% Flip the output correctly. Each column is a subband signal, and should
% be display as the rows.
outsig=outsig.';

switch(plottype)
  case 'image'
    if doclim
      imagesc(xr,yr,outsig,clim);
    else
      imagesc(xr,yr,outsig);
    end;
  case 'contour'
    contour(xr,yr,outsig);
  case 'surf'
    surf(xr,yr,outsig);
end;
set(gca,'YTick',ytickpos);
% Use num2str here explicitly for Octave compatibility.
set(gca,'YTickLabel',num2str(ytick(:)));

axis('xy');
xlabel('Time (s)')
ylabel('Frequency (Hz)')

if nargout>0
  varargout={outsig,fc};
end;


