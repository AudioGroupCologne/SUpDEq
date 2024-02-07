function out = sig_baumgartner2017looming(varargin)
%SIG_BAUMGARTNER2017looming flattens magnitude spectra of HRTFs
%
%   Usage:  Obj = sig_baumgartner2017looming(Obj,C,flow,fhigh)
%           stim = sig_baumgartner2017looming(exp)
%
%   Input parameters:
%     Obj   : reference SOFA object
%     C     : spectral contrast. 1 refers to reference (default), 
%             0 to flat, -1 to flipped spectral magnitude
%     flow  : lower cut-off frequency in Hz. Default is 1 kHz.
%     fhigh : higher cut-off frequency in Hz. Default is 18 kHz.
%     Nfft  : FFT length. Default is 1024.
%     Ntrans: number of frequency bins used to softly transition the 
%             magnitude spectrum. Default is 10.
%     delay : additional delay in ms to allow temporal dispersion in
%             obtained impulse responses. Default is 0.
%
%   Output parameters:
%     Obj   : modified SOFA object
%     stim  : stimulus structure containing 
%             - subject ID (ID)
%             - sampling rate (fs) 
%             - impulse responses (IR referring to contrasts C_IR)
%             - selected azimuths (azi), and
%             - stimuli of the selected experiment (cont for Exp. I; cont and discont for Exp. II; both referring to contrasts specified in C_pair).
%
%   The exp flag may be used to get stimuli and impulse responses:
%     'exp1'    Stimuli from Exp. I.
%     'exp2'    Stimuli from Exp. II.
%
%   Optional flags:
%     'soft'    soft magnitude transition option for use with broadband stimuli
%
%   SIG_BAUMGARTNER2017LOOMING either generates spectrally flattened HRTF
%   representations as used by Baumgartner et al. (2017) in order to induce
%   various degrees of sound externalization, or generates the signals
%   presented to the listeners in one of the two experiments from
%   Baumgartner et al. (2017) as specified by the exp flag.
%
%   References:
%     R. Baumgartner, D. K. Reed, B. Tóth, V. Best, P. Majdak, H. S. Colburn,
%     and B. Shinn-Cunningham. Asymmetries in behavioral and neural responses
%     to spectral cues demonstrate the generality of auditory looming bias.
%     Proceedings of the National Academy of Sciences, 2017. [1]http ]
%     
%     References
%     
%     1. http://www.pnas.org/content/early/2017/08/16/1703247114.abstract
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/signals/sig_baumgartner2017looming.php


% #Author: Robert Baumgartner (2020), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


definput.keyvals.Obj = [];
definput.keyvals.C = 1;
definput.keyvals.flow = 1e3;
definput.keyvals.fhigh = 16e3;
definput.keyvals.Nfft = 2^10;
definput.keyvals.Ntrans = 5; % number of frequency bins used to softly transition the magnitude spectrum
definput.keyvals.delay = 0; % in ms
definput.flags.experiment = {'hrtf','exp1','exp2'};
definput.flags.transition = {'hard','soft'};
definput.flags.source = {'noise','impulse'};

[flags,kv]=ltfatarghelper({'Obj','C','flow','fhigh'},definput,varargin);

%% Spectral contrast manipulation
if flags.do_hrtf
  if isempty(kv.Obj)
    error('Missing reference HRTF.')
  end
  
  out = kv.Obj;
  fs = kv.Obj.Data.SamplingRate;
  ref = shiftdim(kv.Obj.Data.IR,2);

  hrtf = fftreal(ref,kv.Nfft);
  freq = 0:fs/kv.Nfft:fs/2; % frequency vector

  % following processing only done for dedicated frequency range
  idf = freq >= kv.flow & freq <= kv.fhigh; % indices for dedicated frequency range
  Nf = sum(idf); % # frequency bins
  mag = db(abs(hrtf(idf,:,:))); % HRTF magnitudes in dB
  idwf = idf(:) | circshift(idf(:),[1,0]); % include one neighbouring position for evaluation of frequency weighting
  wf = diff(freqtoerb(freq(idwf))); % frequency weighting according to differentiated ERB scale
  wf = repmat(wf(:)/sum(wf),[1,size(mag,2),size(mag,3)]);
  meanmag = repmat(sum(wf.*mag,1),[Nf,1,1]);
  varmag = mag - meanmag;

  ph = angle(hrtf);
  modmag = db(abs(hrtf));
  if flags.do_soft
    ramp = sin(linspace(0,pi/2,kv.Ntrans)).^2;
    fader = [ramp,ones(1,Nf-2*kv.Ntrans),1-ramp]';
    fader = repmat(fader,1,size(modmag,2), size(modmag,3));
    modmag(idf,:,:) = (1-fader).*modmag(idf,:,:) + fader.*(meanmag + kv.C*varmag);
  else
    modmag(idf,:,:) = meanmag + kv.C*varmag;
  end
  mod = ifftreal(10.^(modmag/20).*exp(1i*ph),kv.Nfft);

  mod = circshift(mod,kv.delay/1000*fs,1);
  
  out.Data.IR = shiftdim(mod,1);
end

%% Stimuli from Exp. I % II
if flags.do_exp1 || flags.do_exp2
%   ltfatsetdefaults('dbspl','dboffset',100); % PM: not needed. 
  % Individual data
  data = data_baumgartner2017looming(flags.experiment); % for ID and azimuth
  rawData = struct2cell(data.rawData);
  subj = rawData(1,:);
  azi = rawData(3,:);
  hrtf = data_baumgartner2017looming('hrtf');
  
  % Source signal
  fadeDur = 0.05;
  fs = hrtf(1).Obj.Data.SamplingRate;
  long = noise(fs,1,'white'); % long stimulus for continuous trials
  short = noise(0.5*fs,1,'white'); % short stimulus for discontinuous trials
  
  % Fade in/out
  sinRamp = sin(pi/2*(0:fadeDur*fs-1)/(fadeDur*fs)).^2;
  long = long.*[sinRamp,ones(1,length(long)-2*length(sinRamp)),fliplr(sinRamp)]';
  short = short.*[sinRamp,ones(1,length(short)-2*length(sinRamp)),fliplr(sinRamp)]';
  
  % Band-pass filter
  filtOrder = 4;
  [b_bp,a_bp] = butter(filtOrder,[kv.flow,kv.fhigh]/(fs/2));
  
  % Generate individual stimuli
  C = 0:.5:1; % spectral contrasts
  stim = struct;
  for ss = 1:length(subj)
    
    stim(ss).ID = subj{ss};
    stim(ss).azi = azi{ss};
    
    % Spectral flattening
    stim(ss).fs = fs;
    stim(ss).C_IR = C;
    for ii = 1:length(C)
      Obj = sig_baumgartner2017looming(hrtf(ss).Obj,C(ii),kv.flow,kv.fhigh);
      idpos = Obj.SourcePosition(:,1) == azi{ss} & Obj.SourcePosition(:,2) == 0;
      stim(ss).IR{ii} = squeeze(shiftdim(Obj.Data.IR(idpos,:,:),2));
      stim(ss).short{ii} = SOFAspat(short,Obj,azi{ss},0);
      stim(ss).long{ii} = SOFAspat(long,Obj,azi{ss},0);
    end

    % Band-pass filtering
    for ii = 1:length(C)
      stim(ss).IR{ii} = filter(b_bp,a_bp,stim(ss).IR{ii});
      stim(ss).short{ii} = filter(b_bp,a_bp,stim(ss).short{ii});
      stim(ss).long{ii} = filter(b_bp,a_bp,stim(ss).long{ii});
    end

    % Set level (adjust all stimuli by the same factor)
    targetSPL = 70;
    while any(max(stim(ss).long{ii}(:)) > 1) % avoid clipping
      currentSPL = nan(length(C),2);
      for ii = 1:length(C)
        currentSPL(ii,:) = dbspl(stim(ss).long{ii});
      end
      adjustmentFactor = db2mag(targetSPL-mean(currentSPL(:)));
      for ii = 1:length(C)
        stim(ss).short{ii} = adjustmentFactor*stim(ss).short{ii};
        stim(ss).long{ii} = adjustmentFactor*stim(ss).long{ii};
      end
      targetSPL = targetSPL - 3;
    end
    
  end
  
  % Define combinations
  NC = length(C);
  if flags.do_exp1
    iCpair = unique(nchoosek([1:NC,NC:-1:1], 2),'rows');
  elseif flags.do_exp2
    iCpair = [nchoosek(1:NC,2);nchoosek(NC:-1:1,2)];
  end
  Cpair = C(iCpair);

  % Timings in sec
  isi = 0.1; % inter-stimulus interval
  dur = 1.2; % overall stimulus duration
  xfade = 0.6; % timing of cross-fade
  xfadeDur = 0.01; % duration of cross-fade 

  for ss = 1:length(subj)
    stim(ss).C_pair = Cpair;
    for ip=1:length(Cpair)
      stim(ss).cont{ip} = SpExCue_crossfade(...
        stim(ss).long{iCpair(ip,1)},stim(ss).long{iCpair(ip,2)},...
        fs,dur,xfade,xfadeDur,true);
      stim(ss).discont{ip} = [stim(ss).long{iCpair(ip,1)};zeros(isi*fs,2);stim(ss).long{iCpair(ip,2)}];
    end
  end
    
  out = rmfield(stim,{'long','short'});
  
  if flags.do_exp1
    out = rmfield(out,'discont');
  end
%   ltfatsetdefaults('dbspl','dboffset',93.98); % PM: not needed.
end

end

function varargout = SpExCue_crossfade(sig1,sig2,fs,dur,tcross,durfade,fadeIOflag)
%SpExCue_crossfade creates cross-faded stimulus pair (cos^2 fade)
%   Usage: [sigpair,nM2] = SpExCue_crossfade(sig1,sig2,fs,dur,tcross,durfade,fadeIOflag)
%
%   Input parameters:
%     sig1    : first stimulus (time in first dimension)
%     sig2    : second stimulus (time in first dimension)
%     fs      : sampling rate of signals
%     dur     : overall duration of stimulus pair in sec
%     tcross  : center time of cross-fade in sec
%     durfade : duration of fades (in and out) in sec
%     fadeIOflag : flag to also fade paired signal in and out
%
%   Output parameters:
%     sigpair : cross-faded stimulus pair
%     nM2     : index of center time of cross-fade

if not(exist('fadeIOflag','var')) 
  fadeIOflag = false;
end

n1 = length(sig1); % length of first input signal
n2 = length(sig2); % length of second input signal
ntotal = round(dur*fs); % total length of output signal
ncross = round(tcross*fs); % index of cross-fade center

nfade = 2*round(durfade*fs/2)-1; % closest odd # of time indices for fade
nstop1 = ncross+(nfade-1)/2; % offset index of first input signal
nstart2 = ncross-(nfade-1)/2; % onset index of second input signal
nsig2 = ntotal-nstart2+1; % length of second stimulus part
fadein = sin(0:pi/2/(nfade-1):pi/2);
fadeout = fliplr(fadein);
  
if n1 <= ncross % short signal -> no cross-fade
  
  fadedsig1 = postpad(sig1,ntotal);
  nsig2 = ntotal-ncross+1;
  sig2 = postpad(sig2,nsig2);
  fadedsig2 = cat(1,zeros(ncross-1,size(sig2,2)),sig2); % prepad
  
else % long signal -> cross-fade
  
  fader1 = [ones(1,nstop1-nfade),fadeout];
  fadedsig1 = sig1(1:nstop1,:).*repmat(fader1(:),1,2);
  fadedsig1 = postpad(fadedsig1,ntotal);
    
  fader2 = [fadein,ones(1,nsig2-nfade)];
  fadedsig2 = sig2(n2-nsig2+1:n2,:).*repmat(fader2(:),1,2);
  fadedsig2 = cat(1,zeros(nstart2-1,2),fadedsig2); % prepad
  
end

sigpair = fadedsig1 + fadedsig2;

if fadeIOflag
  fader = [fadein,ones(1,ntotal-2*nfade),fadeout]';
  sigpair = sigpair.*repmat(fader,[1,2]);
end

varargout{1} = sigpair;
if nargout == 2
  varargout{2} = ncross;
end

end


