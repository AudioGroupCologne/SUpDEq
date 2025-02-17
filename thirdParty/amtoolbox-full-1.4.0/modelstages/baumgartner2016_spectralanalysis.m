function varargout = baumgartner2016_spectralanalysis(sig,spl,varargin)
%BAUMGARTNER2016_SPECTRALANALYSIS Spectral analysis
%   Usage:     [mp,fc] = baumgartner2016_spectralanalysis(sig,spl)
%
%   Input parameters:
%     sig     : incoming time-domain signal
%     spl     : sound pressure level (re 20e-6 Pa) in dB. Set to NAN in
%               order to bypass internal level adjustment.
%     type    : flag for target (default) or template
%     name    : identifying string for caching (e.g., 'NH12_baseline') 
%
%   Output parameters:
%     mp      : spectral magintude profile. Dimensions (4-6 optional): 
%               1) frequency, 2) position (polar angle), 3) channel, 4)
%               fiber type, 5) time frame.
%     fc      : center frequencies of auditory filters
%
%   BAUMGARTNER2016_SPECTRALANALYSIS(...) computes temporally integrated
%   spectral magnitude profiles.
%
%   BAUMGARTNER2016_SPECTRALANALYSIS accepts the following optional parameters:
%
%     'tiwin'     Set temporal integration window in seconds. Default is Inf.
%     'ID'        Listener's ID (important for caching).
%     'Condition' Label of experimental condition (also for caching).
%
%   BAUMGARTNER2016_SPECTRALANALYSIS accepts these optional flags:
%
%     'target'    Processing of a target sound (for caching). This is the
%                 default.
%     'template'  Processing of a template sound (for caching).
%     'gammatone' To apply a Gammatone filterbank instead of the Zilany et
%                 al. (2014) model of the auditory periphery.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2016_spectralanalysis.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA M-Signal M-Stats O-Statistics
%   #Author: Robert Baumgartner (2016), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


definput.import={'baumgartner2016'}; % defaults from arg_baumgartner2016
definput.flags.type={'target','template'};

[flags,kv]=ltfatarghelper({},definput,varargin);

% Set cachename
if flags.do_target
  cachenameprefix = 'ireptar_';
  if size(sig,1) > kv.tiwin*kv.fs; cachenameprefix = [cachenameprefix 'tiwin' num2str(kv.tiwin*1e3) 'ms_']; end
  if not(isempty(kv.ID))
    cachenameprefix = [cachenameprefix kv.ID '_'];
  end
  if not(isempty(kv.Condition))
    cachenameprefix = [cachenameprefix kv.Condition '_'];
  end
  %%% equalized conditions (don't run periphery model again on it)
  cachenameprefix = strrep(cachenameprefix,'BB','baseline');
  cachenameprefix = strrep(cachenameprefix,'CL','baseline');
  %%%
else
  cachenameprefix = 'ireptem_';
  if not(isempty(kv.ID))
    cachenameprefix = [cachenameprefix kv.ID '_'];
  end
end
cachenameprefix = [cachenameprefix 'lat' num2str(kv.lat) '_' num2str(spl) 'dB'];

%% Remove pausings (at beginning and end)
Nch = size(sig,3);
if not(isnan(spl))
  idnz = mean(sig(:,:).^2,2) ~=0 | [diff(mean(sig(:,:).^2,2)) ~= 0;1];
  sig = sig(idnz,:,:); 
  % and evaluate broadband ILD
  if Nch == 2
    ILD = dbspl(sig(:,:,2))./dbspl(sig(:,:,1));
    if mean(ILD) < 1
      chdamp = 2;
    else
      ILD = 1./ILD;
      chdamp = 1;
    end
  end
end

%% Gammatone

if flags.do_gammatone
  
  cachename = [cachenameprefix '_gammatone_' num2str(1/kv.space,'%u') 'bpERB'];
  if flags.do_middleear; cachename = [cachename '_middleear']; end
  if flags.do_ihc; cachename = [cachename '_ihc']; end
  [mp,fc] = amt_cache('get',cachename,flags.cachemode);
  if isempty(mp)
    
    % Set level by single factor, not individually for every direction or
    % channel (directions and channels are pooled)
    if not(isnan(spl))
      sigsize = size(sig);
      sig = scaletodbspl(sig(:),spl, 100);
      sig = reshape(sig,sigsize);
    end

    if flags.do_middleear
        miearfilt = middleearfilter(kv.fs);
        sig = lconv(sig,miearfilt(:));
    end

    Lpad = round(5/kv.flow*kv.fs); % add 5 cycles of lowest freq
    sig = postpad(sig,size(sig,1)+Lpad);
    if kv.space == 1
      [mp,fc] = auditoryfilterbank(sig(:,:),kv.fs,...
          'flow',kv.flow,'fhigh',kv.fhigh);
    else
      fc = audspacebw(kv.flow,kv.fhigh,kv.space,'erb');
      [bgt,agt] = gammatone(fc,kv.fs,'complex');
      mp = 2*real(ufilterbankz(bgt,agt,sig(:,:)));  % channel (3rd) dimension resolved!
    end
    Nfc = length(fc);   % # bands

    % IHC transduction
    if flags.do_ihc
      mp = ihcenvelope(mp,kv.fs,'ihc_dau1996');
    end

%     % Set back the channel dimension
%     mp = reshape(mp,[size(mp,1),Nfc,size(sig,2),size(sig,3)]);

    % Averaging over time (RMS)
    if flags.do_target && size(mp,1) > kv.tiwin*kv.fs
      Lframe = round(kv.tiwin*kv.fs); % length of each frame
      Nframes = ceil(size(mp,1)/Lframe); % # frames
      mp = postpad(mp,Lframe*Nframes,0,1);
      mp = reshape(mp,[Lframe,Nframes,Nfc,size(sig,2),Nch]);
      mp = permute(mp,[1,3,4,5,6,2]); % frames as last dimension
      mp = shiftdim(rms(mp));
      time = (0:Lframe:Lframe*Nframes-1)/kv.fs;
    else % integrate over whole duration
      mp = reshape(mp,[size(mp,1),Nfc,size(sig,2),Nch]); % Set back the channel dimension
      mp = shiftdim(rms(mp));
      time = 0;
    end
    
    % Logarithmic transformation (dB) 
    mp = 100 + 20*log10(mp);
    
    amt_cache('set',cachename,mp,fc);
  end
  
  % Limit dynamic range
  if not(isnan(spl))
    mp = min(mp,kv.GT_maxSPL); % maybe +30 because dynamic range was evaluated with broadband noise (40 auditory bands)
    mp = max(mp,kv.GT_minSPL);
  end
  
end


%% Zilany Model

if flags.do_zilany2007 || flags.do_zilany2014
  
  ftd = [0.16,0.23,0.61]; % Liberman (1978)
  ftweights = ftd(kv.fiberTypes)/sum(ftd(kv.fiberTypes));
  for tt = 1:length(kv.fiberTypes)

    ft = kv.fiberTypes(tt);
    
    cachename = [cachenameprefix '_fiberType' num2str(ft)];
    if kv.cihc < 1; cachename = [cachename '_cihc' num2str(kv.cihc)]; end
    if kv.cohc < 1; cachename = [cachename '_cohc' num2str(kv.cohc)]; end
    
    try
      [mp,fc,time] = amt_cache('get',cachename,flags.cachemode);
    catch
      [mp,fc] = amt_cache('get',cachename,flags.cachemode);
    	time = 0:kv.tiwin:(size(mp,5)-1)*kv.tiwin; 
    end
    
    if isempty(mp)

      Nmin = .05*kv.fsmod; % decay time at 700-Hz in response to click
      if length(sig)/kv.fs*kv.fsmod < Nmin
          sig = postpad(sig,ceil(Nmin*kv.fs/kv.fsmod),0,1);
      end
        
      amt_disp(['Compute: ' cachename]);
      Ntar = size(sig,2); % # target angles
      len = ceil(length(sig)/kv.fs*kv.fsmod);
      ANresp = zeros(len,kv.numCF,Ntar,2);
      for ch = 1:Nch
        for ii = 1:Ntar
          
          if Nch > 1 && ch == chdamp;
            spl_mod = ILD(ii)*spl;
          else
            spl_mod = spl;
          end
          
          if flags.do_zilany2007
            [ANout,fc] = zilany2007(spl_mod,sig(:,ii,ch),kv.fs,...
              kv.fsmod,'flow',kv.flow,'fhigh',kv.fhigh,'numCF',kv.numCF);
            ANout = ANout';
          else % zilany2014
            fc = audspace(kv.flow,kv.fhigh,kv.numCF,'erb');
            ANout = zilany2014(scaletodbspl(sig(:,ii,ch),spl_mod),kv.fs,...
              fc,'fiberType',ft,'cohc',kv.cohc,'cihc',kv.cihc);
            ANout=ANout';
          end

          % Compensate for cochlear delay?!?
          
          % Check stimulus onset
          ionset = find(diff(mean(ANout.^2,1)) ~= 0,1,'first');

          ANresp(:,:,ii,ch) = ANout(:,(1:len)+ionset-1)';
          
          amt_disp([num2str(ii+(ch-1)*Ntar) ' of ' num2str(Ntar*Nch) ' done'],'volatile');
        end
        
      end
	  amt_disp();

      % Averaging over time (RMS)
      if flags.do_target && size(ANresp,1) > kv.tiwin*kv.fsmod
        Lframe = kv.tiwin*kv.fsmod; % length of each frame
        Nframes = ceil(size(ANresp,1)/Lframe); % # frames
        ANresp = postpad(ANresp,Lframe*Nframes,0,1);
        ANresp = reshape(ANresp,[Lframe,Nframes,length(fc),size(sig,2),size(sig,3)]);
        ANresp = permute(ANresp,[1,3,4,5,6,2]); % frames as last dimension
        time = (0:Lframe:Lframe*Nframes-1)/kv.fsmod;
      else % integrate over whole duration
%         ANresp = reshape(ANresp,[size(len,length(fc),size(sig,2),size(sig,3)]); % retreive polar dimension if squeezed out
        time = 0;
      end
      mp = shiftdim(mean(ANresp),1);

      amt_cache('set',cachename,mp,fc,time);
    end
    
%     if size(mp,2) ~= size(sig,2) % retreive polar dimension if squeezed out
%         mp = reshape(mp,[size(mp,1),size(sig,2),size(sig,3)]);
%     end
    if tt == 1 % init
      mp_cum = zeros(kv.numCF,size(sig,2),size(sig,3),1,length(time));
      mp_sep = zeros(kv.numCF,size(sig,2),size(sig,3),length(kv.fiberTypes),length(time));
    end
    mp_cum = mp_cum + ftweights(tt)*mp;
    mp_sep(:,:,:,tt,:) = mp;
    
  end
  
  if flags.do_ftcum
    mp = mp_cum;
  else % flags.do_ftopt
    mp = mp_sep;
  end
  
  % adjust frequency range for cached data 
  idf = fc <= kv.fhigh & fc >= kv.flow;
  fc = fc(idf);
  mp = mp(idf,:,:,:,:);

end

% fiber activity gating
% mp = round(mp*0.3);

varargout{1} = mp;
if nargout > 1
  varargout{2} = fc;
  if nargout > 2
    varargout{2} = time;
  end
end

end


