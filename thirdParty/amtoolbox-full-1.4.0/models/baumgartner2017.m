function [E,varargout] = baumgartner2017( target,template,varargin )
%BAUMGARTNER2017 Sound externalization model based on monaural spectral cues
%   Usage:    [E,cues,cueLabels] = baumgartner2017( target,template )
%
%   Input parameters:
%
%     target  : binaural impulse response(s) referring to the directional 
%               transfer function(s) (DFTs) of the target sound(s).
%               Option 1: given in SOFA format -> sagittal plane DTFs will 
%               be extracted internally. 
%               Option 2: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane formatted 
%               according to the following matrix dimensions: 
%               time x direction x channel/ear
%
%     template: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane referring to
%               the perceived lateral angle of the target sound.
%               Options 1 & 2 equivalent to target.
%
%   Output parameters:
%
%     E         : predicted degree of externalization
%
%     cues      : outcomes of individual cues 
%
%     cueLabels : cue labels; cell array with 1st col. denoting acronyms
%                 and 2nd col. for descriptions
%
%   BAUMGARTNER2017(...) is a model for sound externalization.
%   It bases on the comparison of the intra-aural internal representation
%   of the incoming sound with a template and results in a probabilistic
%   prediction of polar angle response.
%
%   BAUMGARTNER2017 accepts the following optional parameters:
%
%     'cueWeights',cW  Set the weights of individual cues to determine the 
%                      final externalization score. 
%
%                      Cue-specific weights (entered as a vector) are ordered as follows: 
%
%                      1 monaural spectral similarity (c.f., Baumgartner et al., 2014). This is the default.
%
%                      2 interaural spectral similarity of ILDs (c.f., Hassager et al., 2016)
%
%                      3 spectral standard deviation of ILDs (c.f., Georganti et al., 2013)
%
%                      4 temporal standard deviation of ILDs (c.f., Catic et al., 2015)
%
%                      5 interaural coherence (c.f., Hassager et al., 2017)
%
%                      6 interaural broadband time-intensity coherence
%
%                      7 difference in sound pressure level
%
%     'fs',fs        Define the sampling rate of the impulse responses. 
%                    Default value is 48000 Hz.
%
%     'S',S          Set the listener-specific sensitivity threshold 
%                    (threshold of the sigmoid link function representing 
%                    the psychometric link between transformation from the
%                    distance metric and similarity index) to S. 
%                    Default value is 1.
%
%     'lat',lat      Set the apparent lateral angle of the target sound to
%                    lat. Default value is 0 degree (median SP).
%
%     'stim',stim    Define the stimulus (source signal without directional
%                    features). As default an impulse is used.
%
%     'fsstim',fss   Define the sampling rate of the stimulus. 
%                    Default value is 48000 Hz.
%
%     'flow',flow    Set the lowest frequency in the filterbank to
%                    flow. Default value is 700 Hz.
%
%     'fhigh',fhigh  Set the highest frequency in the filterbank to
%                    fhigh. Default value is 18000 Hz.
%
%     'space',sp     Set spacing of auditory filter bands (i.e., distance 
%                    between neighbouring bands) to sp in number of
%                    equivalent rectangular bandwidths (ERBs). 
%                    Default value is 1 ERB.
%
%     'do',do        Set the differential order of the spectral gradient 
%                    extraction to do. Default value is 1 and includes  
%                    restriction to positive gradients inspired by cat DCN
%                    functionality.
%
%     'bwcoef',bwc   Set the binaural weighting coefficient bwc.
%                    Default value is 13 degrees.
%
%     'range',c1     Set the range factor of the externalization scores to c1.
%                    Default value is 3.78 from Hassager et al. (2016).
%
%     'offset',c2    Set the offset of the externalization score to c2.
%                    Default value is 1 from Hassager et al. (2016).
%
%     'ILD_JND',L    Set the just noticeable ILD difference to L from the 
%                    internal template. Default value is 1 (dB).
%
%     'ITD_JND',T    Set the just noticeable ITD difference to T from the 
%                    internal template. Default value is 20e-6 (s).
%
%   Requirements:
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
%
%   2) Circular Statistics Toolbox from https://de.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
%
%   References:
%     R. Baumgartner, P. Majdak, H. Colburn, and B. Shinn-Cunningham.
%     Modeling sound externalization based on listener-specific spectral
%     cues. In Acoustics ‘17 Boston: The 3rd Joint Meeting of the Acoustical
%     Society of America and the European Acoustics Association, Boston, MA,
%     Jun 2017.
%     
%
%   See also: exp_baumgartner2017 baumgartner2016_spectralanalysis,
%   baumgartner2016_gradientextraction, baumgartner2014_binauralweighting,
%   baumgartner2021 baumgartner2017_iacc itdestimator data_baumgartner2017
%   sig_li2020
%   demo_baumgartner2017
%   demo_baumgartner2021
%   exp_steidle2019
%   baumgartner2014
%   dietz2011
%   dau1996
%   baumgartner2021
%   baumgartner2016
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/baumgartner2017.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: SOFA M-Stats M-Signal O-Statistics
%   #Author: Robert Baumgartner (2017), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Notification
amt_disp('Note that baumgartner2017 is a preliminary version of baumgartner2021.');
amt_disp('We recommend to use baumgartner2021.');

%% Check input

definput.import={'baumgartner2017','baumgartner2016','baumgartner2014','baumgartner2014_pmv2ppp','localizationerror','amt_cache'};

[flags,kv]=ltfatarghelper(...
  {'fs','S','lat','stim','space','do','flow','fhigh',... %'fsstim'
  'bwcoef','polsamp','rangsamp','mrsmsp','gamma'},definput,varargin);

flags.do_plot = false;

if not(isstruct(target)) && ismatrix(target)
  target = permute(target,[1,3,2]);
%   warning(['Matrix dimensions of target should be: time x direction x channel/ear.' ...
%    'Since 3rd dimension was empty, 2nd dimension was used as channel dimension.'])
end

if not(isstruct(template)) && ismatrix(template)
  template = permute(template,[1,3,2]);
%   warning(['Matrix dimensions of template should be: time x direction x channel/ear.' ... 
%     'Since 3rd dimension was empty, 2nd dimension was used as channel dimension.'])
end

%% Print Settings

if flags.do_print 
  if flags.do_nomrs
    kv.mrsmsp = 0;
  end
  amt_disp(['Settings: PSGE = ' num2str(kv.do,'%1.0f') '; Gamma = ' ...
    num2str(kv.gamma,'%1.0u') '; Epsilon = ' num2str(kv.mrsmsp,'%1.0f') ' deg']);
end


%% Determine lateral angle and extract HRTFs of sagittal plane
 
if isstruct(target) % Targets given in SOFA format
  kv.fs = target.Data.SamplingRate;
  [target,tang] = extractsp( kv.lat,target );
% else
%   fncache = ['latLookup_',template.GLOBAL_ListenerShortName];
%   latLookup = amt_cache('get',fncache,flags.cachemode);
%   if isempty(latLookup)
%     latLookup = itd2angle_lookuptable(template,template.Data.SamplingRate,'dietz2011');
%     amt_cache('set',fncache,latLookup)
%   end
%   tarSig = squeeze(target);
%   kv.lat = wierstorf2013_estimateazimuth(tarSig,latLookup,'fs',kv.fs,'dietz2011','rms_weighting');
%   disp(kv.lat)
end

if isstruct(template) % Template given in SOFA format
  [template,rang] = extractsp( kv.lat,template );
end

% Error handling
% if size(template,2) ~= length(rang)
%   fprintf('\n Error: Second dimension of template and length of polsamp need to be of the same size! \n')
%   return
% end

%% Optional: Middle ear filtering
if flags.do_middleEarFilter
  b=middleearfilter(kv.fs);
  target = filter(b,1,target);
  template = filter(b,1,template);
end

%% Optional: HRTF filtering

dimtar = size(target); % for lconv dim check

if not(isempty(kv.stim))
  target = lconv(target,kv.stim);
end

% check that lconv preserved matrix dimensions (earlier bug in lconv)
if size(target,2) ~= dimtar(2)
  target = reshape(target,[size(target,1),dimtar(2:end)]);
end

% frameLength = round((kv.tempWin*kv.fs));


%% Level difference
SPL = dbspl(target) - dbspl(template);
SPL = max(SPL-kv.ILD_JND,0);
SPL = mean(SPL,3);

%% ITD & IC
[tem.itd,~,tem.iacc] = itdestimator(shiftdim(template,1),'fs',kv.fs,'MaxIACCe');
[tar.itd,~,tar.iacc] = itdestimator(shiftdim(target,1),'fs',kv.fs,'MaxIACCe');
% IACC = tar.iacc/tem.iacc - 1;

tem.ic = baumgartner2017_iacc(squeeze(template),'argimport',flags,kv);
tar.ic = baumgartner2017_iacc(squeeze(target),'argimport',flags,kv);
IC = tar.ic - tem.ic;

%% Filterbank
% [tem.mp,fc] = baumgartner2016_spectralanalysis(template,70,'argimport',flags,kv,'tiwin',0.005,'gammatone','redo');
% tar.mp = baumgartner2016_spectralanalysis(target,70,'argimport',flags,kv,'tiwin',0.005,'gammatone','redo');
[tem.mp,fc] = baumgartner2016_spectralanalysis(template,70,'argimport',flags,kv,'tiwin',kv.tempWin,'gammatone','redo');
tar.mp = baumgartner2016_spectralanalysis(target,70,'argimport',flags,kv,'tiwin',kv.tempWin,'gammatone','redo');
if isscalar(kv.reflectionOnsetTime) % evaluate only direct path (DP)
  idDP = round(kv.reflectionOnsetTime*kv.fs);
  N1ms = round(1e-3*kv.fs); % 1 ms fade out
  taper = [ones(1,idDP-N1ms) , 0.5*(1+cos(linspace(0,pi,N1ms)))];
  temTaper = repmat([taper(:);zeros(size(template,1)-idDP,1)],...
                    [1,size(template,2),size(template,3)]);
  tarTaper = repmat([taper(:);zeros(size(target,1)-idDP,1)],...
                    [1,size(target,2),size(target,3)]);
  [temDP.mp,fc] = baumgartner2016_spectralanalysis(temTaper.*template,70,...
    'argimport',flags,kv,'tiwin',kv.tempWin,'gammatone','redo');
  tarDP.mp = baumgartner2016_spectralanalysis(tarTaper.*target,70,...
    'argimport',flags,kv,'tiwin',kv.tempWin,'gammatone','redo');
end

%% interaural temporal SD of ILDs (Catic et al., 2015)
MinNumTimeFrames = 20;
if size(tem.mp,5) >= MinNumTimeFrames && size(tar.mp,5) >= MinNumTimeFrames
  tem.STild = -diff(tem.mp,1,3); % short-term ILDs
  tar.STild = -diff(tar.mp,1,3);
  ITSD = 1 - mean(std(tar.STild,0,5)./std(tem.STild,0,5));
else
  ITSD = nan;
end

%% Temporal "average" -> max is most appropriate because of RMS and frequency-dependent excitation delay
% if kv.tempWin >= 1
%   tem.mp = max(tem.mp,[],5);
%   tar.mp = max(tar.mp,[],5);
% end
if flags.do_plot
  for ear = 1:2
    subplot(1,2,ear)
    semilogx(fc,squeeze(tar.mp(:,1,ear)))
    xlabel('Frequency (Hz)')
    ylabel('RMS magnitude (dB)')
    hold on
  end
end

%% Spectral cues
if exist('temDP','var')
  tem.psg = baumgartner2016_gradientextraction(temDP.mp,fc,'mgs',1);
  temDP.ild = diff(temDP.mp,1,3);
else
  tem.psg = baumgartner2016_gradientextraction(tem.mp,fc,'mgs',1);
end
tem.ild = diff(tem.mp,1,3);

if exist('tarDP','var')
  tar.psg = baumgartner2016_gradientextraction(tarDP.mp,fc,'mgs',1);
  tarDP.ild = diff(tarDP.mp,1,3);
else
  tar.psg = baumgartner2016_gradientextraction(tar.mp,fc,'mgs',1);
end
tar.ild = diff(tar.mp,1,3);

%% spectral SD of ISLD (Georganti et al., 2013) -> dprime possible
% tar.sdild = std(tar.ild,0,1);
% tem.sdild = std(tem.ild,0,1);
% ref.sdild = std(ref.ild,0,1);
% ISLDspecSD = mean(tar.sdild./tem.sdild,5);

%% Spectral comparison

for iSC = 1:3 % first monaural then interaural
  if iSC == 1 % monaural spectral gradients
    tem.nrep = tem.psg.m;
    tar.nrep = tar.psg.m;
  elseif iSC == 2 % interaural spectral differences
    if exist('temDP','var') && exist('tarDP','var')
      tem.nrep = temDP.ild;
      tar.nrep = tarDP.ild;
    else
      tem.nrep = tem.ild;
      tar.nrep = tar.ild;
    end
  elseif iSC == 3 % spectral SD of interaural differences
    tem.nrep = std(tem.ild,0,1);
    tar.nrep = std(tar.ild,0,1);
    kv.ILD_JND = 0;
  end
  
  % comparison with time average of spectral template
  nrep = {tar.nrep};
  if flags.do_dprime
    nrep{2} = tem.nrep;
  end
  tem.nrep = mean(tem.nrep,5);
  sigma = cell(length(nrep),1);
  for inrep = 1:length(nrep)  
    tmp = repmat(tem.nrep,[1,1,1,1,size(nrep{inrep},5)]);
    if iSC < 3
      delta = abs(tmp-repmat(nrep{inrep},[1,size(tmp,2),1,1,1]));
      delta(delta < kv.ILD_JND) = 0; % limit minimum ILD difference according to JND
    else
      delta = tmp-repmat(nrep{inrep},[1,size(tmp,2),1,1,1]);
    end
    delta = delta./(eps+abs(tmp)); % normalization (Weber fraction)
    sigma{inrep} = mean(delta); % average across frequency bands
    if iSC == 1 % do_intraaural
      sigma{inrep} = baumgartner2014_binauralweighting(sigma{inrep},'argimport',flags,kv);
    end
  end

  % temporal integration
  if length(sigma{1}) == 1
    distmetric = sigma{1};%exp(-kv.S*sigma{1});
  elseif flags.do_dprime % signal detection theory applied to time histograms
  %   figure; histogram(sigma{1}); hold on ; histogram(sigma{2}); legend('target','reference')
    allsigma = [sigma{1}(:);sigma{2}(:)];
    msigma = mean(allsigma);
    sdsigma = std(allsigma);
    mzsigma(1) = mean((sigma{1}-msigma) ./ sdsigma);
    mzsigma(2) = mean((sigma{2}-msigma) ./ sdsigma);
    dprime = max(mzsigma(1)-mzsigma(2),0);
    distmetric = dprime;
  else % temporal weighting according to amount of sensory information available
%     si = exp(-kv.S*sigma{1});
    % figure; plot(squeeze(bsi))
    tweight = mean(mean(abs(tar.nrep)),3); % temporal weighting
    tweight = tweight-min(tweight,[],5); % bounded between 0
    tweight = 2*tweight./max(tweight,[],5); % and 2
    distmetric = sigma{1}.*tweight;
    distmetric = mean(distmetric,5);
  end
  
%   si = distmetric;%exp(-kv.S*distmetric);
  
  if iSC == 1
    MSS = distmetric; % monaural spectral distance
  elseif iSC == 2
    ISS = distmetric; % interaural spectral distance
%     IIC = distmetric; % for ITIC
  elseif iSC == 3
    ISSD = distmetric;
  end
  
end

%% Interaural time-intenstiy coherence (ITIC) -> dprime possible
% IIC = mean(tar.ild)/(mean(tem.ild)+eps);
% if tar.itd == tem.itd && tem.itd == 0
%   TC = 0;
% else
%   TC = (tar.itd-tem.itd)/(eps+tem.itd);
% end
% IC = mean(tar.ild(:))/(eps+mean(tem.ild(:)));
if abs(tar.itd - tem.itd) >= kv.ITD_JND
  ITR = tar.itd/(eps+tem.itd) -1;
else
  ITR = 0;
end
if any(abs(mean(tar.ild) - mean(tem.ild)) >= kv.ILD_JND)
  ILR = mean(tar.ild)./(mean(tem.ild)+eps) -1;
else
  ILR = 0;
end
ITIT = abs( ITR -  ILR );

%% Cue integration/weighting

cues = [MSS; ISS; ISSD; ITSD; IC; ITIT; SPL];
cueLbl = {'MSG','Monaural spectral gradients (c.f., Baumgartner et al., 2014)'; ...
          'ISS','Interaural spectral shape (c.f., Hassager et al., 2016)'; ...
          'ISSD','Interaural spectral standard deviation (c.f., Georganti et al., 2013)'; ...
          'ITSD','Interaural temporal standard deviation (c.f., Catic et al., 2015)'; ...
          'IC','Interaural coherence (c.f., Hassager et al., 2017)'; ...
          'ITIT','Interaural time-intensity trading (ITD vs. ILD)'; ...
          'SPL','Level difference (target - reference)'; ...
         };

if flags.do_intraaural
  kv.cueWeights = 1;
elseif flags.do_interaural
  kv.cueWeights = [0,1];
end

if isscalar(kv.S)
  kv.S = repmat(kv.S,[length(cues),1]);
else
  kv.S = postpad(kv.S(:),length(cues));
end
kv.cueWeights = postpad(kv.cueWeights(:),length(cues))/sum(kv.cueWeights);
si = exp(-cues./kv.S);
bsi = nansum(kv.cueWeights .* si);

E = kv.range*bsi +kv.offset;%max(bsi);%min(1,max(bsi));%geomean(bsi);
if nargout >= 2
  varargout{1} = cues;
  varargout{2} = cueLbl;
end

  
end


