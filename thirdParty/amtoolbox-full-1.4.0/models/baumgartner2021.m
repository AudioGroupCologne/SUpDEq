function [E,varargout] = baumgartner2021( target,template,varargin )
%BAUMGARTNER2021 Sound externalization based on multiple static cues
%   Usage:    [E,cues,cueLabels] = baumgartner2021( target,template )
%
%   Input parameters:
%     target  : binaural impulse response(s) referring to the directional 
%               transfer function(s) (DFTs) of the target sound(s).
%               Option 1: given in SOFA format -> sagittal plane DTFs will 
%               be extracted internally. 
%               Option 2: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane formatted 
%               according to the following matrix dimensions: 
%               time x direction x channel/ear
%     template: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane referring to
%               the perceived lateral angle of the target sound.
%               Options 1 & 2 equivalent to target.
%
%   Output parameters:
%     E         : predicted degree of externalization
%     cues      : outcomes of individual cues 
%     cueLabels : cue labels; cell array with 1st col. denoting acronyms
%                 and 2nd col. for descriptions
%
%   BAUMGARTNER2021(...) is a model framework for auditory
%   externalization perception. It enables to probe the contribution of
%   cue-specific expectation errors and to contrast dynamic versus static  
%   strategies for combining those errors within static listening environments.
%
%   BAUMGARTNER2021 accepts the following optional parameters:
%
%     'cueWeights',cW    Set the weights of individual cues to determine the final externalization score.
%                        Cue-specific weights (entered as a vector) are ordered as follows: 
%
%                        1 monaural spectral similarity (MSS)
%
%                        2 interaural spectral similarity of ILDs (ISS)
%
%                        3 spectral standard deviation of monaural gradients (MSSD)
%
%                        4 spectral standard deviation of ILDs (ISSD)
%
%                        5 interaural broadband time-intensity coherence (ITIT)
%
%                        6 interaural coherence (IC)
%
%                        7 monaural intensity difference (MI)
%
%                        8 temporal standard deviation of ILDs (ITSD). 
%
%                        Default weights are 0.6 for MSS, 0.4 for ISS, and 0 for all others.
%
%     'S',S          Set the cue-specific sensitivity parameter to S. 
%                    1/S represents the slope of sigmoidal mapping
%                    function.
%                    Vector order equivalent to cueWeights.
%                    Default values are determined by the weighted average
%                    sensitivities determined in Baumgartner and Majdak
%                    (2020) - run exp_baumgartner2021('tab2') to show them.
%
%     'lat',lat      Set the apparent lateral angle of the target sound to
%                    lat. Default value is 0 degree (median SP).
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
%   BAUMGARTNER2021 accepts the following flags:
%
%     'LTA'          Looser-takes-all strategy: Model selects minimal 
%                    predicted externalization scores across cues with 
%                    weights larger than zero.
%     'MTA'          Median-takes-all strategy: Model selects median 
%                    predicted externalization scores across cues with 
%                    weights larger than zero.
%     'WTA'          Winner-takes-all strategy: Model selects maximal 
%                    predicted externalization scores across cues with 
%                    weights larger than zero.
%                    
%
%   Requirements:
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
%
%   3) Circular Statistics Toolbox from http://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-
%
%
%   See also: baumgartner2021_mapping,
%   baumgartner2016_spectralanalysis, baumgartner2016_gradientextraction, 
%   baumgartner2014_binauralweighting demo_baumgartner2021 baumgartner2014
%   li2020 baumgartner2016
%
%   References:
%     R. Baumgartner and P. Majdak. Decision making in auditory
%     externalization perception: model predictions for static conditions.
%     Acta Acustica, 5:59, 2021. Publisher: EDP Sciences. [1].html ]
%     
%     References
%     
%     1. https://acta-acustica.edpsciences.org/articles/aacus/abs/2021/01/aacus210038/aacus210038.html
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/baumgartner2021.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: SOFA M-Stats O-Statistics
%   #Author: Robert Baumgartner (2021), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Check input

definput.import={'amt_cache','baumgartner2014','baumgartner2016','baumgartner2021'};
[flags,kv]=ltfatarghelper({'cueWeights','S'},definput,varargin);

flags.do_plot = false;

if not(isstruct(target)) && ismatrix(target)
  target = permute(target,[1,3,2]);
end

if not(isstruct(template)) && ismatrix(template)
  template = permute(template,[1,3,2]);
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
  [target,~] = extractsp( kv.lat,target );
end

if isstruct(template) % Template given in SOFA format
  [template,~] = extractsp( kv.lat,template );
end

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

%% Level difference
MI = dbspl(target) - dbspl(template);
MI = abs(MI);
MI(abs(MI)<kv.ILD_JND) = 0;
MI = MI./dbspl(template);
MI = baumgartner2014_binauralweighting(MI,'argimport',flags,kv); % Eq. 1

%% ITD & IC
[tem.itd,~,tem.iacc] = itdestimator(shiftdim(template,1),'fs',kv.fs,'MaxIACCe');
[tar.itd,~,tar.iacc] = itdestimator(shiftdim(target,1),'fs',kv.fs,'MaxIACCe');

for ipos = 1:size(template,2)
  tem.ic(ipos) = baumgartner2017_iacc(squeeze(template(:,ipos,:)),'argimport',flags,kv);
end
for ipos = 1:size(target,2)
  tar.ic(ipos) = baumgartner2017_iacc(squeeze(target(:,ipos,:)),'argimport',flags,kv);
end
IC = abs(tar.ic-tem.ic) ./ tem.ic; % Eq. 2

%% Filterbank
[tem.mp,fc] = baumgartner2016_spectralanalysis(template,70,'argimport',flags,kv,'tiwin',size(template,1)*kv.fs,'gammatone','redo');
tar.mp = baumgartner2016_spectralanalysis(target,70,'argimport',flags,kv,'tiwin',size(target,1)*kv.fs,'gammatone','redo');

%% interaural temporal SD of ILDs (Catic et al., 2015)
MinNumTimeFrames = 20;
if size(tem.mp,5) >= MinNumTimeFrames && size(tar.mp,5) >= MinNumTimeFrames
  tem.STild = -diff(tem.mp,1,3); % short-term ILDs
  tar.STild = -diff(tar.mp,1,3);
  ITSD = 1 - mean(std(tar.STild,0,5)./std(tem.STild,0,5));
else
  ITSD = nan(size(IC));
end

%% Echo suppression
if isscalar(kv.reflectionOnsetTime) % evaluate only direct path (DP)
  idDP = round(kv.reflectionOnsetTime*kv.fs);
  N1ms = round(1e-3*kv.fs); % 1 ms fade out
  taper = [ones(1,idDP-N1ms) , 0.5*(1+cos(linspace(0,pi,N1ms)))];
  temTaper = repmat([taper(:);zeros(size(template,1)-idDP,1)],...
                    [1,size(template,2),size(template,3)]);
  tarTaper = repmat([taper(:);zeros(size(target,1)-idDP,1)],...
                    [1,size(target,2),size(target,3)]);
  [tem.mp,fc] = baumgartner2016_spectralanalysis(temTaper.*template,70,... % or save to temDP.mp
    'argimport',flags,kv,'tiwin',kv.tempWin,'gammatone','redo');
  tar.mp = baumgartner2016_spectralanalysis(tarTaper.*target,70,...        % or save to tarDP.mp
    'argimport',flags,kv,'tiwin',kv.tempWin,'gammatone','redo');
end

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
tem.psg = baumgartner2016_gradientextraction(tem.mp,fc,'mgs',1,flags.gradients);
tem.ild = diff(tem.mp,1,3);

tar.psg = baumgartner2016_gradientextraction(tar.mp,fc,'mgs',1,flags.gradients);
tar.ild = diff(tar.mp,1,3);


%% Spectral comparison
MSSD = abs(std(tar.psg.m)-std(tem.psg.m))./std(tem.psg.m);
MSSD = baumgartner2014_binauralweighting(MSSD,'argimport',flags,kv);
ISSD = abs(std(tar.ild)-std(tem.ild))./std(tem.ild);

for iSC = 1:2 % first monaural then interaural
  if iSC == 1 % monaural spectral gradients
    tem.nrep = tem.psg.m;
    tar.nrep = tar.psg.m;
  elseif iSC == 2 % interaural spectral differences
    tem.nrep = tem.ild;
    tar.nrep = tar.ild;
  end
  
  % comparison with time average of spectral template
  targetprofile = {tar.nrep};
  if flags.do_dprime
    targetprofile{2} = tem.nrep;
  end
  tem.nrep = mean(tem.nrep,5);
  d_cue = cell(length(targetprofile),1);
  for inrep = 1:length(targetprofile)  
    templateprofile = repmat(tem.nrep,[1,1,1,1,size(targetprofile{inrep},5)]);
    delta = abs(templateprofile-targetprofile{inrep});
    delta(delta < kv.ILD_JND) = 0; % limit minimum ILD difference according to JND
    d_cue{inrep} = mean(delta./(eps+abs(templateprofile))); % Eq. (4)
    if iSC == 1 % do_intraaural
      d_cue{inrep} = baumgartner2014_binauralweighting(d_cue{inrep},'argimport',flags,kv); % Eq. 7
    end
  end

  % temporal integration
  if flags.do_dprime % signal detection theory applied to time histograms
  %   figure; histogram(sigma{1}); hold on ; histogram(sigma{2}); legend('target','reference')
    allsigma = [d_cue{1}(:);d_cue{2}(:)];
    msigma = mean(allsigma);
    sdsigma = std(allsigma);
    mzsigma(1) = mean((d_cue{1}-msigma) ./ sdsigma);
    mzsigma(2) = mean((d_cue{2}-msigma) ./ sdsigma);
    dprime = max(mzsigma(1)-mzsigma(2),0);
    distmetric = dprime;
  elseif size(d_cue{1},5) == 1
    distmetric = d_cue{1};
  else % temporal weighting according to amount of sensory information available
    tweight = mean(mean(abs(tar.nrep)),3); % temporal weighting
    tweight = tweight-min(tweight,[],5); % bounded between 0
    tweight = 2*tweight./max(tweight,[],5); % and 2
    distmetric = d_cue{1}.*tweight;
    distmetric = mean(distmetric,5);
  end
  
  if iSC == 1
    MSS = distmetric; % monaural spectral distance
  elseif iSC == 2
    ISS = distmetric; % interaural spectral distance
  end
  
end

%% Interaural time-intenstiy coherence (ITIC) -> dprime possible
ITR = tar.itd./(eps+tem.itd) -1; % Eq. 5a
ITR(abs(tar.itd - tem.itd) < kv.ITD_JND) = 0;
ILR = mean(tar.ild)./(mean(tem.ild)+eps) -1; % Eq. 5b
ILR(any(abs(mean(tar.ild) - mean(tem.ild)) < kv.ILD_JND,1)) = 0;
ITIT = abs( ITR -  ILR(:) )'; % Eq. 5(a+b)

%% Cue integration/weighting
cues = [MSS; ISS; MSSD; ISSD; ITIT; IC; MI; ITSD];
cueLbl = {'MSS',['Monaural ',flags.gradients,' spectral gradients (c.f., Baumgartner et al., 2014)']; ...
          'ISS','Interaural spectral shape (c.f., Hassager et al., 2016)'; ...
          'MSSD',['Spectral SD of monaural ',flags.gradients,' spectral gradients (c.f., Baumgartner et al., 2014)']; ...
          'ISSD','Spectral SD of interaural spectral differences (c.f., Georganti et al., 2013)'; ...
          'ITIT','Interaural time-intensity trading (ITD vs. ILD)'; ...
          'IC','Interaural coherence (c.f., Hassager et al., 2017)'; ...
          'MI','Monaural intensity difference (target - reference)'; ...
          'ITSD','Interaural temporal standard deviation (c.f., Catic et al., 2015)'; ...
         };


if isscalar(kv.S)
  kv.S = repmat(kv.S,[size(cues,1),1]);
else
  kv.S = postpad(kv.S(:),size(cues,1));
end

E = baumgartner2021_mapping(cues,kv.S,kv.range,kv.offset); % Eq. 8

kv.cueWeights = postpad(kv.cueWeights(:),size(cues,1))/sum(kv.cueWeights);
if flags.do_LTA
  E = nanmin(E(kv.cueWeights>0));
elseif flags.do_MTA 
  E = median(E(kv.cueWeights>0),'omitnan');
elseif flags.do_WTA
  E = nanmax(E(kv.cueWeights>0));   
elseif flags.do_fixedWeights || flags.do_WSM
  E = nansum(kv.cueWeights .* E);
end

if nargout >= 2
  varargout{1} = cues;
  varargout{2} = cueLbl;
end


