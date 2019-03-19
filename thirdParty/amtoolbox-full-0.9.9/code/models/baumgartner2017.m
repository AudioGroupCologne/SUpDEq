function [E,varargout] = baumgartner2017( target,template,varargin )
%BAUMGARTNER2017 Model for sound externalization
%   Usage:    [E,lat] = baumgartner2017( target,template )
%
%   Input parameters:
%     target  : binaural impulse response(s) referring to the directional 
%               transfer function(s) (DFTs) of the target sound(s).
%               Option 1: given in SOFA format -> sagittal plane DTFs will 
%               be extracted internally. 
%               Option 2: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane formated 
%               according to the following matrix dimensions: 
%               time x direction x channel/ear
%     template: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane referring to
%               the perceived lateral angle of the target sound.
%               Options 1 & 2 equivalent to target.
%
%   Output parameters:
%     E       : predicted degree of externalization
%     lat     : lateral angle as predicted by wierstorf2013 model
%
%   BAUMGARTNER2017(...) is a model for sound externalization.
%   It bases on the comparison of the intra-aural internal representation
%   of the incoming sound with a template and results in a probabilistic
%   prediction of polar angle response.
%
%   BAUMGARTNER2017 accepts the following optional parameters:
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
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2017
%
%   3) Circular Statistics Toolbox from http://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-
%
%
%   See also: baumgartner2014_spectralanalysis,
%   baumgartner2014_gradientextraction, baumgartner2014_binauralweighting
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/baumgartner2017.php

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

% AUTHOR: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria

%% Check input

definput.import={'baumgartner2014','baumgartner2014_pmv2ppp','localizationerror','amt_cache'};
definput.keyvals.tempWin = 1; % temporal integration window in sec
definput.flags.normalize = {'regular','normalize'};
definput.flags.cueProcessing = {'intraaural','interaural'};
definput.flags.middleEarFilter = {'','middleEarFilter'};
definput.keyvals.JND = 1.5;
definput.keyvals.c1 = 3.78;
definput.keyvals.c2 = 1;
definput.keyvals.z1 = 0.99;

[flags,kv]=ltfatarghelper(...
  {'fs','S','lat','stim','space','do','flow','fhigh',... %'fsstim'
  'bwcoef','polsamp','rangsamp','mrsmsp','gamma'},definput,varargin);

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
    num2str(kv.gamma,'%1.0u') '; Epsilon = ' num2str(kv.mrsmsp,'%1.0f') ' deg'])
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

%% Middle ear filter
if flags.do_middleEarFilter
  b=middleearfilter(kv.fs);
  target = filter(b,1,target);
  template = filter(b,1,template);
end

%% DTF filtering, Eq.(1)

dimtar = size(target); % for lconv dim check

if not(isempty(kv.stim))
  target = lconv(target,kv.stim);
end

% check that lconv preserved matrix dimensions (earlier bug in lconv)
if size(target,2) ~= dimtar(2)
  target = reshape(target,[size(target,1),dimtar(2:end)]);
end

frameLength = round((kv.tempWin*kv.fs));
Nframes = floor(size(target,1)/frameLength);
if Nframes == 0
  Nframes = 1;
  frameLength = size(target,1);
  target = cat(1,target,zeros(frameLength-size(target,1),size(target,2),size(target,3)));
end
Nang = size(template,2);
bsi = nan(Nframes,Nang);
tem.mp = [];
for iframe = 1:Nframes
  idt = (1:frameLength) + (iframe-1)*frameLength;
  
%% Spectral Analysis, Eq.(2)

[tar.mp,fc] = baumgartner2014_spectralanalysis(target(idt,:,:),'argimport',flags,kv);
if isempty(tem.mp) % integration across whole time range
  tem.mp = baumgartner2014_spectralanalysis(template,'argimport',flags,kv);
end


if flags.do_intraaural
  %% Positive spectral gradient extraction, Eq.(3)
  if kv.do == 1 % DCN inspired feature extraction
    nrep.tem = baumgartner2014_gradientextraction(tem.mp,fc);
    nrep.tar = baumgartner2014_gradientextraction(tar.mp,fc);
  else
    nrep.tem = tem.mp;
    nrep.tar = tar.mp;
  end

  %% Comparison process, Eq.(4)
%   sigma = baumgartner2017comparisonprocess(nrep.tar,nrep.tem); % based on vector product with 0s excluded (nan)
%   sigma = mean(repmat(nrep.tar,1,Nang).*nrep.tem)./mean(nrep.tem.^2);
  dNrep = abs(repmat(nrep.tar,1,Nang)-nrep.tem);
  dNrep(dNrep < kv.JND) = 0;
  sigma = mean(dNrep);

  %% Similarity estimation, Eq.(5)
  si = exp(-kv.S*sigma);
%   si = sigma.^kv.S;
  % si = baumgartner2014_similarityestimation(sigma,'argimport',flags,kv);

  %% Binaural weighting, Eq.(6)
  bsi(iframe,:) = baumgartner2014_binauralweighting(si,'argimport',flags,kv);

  %% Normalize
%   if flags.do_normalize
% %     if not(exist('bsiRef','var'))
%       sigmaRef = baumgartner2017comparisonprocess(nrep.tem,nrep.tem);
%       siRef = sigmaRef.^kv.S;
%       bsiRef = baumgartner2014_binauralweighting(siRef,'argimport',flags,kv);
% %     end
%     bsi(iframe) = bsi(iframe)/bsiRef;
%   end

else % interaural
  
  %% ILDs
  tar.ild = -diff(tar.mp,1,3); % ILD = left - right
  tem.ild = -diff(tem.mp,1,3);
  % figure; plot(fc,tar.ild); hold on; plot(fc,tem.ild); legend('tar','tem')

  %% target-template comparison -> ILD deviation
  dILD = abs(tem.ild-repmat(tar.ild,1,Nang));
  dILD(dILD < kv.JND) = 0; % limit minimum ILD difference according to JND

  %% overall normalized ILD deviation
  dILDnorm(iframe,:) = mean(dILD./abs(tem.ild));
  
  %% Externalization mapping
  bsi(iframe,:) = exp(-kv.S*dILDnorm(iframe,:));

end

%% Scaling
bsi(iframe) = kv.c1*bsi(iframe) +kv.c2;

end
E = max(bsi);%min(1,max(bsi));%geomean(bsi);
if nargout >= 2
  varargout{1} = kv.lat;
end

  
end
