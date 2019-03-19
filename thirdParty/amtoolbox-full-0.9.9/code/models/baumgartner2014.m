function varargout = baumgartner2014( target,template,varargin )
%BAUMGARTNER2014 Model for localization in sagittal planes
%   Usage:    [p,rang] = baumgartner2014( target,template )
%             [p,rang,tang] = baumgartner2014( target,template )
%             [p,rang,tang] = baumgartner2014( target,template,fs,S,lat,stim,fsstim )
%             [err,pred] = baumgartner2014( target,template,errorflag )
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
%     p       : predicted probability mass vectors for response angles 
%               with respect to target positions
%               1st dim: response angle
%               2nd dim: target angle
%     rang    : polar response angles (after regularization of angular 
%               sampling)
%     tang    : polar target angles (usefull if sagittal-plane HRTFs are
%               extracted directly from SOFA object)
%     err     : predicted localization error (acc. to performance measure
%               defined in errorflag*
%     pred    : structure with fields p, rang, tang*
%
%   BAUMGARTNER2014(...) is a model for sound-source localization
%   in sagittal planes (SPs). It bases on the comparison of internal sound 
%   representation with a template and results in a probabilistic
%   prediction of polar angle response.
%
%   BAUMGARTNER2014 accepts the following optional parameters:
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
%     'gamma',G      Set the degree of selectivity 
%                    (slope of the sigmoid link function representing 
%                    the psychometric link between transformation from the
%                    distance metric and similarity index) to G. 
%                    Default value is 6.
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
%     'polsamp',ps   Define the polar-angle sampling of the acoustic data
%                    provided for the current sagittal plane. As default the 
%                    sampling of ARI's HRTFs in the median SP is used, i.e.,
%                    ps = [-30:5:70,80,100,110:5:210] degrees.
%
%     'rangsamp',rs  Define the equi-polar sampling of the response predictions.
%                    The default is rs = 5 degrees.
%
%     'mrsmsp',eps   Set the motoric response scatter eps within the median 
%                    sagittal plane. Default value is 17 degrees.
%
%   BAUMGARTNER2014 accepts the following flags:
%
%     'regular'      Apply spline interpolation in order to regularize the 
%                    angular sampling of the polar response angle. 
%                    This is the default.
%
%     'noregular'    Disable regularization of angular sampling.
%
%     'errorflag'    May be one of the error flags defined in
%                    BAUMGARTNER2014_pmv2ppp or localizationerror.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2014
%
%   3) Circular Statistics Toolbox from http://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-
%
%
%   See also: plot_baumgartner2014, data_baumgartner2014,
%   exp_baumgartner2014, demo_baumgartner2014, baumgartner2014_calibration,
%   baumgartner2014_likelistat, baumgartner2014_pmv2ppp,
%   baumgartner2014_virtualexp, baumgartner2014_spectralanalysis,
%   baumgartner2014_gradientextraction, baumgartner2014_comparisonprocess,
%   baumgartner2014_similarityestimation, baumgartner2014_binauralweighting,
%   baumgartner2014_sensorimotormapping
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791-802, 2014.
%     
%     R. Lyon. All pole models of auditory filtering. In E. R. Lewis, G. R.
%     Long, R. F. Lyon, P. M. Narins, C. R. Steele, and E. Hecht-Poinar,
%     editors, Diversity in auditory mechanics, pages 205-211. World
%     Scientific Publishing, Singapore, 1997.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/baumgartner2014.php

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

definput.import={'baumgartner2014','baumgartner2014_pmv2ppp','localizationerror'};

[flags,kv]=ltfatarghelper(...
  {'fs','S','lat','stim','space','do','flow','fhigh',... %'fsstim'
  'bwcoef','polsamp','rangsamp','mrsmsp','gamma'},definput,varargin);

if not(isstruct(target)) && ndims(target) == 2
  target = permute(target,[1,3,2]);
  warning('Matrix dimensions of target should be: time x direction x channel/ear.')
  warning('Since 3rd dimension was empty, 2nd dimension was used as channel dimension.')
end

%% Print Settings

if flags.do_print 
  if flags.do_nomrs
    kv.mrsmsp = 0;
  end
  amt_disp(['Settings: PSGE = ' num2str(kv.do,'%1.0f') '; Gamma = ' ...
    num2str(kv.gamma,'%1.0u') '; Epsilon = ' num2str(kv.mrsmsp,'%1.0f') ' deg'])
end


%% Extract HRTFs of sagittal plane

if isstruct(target) % Targets given in SOFA format
  kv.fs = target.Data.SamplingRate;
  [target,tang] = extractsp( kv.lat,target );
end

if isstruct(template) % Template given in SOFA format
  [template,kv.polsamp] = extractsp( kv.lat,template );
end

% Error handling
if size(template,2) ~= length(kv.polsamp)
  fprintf('\n Error: Second dimension of template and length of polsamp need to be of the same size! \n')
  return
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

%% Spectral Analysis, Eq.(2)

[ireptar,fc] = baumgartner2014_spectralanalysis(target,'argimport',flags,kv);
ireptem = baumgartner2014_spectralanalysis(template,'argimport',flags,kv);


%% Positive spectral gradient extraction, Eq.(3)

if kv.do == 1 % DCN inspired feature extraction
  nrep.tem = baumgartner2014_gradientextraction(ireptem,fc);
  nrep.tar = baumgartner2014_gradientextraction(ireptar,fc);
else
  nrep.tem = ireptem;
  nrep.tar = ireptar;
end


%% Comparison process, Eq.(4)

sigma = baumgartner2014_comparisonprocess(nrep.tar,nrep.tem);


%% Similarity estimation, Eq.(5)

si = baumgartner2014_similarityestimation(sigma,'argimport',flags,kv);


%% Binaural weighting, Eq.(6)

si = baumgartner2014_binauralweighting(si,'argimport',flags,kv);


%% Sensorimotor mapping, Eq.(7)

[si,rang] = baumgartner2014_sensorimotormapping(si,'argimport',flags,kv);


%% Normalization to PMV, Eq.(8)
p = si ./ repmat(sum(si)+eps,size(si,1),1);


%% Performance measures
if not(isempty(flags.errorflag)) % Simulate virtual experiments
  
  m = baumgartner2014_virtualexp(p,tang,rang,'targetset',kv.exptang);
  err = localizationerror(m,flags.errorflag);
  
elseif not(isempty(flags.ppp)) % Calculate directly via probabilities

  if flags.do_QE_PE_EB
    [err.qe,err.pe,err.pb] = baumgartner2014_pmv2ppp(p,tang,rang,'exptang',kv.exptang);
  else
    err = baumgartner2014_pmv2ppp(p,tang,rang,flags.ppp,'exptang',kv.exptang);
  end
    
end


%% Output
if isempty([flags.errorflag,flags.ppp])
  varargout{1} = p;
  if nargout >= 2
    varargout{2} = rang;
    if nargout >= 3
      try
        varargout{3} = tang;
      catch
        disp('SOFA Object of target DTFs is required to output target angles.')
      end
    end
  end
else
  varargout{1} = err;
  if nargout > 1
    varargout{2} = struct('p',p,'rang',rang,'tang',tang);
  end
end
  
end
