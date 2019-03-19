function varargout = baumgartner2016( target,template,varargin )
%BAUMGARTNER2016 Level-dependent model for localization in sagittal planes
%   Usage:    [p,rang,tang] = baumgartner2016( target,template )
%             [err,pred,m] = baumgartner2016( target,template,errorflag )
%
%   Input parameters:
%     target  : binaural impulse response(s) referring to the directional 
%               transfer function(s) (DFTs) of the target sound(s).
%     template: binaural impulse responses of all available
%               listener-specific DTFs of the sagittal plane referring to
%               the perceived lateral angle of the target sound
%     label:    2-element cell array including labels used for caching. 
%               Provide listener ID (e.g., 'NH12') in label{1}, and
%               target condition information in label{2}.
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
%     m       : item list from virtual experiment. See 
%               help localizationerror for format description.
%
%   BAUMGARTNER2016(...) is a model for sound-source localization
%   in sagittal planes (SPs). It bases on the comparison of internal sound 
%   representation with a template and results in a probabilistic
%   prediction of polar angle response.
%
%   BAUMGARTNER2016 accepts the following optional parameters:
%
%     'ID'           Listeners ID (important for caching).
%
%     'Condition'    Label of experimental condition (also for caching).
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
%     'stim'         Define the stimulus (source signal without directional
%                    features). As default temstim is used.
%
%     'fsstim'       Define the sampling rate of the stimulus. 
%                    Default value is 48000 Hz.
%
%     'temstim'      Define the dummy stimulus used to create the templates.
%                    The default is Gaussian white noise with a duration of 170 ms.
%
%     'SPL',L        Set the SPL of the stimulus to L dB.
%                    Default value is 60dB.
%
%     'SPLtem',Lt    Set the SPL of the templates to a specific SPL of Lt*dB
%                    if Lt is a scalar or define a SPL range between
%                    Lt(1) and Lt(2)*dB if Lt is a two-element vector.
%                    Default range is 40 to 70dB.
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
%     'polsamp',ps   Define the the polar angular sampling of the current
%                    sagittal plane. As default the sampling of ARIs HRTF 
%                    format at the median SP is used, i.e.,
%                    ps = [-30:5:70,80,100,110:5:210] degrees.
%
%     'mrsmsp',mrs   Set the motoric response scatter mrs within the median 
%                    sagittal plane. Default value is 17 degrees.
%
%     'cohc',cohc    OHC scaling factor: 1 denotes normal OHC function (default); 
%                    0 denotes complete OHC dysfunction.
%
%     'cihc',cihc    IHC scaling factor: 1 denotes normal IHC function (default); 
%                    0 denotes complete IHC dysfunction.
%
%     'fiberTypes',fT Types of the fibers based on spontaneous rate (SR) in 
%                     spikes/s: fT=1 for Low SR; fT=2 for Medium SR; 
%                     fT=3 for High SR. Default is fT = 1:3.
%
%   BAUMGARTNER2016 accepts the following flags:
%
%     'zilany2014'   Use the model from Zilany et al. (2009,2014) for spectral 
%                    analysis. This is the default.
%
%     'gammatone'    Use the Gammatone filterbank for spectral analysis. 
%
%     'SPLtemAdapt'  Set SPL of templates acc. to target (*Lt*=*L*). 
%
%     'NHtem'        No adaptation of templates to hearing impairment,
%                    i.e., templates are processed with cohc=cihc=1 and
%                    fT = 1:3.
%
%     'ihc'          Incorporate the transduction model of inner hair 
%                    cells used by Dau et al. (1996).
%
%     'noihc'        Do not incorporate the IHC stage. This is the default.
%
%     'regular'      Apply spline interpolation in order to regularize the 
%                    angular sampling of the polar response angle. 
%                    This is the default.
%
%     'noregular'    Disable regularization of angular sampling.
%
%     'errorflag'    May be one of the error flags defined in
%                    baumgartner2014_pmv2ppp or localizationerror.
%
%     'redoSpectralAnalysis' Flag to redo also spectral analysis based on
%                    zilany2014 model.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2016
%
%
%   See also: data_baumgartner2016,
%   exp_baumgartner2016, baumgartner2016_calibration,
%   baumgartner2014_pmv2ppp,
%   baumgartner2014_virtualexp, localizationerror,
%   baumgartner2016_spectralanalysis
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling the effects of
%     sensorineural hearing loss on auditory localization in the median
%     plane. Trends in Hearing, 20:1-11, 2016.
%     
%     M. S. A. Zilany, I. C. Bruce, and L. H. Carney. Updated parameters and
%     expanded simulation options for a model of the auditory periphery. The
%     Journal of the Acoustical Society of America, 135(1):283-286, Jan.
%     2014.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/baumgartner2016.php

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

%% Check input options 

definput.import={'baumgartner2016'};
posdepnames = {'fs','S','lat','stim','fsstim'};

[flags,kv]=ltfatarghelper(posdepnames,definput,varargin);

% Check default condition label for short stimuli
if length(kv.stim)/kv.fsstim < 0.005 && strcmp(kv.Condition,'Long') % short stimulus
  kv.Condition = '';
end

% Extract sagittal-plane HRTFs
if isstruct(target) % Targets given in SOFA format
  kv.fs = target.Data.SamplingRate;
  [target,tang] = extractsp( kv.lat,target );
end

if isstruct(template) % Template given in SOFA format
  [template,kv.polsamp] = extractsp( kv.lat,template );
end


%% Error handling
if size(template,2) ~= length(kv.polsamp)
  fprintf('\n Error: Second dimension of template and length of polsamp need to be of the same size! \n')
  return
end

%% Stimulus 

% flags.temstim = false;
if isempty(kv.stim) 
    kv.stim = kv.temstim;
%     flags.temstim = true;
    kv.fsstim = kv.fs;
elseif isempty(kv.fsstim) 
    kv.fsstim = kv.fs;
end

if flags.do_headphone% || flags.do_drnl
    hpfilt = headphonefilter(kv.fs);
    kv.stim = lconv(kv.stim,hpfilt(:));
end


%% DTF filtering
if ~isequal(kv.fs,kv.fsstim)
    amt_disp('Sorry, sampling rate of stimulus and HRIRs must be the same!')
    return
end

tmp = lconv(target,kv.stim);
target = reshape(tmp,[size(tmp,1),size(target,2),size(target,3)]); % workaround for lconv bug

tmp = lconv(template,kv.temstim);
template = reshape(tmp,[size(tmp,1),size(template,2),size(template,3)]); % workaround for lconv bug
    
%% Cochlear filter bank -> internal representations

if flags.do_redoSpectralAnalysis
  redoSAflag = 'redo';
else
  redoSAflag = 'normal';
end

% Target profile
[mreptar,fc] = baumgartner2016_spectralanalysis(target,kv.SPL,'target',...
  'argimport',flags,kv,redoSAflag);

% Template profile
if flags.do_NHtem
  kv.cohc = 1;
  kv.cihc = 1;
  kv.fiberTypes = 1:3;
end
if flags.do_SPLtemAdapt
  kv.SPLtem = kv.SPL;
end
if isscalar(kv.SPLtem) % all templates represented at a fixed SPL
  mreptem = baumgartner2016_spectralanalysis(template,kv.SPLtem,'template',...
    'argimport',flags,kv,redoSAflag);
else % average across SPL range
  mreptem = baumgartner2016_spectralanalysis(template,kv.SPLtem(1),'template',...
    'argimport',flags,kv,redoSAflag);
  SPLtemRange = kv.SPLtem(1):10:kv.SPLtem(2);
  for ii = 2:length(SPLtemRange)
    mreptem = mreptem + ...
      baumgartner2016_spectralanalysis(template,SPLtemRange(ii),'template',...
        'argimport',flags,kv,redoSAflag);
  end
  mreptem = mreptem/length(SPLtemRange);
end

%% Positive spectral gradient extraction

if kv.do == 1 && flags.do_psge
  
  % duration-dependence of c2 (inhibitory strength)
  if length(kv.stim)/kv.fs < 5e-3 % halfen if shorter than 5ms
    c2 = kv.psgeshort;
  else
    c2 = 1;
  end
  
  if flags.do_gammatone
    greptar = baumgartner2014_gradientextraction(mreptar,fc,c2);
    greptem = baumgartner2014_gradientextraction(mreptem,fc);
  else % zilany
    greptar = baumgartner2016_gradientextraction(mreptar,fc,'argimport',flags,kv);
    greptem = baumgartner2016_gradientextraction(mreptem,fc,'argimport',flags,kv);
  end
    
else

  amt_disp('Cue extraction different to baumgartner2014','progress')
   
  if kv.do == 1 && flags.do_dcn% DCN inspired feature extraction
    for ch = 1:size(mreptar,3)
      greptem(:,:,ch) = dcn(mreptem(:,:,ch),fc);
      greptar(:,:,ch) = dcn(mreptar(:,:,ch),fc);
    end
  elseif kv.do == 2 % proposed by Zakarauskas & Cynader (1993)
      greptem = diff(mreptem,kv.do);
      greptar = diff(mreptar,kv.do);
  else
      greptem = mreptem;
      greptar = mreptar;
  end
end

%% Comparison process

if flags.do_isd && not(flags.do_intensityweighting) && not(flags.do_diff)
  
  if flags.do_gammatone
    sigma = baumgartner2014_comparisonprocess(greptar,greptem);
  else % zilany
    sigma = baumgartner2016_comparisonprocess(greptar,greptem);
  end

else
  
  amt_disp('Comparison process different to baumgartner2014','progress')
  
  if flags.do_diff
    greptar = diff(greptar);
    greptem = diff(greptem);
  end

  sigma=zeros(size(mreptem,2),size(mreptar,2),size(mreptem,3)); % init
  for ch = 1:size(mreptar,3)
    for it = 1:size(mreptar,2)
      if flags.do_isd
        isd = repmat(greptar(:,it,ch),[1,size(greptem(:,:,ch),2),1]) - greptem(:,:,ch); 
        if kv.do == 0
          sigma(:,it,ch) = sqrt(squeeze(var(isd))); % standard dev. across frequencies (Middlebrooks, 1999)
        else
          if isnan(nanmean(abs(isd)))
            sigma(:,it,ch) = 0;
          else
            if flags.do_intensityweighting
              Ifw = mreptar(:,it,ch)/max(mreptar(:,it,ch));
              isd = isd.*repmat(Ifw,1,size(isd,2));
            end
            sigma(:,it,ch) = nanmean(abs(isd)); % L1-norm across frequencies
          end
        end
      elseif flags.do_corr
        if sum(greptar(:,it,ch)==0) == length(greptar(:,it,ch))
          sigma(:,it,ch) = 0;
        else
          r = corrcoef([greptar(:,it,ch) greptem(:,:,ch)]); % in case of NAN use flag: 'rows','pairwise' (but very slow!!!)
          sigma(:,it,ch) = r(1,2:end);
        end
      end
    end
  end

end


%% OPTIMAL FT SELECTION
Ntang = size(mreptar,2);

if flags.do_ftopt

% Var 1: Winner takes it all
% si = max(si_ft,[],4);

% Var 2: Select for each target the SI distribution of the fiber type
% yielding the minimum internal distance metric (sigma)
if length(kv.fiberTypes) == 3

  sigma_ft = sigma;
  
  sigma = sigma_ft(:,:,:,1,:); % init with ft=1
  N = zeros(Ntang * size(mreptar,3),3);
  N(:,1) = 1;
  for ch = 1:size(mreptar,3)
    for it = 1:Ntang
      for ft = 2:3
        if min(sigma_ft(:,it,ch,ft,:)) < min(sigma(:,it,ch,:,:))
          sigma(:,it,ch,:,:) = sigma_ft(:,it,ch,ft,:);
          N(it+(ch-1)*Ntang,ft) = 1;
        end
      end
    end
  end
  N(:,2) = max(N(:,2) - N(:,3),0);
  N(:,1) = N(:,1) - (N(:,2)+N(:,3));
  N = sum(N);
  fileID = fopen(fullfile(amt_basepath,'modelstages','baumgartner2016_fiberTypeCounter'),'a');
  fprintf(fileID,'%i, %i, %i\n',N(1),N(2),N(3));
  fclose(fileID);

end

end

%% Evaluation of multiple looks
if size(sigma,5) > 1
  sigma = mean(sigma,5);
end

%% Similarity estimation

if flags.do_isd && flags.do_sigmoid
  
  % duration-dependence of gamma (degree of selectivity)
  if length(kv.stim)/kv.fs < 5e-3 % halfen if shorter than 5ms
    kv.gamma = kv.gamma*kv.gammashortfact;
  end
  % duration-dependence of sensitivity
  if length(kv.stim)/kv.fs < 5e-3 % halfen if shorter than 5ms
    kv.S = kv.S*kv.Sshortfact;
  end
    
  if not(flags.do_gammatone)
    sigma = 10*sigma;
  end
  si = baumgartner2014_similarityestimation(sigma,'S',kv.S,'gamma',kv.gamma);

else
  
  amt_disp('Similarity estimation different to baumgartner2014','progress')
  
  if flags.do_zilany2014 && not(flags.do_sigmoid || flags.do_dcn)
    sigma = sigma/10;
  end

  si=zeros(size(sigma)); % init
  for ch = 1:size(mreptar,3)
    for it = 1:size(mreptar,2)

      if flags.do_isd
        if flags.do_sigmoid
          si(:,it,ch) = 1+eps - (1+exp(-kv.gamma*(sigma(:,it,ch)-10*kv.S))).^-1;
        elseif flags.do_exp
          si(:,it,ch) = real(exp(-sigma(:,it,ch)./kv.S));
        elseif flags.do_Gauss
          si(:,it,ch) = real(exp(-(sigma(:,it,ch)./kv.S).^2));
        else
          si(:,it,ch) = interp1([0,kv.SimDL,kv.S,1e10],[1,1,0,0],sigma(:,it,ch),'cubic');
        end
      elseif flags.do_corr
        if flags.do_exp
          si(:,it,ch) = real(exp( (sigma(:,it,ch) - 1)./kv.S ));
        elseif flags.do_Gauss
          si(:,it,ch) = real(exp(-((sigma(:,it,ch) - 1)./kv.S).^2));
        else % power law relationship
          si(:,it,ch) = real(sigma(:,it,ch).^kv.S);
        end
      end

    end
  end
end
  
%% Binaural weighting

si = baumgartner2014_binauralweighting(si,'lat',kv.lat,'bwcoef',kv.bwcoef);


%% Interpolation, prior, sensorimotor mapping
if kv.prior==0
  
  [si,rang] = baumgartner2014_sensorimotormapping(si,...
    'rangsamp',kv.rangsamp,'polsamp',kv.polsamp,'lat',kv.lat,'mrsmsp',kv.mrsmsp,...
    flags.regularization,flags.motoricresponsescatter);
  
else
  
  if flags.do_regular
      rang0 = ceil(min(kv.polsamp)*0.2)*5;    % ceil to 5 deg
      rang = rang0:kv.rangsamp:max(kv.polsamp);
      siint = zeros(length(rang),size(si,2));
      for tt = 1:size(si,2)
          siint(:,tt) = interp1(kv.polsamp,si(:,tt),rang,'spline');
      end
      si = siint;
      si(si<0) = 0; % SIs must be positive (necessary due to spline interp)
  else
      rang = kv.polsamp;
  end

  % Individual prior distribution
  if min(kv.priordist.x) > -90
    kv.priordist.x = [-90;kv.priordist.x(:)];
    kv.priordist.y = [1;kv.priordist.y(:)];
  end
  if max(kv.priordist.x) < 270
    kv.priordist.x = [kv.priordist.x(:);270];
    kv.priordist.y = [kv.priordist.y(:);1];
  end
  priordist = interp1(kv.priordist.x,kv.priordist.y,rang,'linear') .^ kv.prior;
  si = si.*repmat(priordist(:),1,size(si,2));

  % Sensorimotor mapping
  if flags.do_mrs && flags.do_regular && kv.mrsmsp > 0

      angbelow = -90:5:min(rang)-5;
      angabove = max(rang)+5:5:265;
      rang = [angbelow,rang,angabove];
      si = [zeros(length(angbelow),size(si,2)) ; si ; zeros(length(angabove),size(si,2))];

      mrs = kv.mrsmsp/cos(deg2rad(kv.lat)); % direction dependent scatter (derivation: const. length rel. to the circumferences of circles considered as cross sections of a unit sphere)

      x = 0:2*pi/72:2*pi-2*pi/72;
      kappa = 1/deg2rad(mrs)^2; % concentration parameter (~1/sigma^2 of normpdf)
      mrspdf = exp(kappa*cos(x)) / (2*pi*besseli(0,kappa)); % von Mises PDF 
      for tt = 1:size(si,2)
        si(:,tt) = pconv(si(:,tt),mrspdf(:));
      end

  end
end

%% Normalization to PMV
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
          amt_disp('SOFA Object of target DTFs is required to output target angles.')
        end
      end
  end
else
  varargout{1} = err;
  if nargout > 1
    varargout{2} = struct('p',p,'rang',rang,'tang',tang);
    if nargout > 2
      if not(exist('m','var'))
        m = baumgartner2014_virtualexp(p,tang,rang);
      end
      varargout{3} = m;
    end
  end
end
  
end

% function t4 = psge(an,fc,c2)
% %DCN Phenomenological model of dorsal cochlear nucleus (DCN)
% %   Usage:      out = dcn(in)
% %
% %   Input parameters:
% %     an      : spectral profile in dB
% %
% %   Output parameters:
% %     t4      : activity of type IV unit
% 
% %% Parameter Settings
% % c2 = 1; % inhibitory coupling between type II and type IV neurons
% c4 = 1; % coupling between an and type IV neuron
% dilatation = 1; % of tonotopical 1-ERB-spacing between type IV and II neurons
% 
% erb = audfiltbw(fc);
% 
% %% Calculations
% Nb = size(an,1); % # auditory bands
% dt4t2 = round(mean(erb(2:end)./diff(fc))*dilatation); % tonotopical distance between type IV and II neurons
% t4 = zeros(Nb-dt4t2,size(an,2),size(an,3)); % type IV output
% for b = 1:Nb-dt4t2
%   t4(b,:,:) = c4 * an(b+dt4t2,:,:) - c2 * an(b,:,:);
% end
% 
% t4 = (t4 + c2*abs(t4))/2; %disp('only rising edges')
% % t4(t4<0) = nan;
% end

function P = dcn(an,fc)
%DCN Phenomenological model of dorsal cochlear nucleus (DCN)
%   Usage:      P = dcn(an,fc)
%
%   Input parameters:
%     an      : spectral profile in dB
%
%   Output parameters:
%     P       : activity of principal cell

%% Parameter Settings (Table 2, Lomakin & Davis, 2008)
% offsets in octaves
c.WtoI2 = .3;
c.WtoP = .2;
c.I2toP = -.1;
% bandwidth in octaves
bw.ANtoW = 2.5; % effectively .83 due to Gaussian distribution
bw.ANtoI2 = .4;
bw.ANtoP = .4;
bw.WtoI2 = 2.2; % from Fig. 7
bw.WtoP = 2.2; % from Fig. 7
bw.I2toP = .2;
% # projections times strength
Ns.ANtoW = 140*.05;
Ns.ANtoI2 = 48*.55;
Ns.ANtoP = 48*.25;
Ns.WtoI2 = 15*1.4;
Ns.WtoP = 15*.6; % from page 514
Ns.I2toP = 21*2.25;
% non-specific afferents (NSA) leading to spontaneous activity
nsa = 3000; % spikes/s

%% Parameter Settings (Table 1, Reiss & Young, 2005)
% % offsets in octaves
% c.WtoI2 = .3;
% c.WtoP = .05;
% c.I2toP = -.1;
% % bandwidth in octaves
% bw.ANtoI2 = .2;
% bw.ANtoW = 2.0;
% bw.ANtoP = .24;
% bw.WtoI2 = 0.1;%0.05;
% bw.WtoP = 0.1;
% bw.I2toP = .2;
% % # projections times strength
% Ns.ANtoI2 = 23.1;
% Ns.ANtoW = 10.8;
% Ns.ANtoP = 2.75;
% Ns.WtoI2 = 7;
% Ns.WtoP = 0;
% Ns.I2toP = 47.25;
% % non-specific afferents (NSA) leading to spontaneous activity
% nsa = 0; % spikes/s

%% Calculations

Nf = length(fc);
% df = mean(diff(fc));

W = nan(size(an));
for ii = 1:Nf
  nANtoW = fc >= fc(ii)*2^(-bw.ANtoW/2) & fc <= fc(ii)*2^(bw.ANtoW/2);
  win = gausswin(sum(nANtoW)); win = win/sum(win);
  W(ii,:) = Ns.ANtoW * win' * an(nANtoW,:);
% 	W(ii,:) = Ns.ANtoW*mean(an(nANtoW,:));
end
% W = max(W,0);

I2 = nan(size(an));
for ii = 1:Nf
  
  nANtoI2 = fc >= fc(ii)*2^(-bw.ANtoI2/2) & fc <= fc(ii)*2^(bw.ANtoI2/2);
  win = gausswin(sum(nANtoI2)); win = win/sum(win);
  ANtoI2(ii,:) = Ns.ANtoI2 * win' * an(nANtoI2,:);
  
  nWtoI2 = fc >= fc(ii)*2^(c.WtoI2-bw.WtoI2/2) & fc <= fc(ii)*2^(c.WtoI2+bw.WtoI2/2);
  if sum(nWtoI2) == 0; nWtoI2(end) = 1; end % for high fc
  win = gausswin(sum(nWtoI2)); win = win/sum(win);
  WtoI2(ii,:) = Ns.WtoI2 * win' * W(nWtoI2,:);
  
  I2(ii,:) = ANtoI2(ii,:) - WtoI2(ii,:);
% 	I2(ii,:) = Ns.ANtoI2*mean(an(nANtoI2,:)) - Ns.WtoI2*mean(an(nWtoI2,:));
end
% I2 = max(I2,0);

ANtoP = nan(size(an));
WtoP = ANtoP;
I2toP = ANtoP;
P = nan(size(an));
for ii = 1:Nf
  nANtoP = fc >= fc(ii)*2^(-bw.ANtoP/2) & fc <= fc(ii)*2^(bw.ANtoP/2);
  win = gausswin(sum(nANtoP)); win = win/sum(win);
  ANtoP(ii,:) = Ns.ANtoP * win' * an(nANtoP,:);
  
  nWtoP = fc >= fc(ii)*2^(c.WtoP-bw.WtoP/2) & fc <= fc(ii)*2^(c.WtoP+bw.WtoP/2);
  win = gausswin(sum(nWtoP)); win = win/sum(win);
  WtoP(ii,:) = Ns.WtoP * win' * W(nWtoP,:);
  
  nI2toP = fc >= fc(ii)*2^(c.I2toP-bw.I2toP/2) & fc <= fc(ii)*2^(c.I2toP+bw.I2toP/2);
  win = gausswin(sum(nI2toP)); win = win/sum(win);
  I2toP(ii,:) = Ns.I2toP * win' * I2(nI2toP,:);
  
  P(ii,:) = nsa + ANtoP(ii,:) - WtoP(ii,:) - I2toP(ii,:);
% 	P(ii,:) = nsa + Ns.ANtoP*mean(an(nANtoP,:)) - Ns.WtoP*mean(W(nWtoP,:)) - Ns.I2toP*mean(I2(nI2toP,:));
end
% P = P/max(P(:));
% P(P<0.4) = nan;
% P = P(1:90,:);
% P = max(P,0);

end

function v = nanmean (X, varargin) 
% NANMEAN
% v = nanmean (X)

n = sum (~isnan(X), varargin{:});
n(n == 0) = NaN;
X(isnan(X)) = 0;
v = sum (X, varargin{:}) ./ n;

end
