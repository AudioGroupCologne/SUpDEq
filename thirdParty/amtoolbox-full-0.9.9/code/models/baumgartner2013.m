function varargout = baumgartner2013( target,template,varargin )
%BAUMGARTNER2013 Model for sound localization in sagittal planes
%   Usage:    pmv = baumgartner2013( target,template )
%             [pmv,respang,tang] = baumgartner2013(target,template,u,lat,flow,fhigh,stim,fsstim,polsamp)
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
%     template: listener-specific template. 
%               Options 1 & 2 equivalent to target*
%
%   Output parameters:
%     pmv     : predicted probability mass vectors for response angles 
%               with respect to target positions.
%               1st dim: response angle.
%               2nd dim: target angle.
%     respang : polar response angles (after regularization of angular 
%               sampling)
%     tang    : polar target angles in the given sagittal plane
%
%   BAUMGARTNER2013(...) is a model for sound-source localization in
%   sagittal planes (SPs). It bases on the comparison of internal sound 
%   representation with a template and results in a probabilistic
%   prediction of polar angle response.
%
%   BAUMGARTNER2013 accepts the following key/value pairs:
%
%     'fs',fs        Define the sampling rate of the impulse responses. 
%                    Default value is 48000 Hz.
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
%     'lat',lat      Set the perceived lateral angle of the target sound to
%                    lat. Default value is 0 deg (midsagittal plane).
%
%     'u',u          Set the listener-specific uncertainty (standard
%                    deviation of the Gaussian transformation from the
%                    distance metric of the comparison process to the
%                    similarity index) to u. Default value is 2 dB.
%
%     'space',s      Set spacing of auditory filter bands to s numbers of
%                    equivalent rectangular bandwidths (ERBs). 
%                    Default value is 1 ERB.
%
%     'bwsteep',bws  Set the steepness factor bws of the sigmoid function 
%                    applied for binaural weighting of monaural similarity 
%                    indices. Default value is 13 deg.
%
%     'polsamp',ps   Define the the polar angular sampling of the current
%                    SP. As default the sampling of ARI's HRTF format along
%                    the midsagittal plane is used, i.e.,
%                    ps = [-30:5:70,80,100,110:5:210].
%
%   BAUMGARTNER2013 accepts the following flags:
%
%     'gammatone'    Use the Gammatone filterbank for peripheral processing. 
%                    This is the default.
%
%     'langendijk2002_spectralanalysis'        Use a filterbank approximation based on DFT with 
%                    constant relative bandwidth for peripheral processing. 
%                    This was used by Langendijk and Bronkhorst (2002).
%
%     'ihc'          Incorporate the transduction model of inner hair 
%                    cells used by Dau et al. (1996). This is the default.
%
%     'noihc'        Do not incorporate the IHC stage.
%
%     'regular'      Apply spline interpolation in order to regularize the 
%                    angular sampling of the polar response angle. 
%                    This is the default.
%
%     'noregular'    Disable regularization of angular sampling.
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2013
%
%   See also: plot_baumgartner2013, data_baumgartner2013,
%             exp_baumgartner2013
%
%   Demos: demo_baumgartner2013
%
%   References:
%     P. Majdak, R. Baumgartner, and B. Laback. Acoustic and non-acoustic
%     factors in modeling listener-specific performance of sagittal-plane
%     sound localization. Front Psychol, 5(319):pages not available yet,
%     doi:10.3389/fpsyg.2014.00319, 2014.
%     
%     R. Baumgartner. Modeling sagittal-plane sound localization with the
%     application to subband-encoded head related transfer functions.
%     Master's thesis, University of Music and Performing Arts, Graz, June
%     2012.
%     
%     R. Baumgartner, P. Majdak, and B. Laback. Assessment of Sagittal-Plane
%     Sound Localization Performance in Spatial-Audio Applications,
%     chapter 4, pages 93-119. Springer-Verlag GmbH, 2013.
%     
%     T. Dau, D. Pueschel, and A. Kohlrausch. A quantitative model of the
%     effective signal processing in the auditory system. I. Model structure.
%     J. Acoust. Soc. Am., 99(6):3615-3622, 1996a.
%     
%     E. Langendijk and A. Bronkhorst. Contribution of spectral cues to human
%     sound localization. J. Acoust. Soc. Am., 112:1583-1596, 2002.
%     
%     R. Patterson, I. Nimmo-Smith, J. Holdsworth, and P. Rice. An efficient
%     auditory filterbank based on the gammatone function. APU report, 2341,
%     1987.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/baumgartner2013.php

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

%% Notification
amt_disp('Note that baumgartner2013 is a preliminary version of baumgartner2014. If you do not aim to reproduce results from Baumgartner et al. (2013), we recommend to use baumgartner2014 instead.')

%% Check input options 

definput.flags.fbank = {'gammatone','langendijk2002_spectralanalysis','lopezpoveda2001','zilany2007'};
definput.flags.headphonefilter = {'','headphone'};
definput.flags.middleearfilter = {'','middleear'};
definput.flags.ihc = {'ihc','noihc'};
definput.flags.regularization = {'regular','noregular'};

% CP-Falgs:
definput.flags.cp={'std','xcorr'};

definput.keyvals.u=2;           % listener-specific uncertainty in dB
definput.keyvals.lat=0;         % deg
definput.keyvals.flow=700;      % Hz
definput.keyvals.fhigh=18000;   % Hz
definput.keyvals.stim=[];
definput.keyvals.fsstim=[];
definput.keyvals.polsamp=[-30:5:70 80 100 110:5:210];  % polar sampling (for regularization)
definput.keyvals.fs=48000;      % Hz
definput.keyvals.space=1;       % No. of ERBs (Cams)
definput.keyvals.do=0;
definput.keyvals.lvlstim = 40; 	% dBSPL
definput.keyvals.lvltem = 40;  	% dBSPL
definput.keyvals.bwsteep=13;    % steepness in degrees of binaural weighting function

[flags,kv]=ltfatarghelper(...
  {'u','lat','flow','fhigh','stim','fsstim','polsamp',...
  'fs','space','do','lvlstim','lvltem','bwsteep'},definput,varargin);

if isstruct(target) % Targets given in SOFA format
  kv.fs = target.Data.SamplingRate;
  target = extractsp( kv.lat,target );
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
if isempty(kv.stim) 
    kv.stim = [1;0];%[1;zeros(size(target,1),1)];    % impulse
    kv.fsstim = kv.fs;
elseif isempty(kv.fsstim) 
    kv.fsstim = kv.fs;
end

if flags.do_headphone% || flags.do_drnl
    hpfilt = headphonefilter(kv.fs);
    kv.stim = convolve(kv.stim,hpfilt(:));
end

if flags.do_middleear% || flags.do_drnl
    miearfilt = middleearfilter(kv.fs);
    kv.stim = convolve(kv.stim,miearfilt(:));
end


%% DTF filtering
if ~isequal(kv.fs,kv.fsstim)
    error('Sorry, sampling rate of stimulus and HRIRs must be equal!')
end

tmp = convolve(target,kv.stim);
target = reshape(tmp,[size(tmp,1),size(target,2),size(target,3)]);


%% Cochlear filter bank -> internal representations
if flags.do_langendijk2002_spectralanalysis
    
    bpo = kv.space*6; % bands per octave (1 oct. approx. as 6 ERBs)
    ireptem = langendijk2002_spectralanalysis(template,kv.fs,kv.flow,kv.fhigh,bpo);
    ireptar = langendijk2002_spectralanalysis(target,kv.fs,kv.flow,kv.fhigh,bpo);

elseif flags.do_gammatone

    % Determine filterbank
    fc = audspacebw(kv.flow,kv.fhigh,kv.space,'erb');
    Nfc = length(fc);   % # of bands
    [bgt,agt] = gammatone(fc,kv.fs,'classic');
    bgt = real(bgt);  % to ensure octave compatibility
    agt = real(agt);  % to ensure octave compatibility
    
    % Filtering
    ireptar = ufilterbankz(bgt,agt,target(:,:)); % channel (3rd) dimension resolved!
    ireptem = ufilterbankz(bgt,agt,template(:,:)); % resolve 3rd dim
    
    % IHC transduction
    if flags.do_ihc
        ireptar = ihcenvelope(ireptar,kv.fs,'ihc_dau');
        ireptem = ihcenvelope(ireptem,kv.fs,'ihc_dau');
    end
    
    % Set back the channel dimension
    ireptar = reshape(ireptar,[size(target,1),... 
        Nfc,size(target,2),size(target,3)]);
    ireptem = reshape(ireptem,[size(template,1),Nfc,size(template,2),size(template,3)]);
    
    % Averaging over time (RMS)
    ireptar = 20*log10(squeeze(rms(ireptar)));      % in dB!
    ireptem = 20*log10(squeeze(rms(ireptem)));
        
elseif flags.do_drnl   
    
    % Set level
    idnztar = target~=0;    % to ignore pausings
    idnztem = template~=0;  % to ignore pausings
    for ch = 1:size(template,3)
        target(idnztar(:,:,ch)) = setdbspl(target(idnztar(:,:,ch)),kv.lvlstim);
        template(idnztem(:,:,ch)) = setdbspl(template(idnztem(:,:,ch)),kv.lvltem);
    end
    
    % Filtering
    [ireptar,fc] = lopezpoveda2001(target(:,:),kv.fs,'flow',kv.flow,'fhigh',kv.fhigh);  % includes middle ear
    ireptem = lopezpoveda2001(template(:,:),kv.fs,'flow',kv.flow,'fhigh',kv.fhigh);
       
    % IHC transduction
    if flags.do_ihc 
        ireptar = ihcenvelope(ireptar,kv.fs,'ihc_dau');
        ireptem = ihcenvelope(ireptem,kv.fs,'ihc_dau');
    end
    
    % Set back the channel dimension
    ireptar = reshape(ireptar,[size(ireptar,1),... 
        length(fc),size(target,2),size(target,3)]);
    ireptem = reshape(ireptem,[size(ireptem,1),length(fc),size(template,2),size(template,3)]);
    
    % Averaging over time (RMS)
    ireptar = 20*log10(squeeze(rms(ireptar)));
    ireptem = 20*log10(squeeze(rms(ireptem)));

elseif flags.do_zilany2007
    
    fsmod = 100e3;  % Model sampling frequency in Hz
    nf = 200;       % # of AN fibers
  
    fprintf('\n compute internal representation of target set: \n');
    target = [target ; zeros(1e3,size(target,2),size(target,3))]; % concatenate zeros, otherwise comp_zilany2007 complaines about: "reptime should be equal to or longer than the stimulus duration." 
    ireptar = zeros(nf,size(target,2),size(target,3));
    nt = size(target(:,:),2);
    for ii = 1:nt
      [ANout,vfreq] = zilany2007(kv.lvlstim,target(:,ii),kv.fs,...
        fsmod,'flow',kv.flow,'fhigh',kv.fhigh,'nfibers',nf);
      ANout = ANout'-50; % subtract 50 due to spontaneous rate
      ireptar(:,ii) = 20*log10(rms(ANout)); % integrate over time & in db
      fprintf(' %2u of %2u done\n',ii,nt);
    end
  
    fprintf('\n Compute internal representation of template: \n');
    template = [template ; zeros(1e3+1,size(template,2),size(template,3))]; % concatenate zeros, otherwise comp_zilany2007 complaines about: "reptime should be equal to or longer than the stimulus duration." 
    ireptem = zeros(nf,size(template,2),size(template,3));
    nt = size(template(:,:),2);
    for ii = 1:nt
      ANout = zilany2007(kv.lvltem,template(:,ii),kv.fs,...
        fsmod,'flow',kv.flow,'fhigh',kv.fhigh,'nfibers',nf);
      ANout = ANout'-50; % subtract 50 due to spontaneous rate
      ireptem(:,ii) = 20*log10(rms(ANout)); % integrate over time & in dB
      fprintf(' %2u of %2u done\n',ii,nt);
    end
    
end

if size(ireptar,2) ~= size(target,2) % retreive polar dimension if squeezed out
    ireptar = reshape(ireptar,[size(ireptar,1),size(target,2),size(target,3)]);
end


%% Comparison process -> monaural similarity indices (SIs)
si=zeros(size(template,2),size(target,2),size(template,3)); % initialisation
for it = 1:size(target,2)
	si(:,it,:) = langendijk2002_comp(ireptar(:,it,:),ireptem,'s',kv.u, ...
    'argimport',flags,kv);
end


%% Binaural weighting -> binaural SIs
if size(si,3) == 2
    binw = 1./(1+exp(-kv.lat/kv.bwsteep)); % weight of left ear signal with 0 <= binw <= 1
    si = binw * si(:,:,1) + (1-binw) * si(:,:,2);
end


%% Interpolation (regularize polar angular sampling)
if flags.do_regular
    respang0 = ceil(min(kv.polsamp)*0.2)*5;    % ceil to 5°
    respangs = respang0:5:max(kv.polsamp);
    siint = zeros(length(respangs),size(si,2));
    for tt = 1:size(si,2)
        siint(:,tt) = interp1(kv.polsamp,si(:,tt),respangs,'spline');
    end
    si = siint;
    si(si<0) = 0; % SIs must be positive (necessary due to spline interp)
else
    respangs = kv.polsamp;
end


%% Normalization to PMV
pmv = si ./ repmat(sum(si),size(si,1),1);


%% Output
varargout{1} = pmv;
if nargout >= 2
    varargout{2} = respangs;
  if nargout >= 3
    varargout{3} = kv.polsamp; % tang
  end
end
  
end
