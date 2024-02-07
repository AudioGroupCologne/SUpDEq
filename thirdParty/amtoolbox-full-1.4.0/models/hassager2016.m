function Pext = hassager2016(target,template,varargin)
%HASSAGER2016   Sound externalization based on interaural level differences
%
%   Usage: Pext = hassager2016(target,template)
%
%   Input parameters:
%
%     target  : binaural signal for target sound(s).
%               matrix dimensions: time x left/right.
%     template: binaural reference signal.
%               matrix dimensions: time x left/right.
%
%   Output parameters:
%
%     Pext    : predicted degree of externalization
%
%   Optional input parameters:
%
%     'JND',jnd     ILD JND used as minimum limit for evaluation of ILD devations. Default is 1.5 [dB].
%
%     'c1',c1       Gain of ILD-deviation-to-externalization mapping function (decaying exponential) acc. to Eq. (14). Default is 3.78.
%
%     'c2',c2       Intercept of mapping function (Eq. (14)). Default is 1.
%
%     'z1',z1       Scaling parameter of mapping function (Eq. (14)). Default is 0.99.
%
%     'flow',flow   Lowest frequency considered in auditory filterbank. Default is 50 [Hz].
%
%   This model calculates the externalization based on interaural level
%   differences.
%
%   See also: data_hassager2016 sig_hassager2016
%   baumgartner2014
%   sig_li2020
%
%   References:
%     H. G. Hassager, F. Gran, and T. Dau. The role of spectral detail in the
%     binaural transfer function on perceived externalization in a
%     reverberant environment. J. Acoust. Soc. Am., 139(5):2992--3000, 2016.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/hassager2016.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Unknown   
%   #Author: Robert Baumgartner (2019), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.keyvals.JND = 1.5;
definput.keyvals.c1 = 3.78;
definput.keyvals.c2 = 1;
definput.keyvals.z1 = 0.99;
definput.keyvals.fs = 48e3;
definput.keyvals.flow = 50;
definput.keyvals.fhigh = 6000;
definput.keyvals.echoSupStart = 3.8e-3; % start time of echo suppresion
definput.keyvals.echoSupEnd = [10e-3,160e-3]; % fade-out time range of echo suppresion
definput.keyvals.echoSupAtten = 0.01; % suppression attenuation
definput.flags.echoSuppression = {'','echoSuppression'};
definput.flags.middleEarFilter = {'middleEarFilter',''};
[flags,kv]=ltfatarghelper({},definput,varargin);

%% Echo suppression (ES)
if flags.do_echoSuppression
  nSupStart = round(kv.echoSupStart*kv.fs); % start sample of echo suppresion
  nSupEnd = round(kv.echoSupEnd*kv.fs); % fade-out time range of echo suppresion

  NfadeOut = diff(nSupEnd); % # samples of fade out
  ESfadeOut = 0.5*(1-cos(2*pi*(0:NfadeOut)/(2*NfadeOut))); % raised-cosine ramp
  ESfadeOut = ESfadeOut*(1-kv.echoSupAtten)+kv.echoSupAtten; % scaled to supAtten
  ESwin = [kv.supAtten*ones(nSupEnd(1)-nSupStart,1);ESfadeOut(:)];
  
  Ntar = size(target,1);
  if Ntar >= nSupStart
    nEStar = nSupStart:min(nSupEnd(2),Ntar);
    target(nEStar,:,:) = target(nEStar,:,:).*repmat(ESwin,[1,size(target,2),size(target,3)]);
  end
  Ntem = size(template,1);
  if Ntem >= nSupStart
    nEStem = nSupStart:min(nSupEnd(2),Ntem);
    template(nEStem,:,:) = template(nEStem,:,:).*repmat(ESwin,[1,size(template,2),size(template,3)]);
  end
end

%% Middle ear filter
if flags.do_middleEarFilter
  b=middleearfilter(kv.fs);
  target = filter(b,1,target);
  template = filter(b,1,template);
end

%% Excitation patterns
[tar.mp,fc] = baumgartner2014_spectralanalysis(target,'flow',kv.flow,'fhigh',kv.fhigh);
tem.mp = baumgartner2014_spectralanalysis(template,'flow',kv.flow,'fhigh',kv.fhigh);
% figure; plot(fc,tar.mp); hold on; plot(fc,tem.mp); legend('tar-L','tar-R','tem-L','tem-R')

%% ILDs
tar.ild = -diff(tar.mp,1,2); % ILD = left - right
tem.ild = -diff(tem.mp,1,2);
% figure; plot(fc,tar.ild); hold on; plot(fc,tem.ild); legend('tar','tem')

%% target-template comparison -> ILD deviation
dILD = abs(tar.ild-tem.ild);
dILD(dILD < kv.JND) = 0; % limit minimum ILD difference according to JND

%% overall normalized ILD deviation
dILDnorm = mean(dILD./abs(tem.ild));

%% Externalization mapping
Pext = kv.c1*exp(-kv.z1*dILDnorm) +kv.c2;
end


