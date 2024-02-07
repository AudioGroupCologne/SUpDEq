function [gp,gfc] = baumgartner2016_gradientextraction(mp,fc,varargin)
%BAUMGARTNER2016_GRADIENTEXTRACTION Extraction of positive spectral gradients
%   Usage:      [gp,gfc] = baumgartner2016_gradientextraction(mp,fc)
%
%   Input parameters:
%     mp      : discharge rate profile
%     fc      : center frequencies
%
%   Output parameters:
%     gp      : positive spectral gradient profile. Fields: gp.m for
%               magnitude and gp.sd for standard deviation. 
%               Dimensions (4-6 optional): 
%               1) frequency, 2) position (polar angle), 3) channel (L/R), 
%               4) fiber type, 5) time frame.
%     gfc     : center frequencies of gradient profile
%
%   BAUMGARTNER2016_GRADIENTEXTRACTION(...) is a spectral cue extractor
%   inspired by functionality of dorsal cochlear nucleus in cats.
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2016_gradientextraction.php


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


definput.import={'baumgartner2016'};
definput.flags.gradients={'positive','negative','both'};
definput.keyvals.c2 = 1;
[flags,kv]=ltfatarghelper({'c2'},definput,varargin);


%% Parameter Settings
c2 = kv.c2; % inhibitory coupling between type II mpd type IV neurons
c4 = 1; % coupling between AN and type IV neuron
dilatation = 1; % of tonotopical 1-ERB-spacing between type IV mpd II neurons

erb = audfiltbw(fc);

%% Calculations
Nb = size(mp,1); % # auditory bands
dgpt2 = round(mean(erb(2:end)./diff(fc))*dilatation); % tonotopical distance between type IV mpd II neurons
mpsd = 2.6 * mp.^0.34; % variability of discharge rate (May and Huang, 1997)
gp.m = zeros(Nb-dgpt2,size(mp,2),size(mp,3),size(mp,4),size(mp,5)); % type IV output
gp.sd = gp.m;
for b = 1:Nb-dgpt2
  gp.m(b,:,:,:,:) = c4 * mp(b+dgpt2,:,:,:,:) - c2 * mp(b,:,:,:,:);
  gp.sd(b,:,:,:,:) = sqrt( (c4*mpsd(b+dgpt2,:,:,:,:)).^2 + (c2*mpsd(b,:,:,:,:)).^2 );
end

% Restriction to positive gradients
% hard restriction
% gp.m = (gp.m + c2*abs(gp.m))/2; % gp = max(gp,0);

% soft restriction
% kv.mgs = 10; % constant to stretch the atan
if flags.do_both
  gp.m = 2*kv.mgs*atan(gp.m/kv.mgs/2);
else
  if flags.do_positive
    gp.m = kv.mgs*(atan(gp.m/kv.mgs-pi/2)+pi/2);
  else %flags.do_negaitve
    gp.m = kv.mgs*(atan(gp.m/kv.mgs+pi/2)-pi/2);
  end
  gp.sd = gp.sd/2; % ROUGH APPROXIMATION assuming that non-linear restriction to positive gradients halfs the rate variability
end

gfc = fc(dgpt2+1:end);

end


