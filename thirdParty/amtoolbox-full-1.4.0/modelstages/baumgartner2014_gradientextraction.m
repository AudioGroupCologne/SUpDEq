function [gp,gfc] = baumgartner2014_gradientextraction(mp,fc,varargin)
%BAUMGARTNER2014_GRADIENTEXTRACTION Extraction of positive spectral gradients
%   Usage:      [gp,gfc] = baumgartner2014_gradientextraction(mp,fc)
%
%   Input parameters:
%     mp      : spectral magnitude profile in dB
%     fc      : center frequencies
%
%   Output parameters:
%     gp      : positive spectral gradient profile
%     gfc     : center frequencies of gradient profile
%
%
%   BAUMGARTNER2014_GRADIENTEXTRACTION(...) is a spectral cue extractor 
%   inspired by functionality of dorsal cochlear nucleus in cats.
%
%   BAUMGARTNER2014_GRADIENTEXTRACTION accepts the following optional parameters:
%
%     'c2',c2   Inhibitory coupling between type II mpd type IV neurons. 
%               Default is 1.
%     'c4',c4   Inhibitory coupling between AN and type IV neuron. 
%               Default is 1.
%     'spacing',sp  Tonotopical spacing between type IV mpd II neurons in ERBs. 
%                   Default is 1 ERB.
%
%   BAUMGARTNER2014_GRADIENTEXTRACTION accepts the following flags:
%
%     'positive'  Perform positive spectral gradient extraction (default).
%     'negative'  Perform negative spectral gradient extraction.
%     'both'      Perform spectral gradient extraction.
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791--802, 2014.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2014_gradientextraction.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA CircStat M-SIGNAL M-Stats O-Statistics
%   #Author: Robert Baumgartner (2014), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Parameter Settings

definput.keyvals.c2 = 1;
definput.keyvals.c4 = 1;
definput.keyvals.spacing = 1;
definput.flags.gradient = {'positive','negative','both'};

[flags,kv]=ltfatarghelper({'c2','c4','spacing'},definput,varargin);

% if not(exist('c2','var'))
%   c2 = 1; % inhibitory coupling between type II mpd type IV neurons
% end
% c4 = 1; % coupling between AN and type IV neuron
% dilatation = 1; % of tonotopical 1-ERB-spacing between type IV mpd II neurons


%% Calculations
Nb = size(mp,1); % # auditory bands
erb = audfiltbw(fc);
dgpt2 = round(mean(erb(2:end)./diff(fc))*kv.spacing); % tonotopical distance between type IV mpd II neurons
gp = zeros(Nb-dgpt2,size(mp,2),size(mp,3),size(mp,4),size(mp,5)); % type IV output
for b = 1:Nb-dgpt2
  gp(b,:,:,:,:) = kv.c4 * mp(b+dgpt2,:,:,:,:) - kv.c2 * mp(b,:,:,:,:);
end

if flags.do_positive
  gp(gp<0) = 0; % gp = (gp + c2*abs(gp))/2;
elseif flags.do_negative
  gp(gp>0) = 0;
end

% gfc = fc(dgpt2+1:end); % same as AN
gfc = sqrt(fc(1:Nb-dgpt2).*fc(dgpt2+1:end)); % use geometric mean

end


