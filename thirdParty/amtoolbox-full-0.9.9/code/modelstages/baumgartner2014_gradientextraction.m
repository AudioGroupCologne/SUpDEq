function [gp,gfc] = baumgartner2014_gradientextraction(mp,fc,varargin)
%baumgartner2014_gradientextraction - Extraction of positive spectral gradients
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
%   BAUMGARTNER2014_GRADIENTEXTRACTION(...) is a spectral cue extractor
%    inspired by functionality of dorsal cochlear nucleus in cats.
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
%     Acoustical Society of America, 136(2):791-802, 2014.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2014_gradientextraction.php

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

% AUTHOR: Robert Baumgartner

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
