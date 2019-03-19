function [gp,gfc] = baumgartner2016_gradientextraction(mp,fc,varargin)
%baumgartner2016_gradientextraction - Extraction of positive spectral gradients
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
%    inspired by functionality of dorsal cochlear nucleus in cats.
%
%   References:
%     R. Baumgartner, P. Majdak, and B. Laback. Modeling sound-source
%     localization in sagittal planes for human listeners. The Journal of the
%     Acoustical Society of America, 136(2):791-802, 2014.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2016_gradientextraction.php

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

definput.import={'baumgartner2016'};
definput.keyvals.c2 = 1;
[flags,kv]=ltfatarghelper({'c2'},definput,varargin);

% if isempty(varargin) || not(isempty(varargin)) && not(isstruct(varargin{1}))
%   definput.import={'baumgartner2016','amt_cache'};
%   definput.keyvals.c2 = 1;
%   [flags,kv]=ltfatarghelper({'c2'},definput,varargin);
% else % kv and flags directly transfered
%   kv = varargin{1};
%   flags = varargin{2};
%   kv.c2 = 1;
% end

%% Parameter Settings
% if not(exist('c2','var'))
  c2 = kv.c2; % inhibitory coupling between type II mpd type IV neurons
% end
c4 = 1; % coupling between AN and type IV neuron
dilatation = 1; % of tonotopical 1-ERB-spacing between type IV mpd II neurons

erb = audfiltbw(fc);

%% Calculations
Nb = size(mp,1); % # auditory bands
dgpt2 = round(mean(erb(2:end)./diff(fc))*dilatation); % tonotopical distance between type IV mpd II neurons
mpm = mp;
mpsd = 2.6 * mpm.^0.34; % variability of discharge rate (May and Huang, 1997)
gp.m = zeros(Nb-dgpt2,size(mp,2),size(mp,3),size(mp,4),size(mp,5)); % type IV output
gp.sd = gp.m;
for b = 1:Nb-dgpt2
  gp.m(b,:,:,:,:) = c4 * mpm(b+dgpt2,:,:,:,:) - c2 * mpm(b,:,:,:,:);
  gp.sd(b,:,:,:,:) = sqrt( (c4*mpsd(b+dgpt2,:,:,:,:)).^2 + (c2*mpsd(b,:,:,:,:)).^2 );
end

% Restriction to positive gradients
% hard restriction
% gp.m = (gp.m + c2*abs(gp.m))/2; % gp = max(gp,0);

% soft restriction
% kv.mgs = 10; % constant to stretch the atan
gp.m = kv.mgs*(atan(gp.m/kv.mgs-pi/2)+pi/2);
gp.sd = gp.sd/2; % ROUGH APPROXIMATION assuming that non-linear restriction to positive gradients halfs the rate variability

gfc = fc(dgpt2+1:end);

% if nargout > 1

end
