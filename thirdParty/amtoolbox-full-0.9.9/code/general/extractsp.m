function varargout = extractsp(lat,varargin)
%EXTRACTSP Sagittal plane (SP) HRTFs from measurement data
%   Usage:    [sphrtfs,polangs] = extractsp( lat,hM,pos )
%             [sphrtfs,polangs] = extractsp( lat,Obj )
%
%   Input parameters:
%     lat     : lateral angle of the SP
%     Obj     : HRIR Data in SOFA format
%     hM      : matrix containing head-related impulse responses in ARI
%               Format
%               Dimensions: time,position,channel 
%               (for more details see doc: HRTF format description)
%     pos     : source-position matrix referring to 2nd dimension of hM and 
%               formated acc. to meta.pos (ARI format).
%               6th col: lateral angle. 7th col: polar angle
%
%   Output parameters:
%     sphrtfs : all available HRTFs in the current SP, sorted acc. to
%               ascending polar angle
%     polangs : corresponding polar angles
%
%   EXTRACTSP(...) extracts all HRTFs available for a specific SP or
%   lateral angle. In order to result in a representative HRTF template,
%   demands are made on:
%
%     1) lateral tolerance to be as small as possible, but min. 2 and
%        max. 5.
% 
%     2) polar angles to range from max. -30 deg. to min. 210 deg.,
%
%     3) gaps of polar angles to be max. 30 deg. large.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/extractsp.php

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

if nargin == 2 || isstruct(varargin{1}) % SOFA input
  
  Obj = varargin{1};
	hM = permute(double(Obj.Data.IR),[3 1 2]);

  pos(:,1)=bsxfun(@times,Obj.SourcePosition(:,1),ones(Obj.API.M,1));
  pos(:,2)=bsxfun(@times,Obj.SourcePosition(:,2),ones(Obj.API.M,1));
  [pos(:,6), pos(:,7)]=sph2hor(pos(:,1),pos(:,2));
  
elseif nargin == 3 % aRI Format
  
  hM = varargin{1};
  pos = varargin{2};
  
else return
end
  

%% Extract SP

dlat = 1;      % initial lateral tolerance (+/-) in deg
pol = [0,0];   % initial polar angles
dx = 0.01;

while isempty(pol) || ...
        (min(pol) > -30+dx || max(pol) < 210-dx ... % ensure that important polar range is included
        || max(diff(pol))>30)...            % and gaps are <= 30deg
        && dlat <= 5                        % but keep tolerance smaller than 5 deg

  idx=find( pos(:,6)>=-(dlat+dx)/2+lat & pos(:,6)<=(dlat+dx)/2+lat );
  latActual = pos(idx,6); % actual lateral angles
  [pol,polidx]=unique(real(pos(idx,7)));   % sorted polar angles
  latActual = latActual(polidx);
  sphrtfs=double(hM(:,idx,:));  % unsorted DTFs of SP
  sphrtfs=sphrtfs(:,polidx,:);          % sorted DTFs of SP

  dlat = dlat + 1;  % increase dlat

end

% assure polar-angle range to be within [-90,270[ (may occur due to
% numerical inaccuracy)
if sum(pol >= 270) || sum(pol < -90)
  pol = mod(pol+90,360)-90;
  [pol,idsort] = sort(pol);
  sphrtfs = sphrtfs(:,idsort,:);
end

varargout{1} = sphrtfs;
if nargout >= 2
  varargout{2} = pol;
  if nargout == 3
    varargout{3} = latActual;
  end
end

end
