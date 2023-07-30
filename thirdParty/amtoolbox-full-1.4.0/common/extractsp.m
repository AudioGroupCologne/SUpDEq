function varargout = extractsp(lat,varargin)
%EXTRACTSP Sagittal plane (SP) HRTFs from measurement data
%   Usage:    [sphrtfs,polangs] = extractsp( lat,hM,pos )
%             [sphrtfs,polangs] = extractsp( lat,Obj )
%             [sphrtfs,polangs,latangs] = extractsp( lat,hM,pos )
%             [sphrtfs,polangs,latangs,idx] = extractsp( lat,hM,pos )
%             [sphrtfs,polangs,latangs] = extractsp( lat,Obj )
%             [sphrtfs,polangs,latangs,idx] = extractsp( lat,Obj )
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
%     polangs : corresponding polar angles (deg)
%     latangs : corresponding actual lateral angles (deg)
%     idx     : index to the considered positions
%
%   EXTRACTSP(...) extracts all HRTFs available for a specific SP or
%   lateral angle. In order to result in a representative HRTF template,
%   demands are made on:
%
%   1) lateral tolerance to be as small as possible, but min. 2 and max. 5.
% 
%   2) polar angles to range from max. -30 deg. to min. 210 deg.,
%
%   3) gaps of polar angles to be max. 30 deg. large.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/extractsp.php


%   #Author: Robert Baumgartner, Acoustics Research Institute, Vienna, Austria
%   #Author: Piotr Majdak (2019)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

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

% assure polar-angle range to be within [-90,270[ (may occur due to
% numerical inaccuracy)
if sum(pos(:,7) >= 270) || sum(pos(:,7) < -90)
  pos(:,7) = mod(pos(:,7)+90,360)-90;
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
  [~,polidx]=unique(real(pos(idx,7)));  % sort polar angles
  pol=pos(idx(polidx),7);               % sorted polar angles
  latActual = pos(idx(polidx),6);       % actual lateral angles
  sphrtfs=double(hM(:,idx(polidx),:));  % sorted DTFs of SP

  dlat = dlat + 1;  % increase dlat

end


varargout{1} = sphrtfs;
if nargout >= 2
  varargout{2} = pol;
  if nargout >= 3
    varargout{3} = latActual;
    if nargout == 4
        varargout{4}=idx(polidx);   % index to the considered positions
    end
  end
end

end


