function E = baumgartner2021_mapping(d,S_cue,Erange,Eoffset,varargin)
%BAUMGARTNER2021_MAPPING externalization mapping function of baumgartner2021 model
%   Usage:    E = baumgartner2021_mapping(d,S_cue,Erange,Eoffset)
%
%   Input parameters:
%     d:        deviation metric
%
%     S_cue:    cue-specific sensitivity parameter (inverse slope) 
%
%     Erange:   range of rating scale
%
%     Eoffset:  scale offset
%
%   Output parameters:
%     E:        externalization score within rating scale
%
%   Description: 
%   ------------
%
%   BAUMGARTNER2021_MAPPING(...) represents a sigmoidal mapping function 
%   scaled by?Erange, shifted by?Eoffset, and slope-controlled by S_cue 
%   used to map deviation metrics?to externalization ratings within the
%   baumgartner2021 model.
%
%   BAUMGARTNER2021_MAPPING accepts the following flags:
%
%     'single'       Single-sided externalization ratings with respect to 
%                    reference sound. This is the default.
%
%     'two'          Two-sided externalization ratings.
%
%   See also: baumgartner2021 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2021_mapping.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: SOFA M-Stats O-Statistics
%   #Author: Robert Baumgartner (2021), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% AUTHOR : Robert Baumgartner

definput.flags.sided={'single','two'};
[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_single % sided
  Erange = Erange*2;
end

E = Erange./(1+exp(d./S_cue))+Eoffset;


