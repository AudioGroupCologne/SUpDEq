function [directionCues cueEnergies] = takanen2013_cueconsistency(directionCues, cueEnergies,fc)
%takanen2013_cueconsistency Check the consistency before the cues are combined                 
%   Usage: [directionCues cueEnergies] = takanen2013_cueconsistency(directionCues, cueEnergies,fc)
%
%   Input parameters:
%        directionCues : sturcture containing the directional cues obtained
%                        from the model of the MSO, LSO, and wideband MSO
%                        of the two hemispheres%                    
%        cueEnergies   : sturcture containing the energies of the directional
%                        cues 
%
%        fc            : characteristic frequencies
%
%   Output parameters:
%        directionCues : sturcture containing the directional cues obtained
%                        from the model of the MSO, LSO, and wideband MSO
%                        of the two hemispheres%                    
%        cueEnergies   : sturcture containing the energies of the directional
%                        cues 
%
%   The role of the interaural time shifts, envelope time shifts, and level
%   differences in localization are different in different frequency
%   ranges. In order to emulate this when the directional cues are
%   combined, the consistency of the cues is checked and the energies of
%   the different cues are affected. More detailed description about the
%   process can be found in Takanen, Santala, Pulkk 2013 (Sec. 3.2.4)
%
%   See also: takanen2013
%
%   References:
%     M. Takanen, O. Santala, and V. Pulkki. Visualization of functional
%     count-comparison-based binaural auditory model output. Hearing
%     research, 309:147-163, 2014. PMID: 24513586.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/takanen2013_cueconsistency.php

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

%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland


%% ------ Computation ----------------------------------------------------


% 1) the given directional cues of the two hemispheres should not both
% point to the opposite hemispheres, indicating a failure to detect the
% direction at that time instant.

%first the MSO model
cueEnergies.leftMso(((directionCues.leftMso<0).*(directionCues.rightMso<0))==1) =0;
cueEnergies.rightMso(((directionCues.rightMso<0).*(directionCues.leftMso<0))==1) =0;

directionCues.leftMso((directionCues.leftMso<0).*(directionCues.rightMso>...
    abs(directionCues.leftMso))==1)=0;
directionCues.rightMso((directionCues.rightMso<0).*(directionCues.leftMso>...
    abs(directionCues.rightMso))==1)=0;

%then the same for LSO
cueEnergies.leftLso(((directionCues.leftLso<0).*(directionCues.rightLso<0))==1) =0;
cueEnergies.rightLso(((directionCues.rightLso<0).*(directionCues.leftLso<0))==1) =0;

directionCues.leftLso((directionCues.leftLso<0).*(directionCues.rightLso>...
    abs(directionCues.leftLso))==1)=0;
directionCues.rightLso((directionCues.rightLso<0).*(directionCues.leftLso>...
    abs(directionCues.rightLso))==1)=0;

%and for wideband MSO
cueEnergies.leftWbMso(((directionCues.leftWbMso<0).*(directionCues.rightWbMso<0))==1) =0;
cueEnergies.rightWbMso(((directionCues.rightWbMso<0).*(directionCues.leftWbMso<0))==1) =0;

directionCues.leftWbMso((directionCues.leftWbMso<0).*(directionCues.rightWbMso>...
    abs(directionCues.leftWbMso))==1)=0;
directionCues.rightWbMso((directionCues.rightWbMso<0).*(directionCues.leftWbMso>...
    abs(directionCues.rightWbMso))==1)=0;

% 2) Wideband MSO cues are considered as valid in the low frequencis only if
%   the cue points more to the side than the narrowband MSO and LSO cues
limit = find(fc<1000,1,'last');
ind = (((directionCues.leftWbMso<directionCues.leftMso)+(directionCues.leftWbMso<directionCues.leftLso))...
    +(directionCues.leftWbMso<30))>0;
cueEnergies.leftWbMso(ind(:,1:limit))=0;
ind = (((directionCues.rightWbMso<directionCues.rightMso)+(directionCues.rightWbMso<directionCues.rightLso))...
    +(directionCues.rightWbMso<30))>0;
cueEnergies.rightWbMso(ind(:,1:limit))=0;

