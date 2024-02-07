function obj = hohmann2002_clearstate(obj)
%HOHMANN2002_CLEARSTATE  Reset states of hohmann2002 filters
%   Usage: filter = hohmann2002_clearstate(filter)
%          fb = hohmann2002_clearstate(fb)
%
%   filter = HOHMANN2002_CLEARSTATE(filter) resets the states of the
%   filter created by HOHMANN2002_FILTER
%
%   fb = HOHMANN2002_CLEARSTATE(fb) resets the states of the
%   filterbank fb created by HOHMANN2002
% 
%   delay = HOHMANN2002_CLEARSTATE(delay) resets the states of the
%   delay created by HOHMANN2002_DELAY
%
%   synth = HOHMANN2002_CLEARSTATE(synth) resets the states of the
%   sinthesis filterbank created by HOHMANN2002_SYNTH
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/hohmann2002_clearstate.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified  
%   #Author   : Universitaet Oldenburg, tp (2002 - 2007)
%   #Author   : Piotr Majdak (2016)
%   Adapted from function gfb_*_clear_state

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if ~isfield(obj,'type'), error('Type of the object missing'); end
switch(obj.type)
  case 'gfb_Filter'
    obj.state = zeros(1, obj.gamma_order);
  case 'gfb_analyzer'
    for band = 1:length(obj.center_frequencies_hz)
        obj.filters(1, band).state = zeros(1, obj.filters(1, band).gamma_order);
    end
  case 'gfb_Delay'
    obj.memory = zeros(size(obj.memory));
  case 'gfb_Synthesizer'
    obj.delay.memory = zeros(size(obj.delay.memory));
  otherwise
    error('Unknown type of HOHMANN2002 filter object');    
end


