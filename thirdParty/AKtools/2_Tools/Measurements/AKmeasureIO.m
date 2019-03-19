% This script is called by AKmeasureDemo from AKtools
% See AKmeasureDemo.m for examples

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License.

if strcmpi(s.audio_engine, 'playrec')
    
    if ~isfield(s, 'playDevice')
        s.playDevice = false;
    end
    if ~isfield(s, 'recordDevice')
        s.recordDevice = false;
    end
    if ~isfield(s, 'bufferSize')
        s.bufferSize = false;
    end
    
    [s.playDevice, s.recordDevice, s.bufferSize, s.fs] = AKplayrecSetup(s.playDevice, s.recordDevice, s.bufferSize, s.fs);
    
    
elseif strcmpi(s.audio_engine, 'pa_wavplay')
    
    if ~isfield(s, 'playDevice')
        s.playDevice = false;
    end
    if ~isfield(s, 'recordDevice')
        s.recordDevice = false;
    end
    
    [s.playDevice, s.recordDevice] = AKpa_wavplaySetup(s.playDevice, s.recordDevice);
    
end

if s.playDevice < 0 || s.recordDevice < 0
    error('AKmeasureIO:device','No play or record device selected.');
end

clear h
