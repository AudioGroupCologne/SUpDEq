% [rec, startTime] = AKplayrec(playrecMode, playBuffer, playChanList, recChanList, recDuration)
% 
% manages audio input/output using playrec externals compiled from
% http://www.playrec.co.uk
%
% you must call AKplayrecSetup before using AKplayrec
% See AKmeasureDemo.m for examples
% 
% play audio:
% AKplayrec('play', playBuffer, playChanList);
%
% record audio:
% rec = AKplayrec('rec', [], [], recChanList, recDuration);
%
% play and record audio
% rec = AKplayrec(playrecMode, playBuffer, playChanList, recChanList, recDuration);
%
% I N P U T:
% playMode     - 'play'    for audio playback
%                'rec'     for audio recording
%                'playrec' for audio playback and recording
% playBuffer   - audio time signal for playback of size [N C] where N is
%                the duration in samples and C the number of channels. If C
%                does not match the number of channels provided in
%                playChanList, the first channel is copied to all channels
%                given in playChanList
% playChanList - vector specifying the channels for audio playback
%                e.g. [1 3] to play audio to the first and third channel
% recChanList  - vector specifiying the channels for audio recording
% recDuration  - duration for recording audio in samples. If recDuration is
%                not passed for the 'playrec' mode, it is set to the length
%                of the playBuffer
%
% O U T P U T
% rec          - recorded samples of size [M C], where M equals the
%                recDuration and C the number of recording channels
%                specified in recChanList
% startTime    - the the playback/recording started. This is obtained from
%                matlabs 'now' command right after starting playrec. Use
%                datestr(startTime,'YYYY-mm-dd hh:MM:SS.FFF') to see the
%                start time including full data and millisenconds
%
% 10/2016  - fabian.brinkmann@tu-berlin.de
%

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

function [rec, startTime] = AKplayrec(playrecMode, playBuffer, playChanList, recChanList, recDuration)

switch lower(playrecMode)
    case 'play'
        
        % check inputBuffer size
        if size(playBuffer, 2) ~= numel(playChanList)
            playBuffer = repmat(playBuffer(:,1), [1 numel(playChanList)]);
        end
        
        % start playback
        pageNum = playrec(playrecMode, playBuffer, playChanList);
        
        % get starting time
        startTime = now;
        
        % reset missed samples count
        playrec('resetSkippedSampleCount')
        
        % check for missed samples
        while ~playrec('isFinished', pageNum)
            missedSamples = playrec('getSkippedSampleCount');
            pause(.1)
        end
        
        % no output in playback only mode
        rec = [];
        
        % display results
        if missedSamples
            fprintf('\n%u samples missed during playback! Try to change the page size or re-run AKplayrecSetup if this occurs again.\n\n', missedSamples)
        else
            fprintf('\nplayback successful\n\n')
        end
        
    case 'rec'
        
        % start recording and reset missedSample counter
        pageNum = playrec(playrecMode, recDuration, recChanList);
        
        % get starting time
        startTime = now;
        
        % reset missed samples count
        playrec('resetSkippedSampleCount')
        
        % check for missed samples
        while ~playrec('isFinished', pageNum)
            missedSamples = playrec('getSkippedSampleCount');
            pause(.1)
        end
        
        % get recorded audio
        rec = playrec('getRec', pageNum);
        
        % display results
        if missedSamples
            fprintf('\n%u samples missed during recording! Try to change the page size or re-run AKplayrecSetup if this occurs again.\n\n', missedSamples)
        else
            fprintf('\nrecording successful\n\n')
        end
        
    case 'playrec'
        
        % check record duration
        if ~exist('recDuration', 'var')
            recDuration = -1;
        end
        
        % check inputBuffer size
        if size(playBuffer, 2) ~= numel(playChanList)
            playBuffer = repmat(playBuffer(:,1), [1 numel(playChanList)]);
        end
        
        % start playback and recording
        pageNum = playrec(playrecMode, playBuffer, playChanList, recDuration, recChanList);
        
        % get starting time
        startTime = now;
        
        % reset missed samples count
        playrec('resetSkippedSampleCount')
        
        while ~playrec('isFinished', pageNum)
            missedSamples = playrec('getSkippedSampleCount');
            pause(.1)
        end
        
        rec = playrec('getRec', pageNum);
        
        if missedSamples
            fprintf('\n%u samples missed during playback and recording! Try to change the page size or re-run AKplayrecSetup if this occurs again.\n\n', missedSamples)
        else
            fprintf('\nplayback & recording successful\n\n')
        end
        
end

% clear output if not wanted
if ~nargout
    clear rec
end
