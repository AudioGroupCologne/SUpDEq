% IR_cut = AKcutIRmixingTime(IR, fs, stop_time, onset_threshold_dB, peak_secure_margin)
% Cuts the impulse response (IR) from onset position to stop_time for the
% estimation of the perceptual mixing time from data beased predictors
%
% see AKperceptualMixingTimeDemo.m for examples
%
% I N P U T:
% IR                 - impulse response of size [N C] where N is the number
%                      of samples, and C the number of channels
% fs                 - sampling rate in Hz
% stop_time          - stopping time in ms, if = 0 apply no shortening
% onset_threshold_dB - peak criterion (onset_threshold*maximum)
%
% O U T P U T:
% IR_cut             - IR from onsetposition to stop_time
%
%
% A. Lindau, L. Kosanke, 2011
% alexander.lindau@tu-berlin.de
% audio communication group
% Technical University of Berlin

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
%------------------------------------------------------------------------
function IR_cut = AKcutIRmixingTime(IR, fs, stop_time, onset_threshold_dB, peak_secure_margin)

% calculate linear value of onset threshold from dB value
onset_threshold=10^(onset_threshold_dB/20);

% for all IR-channels, find position where onset threshold value is reached
go_on=1;
k=0;

for i=1:size(IR,2)
    MAX=max(abs(IR(:,i)));
    % for full lenght of channel, find peak position
    while go_on
        k=k+1;
        
        % speichere beginn der IR je Kanal in "del" ab
        if abs(IR(k,i)) > MAX*(onset_threshold)
            del(i)=k; %#ok<AGROW>
            go_on=0;
        end
    end
    go_on=1;k=0;
end

% convert stop_time in [ms] to samples and shorten
if stop_time ~= 0
    
    stop_time = floor(stop_time/1000*fs); % IR length from peak to stop_time in samples
    
    if min(del) <= peak_secure_margin
        peak_secure_margin = 0; % ignore peak_secure_margin
    end
    
    % allocate space
    IR_cut = zeros((min(del)-peak_secure_margin+stop_time) - (min(del)-peak_secure_margin) + 1, size(IR,2));
    
    for j=1:size(IR,2)
        IR_cut(:,j) = IR(min(del)-peak_secure_margin:min(del)-peak_secure_margin+stop_time,j);
    end
 
else % ... no shortening
    
    if min(del) <= peak_secure_margin
        peak_secure_margin = 0; % ignore peak_secure_margin
    end
    
    stop_time = length(IR)-(min(del)-peak_secure_margin); % IR length from peak to end in samples
    
    % allocate space
    IR_cut = zeros((min(del)-peak_secure_margin+stop_time) - (min(del)-peak_secure_margin) + 1, size(IR,2));
       
    for j=1:size(IR,2)
        IR_cut(:,j) = IR(min(del)-peak_secure_margin:min(del)-peak_secure_margin+stop_time,j);
    end
    
end