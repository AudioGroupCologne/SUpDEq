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

fprintf('\n********************* reference measurement *********************\n')

% record and get time data
clip            = 1;
Nclip_reduction = 0;
impulseOffset   = 64;

s.clipReductionReference = 6;

% generate output buffer
if strcmpi(s.referenceType, 'latency')
    referenceBuffer = AKdirac(s.fs, 1, impulseOffset);
elseif strcmpi(s.referenceType, 'complex')
    referenceBuffer = s.sweep;
else
    error('s.referenceType not valid. Must be false, ''latency'' or ''complex''')
end

if strcmpi(s.audio_engine, 'pa_wavplay')
    referenceBuffer = [zeros(size(referenceBuffer, 1), s.referenceChOut-1) referenceBuffer];
end
    
while clip == 1
    % playback and record, scale according to level and clip reduction

    if strcmpi(s.audio_engine, 'playrec')

        tmp = AKplayrec('playrec', referenceBuffer*10^(s.referenceLevelOut_dB/20)*10^(-Nclip_reduction*s.clipReductionReference/20), s.referenceChOut, s.referenceChIn);
            
    elseif strcmpi(s.audio_engine, 'pa_wavplay')
            
        % record the sweep
        %  pa_wavplayrecord(playbuffer,[playdevice],[samplerate], [recnsamples], [recfirstchannel], [reclastchannel], [recdevice], [devicetype])
        tmp = pa_wavplayrecord(referenceBuffer*10^(s.referenceLevelOut_dB/20)*10^(-Nclip_reduction*s.clipReductionReference/20),...
                               s.playDevice, s.fs, 0, s.referenceChIn, s.referenceChIn, s.recordDevice);
            
    end
        
    % check for clipping in case of clipping reduce level
    if max(abs(tmp(:))) >= .9999
        disp('Clipping: Auto decrease level...')
        Nclip_reduction = Nclip_reduction + 1;
        clip = 1;
    else
        clip =  0;
    end
end

% get time
dateStr = [datestr(now, 'yyyy.mm.dd HH-MM-SS') 'h'];
   
% remove referenceLevelOut
tmp = tmp * 10^(-s.referenceLevelOut_dB/20);  

if strcmpi(s.referenceType, 'complex')
    % save for deconvolution
    reference = tmp;
    
    % get IR from reference measurement
    tmp = AKdeconv(tmp, s.sweep, s.fs, 'x_inv_dyn', x.dynamic); 
elseif strcmpi(s.referenceType, 'latency')
    tmp = circshift(tmp, -impulseOffset);
    reference = tmp;
end

% truncate to half a second
tmp = tmp(1:ceil(.5*s.fs));
   
% find indice of the largest element in IR
[~,I] = max(tmp);
% get system latency
s.referenceLatency = I -1;
% display detected latency
disp(['Latency: ' num2str(s.referenceLatency) ' samples'])

% save reference measurement
if data.ref
    save(fullfile(data.dir, 'Reference', ['reference_' dateStr '_' s.referenceType '.mat']), 'reference')
    disp(['Reference saved to ''' ['reference_' dateStr '_' s.referenceType '.mat''']])
end

% plot reference
AKf(30,15)
set(gcf, 'Name', 'Reference measurement (not inverted)', 'NumberTitle', 'off')
AKpMulti(tmp, '1a', 'fs', s.fs)
if data.plot
    saveas(gcf, fullfile(data.dir, 'Plots', ['reference_' dateStr '_' s.referenceType '.pdf']))
    disp(['Plot saved to ' ['''reference_' dateStr '_' s.referenceType '.pdf''']])
end
if data.close
    close(gcf);
end
    
clear I Nclip_reduction referenceBuffer nAVG tmp impulseOffset dateStr clip

