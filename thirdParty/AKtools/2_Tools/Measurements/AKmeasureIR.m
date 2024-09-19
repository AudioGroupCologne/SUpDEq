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

Nin  = numel(s.chIn);

if strcmpi(s.chMode, 'single')
    Nout = numel(s.chOut);
else
    Nout = 1;
end

% allocate memory for raw and final data
raw = zeros(2^x.NFFT, Nin, Nout);
ir  = raw;

for Nsrc = 1:Nout
    
    s.measurementTime = datestr(now);
    fprintf('\n********************* %s *********************\n', s.measurementTime)
    
    % get current source - output channel(s)
    if strcmpi(s.chMode, 'single')
        curSrc = s.chOut(Nsrc);
    else
        curSrc = s.chOut;
    end
    
    % record and get time data
    nAVG = 0;
    Nclip_reduction = 0;
    
    while nAVG < s.average
        % playback and record, scale according to level and clip reduction
        if strcmpi(s.audio_engine, 'playrec')
            
            tmp = AKplayrec('playrec', s.sweep*10^(s.levelOut_dB/20)*10^(-Nclip_reduction*s.clipReduction/20), curSrc, s.chIn);
            
        elseif strcmpi(s.audio_engine, 'pa_wavplay')
            
            % generate output buffer
            wavBuffer           = zeros(size(s.sweep,1), max(curSrc));
            wavBuffer(:,curSrc) = repmat(s.sweep, 1, numel(curSrc));
            
            % record the sweep
            % pa_wavplayrecord(playbuffer,[playdevice],[samplerate], [recnsamples], [recfirstchannel], [reclastchannel], [recdevice], [devicetype])
            tmp = pa_wavplayrecord(wavBuffer*10^(s.levelOut_dB/20)*10^(-Nclip_reduction*s.clipReduction/20),...
                       s.playDevice, s.fs, 0, min(s.chIn), max(s.chIn), s.recordDevice);
                   
            % get rid of extra channels, and resort
            if size(tmp, 2) > numel(s.chIn)
                tmp = tmp(:, s.chIn - min(s.chIn) + 1);
            end
            
        end
        
        % initialize/clear output data
        if nAVG == 0
            rec_raw = zeros(size(tmp));
        end
        
        % average
        rec_raw = rec_raw + tmp;
        
        % check for clipping in case of clipping reduce level and reset
        % average counter. increase average counter elsewise
        if max(abs(tmp(:))) >= .9999
            disp('Clipping: Auto decrease level...')
            Nclip_reduction = Nclip_reduction + 1;
            nAVG = 0;
        else
            nAVG = nAVG + 1;
        end
    end
    
    for Nrec = 1:Nin
        if Nrec == 1
            fprintf('Peak level:   ch. %2.i ',s.chOut(Nsrc));
        else
            fprintf('                     ');
        end
        fprintf('-> ch. %2.i    %7.2f  [dBFS]\n',...
            s.chIn(Nrec),db(max(abs(rec_raw(:,Nrec)))));
    end
    for Nrec = 1:Nin
        if Nrec == 1
            fprintf('RMS level:    ch. %2.i ',s.chOut(Nsrc));
        else
            fprintf('                     ');
        end
        fprintf('-> ch. %2.i    %7.2f  [dBFS]\n',...
            s.chIn(Nrec),db(rms(rec_raw(:,Nrec))));
    end
    % scale according to averages and clip reduction
    rec_raw = rec_raw / nAVG * 10^(Nclip_reduction*s.clipReduction/20);
    
    % save to output mat-file
    raw(:, :, Nsrc) = rec_raw;
    
    % deconvolution
    if strcmpi(s.referenceType, 'complex')
        ir(:,:,Nsrc) = AKdeconv(raw(:,:,Nsrc), reference, s.fs, 'x_inv_dyn', x.dynamic, 'deconv_type', s.deconvType);
    elseif strcmpi(s.referenceType, 'latency') || s.referenceLatency
        ir(:,:,Nsrc) = AKdeconv(raw(:,:,Nsrc), s.sweep, s.fs, 'x_inv_dyn', x.dynamic, 'deconv_type', s.deconvType, 'do_circshift', -s.referenceLatency);
    else
        ir(:,:,Nsrc) = AKdeconv(raw(:,:,Nsrc), s.sweep, s.fs, 'x_inv_dyn', x.dynamic, 'deconv_type', s.deconvType);
    end
    
    % apply level calibration
    if s.calibrate
        ir(:,:,Nsrc) = ir(:,:,Nsrc) / s.calibrate_amplitude_per_Pascal;
    end
    
    % apply subsonic filter
    if s.subSonic
        ir(:,:,Nsrc) = AKfilter(ir(:,:,Nsrc), 'hp', 20, zeros(1,Nin), s.fs, 4, 'butter');
    end
    
end

% cut
if s.N && s.N > size(ir,1)
    ir = ir(1:s.N,:,:);
end

% plot data
clear h

% set the physical data unit for plotting
if s.calibrate
    unit = 'Pa';
else
    unit = false;
end

for Nsrc = 1:Nout
    h(Nsrc).f = AKf(20,20); %#ok<SAGROW>
    set(h(Nsrc).f, 'CloseRequestFcn', 'disp(''Figures can be closed after saving'')') % prevent manually closing the figure before saving
    set(h(Nsrc).f, 'name', ['Output channel: ' num2str(s.chOut(Nsrc)) ' (leave open for saving)'], 'NumberTitle','off')
    subplot(2,1,1)
    AKp(ir(:,:,Nsrc), 'et2d', 'c', 'cyc', 'fs', s.fs, 'unit', unit)
    subplot(2,1,2)
    AKp(ir(:,:,Nsrc), 'm2d', 'c', 'cyc', 'fs', s.fs, 'unit', unit)
end

clear Nin Nout Nsrc rec_raw nAVG Nclip_reduction tmp wavBuffer curSrc sweep unit
