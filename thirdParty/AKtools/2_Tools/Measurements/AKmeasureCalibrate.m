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

fprintf('\n********************* input level calibration *********************\n')

if islogical(s.calibrate) % -------------------- measured level calibration
    
    % check number of input channels
    if numel(s.calibrateChIn) ~= 1
        error('AKcalibrate:InputChannel', 'You can only use a single input channel for calibration')
    end
    
    T = 5*s.fs;
    
    input('Apply calibrator to microphone and press Enter')
    
    % record calibration signal
    if strcmpi(s.audio_engine, 'playrec')
        
        calibration = AKplayrec('rec', [], [], s.calibrateChIn, T);
        
    elseif strcmpi(s.audio_engine, 'pa_wavplay')
        
        % record the sweep
        %  pa_wavplayrecord(playbuffer,[playdevice],[samplerate], [recnsamples], [recfirstchannel], [reclastchannel], [recdevice], [devicetype])
        calibration = pa_wavrecord(s.calibrateChIn, s.calibrateChIn, T, s.fs, s.recordDevice);
        
    end
    
    % check for clipping in case of clipping reduce level
    if max(abs(calibration(:))) >= .9999
        error('AKcalibrate:Clipping', 'Input signal clipped. Please decrease the calibration level or the microphone pre-amplification gain')
    end
    
    % get time
    dateStr = [datestr(now, 'yyyy.mm.dd HH-MM-SS') 'h'];
    
    % detect input level
    s.calibrateLevelIn = rms(calibration) * sqrt(2);
    
    % get correction factor
    s.calibrate_amplitude_per_Pascal = s.calibrateLevelIn / ( 2e-5 * 10^(s.calibrationLevel/20) );
    s.calibrate_dBFs_per_Pascal      = 20*log10(s.calibrate_amplitude_per_Pascal);
    
    % save calibration measurement
    if data.calib
        save(fullfile(data.dir, 'Reference', ['calibration_' dateStr '.mat']), 'calibration')
        disp(['Recorded calibration signal saved to ''' ['calibration_' dateStr '.mat''']])
    end
    
    AKf(30,15)
    set(gcf, 'Name', 'Level Calibration', 'NumberTitle', 'off')
    subplot(2,1,1)
        AKp(calibration / s.calibrate_amplitude_per_Pascal, 't2d', 'fs', s.fs, 'unit', 'pa', 'xu', 's')
        hold on
        dashline([0 T/s.fs],  [2e-5 * 10^(s.calibrationLevel/20) 2e-5 * 10^(s.calibrationLevel/20)], 2, 2, 2, 2, 'r', 'linewidth', 2)
        dashline([0 T/s.fs], -[2e-5 * 10^(s.calibrationLevel/20) 2e-5 * 10^(s.calibrationLevel/20)], 2, 2, 2, 2, 'r', 'linewidth', 2)
        legend('time signal', 'desired level', 'location', 'SouthEast')
        title(['time signal - sensitivity: ' num2str(s.calibrate_amplitude_per_Pascal) ' amplitude per Pascal'])
    subplot(2,1,2)
        AKp(calibration / s.calibrate_amplitude_per_Pascal, 'et2d', 'fs', s.fs, 'unit', 'pa', 'dr', 20, 'xu', 's')
        hold on
        dashline([0 T/s.fs],  [s.calibrationLevel s.calibrationLevel], 2, 2, 2, 2, 'r', 'linewidth', 2)
        legend('time signal', 'desired level', 'location', 'SouthEast')
        title(['time signal - sensitivity: ' num2str(s.calibrate_dBFs_per_Pascal) ' dB(Fs) per Pascal'])
            
    if data.plot
        saveas(gcf, fullfile(data.dir, 'Plots', ['calibration_' dateStr '.jpg']))
        disp(['Plot saved to ' ['''calibration_' dateStr '.jpg''']])
    end
    if data.close
        close
    end
    
else % ------------------------------------------- manual level calibration
    
    s.calibrate_dBFs_per_Pascal      =  s.calibrate;
    s.calibrate_amplitude_per_Pascal = 10^(s.calibrate/20);
    
    fprintf('\nManually calibrated to a sensitivity of %2.2f dB Fs per 1 Pascal\n', s.calibrate_dBFs_per_Pascal)
    
end

% --------------------------- compensate influence of reference measurement
if strcmpi(s.referenceType, 'complex')
    % get level at calibration frequency
    tmp = AKdeconv(reference, s.sweep, s.fs, 'x_inv_dyn', x.dynamic);
    
    if size(tmp,1) < s.fs
        tmp(end+1:s.fs) = 0;
    end
    
    tmp = abs(fft(tmp));
    tmp = tmp( round( s.calibrationFrequency / (s.fs/size(tmp,1)) ) + 1 );
    
    s.calibrate_amplitude_per_Pascal = s.calibrate_amplitude_per_Pascal / tmp;
    s.calibrate_dBFs_per_Pascal      = 20*log10(s.calibrate_amplitude_per_Pascal);
    
    fprintf('\nInput calibration finished\n')
    fprintf('sensitivity of %2.2f dB Fs per 1 Pascal detected\n', s.calibrate_dBFs_per_Pascal)
    fprintf('Corrected for the influence of the reference measurent (%2.2f dB @ %i Hz)\n', 20*log10(tmp), s.calibrationFrequency)
else
    fprintf('\nInput calibration finished\n')
    fprintf('sensitivity of %2.2f dB Fs per 1 Pascal detected\n', s.calibrate_dBFs_per_Pascal)
end

fprintf('WARNING: RE-CALIBRATE AFTER CHANGING ANY SETTINGS IN THE INPUT CHAIN!\n\n')


clear dateStr calibration tmp ans