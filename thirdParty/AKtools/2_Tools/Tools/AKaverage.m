% data = AKaverage(data, averageMode, phaseMode, phaseParameter, averageWeights)
% Can be used to average impulse responses in different ways. Sometimes it
% might be helpfull to align the data before -> see AKalign() for this!
%
% e.g. use
% AKaverage(data, 'time') simply averages the impulse responses, while
% AKaverage(data, 'abs') averages the magnitude spectra, and
% AKaverage(data, 'power') the power spectra
%
%
% I N P U T:
% data           - single or multi channel impulse responses of size
%                  [N M C] where N is the number of samples, M the number
%                  of measurements, and C the number of channels.
% averageMode    - a string the defines the average mode
%                  'time'   : averages the impulse responses
%                             (this is the default)
%                  'complex': averages the complex spectra
%                  'abs'    : averages the magnitude spectra
%                  'power'  : averages the power spectra, i.e.
%                             abs(fft())^2. Note that the square root is
%                             taken after averaging!
%                  'db'     : averages the log. magnitude spectra
% phaseMode      - boolean, or string that specifies how the phase is
%                  averaged in case the averageMode is 'abs', 'power' or
%                  'db'. phaseMode requires additional parameters
%                  (phaseParameter, see below)
%                  false    : ignore the phase. This is the default and
%                             results in zero phase impulse responses
%                  'copy'   : copies the phase of one impulse response
%                             inside data to the averaged impulse response
%                  'unwrap' : averages across the unwrapped phase. Note
%                             that the squared phase is averaged, and the
%                             square root is taken afterwards if
%                             averageMode='power'
%                  'min'    : calculated the minimum phase for the averaged
%                             impulse response
%                  'lin'    : calculated the linear phase for the averaged
%                             impulse response
%                  'zero'   : see false
% phaseParameter - boolean, scalar, or two element vector that sets the
%                  parameters for the following phaseMode
%                  'copy'   : [m] integer that defines which measurement is
%                             used for copying the phase (default = 1)
%                  'unwrap' : [fs fAlign] averaging the unwrapped phase
%                             might produce errors if the phase is noisy at
%                             low frequencies. In this case the phase
%                             values can be aligned at the frequency fAlign
%                             in Hz within the range of 2pi. fs specifies
%                             the sampling frequency. The default=false
%                             does not align the phase before averaging
%                  'min'    : [fs NTTF_double] See AKphaseManipulation
%                             (default = [44100 2])
%                  'lin'    : [fs] See AKphaseManipulation
%                             (default = 44100);
% averageWeights - M element vector that gives weights for averaging the
%                  data (M = number of measuremts). Before applying, the
%                  weights are normalized for their sum to equal one.
%                  The default=false will equally average across
%                  measurements.
%
% O U T P U T:
% data - averaged impulse responses of size [N C], where N is the number of
%        samples, and C the number of channels from the input data
%
% 
% 12/2016 - fabian.brinkmann@tu-berlin.de

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
function data = AKaverage(data, averageMode, phaseMode, phaseParameter, averageWeights)

% --------------------------------------------------- set default parameter 
if ~exist('averageMode', 'var')
    averageMode = 'time';
end
if ~exist('phaseMode', 'var')
    phaseMode = false;
end
if ~exist('phaseParameter', 'var')
    phaseParameter = false;
end
if ~exist('averageWeights', 'var')
    averageWeights = false;
end


% --------------------------------- handle input with more than one channel
if size(data,3)>1
    % allocate space for output
    dataOut = zeros(size(data,1), size(data,3));
    
    % get average for each page of data
    for cc = 1:size(data,3)
        dataOut(:,cc) = AKaverage(data(:,:,cc), averageMode, phaseMode, phaseParameter, averageWeights);
    end
    
    % save output
    data = dataOut;
    return
end


% ------------------------------------------------------------- check input
if any(averageWeights)
    % check dimension of the weights
    if numel(averageWeights) ~= size(data,2)
        error('AKaverage:Input', 'averageWeights must have as many elements as data has columns.')
    end
    
    % normalize and correct dimension
    averageWeights = averageWeights ./ sum(averageWeights);
    averageWeights = reshape(averageWeights, [1 numel(averageWeights)]);
end

% -------------------------------------------------------- average the data
% (use mean for averageing if no weights should be applied, or apply
% weights and sum the data otherwise) 

switch lower(averageMode)
    case 'time'
        if any(averageWeights)
            data = sum(AKm(data, averageWeights, '*'), 2);
        else
            data = mean(data, 2);
        end
        
    case 'complex'
        if any(averageWeights)
            data = sum(AKm(fft(data), averageWeights, '*'), 2);
        else
            data = mean(fft(data), 2);
        end
        data = ifft(data, 'symmetric');
        
    case {'abs','db','power'}
        
        % get the averaged magnitude spectrum
        if strcmpi(averageMode, 'abs')
            if any(averageWeights)
                data_mag = abs(fft(data));
                data_mag = sum(AKm(data_mag, averageWeights, '*'), 2);
            else
                data_mag = mean(abs(fft(data)), 2);
            end
        elseif strcmpi(averageMode, 'power')
            if any(averageWeights)
                data_mag = abs(fft(data)).^2;
                data_mag = sum(AKm(data_mag, averageWeights, '*'), 2);
                data_mag = sqrt(data_mag);
            else
                data_mag = sqrt(mean(abs(fft(data)).^2, 2));
            end
        elseif strcmpi(averageMode, 'db')
            if any(averageWeights)
                data_mag = db(fft(data));
                data_mag = sum(AKm(data_mag, averageWeights, '*'), 2);
                data_mag = 10.^(data_mag./20);
            else
                data_mag = mean(db(fft(data)), 2);
                data_mag = 10.^(data_mag./20);
            end
        else
            error('AKaverage:Input', 'averageMode must be ''time'', ''complex'', ''abs'', or ''db''')
        end
        
        % get the phase response
        if strcmpi(phaseMode, 'copy')
            % default parameter
            if ~exist('phaseParameter', 'var')
                phaseParameter = 1;
            end
            
            % copy phase
            data_ang = angle(fft(data(:,phaseParameter)));
            
            % get the complex spectrum, and impulse response
            data = data_mag .* exp(1i .* data_ang);
            data = ifft(data, 'symmetric');
            
        elseif strcmpi(phaseMode, 'unwrap')
            % default parameter
            if ~exist('phaseParameter', 'var')
                phaseParameter = false;
            end
            
            % get unwrapped phase
            [data_ang, isEven] = AKboth2singleSidedSpectrum( unwrap( angle( fft( data ) ) ) );
            
            % align values at fAlign to the first channel
            if all(phaseParameter)
                fAlign        = phaseParameter(2);
                fs            = phaseParameter(1);
                fID           = round(fAlign/(fs/size(data,1))) + 1;
                phaseMismatch = data_ang(fID,1) - data_ang(fID,:);
                phaseMismatch = round(phaseMismatch ./ (2*pi)) * 2*pi;
                data_ang      = AKm(data_ang, phaseMismatch, '+');
            end
            
            % average phase responses
            if strcmpi(averageMode, 'power')
                if any(averageWeights)
                    data_ang = sum(AKm(data_ang.^2, averageWeights, '*'), 2);
                    data_ang = -sqrt(data_ang);
                else
                    data_ang = -sqrt(mean(data_ang.^2, 2));
                end
            else
                if any(averageWeights)
                    data_ang = sum(AKm(data_ang, averageWeights, '*'), 2);
                else
                    data_ang = mean(data_ang, 2);
                end
            end
            
            % get both sided spectrum
            data_ang = AKsingle2bothSidedSpectrum(data_ang, isEven);
            
            % get the complex spectrum, and impulse response
            data = data_mag .* exp(1i .* data_ang);
            data = ifft(data, 'symmetric');
            
        elseif any(strcmpi(phaseMode, {'min' 'lin'}))
            % default parameter
            if ~exist('phaseParameter', 'var')
                phaseParameter = [44100 2];
            elseif numel(phaseParameter) == 1
                phaseParameter = [phaseParameter 0];
            end
            
            data = AKphaseManipulation(ifft(data_mag, 'symmetric'), phaseParameter(1), phaseMode, phaseParameter(2), 0);
            
        elseif (islogical(phaseMode) && ~phaseMode) ...
                || strcmpi(phaseMode, 'zero')
            % get zero phase impulse response
            data = ifft(data_mag, 'symmetric');
            
        else
            error('AKaverage:Input', 'phaseMode must be ''copy'', ''unwrap'', ''min'', ''lin'', ''zero'' or false')
        end
        
        
    otherwise
        error('AKaverage:Input', 'averageMode must be ''time'', ''complex'', ''abs'', ''power'', or ''db''')
end
