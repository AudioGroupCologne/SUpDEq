% [dataNormalized, values, gains] = AKnormalize(data, type, operation, channelwise, value, f, frac, fs, doPlot)
%
% can be used to normalize the max, mean, or rms of data to any value
% See AKnormalizeDemo.m for example usage
%
%
% I N P U T
% data         - time signal of size [N M C], where N is the number of
%                samples, M the number of measurements, and C the number of
%                channels
% type         - 'time' : normalizes the time signal to 'value'
%                'abs'  : normalizes the magnitude spectrum to 'value'
%                'dB'   : normalizes the dB magniute spectrum to 'value'
%                (default = 'time')
% operation    - 'max'  : finds the maximum of the absolute data and
%                         normalizes it to 'value' if type='time', OR
%                         finds the maximum of data (without the absolute)
%                         and normalizes it to 'value'
%                'mean' : normalizes the mean to 'value'
%                'rms'  : normalizes the route mean square to 'value'
%                (default is 'max')
% channelwise  - 'each' : normalize each measurement of a channel
%                         separately
%                'max'  : normalize the maximum across measurements per
%                         channel
%                'min'  : normalize the minimum across measurements per
%                         channle
%                'mean' : normalize the mean across measurements per
%                         channel
%                'rms'  : normalize the route mean square across measure-
%                         ments per channle.
%                (default is 'max')
% value        - normalizes to 'value' which can be a scalar of a vector
%                with C elements (C = number of channels). In this case
%                each channel is normalized to the corresponding value
%                (default is 0 dB for type='dB', and 1 otherwise)
% f            - two element vector specifiying upper and lower frequency 
%                bound [Hz] for normalization if type={'abs' 'dB'}, OR
%                scalar specifiying the center frequency [Hz] for
%                normalization (see 'frac', default is [0 fs/2])
% frac         - bandwith for normalization in fractional octaves, if 'f'
%                is a scalar, e.g. 3 will normalize within a third octave
%                around 'f'
%                (default = false)
% fs           - sampling frequency in Hz (default = 44100)
% doPlot       - show a plot of the results (default = false)
%
% O U T P U T
% dataNormalized - the normalized input data
% values         - values of each channel before normalization
% gains          - the gain that was applied for normalization
%
% 11/2016  -  fabian.brinkmann@tu-berlin.de
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
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function [dataNormalized, values, gains] = AKnormalize(data, type, operation, channelwise, value, f, frac, fs, doPlot)

% -------------------------------------------------- set default parameters
if ~exist('type', 'var')
    type = 'time';
end
if ~exist('operation', 'var')
    operation = 'max';
end
if ~exist('channelwise', 'var')
    channelwise = 'max';
end
if ~exist('value', 'var')
    if strcmpi(type, 'db')
        value = 0;
    else
        value = 1;
    end
end
if ~exist('doPlot', 'var')
    doPlot = false;
end
if ~exist('fs', 'var')
    fs = 44100;
end
if ~exist('frac', 'var')
    frac = false;
end
if ~exist('f', 'var')
    f = [0 fs/2];
end

% -------------------------------------------------------------- copy input
if doPlot
    dataRaw = data;
end
dataNormalized = data;

% ------------------------------------------------------------- check input
if size(value,1)>size(value,2)
    value = value';
end
if numel(value)~=1 && numel(value)~=size(data,2)
    error('AKnormalize:Input', '''value'' must be a scalar or a vector with as many elements as ''data'' has channels.')
end
if any(value==0) && any(strcmpi(type, {'time', 'abs'}))
    error('AKnormalize:Input', 'A ''value'' of 0 will results in division by 0 and ist not allowed for ''time'' and ''abs'' normalization')
end

% ----------------------------------------------- handle multi-channel data
if size(data,3) > 1
    % allocate space
    dataNormalized = zeros(size(data));
    values         = zeros(size(data,3), size(data,2));
    gains          = values;
    
    for cc = 1:size(data,3)
        [dataNormalized(:,:,cc), values(cc,:), gains(cc,:)] = AKnormalize(data(:,:,cc), type, operation, channelwise, value, f, frac, fs, doPlot);
    end
    return
end

% ------------------------------------ transform data to the desired domain
switch lower(type)
    case 'time'
        
    case 'abs'
        data = AKboth2singleSidedSpectrum( abs( fft(data) ) );
    case 'db'
        data = AKboth2singleSidedSpectrum( db( fft(data) ) );
    otherwise
        error('AKnormalize:Input', ['invalid ''type'' passed (' type ')'])
end

% -------------------------------------------- get bounds for normalization
switch lower(type)
    case 'time'
        ID = [1 size(dataNormalized,1)];
    otherwise
        
        % get the limits for searching the maxima
        if numel(f)==1
            % frequency spacing
            df = fs/size(dataNormalized,1);
            
            % get frequency limits
            f = f/2^(1/(2*frac));
            f = [f f*2^(1/frac)];
            
            % get corresponding indicees
            ID = round(f/df)+1;
        else
            % frequency spacing
            df = fs/size(dataNormalized,1);
            
            % get corresponding indicees
            ID = round(f/df)+1;
        end
end

% ----------------------------------------- get the value for normalization
switch lower(operation)
    case 'max'
        if strcmpi(type, 'time')
            values = max(abs(data(ID(1):ID(2),:)));
        else
            values = max(data(ID(1):ID(2),:));
        end
    case 'mean'
        values = mean(data(ID(1):ID(2),:));
    case 'rms'
        values = rms(data(ID(1):ID(2),:));
    otherwise
        error('AKnormalize:Input', ['operation must be ''abs'', ''mean'', or ''rms'' but is ''' operation ''''])
end

if strcmpi(type, 'db')
    valuesLin = 10.^(values/20);
    valueLin  = 10.^(value/20);
else
    valuesLin = values;
    valueLin  = value;
end

% ------------------------------------------------- apply the normalization
if strcmpi(channelwise, 'each')
    % normalize to 1
    dataNormalized = AKm(dataNormalized, valuesLin, '/');
elseif any(strcmpi(channelwise, {'max' 'min' 'mean' 'rms'}))
    % normalize to 1
    dataNormalized = AKm(dataNormalized, eval([channelwise '(valuesLin)']), '/');
else
     error('AKnormalize:Input', ['channelwise must be ''each'', ''max'', ''min'', ''mean'', or ''rms'' but is ''' channelwise ''''])
end

% normalize to value
dataNormalized = AKm(dataNormalized, valueLin, '*');

% save the applied gains
if strcmpi(channelwise, 'each')
    gains = 1 ./ valuesLin;
else
    gains = 1 ./ eval([channelwise '(valuesLin)']);
end

gains = AKm(gains, valueLin, '*');

% -------------------------------------------------------------------- plot
if doPlot
    if strcmpi(type, 'time')
        plotMode = 't2d';
    else
        plotMode = 'm2d';
    end
    
    AKf
    subplot(1,2,1)
        AKp(dataRaw, plotMode)
        title('data (not normalized)')
    subplot(1,2,2)
        AKp(dataNormalized, plotMode)
        if strcmpi(type, 'time')
            title('data (normalized)')
        else
            title('data (normalized, red lines show normalization range)')
        end
        
    if ~strcmpi(type, 'time')
        hold on
        plot([f(1) f(1)], get(gca, 'yLim'), '--', 'color', AKcolors('r'))
        plot([f(2) f(2)], get(gca, 'yLim'), '--', 'color', AKcolors('r'))
    end
    
end

% ------------------------------------------------- check if we have output
if ~nargout
    clear dataNormalized values
end
