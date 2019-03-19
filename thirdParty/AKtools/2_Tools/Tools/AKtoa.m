% [toa, toaSmoothed] = AKtoa(data, toaMethod, toaParameter, upSample, smoothing, doPlot)
% estimated the time of arrival (TOA) in impulse responses either using a
% threshold based onset detection or by cross-correlation between the input
% IRs and their corresponding minimum phase version. To achieve sub-sample
% quantization of the TAO estimated, data can be up-sapled.
%
% e.g. use
% AKtoa(data, 'onset'), or
% AKtoa(data, 'cc'),
% to estimate the TOAs for the 10 times up-sampled impulse responses.
%
% For more examples see AKtoaDemo.m
%
% Note that cross-correlation requeires longer calculation times but might
% work better for anechoic impulse responses, while onset detection is more
% robust for reverberant cases. In any case, passing only a the part of the
% impulse responses that contains the onset will decrease computation time
% and often improve the results.
%
%
% I N P U T:
% data         - single or multi channel impulse responses of size [N M C]
%                where N is the number of samples, M the number of measure-
%                ments, and C the number of channels.
% toaMethod    - 'onset': uses onset detection from AKonsetDetekt.m for
%                         estimating the TOA
%              - 'cc'   : uses cross correlation beteen data and the
%                         minimum-phase version of data to estimate the TOA
%                         by finding the point of maximum cross-correlation
%                         with AKcrossCorr.m
% toaParameter - a cell that specifies the parameters for TOA estimation in
%                dependency of the toaMethod
%                'onset': {onset_threshold threshMode} see AKonsetDetect
%                         for documentation.
%                'cc'   : {fs corrRange absolute} fs gives the sampling
%                         rate in Hz. See AKcrossCorr for documentation of
%                         remaining parameters.
%                false sets the defualt parameters which are {-6 'rel'}
%                for TOA estimation by onset detection, and 
%                {44100 [0 round(size(data,1)/4)] false} if using cross-
%                correlation. Note that negative values for the corrRange
%                don't make sense if using a minimum-phase response for
%                cross correlation. They are thus clipped to zero.
% upSample     - positve integer that specifies the up-sampling factor.
%                The default false will not apply up-sampling.
% smoothing    - specify a value between 0 and 1 to apply smoothing splines
%                to the TOA estimates (1 no smoothing). The default false
%                will not smooth the estimated TOAs.
% doPlot       - true or false to show or omit plotting
%
% O U T P U T:
% toa          - estimated TOAs in samples with a quantization of
%                1/upSample. toa has size [M C] (see INPUT 'data')
%
%
% 12/2016  -  fabian.brinkmann@tu-berlin.de

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
function [toa, toaSmoothed] = AKtoa(data, toaMethod, toaParameter, upSample, smoothing, doPlot)

% ------------------------------------------------ 1. set default parameter
if ~exist('toaMethod', 'var')
    toaMethod = 'onset';
end
if ~exist('toaParameter', 'var')
    if strcmpi(toaMethod, 'onset')
        toaParameter = {6 'rel'};
    elseif strcmpi(toaMethod, 'cc')
        toaParameter = {44100 round(size(data,1)/4) false};
    end
elseif ~iscell(toaParameter)
    if strcmpi(toaMethod, 'onset')
        toaParameter = {6 'rel'};
    elseif strcmpi(toaMethod, 'cc')
        toaParameter = {44100 round(size(data,1)/4) false};
    end
end
if ~exist('upSample', 'var')
    upSample = false;
end
if ~exist('smoothing', 'var')
    smoothing = false;
end
if ~exist('doPlot', 'var')
    doPlot = false;
end

% -------------------------------------------------------- 2. estimate TOAs
% allocate space for output
toa = zeros(size(data,3), size(data,2));

% loop across channels
for cc = 1:size(data,3)
    
    d = data(:,:,cc);
    
    if strcmpi(toaMethod, 'onset')
        % get parameter
        threshold  = toaParameter{1};
        threshMode = toaParameter{2};
        % get TOAs by onset detection
        toa(cc,:) = AKonsetDetect(d, upSample, threshold, threshMode);
        
    elseif strcmpi(toaMethod, 'cc')
        % get, and check parameter
        fs        = toaParameter{1};
        corrRange = toaParameter{2};
        absolute  = toaParameter{3};
        
        if numel(corrRange) == 1
            corrRange = [0 corrRange]; %#ok<AGROW>
        elseif corrRange(1) < 0
            corrRange = [0 corrRange(2)];
        end
        % generate minimum phase IRs
        dMin = AKphaseManipulation(d, fs, 'min', 1, 0);
        % get TOAs by cross correlation between original and min-phase IRs
        toa(cc,:) = AKcrossCorr(d, dMin, corrRange, upSample, absolute);
    else
        error('AKtoa:Input', 'toaMethod must be either ''onset'', or ''correlation''')
    end
    
end

% ------------------------------------------------ 3. smooth estimated TOAs
toaSmoothed = nan(size(toa));

if isnumeric(smoothing)
    for cc = 1:size(data,3)
        
        smSpline          = fit((1:size(toa,2))', toa(cc,:)', 'smoothingspline', 'SmoothingParam', smoothing);
        toaSmoothed(cc,:) = feval(smSpline, 1:size(toa,2));
        
    end
end

% ----------------------------------------------------- 4. plot the reslult
if doPlot
    
    AKf(15*size(data,3),10)
    
    for cc = 1:size(data, 3)
        % get the range for ploting
        xlims = [round(min(toa(cc,:)))-50 round(max(toa(cc,:)))+50];
        xlims(1) = max(xlims(1), 0);
        xlims(2) = min(xlims(2), size(data,1));
        
        subplot(1,size(data,3),cc)
        % plot IRs
        AKp(data(1:xlims(2),:,cc), 'et3d', 'dr', [-40 0], 'cm', 'gray', 'xu', 'n', 'cb', false, 'norm_d', 0)
        xlabel('measurements')
        title(['Channel ' num2str(cc) ': energy time curve (gray) and TOA'])
        % plot detected TOAs
        hold on
        stairs(.5:size(data,2)+.5, [toa(cc,:) toa(cc,end)], 'color', AKcolors('r'), 'linewidth', 1.2)
        if smoothing
            plot3(1:size(data,2), toaSmoothed(cc,:), 1.5*ones(size(data,2), 1), 'color', AKcolors('g'), 'linewidth', 1.2)
        end
        
        ylim(xlims)
        
        if cc == size(data, 3) && smoothing
            legend('TOA', 'TOA smoothed', 'location', 'best')
        end
        
    end
end