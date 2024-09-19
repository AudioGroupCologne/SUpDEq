% [ccArgMax, ccSign, ccMax] = AKcrossCorr(x, y, corrRange, upSample, absolute)
% finds the maximum cross correlation between x, and y with the possibility
% of up-sampling for sub-sample accuracy.
%
% For example usage see AKtoaDemo.m and AKtoa.m
%
% I N P U T:
% x, y      - single or multi channel time signals of size [N M C], where N
%             is the number of samples, M the number of measurements, and C
%             the number of channels. If only one signal is multi channel,
%             the other is repeated accordingly. If x, and y have a
%             different number of samples, the longer one is truncated
%             accordingly.
% corrRange - one or two element vector that sets the search range in
%             samples to find the maximum cross correlation. E.g. [-10 20]
%             will shift y for -10*upSample samples to the back and
%             20*upSample samples to the front, and [10] will shift y
%             +/-10*upSample samples to the front and back.
%             Note that shifting is not done by a circshift, but by cutting
%             parts of x and y. This means that only a small number of
%             samples is used for cross-correlation if the search range is
%             too large. Because this can cause errornous results, make
%             sure to adjust the search range to your needs.
%             The default is round(size(x,1)/4) samples.
% upSample  - positve integer that specifies the up-sampling factor.
%             Default is 10. Note that the corrRange is extended
%             accordingly, i.e. a corrRange of 36, and an up-sampling of 10
%             will results in 360 samples shift in each direction in the
%             up-sampled domain. Passing 0, 1, or false will omit
%             upsampling.
% absolute  - true to look for the absolute maximum, i.e. ignore the sign
%             of the cross-correlation, or false to look for the largest
%             positive value (default = false)
%
% O U T P U T:
% ccArgMax  - shift of y relative to x that resulted in the largest
%             cross-correlation, i.e. a value of 1.5 means that the largest
%             correlation occured if y was delayed by 1.5 samples.
% ccSign    - sign of the cross-correlation at ccArgMax.
% ccMax     - maximum value of the cross correlation.
%
% 11/2016 - fabian.brinkmann@tu-berlin.de

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
function [ccArgMax, ccSign, ccMax] = AKcrossCorr(x, y, corrRange, upSample, absolute)

%% ---------------------------------------------- 1. check size of x, and y
if size(x,1) < size(y,1)
    y = y(1:size(x,1));
elseif size(x,1) > size(y,1)
    x = x(1:size(y,1));
end

if size(x,2) == 1 && size(y,2) ~= 1
    x = repmat(x, [1 size(y,2) 1]);
elseif size(y,2) == 1 && size(x,2) ~= 1
    y = repmat(y, [1 size(x,2) 1]);
end

if size(x,3) == 1 && size(y,3) ~= 1
    x = repmat(x, [1 1 size(y,3)]);
elseif size(y,3) == 1 && size(x,3) ~= 1
    y = repmat(y, [1 1 size(x,3)]);
end

%% --------------------------------------------- 2. defualt input parameter
if ~exist('absolute', 'var')
    absolute = false;
end
if ~exist('upSample', 'var')
    upSample = 10;
elseif ~upSample
    upSample = 1;
end
if ~exist('corrRange', 'var')
    corrRange = round(size(x,1)/4);
end

if numel(corrRange) == 1
    corrRange = [-corrRange corrRange];
end

corrRange = sort(corrRange);

%% ------------------------------------------- 3. handle multi channel data

if size(x,3) > 1
    
    % allocate space for output data
    ccArgMax = zeros(size(x,3), size(x,2));
    ccSign   = ccArgMax;
    
    % loop over channels
    for cc = 1:size(x,3)
        [ccArgMax(cc,:), ccSign(cc,:)] = AKcrossCorr(x(:,:,cc), y(:,:,cc), corrRange, upSample, absolute);
    end
    
    return
end


%% -------------------------------------------------- 4. up sample x, and y
if upSample > 1
    
    % allocate
    xUp = zeros(size(x,1)*upSample, size(x,2));
    yUp = xUp;
    
    % upsample
    for mm = 1:size(x,2)
        if upSample > 1
            xUp(:,mm) = interp(x(:,mm), upSample);
            yUp(:,mm) = interp(y(:,mm), upSample);
        else
            xUp(:,mm) = x(:,mm);
            yUp(:,mm) = y(:,mm);
        end
    end
    
    % delete emporary variables
    x = xUp;
    y = yUp;
    
    clear xUp yUp cc
    
end

% ------------------------------------------------------ 4. cross correlate

% amount of shift that is applied to y to cross correlate it to x
shift = corrRange(1)*upSample:corrRange(2)*upSample;

% allocate space for raw cross correlation values
crossCorr = zeros(numel(shift), size(x,2));

% calcualte cross correlation
for mm = 1:numel(shift)
    
    % un-normalized cross correlation
    if shift(mm) < 0
        % shift y to the left
        crossCorr(mm,:) = sum( x(1:end-abs(shift(mm)), :) .* y(abs(shift(mm))+1:end, :) );
    else
        % shift x to the left
        crossCorr(mm,:) = sum( x(shift(mm)+1:end, :)      .* y(1:end-shift(mm), :) );
    end
    
    % normalize it to the number of elements
    crossCorr(mm,:) = crossCorr(mm,:) / (size(x,1)-abs(shift(mm)));
    
end

% normalize it to the standard deviations
crossCorr = crossCorr ./ (std(x) .* std(y));

% find maximum
if absolute
    [ccSign, ccID] = max(abs(crossCorr));
else
    [ccSign, ccID] = max(crossCorr);
end

ccMax = crossCorr(ccID);

% restore the sign
for mm = 1:size(x,2)
    ccSign(mm) = sign(crossCorr(ccID(mm), mm));
end

% get the shift between x and y
ccArgMax = shift(ccID) / upSample;