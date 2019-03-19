% AKfrequencyTicks(unit, h, f)
% sets the ticks and tick-labels on an x-axis that shows a frequency plot
%
% I N P U T
% unit - the desired unit 'Hz' (default) or 'kHz'
% h    - handle to the axis (default = gca)
% f    - a vector of frequencies in Hz that should have a text label
%        by default this is [1 100 1e3 1e3 1e5].
%
% 10/2017 - fabian.brinkmann@tu-berlin.de

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
function AKfrequencyTicks(unit, h, f)

if ~exist('unit', 'var')
    unit = 'Hz';
end

if ~exist('h', 'var')
    h = gca;
end

% set the x ticks
xTicks    = [1:10 20:10:100 200:100:1000 2000:1000:10000 20000:10000:100000];

% set the x labels
if strcmpi(unit, 'Hz')
    xLabel = {     1 '' '' '' '' '' '' '' '' ...
                  10 '' '' '' '' '' '' '' '' ...
                 100 '' '' '' '' '' '' '' '' ...
                '1k' '' '' '' '' '' '' '' '' ...
               '10k' '' '' '' '' '' '' '' '' ...
              '100k'};
else
    xLabel = { 0.001 '' '' '' '' '' '' '' '' ...
                0.01 '' '' '' '' '' '' '' '' ...
                 0.1 '' '' '' '' '' '' '' '' ...
                   1 '' '' '' '' '' '' '' '' ...
                  10 '' '' '' '' '' '' '' '' ...
                 100};
end

% write labels for additional frequencies
if exist('f', 'var')
    for nn = 1:numel(f)
        if f(nn) >= 1
            id = find(xTicks == f(nn));
            if id
                if strcmpi(unit, 'Hz')
                    if f(nn) < 1e3
                        xLabel{id} = f(nn);
                    else
                        xLabel{id} = [num2str(f(nn)/1000) 'k'];
                    end
                else
                    xLabel{id} = f(nn)/1000;
                end
            else
                id = find(xTicks < f(nn), 1, 'last');
                xTicks = [xTicks(1:id) f(nn) xTicks(id+1:end)];
                if strcmpi(unit, 'Hz')
                    if f(nn) < 1e3
                        xLabel = [xLabel(1:id) f(nn) xLabel(id+1:end)];
                    else
                        xLabel = [xLabel(1:id) [num2str(f(nn)/1000) 'k'] xLabel(id+1:end)];
                    end
                else
                    xLabel = [xLabel(1:id) f(nn)/1000 xLabel(id+1:end)];
                end
            end
        end
    end
end

% set last value on x-axis
xLims = get(h, 'XLim');
xMax  = xLims(2) - rem(xLims(2), 10e3);

if xMax
    id = xTicks <= xMax;
    xTicks = xTicks(id);
    xLabel = xLabel(id);
end

if xMax > 10e3
    if strcmpi(unit, 'Hz')
        xLabel{end} = [num2str(xMax/1000) 'k'];
    else
        xLabel{end} = num2str(xMax/1000);
    end
end

% set the axis
set(h, 'XTick', xTicks, 'XTickLabel', xLabel, 'XScale', 'log', 'XLim', xLims)