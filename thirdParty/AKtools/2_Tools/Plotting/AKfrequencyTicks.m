% [fTick, fLabel] = AKfrequencyTicks(unit, h, f)
% sets the ticks and tick-labels on an x-axis that shows a frequency plot
%
% I N P U T
% unit - the desired unit 'Hz' (default) or 'kHz'
% h    - if you want to set the frequency axis for a plot, pass the handle
%        to the axis (e.g. gca), otherwise pass false
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
function [fTick, fLabel] = AKfrequencyTicks(unit, h, f)

if ~exist('unit', 'var')
    unit = 'Hz';
end

if ~exist('h', 'var')
    h = false;
end

% set the x ticks
fTick    = [1:10 20:10:100 200:100:1000 2000:1000:10000 20000:10000:100000];

% set the x labels
if strcmpi(unit, 'Hz')
    fLabel = {     1 '' '' '' '' '' '' '' '' ...
                  10 '' '' '' '' '' '' '' '' ...
                 100 '' '' '' '' '' '' '' '' ...
                '1k' '' '' '' '' '' '' '' '' ...
               '10k' '' '' '' '' '' '' '' '' ...
              '100k'};
else
    fLabel = { 0.001 '' '' '' '' '' '' '' '' ...
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
            id = find(fTick == f(nn));
            if id
                if strcmpi(unit, 'Hz')
                    if f(nn) < 1e3
                        fLabel{id} = f(nn);
                    else
                        fLabel{id} = [num2str(f(nn)/1000) 'k'];
                    end
                else
                    fLabel{id} = f(nn)/1000;
                end
            else
                id = find(fTick < f(nn), 1, 'last');
                fTick = [fTick(1:id) f(nn) fTick(id+1:end)];
                if strcmpi(unit, 'Hz')
                    if f(nn) < 1e3
                        fLabel = [fLabel(1:id) f(nn) fLabel(id+1:end)];
                    else
                        fLabel = [fLabel(1:id) [num2str(f(nn)/1000) 'k'] fLabel(id+1:end)];
                    end
                else
                    fLabel = [fLabel(1:id) f(nn)/1000 fLabel(id+1:end)];
                end
            end
        end
    end
end

% set the axis
if strcmpi(class(h), 'matlab.graphics.axis.Axes')
    set(h, 'XTick', fTick, 'XTickLabel', fLabel, 'XScale', 'log')
end

% check if output should be returned
if ~nargout
    clear
end