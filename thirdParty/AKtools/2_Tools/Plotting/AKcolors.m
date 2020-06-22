% colors = AKcolors(input)
%
% generates the color space used by AKp.m from the Matlab default and some
% additional colors and outputs it.
%
% To see all colors call
% AKcolors(true)
%
% To get all colors call
% colors = AKcolors
%
% To get the RGB values for a single color call
% AKcolors(char)
% where char is one of characters contained in colors.string as obtained
% from colors = AKcolors
%
% 10/2016  - fabian.brinkmann@tu-berlin.de

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
function colors = AKcolors(input)

if ~nargin
    input = false;
end

colors.names = {
                1  'black'      'k'
                2  'white'      'w'
                3  'blue'       'b'
                4  'red'        'r'
                5  'yellow'     'y'
                6  'purple'     'p'
                7  'green'      'g'
                8  'light blue' 'l'
                9  'dark red'   'd'
                10 'magenta'    'm'
                11 'cyan'       'c'
                12 'orange'     'o'
                };

colors.rgb = [
                0         0         0       % 1.  black      'k'
                1         1         1       % 2.  white      'w'
                0         0.4470    0.7410  % 3.  blue       'b'
                0.8500    0.2250    0.0980  % 4.  red        'r'
                0.9290    0.6940    0.1250  % 5.  yellow     'y'
                0.4940    0.1840    0.5560  % 6.  purple     'p'
                0.2660    0.6740    0.1880  % 7.  green      'g'
                0.3010    0.7450    0.9330  % 8.  light blue 'l'
                0.6350    0.0780    0.1840  % 9.  wine red   'w'
                .8        .2        .8      % 10. magenta    'm'
                .1        .8        .8      % 11. cyan       'c'
                1         .54       0       % 12. orange     'o'
              ];
          
colors.string = 'kwbrypgldmco';

if ischar(input)
    colors = colors.rgb(strfind(colors.string, input),:);
end

if input && ~ischar(input)
    
    AKf(20,20)
    set(gcf, 'numberTitle', 'off', 'name', 'AK colors', 'toolBar', 'none', 'menuBar', 'none')
    data = repmat(linspace(1, 0, numel(colors.string)), [2 1]);
    AKp(data, 't2d', 'c', colors.string, 'xu', 'n', 'dr', [0 1.05], 'lw', 5)
    grid off; box off; axis off
    set(gcf, 'color', [.8 .8 .8]);
    title ''
    
    for nn = 1:numel(colors.string)
        text(1.01, data(1,nn)+.03, [num2str(colors.names{nn,1}) ': ' colors.names{nn,2} ' (''' colors.names{nn,3} '''), RGB: ' ...
                              num2str(colors.rgb(nn,1)) '|' num2str(colors.rgb(nn,2)) '|' num2str(colors.rgb(nn,3))])
    end
    
end

if ~nargout
    clear colors
end