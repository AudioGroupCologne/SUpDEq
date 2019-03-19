% AKcm = AKcolormaps(AKcm, N)
% creates the N step colormap specified by the string AKcm
% see AKplotDemo.m for examples on how to set colomaps in a plot
%
% I N P U T
% AKcm - string specifiying the desired colormap.
%        this can be one of the custom colormaps listed below, or a map
%        from the color brewer set - call AKcolormaps() for a list of
%        brewer color maps.
%
%        Custom color maps
%        'AKbwr': going from blue to white to red (Hagen Wierstorf style)
%        'AKgyr': going from green to yellow to red (The danger scale)
%        'AKwr' : going from white to red
%        'LGBT' : if you like the look of multicolored maps, but 'yet'
%                 isn't your type you can try 'LGBT', 'LGBTflag', or
%                 'LGBTcat'
%        
% N    - Number of steps / colors in the colormap (default is 128)
%
% O U T P U T
% AKcm - colormap, that is a [N x 3] matrix with rgb values
%
% 06/2016 - fabian.brinkmann@tu-berlin.de, initial dev.
% 10/2016 - helmholz@campus.tu-berlin.de, added brewer colormap support

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
function AKcm = AKcolormaps(AKcm, N)

if ~nargin
    AKcolormapsBrewer
else
    
    if ~exist('N', 'var')
        N = 128;
    end
    
    switch upper(AKcm)
        case 'AKBWR'
            AKcm = makeColorMap([.2 .2 1], [1 1 1], [1 .2 .2], N);
        case 'AKWR'
            AKcm = makeColorMap([1 1 1], [1 .2 .2], N);
        case 'AKGYR'
            AKcm = makeColorMap([18 159 73]/255, [239 238 36]/255, [238 36 36]/255, N);
        case {'LGBT' 'LGBTFLAG' 'LGBTCAT'}
            LGBT = [220  38  36   % red
                    251  91  19   % orange
                    241 221  45   % yellow
                     95 185  95   % green
                     42  57 187   % blue
                    162  42 164]; % ourple
            LGBT = flipud(LGBT) / 255;
            
            x    = linspace(1, N, size(LGBT,1));
            
            if strfind(upper(AKcm), 'FLAG') %#ok<*STRIFCND> contains is recommended and was introduced in 2016b
                AKcm = LGBT;
            elseif strfind(upper(AKcm), 'CAT')
                AKcm = repmat(LGBT, [ceil(N/size(LGBT,1)) 1]);
                AKcm = AKcm(1:N, :);
            else
                AKcm = interp1(x, LGBT, 1:N, 'linear');
            end
            
        otherwise
            try
                brewer_struct = AKcolormapsBrewer(AKcm);
                AKcm = cbrewer(brewer_struct, AKcm, N, 'PCHIP');
            catch
                error(['AKcolormap ' AKcm ' not defined'])
            end
    end
end
