% AKpMulti(data, plot_type, varargin)
% can be used for showing multiple plots of the data in a single figure.
% AKpMulti parses the plot_type and passes the input to AKp, for example
%
% AKpMulti(AKnoise, {'t2d' 'm2d'})
% plots the time signal and spectrum of a pink noise. AKpMulti can also be
% called using AKp, e.g.
% AKp(AKnoise, {'t2d' 'm2d'})
%
% See AKplotDemo.m for examples
%
%
% I N P U T
% data      - see help AKp
% plot_type - cell array with different plot_type as documented in AKp, or
%             string (see predefined plot types for more information)
%             The string in each cell determine the current plot type, and
%             the size of the cell array determines the subplot layout,
%             e.g.
%             {'t2d'; 'm2d'}
%             plot the impulse response and spectrum in a figure with 2
%             vertically arranges subplots
%             {'t2d', 'm2d'}
%             plot the impulse response and spectrum in a figure with 2
%             horizontally arranges subplots
%             {'t2d', 'm2d'; 'p2d' []}
%             plots a figure with a 2x2 subplot matrix, where the last
%             subplot is empty
%             
% varargin  - is passed to AKp. See AKp for documentation
%
% P R E D E F I N E D  P L O T  T Y P E S
%  '1a' : {'t2d' 'et2d'; 'm2d' 'p2d'}
%  '1b' : {'t3d' 'et3d'; 'm3d' 'p3d'}
%  '2a' : {'t2d'; 'm2d'}
%  '2b' : {'t3d'; 'm3d'}
%
% O U T P U T
% dataOut   - struct holding the output data from AKp. See AKp for
%             documentation
%
% 10/2016  -  fabian.brinkmann@tu-berlin.de

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
% limitations under  the License.
function dataOut = AKpMulti(data, plot_type, varargin)

if ~exist('plot_type', 'var')
    plot_type = '1a';
end

if iscell(plot_type) % <- user defined plots
    
    % make figure if none exists
    if isempty(findall(0,'Type','Figure'))
        AKf
    end
    
    % get dimensions
    [N, M] = size(plot_type);
    
    % run the plots
    for nn = 1:N
        for mm = 1:M
            
            % check if current subplot is empty
            if ~isempty(plot_type{nn,mm})
                subplot(N, M, (nn-1)*M+mm)
                dOut = AKp(data, plot_type{nn,mm}, varargin{:});
                
                % copy output data of single AKp calls to output data of
                % AKpMulti
                f = fieldnames(dOut);
                for ff = 1:numel(f)
                    dataOut((nn-1)*M+mm, 1).(f{ff}) = dOut.(f{ff});
                end
            end
            
        end
    end
    
else  % <- pre defined plots
    
    % re-define plot_type
    switch lower(plot_type)
        case '1a'
            plot_type = {'t2d' 'et2d'; 'm2d' 'p2d'};
        case '1b'
            plot_type = {'t3d' 'et3d'; 'm3d' 'p3d'};
        case '2a'
            plot_type = {'t2d'; 'm2d'};
        case '2b'
            plot_type = {'t3d'; 'm3d'};
    end
    
    % run AKpMulti again
    dataOut = AKpMulti(data, plot_type, varargin{:});
    
end

% clear output data if not wanted
if ~nargout
    clear dataOut
end
