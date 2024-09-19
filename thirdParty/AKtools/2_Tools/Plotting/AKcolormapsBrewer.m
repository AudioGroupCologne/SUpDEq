% Only to be called from AKcolormaps. See AKcolomaps for documentation.
% See AKplotDemo.m for examples
% 
% 10/2016 - helmholz@campus.tu-berlin.de

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
function fieldName = AKcolormapsBrewer(color)

if nargin
    
    % check if the brewer colormap exists
    fieldName = '';
    load('colorbrewer.mat');
    
    fields = fieldnames(colorbrewer);
    for f = 1:length(fields)
        if isfield(colorbrewer.(fields{f}),color)
            fieldName = fields{f};
            return; % found name in colorbrewer.(fields{f})
        end
    end
    
    if ~nargout
        clear fieldName
    end
    
else
    
    % show a plot of all brewer colormaps
    if isempty(findall(0,'Type','Figure'))
        AKf
    end
    set(gcf, 'numberTitle', 'off', 'name', 'Brewer colormaps (http://colorbrewer.org/)', 'toolBar', 'none', 'menuBar', 'none')
    imagesc(imread(which('cbrewer_preview.jpg')))
    axis off
    
    if nargout
        fieldName = [];
    end
    
end
