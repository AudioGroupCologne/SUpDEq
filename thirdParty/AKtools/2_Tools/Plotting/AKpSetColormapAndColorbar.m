% [dr,c_obj] = AKpSetColormapAndColorbar(cm, cb, cr, dr, cl)
%
% sets properties of the colormap and colorbar
%
% AKpSetColormapAndColorbar('RdBu', 'South', 1, [0 5], 'dB')
% will create a colormap from Red to Blue with a colorbar at the bottom of
% the plot. The range will be restricted to 0-5, with a resolution of 1,
% i.e 5 colors will be used. The colorbar is labeld with 'dB'.
%
% See AKplotDemo.m for examples
%
% I N P U T
% cm            : (a) string specifying the colormap (cf. doc colormap)
% (RdBu_flip)         append '_flip' to invert the colormap
%                     (e.g. 'gray_flip' wil have white assigned to the
%                      smallest and black assigned to the largest value)
%                 (b) string specifying AKcolormaps (see AKcolormaps.m)
%                     append '_flip' to invert the colormap
%                 (c) colormap as [N x 3] matrix with values between zero
%                     and one, where each row represents a color and N the
%                     number of steps. cr will be obsolete in this case
%                 (default is RdBu_flip, from the color brewer set)
% cb            : location of colorbar. 0 for not showing the colorbar.
%                 (default is 'EastOutside', see doc colorbar)
% cr            : scalar value specifying the colormap resolution. The 
%                 minimum value defined by 'dr' will be adjusted if cr does
%                 not perfectly fit the range.
%                 (By default, a colormap with 128 steps is created)
% dr            : Dynamic range of the colormap specified by a two element
%                 vector [cLow cHigh].
%                 (By default the current range is used)
% cl            : string giving the label of the colorbar (default = false)
%
% O U T P U T
% dr            : dynamic range of the colormap (might have been changed)
% c_obj         : handle to the colorbar
%
% fabian.brinkmann@tu-berlin.de


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
function [dr,c_obj] = AKpSetColormapAndColorbar(cm, cb, cr, dr, cl)

if ~exist('cm', 'var')
    cm = 'RdBu_flip';
end
if ~exist('cb', 'var')
    cb = 'EastOutside';
end
if ~exist('cr', 'var')
    cr = false;
end
if ~exist('dr', 'var')
    dr = get(gca, 'clim');
end
if ~exist('cl', 'var')
    cl = false;
end

% check if we flip the colormap
if ischar(cm)
    if strfind(cm, '_flip')
        cm = cm(1:end-5);
        cm_inv = true;
    else
        cm_inv = false;
    end
else
    cm_inv = false;
end

% set colormap
if strcmpi(cm(1:2), 'AK') || strcmpi(cm(1:4), 'LGBT') || ~isempty(AKcolormapsBrewer(cm))
    if islogical(cr)
        
        % fixed resolution of 128 steps
        AKcm = AKcolormaps(cm);
        c_obj = colormap(AKcm);
        
    else
        % resolution according to cr
        num_steps = dr(2):-(abs(cr)):dr(1);
        % adjust zmin, to make sure that cr fits into the data range
        if num_steps(end) > dr(1)
            dr(1) = num_steps(end) - abs(cr);
        end
        clear num_steps
        % still rounding is needed to ensure that an integer number is
        % passed to colormap (numerical inaccuracies if zmin and zmax
        % are not integer values)
        AKcm = AKcolormaps(cm, round(abs(dr(2)-dr(1))/cr));
        c_obj = colormap(AKcm);
    end
    
elseif ischar(cm)
    try
        if islogical(cr)
            % fixed resolution of 128 steps
            eval(['c_obj = colormap(' cm '(128));'])
        else
            % resolution according to cr
            num_steps = dr(2):-(abs(cr)):dr(1);
            % adjust zmin, to make sure that cr fits into the data range
            if num_steps(end) > dr(1)
                dr(1) = num_steps(end) - abs(cr);
            end
            clear num_steps
            % still rounding is needed to ensure that an integer number is
            % passed to colormap (numerical inaccuracies if zmin and zmax
            % are not integer values)
            eval(['obj = colormap(' cm '(round(abs(dr(2)-dr(1))/cr)));'])
        end
    catch
        error(['AKcolormap ' cm ' not defined'])
    end
else
    c_obj = colormap(cm);
end

% flip the colormap
if cm_inv
    c_obj = colormap(flipud(colormap));
end

% make a colorbar
if cb ~= 0
    c_obj = colorbar('location', cb);
    
    % label it
    if iscell(cl)
        ylabel(c_obj, cl)
    elseif cl ~= 0
        ylabel(c_obj, cl)
    end
    
end

% limit the colorrange
set(gca, 'CLim', [dr(1) dr(2)]);

% clear output if not required
if ~nargout
    clear dr c_obj
end
