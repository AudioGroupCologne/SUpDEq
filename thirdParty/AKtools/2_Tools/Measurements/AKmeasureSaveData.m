% This script is called by AKmeasureDemo.m from AKtools
% % See AKmeasureDemo.m for examples

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

if data.raw || data.IR
    fprintf('\n');
    
    
    %% get filename
    g.f = AKf(20,4);
    set(g.f, 'units', 'normalized')
    posTmp = get(g.f, 'Position');
    posTmp(2) = .075;
    set(g.f, 'Position', posTmp)
    clear posTmp
    set(g.f, 'DockControls', 'off', 'MenuBar', 'none', 'name', 'Enter file name (files will be saved to data.dir)', 'NumberTitle','off')
    if isfield(s, 'name')
        g.e = uicontrol('Style','edit', 'units', 'normalized', 'position', [0.05 0.5 0.9 .4], 'fontsize', 16, 'string', s.name);
    else
        g.e = uicontrol('Style','edit', 'units', 'normalized', 'position', [0.05 0.5 0.9 .4], 'fontsize', 16);
    end
    uicontrol('Style','pushbutton', 'string', 'save',    'units', 'normalized', 'position', [0.05 0.1 .45 .4], 'fontsize', 16, 'callback', 'if ~isempty(get(g.e, ''string'')); s.name=get(g.e, ''string''); Continue = true; close gcf; else set(g.e, ''string'', ''enter filename here''); end');
    uicontrol('Style','pushbutton', 'string', 'discard', 'units', 'normalized', 'position', [0.50 0.1 .45 .4], 'fontsize', 16, 'callback', 'Continue=false; close gcf;');
    %%
    uiwait(g.f)
    pause(.1)
    if Continue
        % check if there is data to override
        if exist(fullfile(data.dir, 'Data', [s.name '_setup.mat']), 'file') || ...
                exist(fullfile(data.dir, 'Data', [s.name '_raw.mat']), 'file')   || ...
                exist(fullfile(data.dir, 'Data', [s.name '_ir.mat']), 'file')
            
            g.f = AKf(10,2);
            set(g.f, 'DockControls', 'off', 'MenuBar', 'none', 'name', 'file exists - override?', 'NumberTitle','off')
            uicontrol('Style','pushbutton', 'string', 'cancel', 'units', 'normalized', 'position', [0.05 0.1 .4 .8], 'fontsize', 16, 'callback', 'saveData = false; close gcf');
            uicontrol('Style','pushbutton', 'string', 'save', 'units', 'normalized', 'position', [0.55 0.1 .4 .8], 'fontsize', 16, 'callback', 'saveData = true; close gcf');
            
            uiwait(g.f);
            pause(.1);
        else
            saveData = true;
        end
        
        % save the data
        if saveData
            % remove the sweep from the save data
            sweep = s.sweep;
            s     = rmfield(s, 'sweep');
            
            if data.meta
                file = fullfile(data.dir, 'Data', [s.name '_setup.mat']);
                save(file, 's', 'm', 'x', 'data');
                disp('Setup, meta data, and excitation signal saved to ...');
                fprintf('% 64s\n', ['''',file,'''']);
            end
            
            if data.raw
                file = fullfile(data.dir, 'Data', [s.name '_setup.mat']);
                save(file, 'raw', 'sweep');
                disp('Excitation signal and recorded sweep saved to ...');
                fprintf('% 64s\n', ['''',file,'''']);
            end
            
            if data.IR
                file = fullfile(data.dir, 'Data', [s.name '_ir.mat']);
                save(file, 'ir');
                disp('Impulse response(s) saved to ...');
                fprintf('% 64s\n', ['''',file,'''']);
            end
            
            if data.plot && exist('h', 'var')
                for nn = 1:numel(h)
                    if ishandle(h(nn).f)
                        file = fullfile(data.dir, 'Plots', ...
                            sprintf('%s_ch%d.pdf', s.name, s.chOut(nn)));
                        saveas(h(nn).f, file);
                        disp('Plot saved to ...');
                        fprintf('% 64s\n', ['''',file,'''']);
                    end
                end
            end
            
            s.sweep = sweep;
            clear sweep file;

        end  
    end
    
    % make figures closeable again
    for nn = 1:numel(h)
        if ishandle(h(nn).f)
            set(h(nn).f, 'CloseRequestFcn', 'closereq')
            
			% close figures
            if data.close
                close(h(nn).f);
            end
        end
    end
end

clear h g nn saveData Continue;
