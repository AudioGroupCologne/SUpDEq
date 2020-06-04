%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% script supdeq_start
%
% Running this script will add all required paths and install SOFA, 
% if not previously installed on your system. 
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

%%
% Start script from AKtools, copied and edited with permission from 
% Fabian Brinkmann, Audio Communication Group, TU Berlin

% -------- check if third party tools are already installed outside supdeq
installed.AK        = false;
installed.auditory  = false;
installed.SOFA      = false;
installed.SOFiA     = false;

if exist('AKsht', 'file')
    if ~strfind(which('AKsht'), 'SUpDEq')
        installed.AK = true;
    end
end
if exist('MakeERBFilters', 'file')
    if ~strfind(which('MakeERBFilters'), 'SUpDEq')
        installed.auditory = true;
    end
end
if exist('SOFAload', 'file')
    if ~strfind(which('SOFAload'), 'SUpDEq')
        installed.SOFA = true;
    end
end
if exist('sofia_itc', 'file')
    if ~strfind(which('sofia_itc'), 'SUpDEq')
        installed.SOFiA = true;
    end
end

% ----------------------------------------------------------- add all paths
addpath(pwd, '-end')
addpath(genpath(fullfile(pwd,'evaluation')), '-end')
addpath(genpath(fullfile(pwd,'materials')), '-end')
addpath(genpath(fullfile(pwd,'thirdParty')), '-end')

% add SOFA correctly using the start script (try if the silent option works)
try
    run(fullfile(pwd, 'thirdParty', 'sofa-api-mo-1.0.2', 'API_MO', 'SOFAstart(''silent'')'))
catch %#ok<CTCH>
    run(fullfile(pwd, 'thirdParty', 'sofa-api-mo-1.0.2', 'API_MO', 'SOFAstart'))
end

% - remove third party tools if they were already installed outside supdeq
if installed.AK
    rmpath(genpath(fullfile(pwd, 'thirdParty', 'AKtools')))
end
if installed.auditory
    rmpath(genpath(fullfile(pwd, 'thirdParty', 'amtoolbox-full-0.9.9')))
end
if installed.SOFA
    rmpath(genpath(fullfile(pwd, 'thirdParty', 'sofa-api-mo-1.0.2')))
end
if installed.SOFiA
    rmpath(genpath(fullfile(pwd, 'thirdParty', 'SOFiA R13_MIT-License')))
end

% --------------------------------------------------------- show the output
fprintf('\n\n')
disp('*******************************************************************')
disp('***** SUpDEq - Spatial Upsampling by Directional Equalization *****')
disp('*******************************************************************')
disp('**                                                               **')
disp('**        (C) Christoph Pörschmann and Johannes M. Arend         **')
disp('**             TH Köln - University of Applied Sciences          **')
disp('**            Institute of Communications Engineering            **')
disp('**      Department of Acoustics and Audio Signal Processing      **')
disp('**          Betzdorfer Str. 2, D-50679 Cologne, Germany          **')
disp('**                                                               **')
disp('**           Run supdeq_demo now to test the toolbox...          **')
disp('**                                                               **')
disp('*******************************************************************')
disp('***** SUpDEq - Spatial Upsampling by Directional Equalization *****')
disp('*******************************************************************')
fprintf('\n\n')

clear installed

% ----------------------------------------------- check matlab dependencies
toolboxes = {'Signal Processing Toolbox';
             'Signal_Toolbox'};

for tt = 1:size(toolboxes, 2)
    if ~license('test', toolboxes{2,tt})
        warning('SUpDEq:MatlabToolboxes', ['SUpDEq needs the ' toolboxes{1,tt} ' which is currently not installed. This might limit the functionality of SUpDEq.'])
    end
end

clear tt toolboxes

