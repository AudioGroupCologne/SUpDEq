% This script installes AKtools by adding the folders to the Matlab search
% path. Change the working directory to the folder that contains this file,
% and execute it to instll AKtools

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

% -------- check if third party tools are already installed outside AKtools
installed.playrec  = false;
installed.pa_wav   = false;
installed.SOFA     = false;
installed.auditory = false;

if exist('playrec', 'file')
    if ~strfind(which('playrec'), 'AKtools') %#ok<*STRIFCND>
        installed.playrec = true;
    end
end
if exist('pa_wavplayrecord', 'file')
    if ~strfind(which('pa_wavplayrecord'), 'AKtools')
        installed.pa_wav = true;
    end
end
if exist('SOFAload', 'file')
    if ~strfind(which('SOFAload'), 'AKtools')
        installed.SOFA = true;
    end
end
if exist('MakeERBFilters', 'file')
    if ~strfind(which('MakeERBFilters'), 'AKtools')
        installed.auditory = true;
    end
end

% ----------------------------------------------------------- add all paths
addpath(pwd, '-end')
addpath(genpath(fullfile(pwd, '0_DemoData')), '-end')
addpath(genpath(fullfile(pwd,'1_DemoScripts')), '-end')
addpath(genpath(fullfile(pwd,'2_Tools')), '-end')

% add SOFA correctly using the start script (try if the silent option works)
rmpath(genpath(fullfile(pwd,'2_Tools', 'ThirdParty', 'SOFA_API_MO_1.0.2')))
if exist('SOFAgetVersion', 'file') == 2
    disp('Using installed SOFA API for Matlab located at')
    disp(which('SOFAgetVersion'))
else
    disp('Installing SOFA API for Matlab located from')
    disp(fullfile(pwd, '2_Tools', 'ThirdParty', 'SOFA_API_MO'))
    try
        run(fullfile(pwd, '2_Tools', 'ThirdParty', 'SOFA_API_MO', 'SOFAstart(''silent'')'))
    catch %#ok<CTCH>
        run(fullfile(pwd, '2_Tools', 'ThirdParty', 'SOFA_API_MO', 'SOFAstart.m'))
    end
end

% - remove third party tools if they were already installed outside AKtools
if installed.playrec
    rmpath(genpath(fullfile(pwd, '2_Tools', 'ThirdParty', 'playrec')))
end
if installed.pa_wav
    rmpath(genpath(fullfile(pwd, '2_Tools', 'ThirdParty', 'portaudio_wavplay')))
end
if installed.SOFA
    rmpath(genpath(fullfile(pwd, '2_Tools', 'ThirdParty', 'SOFA_API_MO')))
end
if installed.auditory
    rmpath(genpath(fullfile(pwd, '2_Tools', 'ThirdParty', 'auditoryToolbox')))
end

% --------------------------------------------------------- show the output
fprintf('\n\n')
disp('*******************************************************************')
disp('*********************   AKtools for Matlab   **********************')
disp('*******************************************************************')
disp('**                                                               **')
disp('**  (C) Fabian Brinkmann, Audio Communication Group, TU Berlin   **')
disp('**     The licence can be found in the AKtools main folder       **')
disp('**                                                               **')
disp('**      If you want to keep AKtools in your path, click <a href="matlab:savepath; disp(''AKtools were saved to Matlab search path'')">here</a>     **')
disp('**         If you want to remove it later run AKtoolsStop        **')
disp('**                                                               **')
disp('**          Some demos require the ''FABIAN head-realted          **')
disp('**                   transfer function data set''                 **')
disp('**         <a href="https://dx.doi.org/10.14279/depositonce-5718.2">https://dx.doi.org/10.14279/depositonce-5718.2</a>        **')
disp('**                                                               **')
disp('**     Audio playback and recording on Windows needs ASIO4all    **')
disp('**                 <a href="http://www.asio4all.com">http://www.asio4all.com</a>                       **')
disp('**                                                               **')
disp('**                                                               **')
disp('*******************************************************************')
disp('*********************   AKtools for Matlab   **********************')
disp('*******************************************************************')
fprintf('\n\n')

clear installed

% ----------------------------------------------- check matlab dependencies
toolboxes = {'Signal Processing Toolbox' 'Statistics and Machine Learning Toolbox';
             'Signal_Toolbox'            'Statistics_Toolbox'};

for tt = 1:size(toolboxes, 2)
    if ~license('test', toolboxes{2,tt})
        warning('AKtools:MatlabToolboxes', ['AKtools needs the ' toolboxes{1,tt} ' which is currently not installed. This might limit the functionality of AKtools.'])
    end
end

clear tt toolboxes

% ----------------------------------------- add context menues to startup.m
startup_file = 'startup.m';
if exist(startup_file, 'file')
    % open startup.m and check if the context menues are already added
    startup_file = which(startup_file);
    fid          = fopen(startup_file, 'r');
    startup_line = fscanf(fid, '%s'); % doesn't contain whitespaces or \n
    if fclose(fid)
        warning('AKtoolsStart:startup','startup_file could not be closed');
    end
    
    if contains(startup_line, 'addedByAKtools')
        add_context = false;
    else
        add_context = true;
    end
    
else
    % find or generate a directory to create a startup file
    if isempty(userpath)
        userpath('reset');
    end
    if ~exist(userpath, 'dir')
        mkdir(userpath);
    end
    
    startup_file = fullfile(userpath, startup_file);
    add_context  = true;
end

if add_context
    try
        % add context menues to the end of the startup files
        fid         = fopen(startup_file, 'a');
        fprintf(fid, '\n');
        fprintf(fid, '%%----------------------------------------------- start of addedByAKtools\n');
        fprintf(fid, 'whatClasses = {''double'', ''java.lang.Object''};\n');
        fprintf(fid, 'menuName    = ''AKtools context-menu'';\n');
        fprintf(fid, 'menuItems   = {''AKp(2D-plot)'', ''AKp(3D-plot)''};\n');
        fprintf(fid, 'menuActions = {''AKp($1, ''''1a'''')'', ''AKp($1, ''''1b'''')''};\n');
        fprintf(fid, 'com.mathworks.mlwidgets.workspace.MatlabCustomClassRegistry.registerClassCallbacks(whatClasses,menuName,menuItems,menuActions);\n');
        fprintf(fid, 'clear whatClasses menuName menuItems menuActions\n');
        fprintf(fid, '%% ------------------------------------------------- end of addedByAKtools\n');
        fclose(fid);
        
        % run everything to add the context menues right away
        whatClasses = {'double', 'java.lang.Object'};
        menuName    = 'AKtools context-menu';
        menuItems   = {'AKp(2D-plot)', 'AKp(3D-plot)'};
        menuActions = {'AKp($1)', 'AKp($1, ''''1b'''')'};
        com.mathworks.mlwidgets.workspace.MatlabCustomClassRegistry.registerClassCallbacks(whatClasses,menuName,menuItems,menuActions);
        clear whatClasses menuName menuItems menuActions
        
        % let the user now what AKtools just did
        fprintf('AKtools added content to %s:\n', startup_file)
        disp('Context menues for plotting variables from the workspace are now available.')
        disp('Plotting can now be done by right clicking a time domain variable in the')
        disp('workspace and select ''AKp(2D-plot)'', or ''AKp(3D-plot)''')
    catch
        % we are a bit sad that this did not work - but we wont tell anyone
    end
end

clear startup_file startup_line fid add_context
