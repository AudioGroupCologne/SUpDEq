function data = data_majdak2013(varargin)
%DATA_MAJDAK2013 Listener specific localization in saggital planes
%   Usage: data = data_majdak2013(condition)
% 
%   Output parameters:
%     data.subject_str : listener ID (string)
%     data.pos         : experimental data matrix conaining 9 colums
%
%                        col 1: target azimuth
%
%                        col 2: target elevation
%
%                        col 3: response azimuth
%
%                        col 4: response elevation
%
%                        col 5: target lateral angle
%
%                        col 6: target polar angle
%
%                        col 7: response lateral angle
%
%                        col 8: response polar angle
%
%                        col 9: response polar angle adjusted (mirrored if confusion)
%     data.session     : session ID
%     data.session_str : session ID (string)
%     data.trial       : trial ID (within each session)
%     data.is_control  : listener group (control vs tested condition)
% 
%
%   DATA_MAJDAK2013(condition) returns listener-specific experimental data
%   from Majdak et al.  (2013) testing localization performance in sagittal
%   planes for low-pass filtered and spectrally warped DTFs.
%
%   The condition flag may be one of:
%
%     'all'     data for all conditions. This is the default.
%     'control' only data of the control group (listeners trained with
%               low-pass filtered (at 8.5kHz) DTFs).
%     'target'  only data of the target group (listeners trained with
%               spectrally warped (2.8-16kHz warped to 2.8-8.5kHz) DTFs).
%     'fig6'    the data necessary to reproduce Figure 6 from the
%               associated publication
%
%   References:
%     P. Majdak, T. Walder, and B. Laback. Effect of long-term training on
%     sound localization performance with spectrally warped and band-limited
%     head-related transfer functions. J. Acoust. Soc. Am., 134:2148--2159,
%     2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_majdak2013.php


%   #Author: Robert Baumgartner (2015)
%   #Author: David Poirier-Q. (2022)
%   #Author: Clara Hollomey (2023)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Check input options

% Define input flags
definput.flags.testcondition = { 'all', 'control', 'target', 'fig6' };
definput.flags.condition = {'BB','LP','W', 'fig6'};

% Parse input options
[flags, ~] = ltfatarghelper({},definput,varargin);


%% Extract data

% load data
if flags.do_fig6
    data = amt_load('majdak2013', 'data_fig6.mat');
else
    data = amt_load('majdak2013', 'data.mat');
end
%data = data.data;

% apply filter if need be
%if( ~strcmp(flags.condition, 'all')  )
if ~flags.do_all && ~flags.do_fig6
    
    % select data
    selVect = data.is_control == contains(flags.condition, 'control');
    
    % apply filter across struct fields
    fieldNames = fieldnames(data);
    for iField = 1:length(fieldNames)
        data.(fieldNames{iField}) = data.(fieldNames{iField})(selVect, :);
    end

end

if ~flags.do_fig6
    x = data;
    C = find(ismember(x.condition,flags.condition));
    for ll = 1:length(x.subject)

      data(ll).mtx = x.subject(ll).expData{C}(:,1:8);
      data(ll).id = x.subject(ll).id;

    end
else
    data = data.data;
end
% todo:
% - create exp_majdak script
% - integrate sub functions (exponential, tick bar)
% 
% majdak todo:
% - data need to go to auxdata 
% - the access to the data needs to be integrated with the already existing data_majdak2013 - but this will depend on the way you need these data (an those I sent you previously)
% - GetData needs to be integrated (into the data_majdak2013 function, I guess)
% - the format of exp_ function need to be adapted (see other exp_ files as an example). 

