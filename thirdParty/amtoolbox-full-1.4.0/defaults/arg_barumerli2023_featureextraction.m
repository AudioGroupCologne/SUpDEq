function definput=arg_barumerli2023_featureextraction(definput)
    definput.flags.compute = {'all', 'template', 'target'};
    definput.flags.feature_monaural = {'pge', 'dtf', 'monaural_none', 'reijniers'};
    definput.flags.source = {'source_broadband', 'source'};
    
    % filterbank options
    definput.keyvals.fs = 48e3;
    definput.keyvals.flow = 700;
    definput.keyvals.fhigh = 18e3;
    
    definput.keyvals.space = 1; % filterbank spacing
    definput.keyvals.monoaural_bw = [0 Inf]*1e3;
    
    % target directions
    definput.keyvals.targ_az = [];
    definput.keyvals.targ_el = [];

    % sound source (time domain)
    % (default) broad band sound source at 0dB
    definput.keyvals.source_ir = [];
    definput.keyvals.source_fs = 0;


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_barumerli2023_featureextraction.php



