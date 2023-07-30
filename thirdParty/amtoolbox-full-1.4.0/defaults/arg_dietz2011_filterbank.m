function definput=arg_dietz2011_filterbank(definput)
 
  % Parameters for filtering the haircell output
  definput.keyvals.filter_order = 2;            % used for both env and fine
  definput.keyvals.filter_attenuation_db = 10;  % used for both env and fine
  % Finestructure filter
  definput.keyvals.fine_filter_finesse = 3;
  % Envelope/modulation filter
  definput.keyvals.mod_center_frequency_hz = 135;
  definput.keyvals.mod_filter_finesse = 8; % => bandwidth: 16.9 Hz
  % ILD filter
  definput.keyvals.level_filter_order = 2;
  definput.keyvals.level_filter_cutoff_hz = 30;


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_dietz2011_filterbank.php



