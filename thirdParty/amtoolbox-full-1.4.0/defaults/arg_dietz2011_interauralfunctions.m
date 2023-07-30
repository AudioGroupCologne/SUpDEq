function definput=arg_dietz2011_interauralfunctions(definput)

  definput.keyvals.signal_level_dB_SPL = 70;
  definput.keyvals.tau_cycles  = 5; % see Fig. 3c in Dietz (2011)
  definput.keyvals.compression_power = 0.4;

  % ask for simulating temporal resolution of binaural processor
  % this returns the *_lp values
  definput.flags.lowpass = {'lowpass','nolowpass'};


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_dietz2011_interauralfunctions.php



