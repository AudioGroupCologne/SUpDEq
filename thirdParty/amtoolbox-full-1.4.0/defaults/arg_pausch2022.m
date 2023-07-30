function definput = arg_pausch2022(definput)

definput.flags.type_stf = {'hartf_front','hartf_rear','hrtf'};
definput.flags.type_mod = {'pausch','kuhn','woodworth','woodworth_ext'};
definput.flags.plot = {'no_plot','plot'};

definput.keyvals.x5 = [];
definput.keyvals.d9 = [];
definput.keyvals.d10 = [];
definput.keyvals.Theta3 = [];
definput.keyvals.azi_min = 0;
definput.keyvals.azi_max = 180;
definput.keyvals.azi_res = 2.5;
definput.keyvals.c = 343;
definput.keyvals.bp_fc_low = 500;
definput.keyvals.bp_fc_high = 1500;

definput.keyvals.col_meas = [0.7, 0.7, 0.7];
definput.keyvals.col_mod1 = [0.4, 0.4, 0.4];
definput.keyvals.col_mod2plus = [0, 0, 0];
definput.keyvals.col_mod3 = [0, 84, 159]/255;

definput.keyvals.col_shad_meas = [0.8, 0.8, 0.8];
definput.keyvals.col_shad_mod1 = [0.7, 0.7, 0.7];
definput.keyvals.col_shad_mod2plus = [0.4, 0.4, 0.4];
definput.keyvals.col_shad_mod3 = [142, 186, 229]/255;

definput.keyvals.alpha = 0.2;




% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_pausch2022.php



