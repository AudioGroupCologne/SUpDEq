function definput = arg_modfilterbank(definput)
% ARG_MODFILTERBANK
%
%   #License: GPL
%   #Author: Alejandro Osses (2020): Initial version
%   #Author: Clara Hollomey (2021): Adapted to AMT
%   #Author: Piotr Majdak (2021): Adapted to AMT 1.0
%
%   References:
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band
%     carriers. J. Acoust. Soc. Am., 102:2892--2905, 1997a.
%     
%     M. Jepsen, S. Ewert, and T. Dau. A computational model of human
%     auditory signal processing and perception. J. Acoust. Soc. Am.,
%     124(1):422--438, 2008.
%     
%     J. Verhey, T. Dau, and B. Kollmeier. Within-channel cues in
%     comodulation masking release (cmr): experiments and model predictions
%     using a modulation-filterbank model. J. Acoust. Soc. Am.,
%     106:2733--2745, 1999.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_modfilterbank.php


% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.flags.mfb = {'mfb','no_mfb'};
definput.flags.modfilterlimit   = {'mfc_upper_limit','no_mfc_upper_limit'};
definput.flags.modfilter_150Hz_LP = {'LP_150_Hz','no_LP_150_Hz','LP_150_Hz_att'};
definput.flags.phase_insens = {'phase_insens_hilbert','no_phase_insens'};
% Attenuation factor applied to mod filters above 10 Hz (only applied if phase_insens_hilbert is on):
definput.flags.att_factor       = {'att_factor','no_att_factor'}; 

definput.keyvals.Q_mfb = 2;
definput.keyvals.mfc_upper_limit_max = 1000; % Hz, maximum upper limit

definput.groups.mfb_dau1997 = { 'no_mfc_upper_limit',...
                                'no_LP_150_Hz', ...
                                'no_att_factor'};

definput.groups.mfb_verhey1999 = {'mfc_upper_limit',...
                            'no_LP_150_Hz', ...
                            'no_att_factor'};

definput.groups.mfb_jepsen2008 = {'mfc_upper_limit',...
                            'LP_150_Hz', ...
                            'att_factor'};

definput.groups.mfb_king2019 = {'mfc_upper_limit',...
                            'no_LP_150_Hz', ...
                            'no_att_factor'}; 
                        
definput.groups.mfb_osses2021_att_gain = {'mfc_upper_limit',...
                            'LP_150_Hz_att', ...
                            'att_factor'}; % corresponds to the "att. gain" condition from Appendix C in Osses and Kohlrausch (2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


