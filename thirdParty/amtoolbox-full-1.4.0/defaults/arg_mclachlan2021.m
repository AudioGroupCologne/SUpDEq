function definput=arg_mclachlan2021(definput)
%ARG_MCLACHLAN
%   #License: GPL
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_mclachlan2021.php

definput.keyvals.targ_az = [];          % azimuth target angle in degrees
definput.keyvals.targ_el = [];          % elevation target angle in degrees
definput.keyvals.num_exp = 500; % # of virtual experimental runs
definput.keyvals.SHorder = 5; % spherical harmonics order

% gammatone filterbank parameters
definput.keyvals.fb_ch = 30;
definput.keyvals.fb_low = 300;
definput.keyvals.fb_high = 15e3;
definput.keyvals.source_ir = 0; % broad band sound source at 0dB
definput.keyvals.SNR = 75;      % signal to noise ratio
definput.keyvals.dt = 0.005;    % time step duration in seconds
definput.keyvals.rot_type = 'yaw'; % rotation type ('yaw','pitch','roll')
definput.keyvals.rot_size = 0.1; % rotation amount in degrees
definput.keyvals.stim_dur = 0.1; % stimulus duration in seconds

% parameters of the model computed in the supplementary material
definput.keyvals.sig_itd0 = 0.569;  % standard deviation on initial look
definput.keyvals.sig_itdi = 2*0.569;  % standard deviation on subsequent looks
definput.keyvals.sig_I = 3.5;   % standard deviation on internal noise
definput.keyvals.sig_S = 3.5;   % standard deviation on variation of source spectrum
definput.keyvals.sig = 5; 

definput.groups.group_mclachlan2021={'rot_size',10,'num_exp',10,'sig_itdi',0.8}; % force to use these keyvalues (or flags) even if they are previously defined
definput.groups.fig6={'rot_size',0,'num_exp',10,'sig_itdi',0.8};
definput.groups.fig7a={'rot_size',10,'num_exp',10,'sig_itdi',0.8};
definput.groups.fig7b={'rot_type','pitch','rot_size',10,'num_exp',10,'sig_itdi',0.8};
definput.groups.fig7c={'rot_type','roll','rot_size',10,'num_exp',10,'sig_itdi',0.8};
definput.groups.fig1={'rot_size',0,'num_exp',50,'sig_itd0',1,'sig_itdi',2,'sig_I',10};
definput.groups.group_AMTexample1={'targ_az',30,'targ_el',30,'rot_size',10,'num_exp',1,'sig_itdi',1000,'sig_I',1000};
definput.groups.group_AMTexample2={'targ_az',30,'targ_el',30,'rot_size',10,'num_exp',10,'sig_itdi',1000,'sig_I',6,'sig_S',6};
definput.groups.group_AMTexample3={'targ_az',30,'targ_el',30,'rot_size',10,'num_exp',1,'sig_I',1000};

