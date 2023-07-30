% generate mpar
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_eurich2022.php


function mpar = eurich2022_parameters

mpar.fs                             = 48000;
mpar.GT_filters_per_ERBaud          = [1]; % filters per ERB (narrower spacing causes very long computation times)
mpar.GT_bandwidth_factor            = 1; % "standard" filter bandwidth
mpar.GT_lowest_center_frequency     = 67; %Hz
mpar.GT_fix_center_frequency        = 500; % one filter will be centered here --> fc
mpar.GT_highest_center_frequency    = 1000; % just central channel for now
mpar.GT_filterorder                 = 4;
mpar.interference_sigma             = [0.5];
mpar.iKernelThresh                  = 1e-3; % treshold above which a value of the Gaussian filter window is used
mpar.rho_max                        = [0.9]; 
mpar.monaural_internal_noise_sigma  = [0.4];%0.2
mpar.binaural_internal_noise_sigma  = [0.2]; %5


