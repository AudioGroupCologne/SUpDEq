function g = lyon2011_stageg(CAR_coeffs, relative_undamping)
% LYON2011_STAGEG obtain unity gain
%
%   Usage:
%     g = lyon2011_stageg(CAR_coeffs, relative_undamping)
%
%   Input parameters:
%     CAR_coeffs         : struct containing the CARFAC coefficients
%     relative_undamping : relative undamping
%
%   Output parameters:
%     g                  : gain
%
%   Returns the stage gain g needed to get unity gain at DC
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lyon2011_stageg.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #License: Apache2
%   #Author: Richard F. Lyon (2013): original implementation (https://github.com/google/carfac)
%   #Author: Amin Saremi (2016): adaptations for the AMT
%   #Author: Clara Hollomey (2021): integration in the AMT 1.0

%   This file is licensed unter the Apache License Version 2.0 which details can 
%   be found in the AMT directory "licences" and at 
%   <http://www.apache.org/licenses/LICENSE-2.0>. 
%   You must not use this file except in compliance with the Apache License 
%   Version 2.0. Unless required by applicable law or agreed to in writing, this 
%   file is distributed on an "as is" basis, without warranties or conditions 
%   of any kind, either express or implied.




r1 = CAR_coeffs.r1_coeffs;  % at max damping
a0 = CAR_coeffs.a0_coeffs;
c0 = CAR_coeffs.c0_coeffs;
h  = CAR_coeffs.h_coeffs;
zr = CAR_coeffs.zr_coeffs;
r  = r1 + zr .* relative_undamping;
g  = (1 - 2*r.*a0 + r.^2) ./ (1 - 2*r.*a0 + h.*r.*c0 + r.^2);


