function stage_state = lyon2011_spatialsmooth(coeffs, stage_state)
%LYON2011_SPATIALSMOOTH spatial smoothing using FIR coefficients
%
%   Usage: stage_state = lyon2011_spatialsmooth(coeffs, stage_state)
%
%   Input parameters:
%     coeffs      : struct containing coeffs from AGC stage
%     stage_state : smoothed coefficients
%
%   Output parameters:
%     stage_state : smoothed state array
%
%   ´LYON2011_SPATIALSMOOTH´ performs spatial smoothing. This 
%   file is part of an implementation of Lyon's cochlear model:
%   "Cascade of Asymmetric Resonators with Fast-Acting Compression"
%
%   See also: lyon2011
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lyon2011_spatialsmooth.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #License: Apache2
%   #Author: Richard F. Lyon (2013): original implementation (https://github.com/google/carfac)
%   #Author: Amin Saremi (2016): adaptations for the AMT
%   #Author: Clara Hollomey (2021): integration in the AMT 1.0
%   #Author: Richard Lyon (2022): bug fixes for AMT
%   #Author: Mihajlo Velimirovic (2022): implementation of the option ihc_potential

% This file is licensed unter the Apache License Version 2.0 which details can 
% be found in the AMT directory "licences" and at 
% <http://www.apache.org/licenses/LICENSE-2.0>. 
% You must not use this file except in compliance with the Apache License 
% Version 2.0. Unless required by applicable law or agreed to in writing, this 
% file is distributed on an "as is" basis, without warranties or conditions 
% of any kind, either express or implied.

n_iterations = coeffs.AGC_spatial_iterations;

FIR_coeffs = coeffs.AGC_spatial_FIR;
switch coeffs.AGC_spatial_n_taps
  case 3
    for iter = 1:n_iterations
      stage_state = ...
        FIR_coeffs(1) * stage_state([1, 1:(end-1)], :) + ...
        FIR_coeffs(2) * stage_state + ...
        FIR_coeffs(3) * stage_state([2:end, end], :);
    end
  case 5  % 5-tap smoother duplicates first and last coeffs:
    for iter = 1:n_iterations
      stage_state = ...
        FIR_coeffs(1) * (stage_state([1, 2, 1:(end-2)], :) + ...
        stage_state([1, 1:(end-1)], :)) + ...
        FIR_coeffs(2) *  stage_state + ...
        FIR_coeffs(3) * (stage_state([2:end, end], :) + ...
        stage_state([3:end, end, end-1], :));
    end
  otherwise
    error('Bad AGC_spatial_n_taps in lyon2011_spatialsmooth');
end

