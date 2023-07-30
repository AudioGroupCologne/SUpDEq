function conductance = lyon2011_detect(x_in)
%LYON2011_DETECT calculates conductance using a sigmoidal detection nonlinearity
%
%   Usage: conductance = lyon2011_detect(x_in)
%
%   Input parameters:
%     x_in : input signal
%
%   Output parameters:
%     conductance : conductance
%
%   An IHC-like sigmoidal detection nonlinearity for the CARFAC.
%   Resulting conductance is in about [0...1.3405]
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lyon2011_detect.php


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




a = 0.175;   % offset of low-end tail into neg x territory
% this parameter is adjusted for the book, to make the 20% DC
% response threshold at 0.1

set = x_in > -a;
z = x_in(set) + a;

% zero is the final answer for many points:
conductance = zeros(size(x_in));
conductance(set) = z.^3 ./ (z.^3 + z.^2 + 0.1);


%% other things I tried:
%
% % zero is the final answer for many points:
% conductance = zeros(size(x_in));
%
% order = 4;  % 3 is a little cheaper; 4 has continuous second deriv.
%
% % thresholds and terms involving just a, b, s are scalar ops; x are vectors
%
% switch order
%   case 3
%     a = 0.15;  % offset of low-end tail into neg x territory
%     b = 1; % 0.44;   % width of poly segment
%     slope = 0.7;
%
%     threshold1 = -a;
%     threshold2 = b - a;
%
%     set2 = x_in > threshold2;
%     set1 = x_in > threshold1 & ~set2;
%
%     s = slope/(2*b - 3/2*b^2);  % factor to make slope at breakpoint
%     t = s * (b^2 - (b^3) / 2);
%
%     x = x_in(set1) - threshold1;
%     conductance(set1) = s * x .* (x - x .* x / 2);  % x.^2 - 0.5x.^3
%
%     x = x_in(set2) - threshold2;
%     conductance(set2) = t + slope * x ./ (1 + x);
%
%
%   case 4
%     a = 0.24;  % offset of low-end tail into neg x territory
%     b = 0.57;   % width of poly segment; 0.5 to end at zero curvature,
%     a = 0.18;  % offset of low-end tail into neg x territory
%     b = 0.57;   % width of poly segment; 0.5 to end at zero curvature,
%     % 0.57 to approx. match curvature of the upper segment.
%     threshold1 = -a;
%     threshold2 = b - a;
%
%
%     set2 = x_in > threshold2;
%     set1 = x_in > threshold1 & ~set2;
%
%     s = 1/(3*b^2 - 4*b^3);  % factor to make slope 1 at breakpoint
%     t = s * (b^3 - b^4);
%
%     x = x_in(set1) - threshold1;
%     conductance(set1) = s * x .* x .* (x - x .* x);  % x.^3 - x.^4
%
%     x = x_in(set2) - threshold2;
%     conductance(set2) = t + x ./ (1 + x);
%
% end
%


