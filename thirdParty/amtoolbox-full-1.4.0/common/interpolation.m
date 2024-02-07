function out = interpolation(xx, yy, x, method)
%INTERPOLATION performs cubic, pchip, or linear interpolation
%
%   Usage: out = interpolation(xx, yy, x, method);
%
%
%   INTERPOLATION finds the values of the
%   underlying function yy=F(xx) at x
%   via pchip, cubic, or linear interpolation
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/interpolation.php


%   #Author: Clara Hollomey (2022): adaptations for AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if nargin < 4
    method = 'pchip';
end

if (length(xx) > 1)
    switch method
        case 'cubic'
          out = interp1(xx, yy, x, method);    
        case 'pchip'
          out = interp1(xx, yy, x, method);
        case 'linear'
          xx = [-10^9 xx 10^9];
          yy = [yy(1) yy yy(end)]; 
          out = interp1(xx, yy, x, method);   
    end 
else
    out = yy * ones(size(x));
end


