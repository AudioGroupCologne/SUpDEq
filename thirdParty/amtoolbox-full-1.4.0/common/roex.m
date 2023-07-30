function W = roex(fc, fs)
%ROEX  roex fit as proposed in patterson1982
%   Usage: W = roex(fc,fs,'type');
%
%   Input parameters:
%      fc    : center frequency in Hz.
%      fs    : sampling rate in Hz.
%      type  : 'p', 'pr', or 'pwt'.
%
%   Output parameters:
%      W     :  vector containing the roex-fitted values.
%
%   ROEX(fc, fs, varagin) computes the roex fit to
%   notched-noise masking data as proposed in patterson1982.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/roex.php


%   #Author: Clara Hollomey (2021): adaptations for AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

g = (0.001: 0.001 : 0.8) *fs/fc; 
% corresponds to deltaF/f0, with boundaries as in patterson1982

r = 0.0001; %approximates shallow tail outside of passband
p = 25; %reflects broadening of filter passband
    
%summation of roex approximation and integration tail
W = (1-r)*(1+p*g).*exp(-p*g) + r +...
   -(1-r)*p^(-1)*(2+p*g).*exp(-p*g)+r*g;

%build symmetric filter 
W = [flip(W), W];

