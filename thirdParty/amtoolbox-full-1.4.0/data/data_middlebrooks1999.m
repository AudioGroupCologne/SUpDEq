function data = data_middlebrooks1999
%DATA_MIDDLEBROOKS1999 Statistics about non-individualized HRTFs
%   Usage: data = data_middlebrooks1999
%
%   DATA_MIDDLEBROOKS1999 returns statistics summary from Fig. 13
%   (Middlebrooks, 1999b) showing the effect of non-individualized HRTFs.
%
%   Statistics of those parameters are stored as .mean and .quantiles*
%   representing the arithmetic mean and {0,5,25,50,75,95,100} quantiles,
%   respectively.
%
%   The data struct comprises the following fields:
%
%     'le_own'    local lateral RMS error (LE) when localizing with own
%                 HRTFs
%     'le_other'  LE when localizing with others' HRTFs
%     'lb_own'    magnitude of lateral bias (LB) when localizing with own
%                 HRTFs; upper-rear quadrant excluded from analysis
%     'lb_other'  LB when localizing with others' HRTFs
%     'qe_own'    quadrant error rate (QE) when localizing with own
%                 HRTFs
%     'qe_other'  QE when localizing with others' HRTFs
%     'pe_own'    local polar RMS error (PE) when localizing with own
%                 HRTFs
%     'pe_other'  PE when localizing with others' HRTFs
%     'pb_own'    magnitude of polar bias (PB) when localizing with own
%                 HRTFs; upper-rear quadrant excluded from analysis
%     'pb_other'  PB when localizing with others' HRTFs
%
%
%   References:
%     J. C. Middlebrooks. Virtual localization improved by scaling
%     nonindividualized external-ear transfer functions in frequency. J.
%     Acoust. Soc. Am., 106:1493--1510, 1999.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_middlebrooks1999.php


%   #Author: Robert Baumgartner
%   #Author: Roberto Barumerli

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% Quantiles: {0,5,25,50,75,95,100}%

% LE
data.le_own.quantiles = [10.93, 12.12, 13.07, 14.19, 15.61, 17.06, 19.66];
data.le_own.mean = 14.54;
data.le_other.quantiles = [12.25, 14.70, 16.06, 16.88, 19.46, 22.55, 25.45];
data.le_other.mean = 17.04;

% LB
data.lb_own.quantiles = [-0.01, 0.15, 2.43, 3.80, 5.81, 7.11, 8.39];
data.lb_own.mean = 3.85;
data.lb_other.quantiles = [0.12, 0.42, 1.71, 4.23, 11.82, 13.39, 14.17];
data.lb_other.mean = 6.01;

% QE
data.qe_own.quantiles = [0,0,1,3.5,5,13,17];
data.qe_own.mean = 4.5;
data.qe_other.quantiles = [7.5,8,12.5,19,27.5,38,39];
data.qe_other.mean = 21;

% PE
data.pe_own.quantiles = [21,23,25,27,30,34,36];
data.pe_own.mean = 28;
data.pe_other.quantiles = [23,33,38,42,48,54,55];
data.pe_other.mean = 42.5;

% EB
data.pb_own.quantiles = [1,2.5,6,10,13,20,25.5];
data.pb_own.mean = 10;
data.pb_other.quantiles = [0.5,2,7,18,29,42,52.5];
data.pb_other.mean = 19;
end


