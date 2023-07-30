function data = data_klingel2022
%DATA_KLINGEL2022 Data from of Klingel & Laback (2022)
%   Usage: data = data_klingel2022
%
%   Output parameters:
%     data   : structure containing the data
%
%   Lateralization data with spatially inconsistent ITD and ILD cues
%
%   The field data_exp1 comprises the raw data of experiment 1:
%
%     'dimension1'   trial
%
%     'dimension2'   frequency band (-1 = multi, 1 = low,
%                          2 = mid-low, 3 = mid-high, 4 = high)
%                          ITD azimuth (deg)
%                          ILD azimuth (deg)
%                          response azimuth (deg)
%                          response time (ms)
%                          item repetition
%
%     'dimension3'    testing time (pretest, posttest)
%
%     'dimension4'   subject (odd numbers belong to ITD group,
%                          even numbers belong to ILD group)
%
%   The field data_exp2 comprises the raw data of experiment 2:
%
%     'dimension1'   trial
%
%     'dimension2'   frequency band (2 = mid-low, 3 = mid-high)
%                    ITD azimuth (deg) 
%                    ILD azimuth (deg)
%                    response azimuth (deg)
%                    response time (ms)
%                    item repitition
%
%     'dimension3'    testing time (pretest, posttest)
%
%     'dimension4'   subject (odd numbers belong to ITD group,
%                          even numbers belong to ILD group)
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_klingel2022.php


%   #Author: Maike Klingel (2022)
%   #Author: Clara Hollomey (2022): adaptations for AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


x = amt_load('klingel2022', 'data_exp1.mat');
y = amt_load('klingel2022', 'data_exp2.mat');
data.exp1 = x.data_exp1;
data.exp2 = y.data_exp2;


