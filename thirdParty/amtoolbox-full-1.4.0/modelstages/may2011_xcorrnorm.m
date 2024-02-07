function [xcorr,lags] = may2011_xcorrnorm(inL,inR,maxLag,bDetrend,bNorm)
%MAY2011_XCORRNORM   Normalized time-domain cross-correlation function
%
%   Usage:
%     [XCORR,LAGS] = may2011_xcorrnorm(INL,INR)
%     [XCORR,LAGS] = may2011_xcorrnorm(INL,INR,MAXLAG,bDETREND,bNORM)
%   
%   Input parameters:
%     INL      : left input arranged as  [nSamples x nChannels]
%     INR      : right input arranged as [nSamples x nChannels]
%     MAXLAG   : computation is performned over the lag range -MAXLAG:MAXLAG
%                (default, MAXLAG = nSamples-1)
%     bDETREND : substract mean     (default, bDETREND = true)
%     bNORM    : normalization flag (default, bNORM    = true)
%
%   Output parameters:
%     XCORR    : cross-correlation function [nSamples x nChannels] 
%     LAGS     : time lags of cross-correlation function [2*MAXLAG+1 x 1]
%
%   MAY2011_XCORRNORM calculates the normalized cross-correlation
%   function. It outputs both, the function itself, and the time lags used.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/may2011_xcorrnorm.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Author: Tobias May (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!');
end

% Set default values
if nargin < 3 || isempty(maxLag);   maxLag   = size(inL,1)-1; end
if nargin < 4 || isempty(bDetrend); bDetrend = true;          end
if nargin < 5 || isempty(bNorm);    bNorm    = true;          end
    
% MEX processing
[xcorr,lags] = comp_may2011_xcorrnorm(inL,inR,maxLag,bDetrend,bNorm);

%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************


