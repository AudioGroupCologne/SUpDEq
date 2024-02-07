function out = moore2016_cochlea( s, filenameFilter )
%MOORE2016_COCHLEA performs outer and middle ear filtering
%
%   Input parameters:
%     s : input signal
%     filenameFilter : name of the file containing the filter coefficients
%
%   Output parameters:
%     out : filtered signal
%
%   This code corresponds to the outer and middle ear filtering in the binaural loudness model moore2016
%   in the version for TVL 2016 based on ANSI S3.4-2007 and Moore & Glasberg (2007).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/moore2016_cochlea.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author: Josef Schlittenlacher (2018): original code
%   #Author: Clara Hollomey (2021): integration in the AMT

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


    %load([filenameFilter]);
    data = amt_load('moore2016', 'ff_32000.mat');
    vecCoefficients = data.vecCoefficients;
    out(:,1) = conv( vecCoefficients, s(:,1) );
    if ( min( size(s) ) > 1 )
        out(:,2) = conv( vecCoefficients, s(:,2) );
    else
        out(:,2) = out(:,1);
    end
    out = out((1025):(end-1024),:);
%     out = [zeros(1024,2); s; zeros(1024,2)];  % without filter
end


