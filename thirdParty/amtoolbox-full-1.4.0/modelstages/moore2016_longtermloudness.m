function Nlt = moore2016_longtermloudness( Nst )
%MOORE2016_LONGTERMLOUDNESS calculates the long term loudness
%
%   Usage: Nlt = moore2016_longtermloudness( Nst )
%
%   Input parameters:
%     Nst : short term loudness
%
%   Output parameters:
%     Nlt : long term loudness
%
%   This code calculates the long term loudness by assemling successive short term
%   loudness frames and applying the corresponding AGC stage. It is part of the binaural 
%   loudness model moore2016 in the version for TVL 2016 based on ANSI S3.4-2007 and Moore & Glasberg (2007).
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/moore2016_longtermloudness.php


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

aA = 0.01;
aR = 0.00133;

Nlt = zeros( size(Nst) );

Nlt(1) = moore2016_agcnextframe( 0, Nst(1), aA, aR );
for i=2:length(Nst)
    Nlt(i) = moore2016_agcnextframe( Nlt(i-1), Nst(i), aA, aR );
end


