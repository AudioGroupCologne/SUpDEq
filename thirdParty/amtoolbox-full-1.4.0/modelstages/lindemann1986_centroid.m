function d = lindemann1986_centroid(cc)
%LINDEMANN1986_CENTROID Calculates the centroid for a cross-correlation
%   Usage: d = lindemann1986_centroid(cc)
%
%   Input parameters:
%       cc  : Lindemann cross-correlation. Dim: 1 x delay line length
%
%   Output parameters:
%       d   : lindemann1986_centroid in the range -1..1~ms
%
%   LINDEMANN1986_CENTROID(cc) calculates the centroid for a given
%   cross-correlation from the Lindemann model.
%
%   The centroid is computed by (see Lindemann (1986a), page 1613, eq. 22):
%
%              M                  M
%       d = ( sum m*Psi(m) ) / ( sum Psi(m) )
%             m=-M               m=-M
%
%   where M is half the length of the delay line -M,...,M.
%
%   See also: lindemann1986
%
%   References:
%     W. Lindemann. Extension of a binaural cross-correlation model by
%     contralateral inhibition. I. Simulation of lateralization for
%     stationary signals. J. Acoust. Soc. Am., 80:1608--1622, 1986.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/lindemann1986_centroid.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Author: Hagen Wierstorf (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% ------ Checking of input  parameters -----------------------------------

narginchk(1,1);

if ~isnumeric(cc) || ~isvector(cc)
    error('%s: cc has to be a numeric vector signal!',upper(mfilename));
end
% Ensure size(cc) = delay line length x 1
if size(cc,1)==1
    cc = cc';
end


% ------ Computation -----------------------------------------------------
% Calculate the length of the delay line as -M:M
m = linspace(-1,1,length(cc))';
% Calculate the centroid using the -M:M delay line
d = sum(m.*cc) / sum(cc);



