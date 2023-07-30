function inoutsig = scaletodbspl(inoutsig,spl,varargin)
%SCALETODBSPL  Scale a signal to have a specific SPL (in dB)
%
%   Usage: 
%     outsig = scaletodbspl(insig,spl);
%     outsig = scaletodbspl(insig,spl,dboffset);
%     outsig = scaletodbspl(insig,spl,'ac');
%     outsig = scaletodbspl(spl);
%     outsig = scaletodbspl(spl,[],dboffset);
%
%   SCALETODBSPL(insig,spl) scales the signal insig such that its 
%   sound pressure level (SPL) corresponds to spl (in dB). 
%   For the scaling, the AMT default level convention is used, which,
%   unless modified by the user, is SPL of 93.98 dB for an RMS of 1.
%
%   SCALETODBSPL(insig,spl,dboffset) scales the signal insig such 
%   that its SPL corresponds to spl (in dB). For the scaling, the 
%   user's level convention provided by dboffset defining the SPL 
%   for an RMS of 1.
%
%   SCALETODBSPL(spl) returns the scaling constant required to scale a signal
%   having RMS of 1 to SPL of spl (in dB), corresponding to the AMT default 
%   level convention. 
%
%   SCALETODBSPL(0) returns the reference ref used to calculate the SPL 
%   in the formula SPL = 20*log_10 (signal/ref)
%
%   SCALETODBSPL(spl,[],dboffset) returns the scaling constant required 
%   to scale a signal having RMS of 1 to SPL of spl (in dB), corresponding 
%   to the user's level convention provided by dboffset defining the SPL 
%   for an RMS of 1.
%
%   If the input is a matrix, it is assumed that each column is a signal.
%
%   SCALETODBSPL(insig,spl,'ac') does the same, but considers only the AC
%   component of the signal (i.e. the mean is removed).
%
%   References:
%     B. C. J. Moore. An Introduction to the Psychology of Hearing. Academic
%     Press, 5th edition, 2003.
%     
%
%   See also: dbspl
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/scaletodbspl.php


%   #Author: Peter L. Soendergaard (2009)
%   #Author: Piotr Majdak (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% ------ Checking of input parameters ---------

definput.flags.mean={'noac','ac'};
definput.keyvals.dboffset=-20*log10(20e-6); % default: 20 �Pa = 0 dB
[flags,kv]=ltfatarghelper({'dboffset'},definput,varargin);

narginchk(1,5);

if ~isnumeric(inoutsig)
  error('%s: insig must be numeric.',upper(mfilename));
end;

% In the code below, "scaletodbspl" obtains the reference level from "dbspl"
% by calling "dbspl(1)", which will return only the offset measured in dB.

if nargin==1 || isempty(spl)
  % Special mode, only the level has been given
  spl=inoutsig;
  
  if ~isscalar(spl) 
    error('%s: spl must be a scalar.',upper(mfilename));
  end;

  inoutsig=gaindb(1,spl-dbspl(1,'dboffset',kv.dboffset));
  return;
end;


if ~isnumeric(spl) || ~isscalar(spl) 
  error('%s: spl must be a scalar.',upper(mfilename));
end;

% if (nargin<3) || (~ischar(options))
%   options='';
% end;


% ------ Computation --------------------------
if isinf(spl), 
  inoutsig=zeros(size(inoutsig)); 
else
  if isvector(inoutsig)
    inoutsig = gaindb(inoutsig/rms(inoutsig,flags.mean),...
      spl-dbspl(1,'dboffset',kv.dboffset));
  else
    % If we have a matrix, set the level for every column.
    for ii=1:size(inoutsig,2);
      inoutsig(:,ii) = gaindb(inoutsig(:,ii)/rms(inoutsig(:,ii),flags.mean),...
        spl-dbspl(1,'dboffset',kv.dboffset));
    end;
  end;
end


