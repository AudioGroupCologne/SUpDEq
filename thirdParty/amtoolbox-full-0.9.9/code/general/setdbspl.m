function inoutsig = setdbspl(inoutsig,lvl,varargin);
%SETDBSPL  Set level of signal in dB
%   Usage: outsig = setdbspl(insig,lvl);
%          outsig = setdbspl(insig,lvl,'ac');
%
%   SETDBSPL(insig,lvl) sets the SPL (sound pressure level) of the signal
%   insig to lvl dB, using the convention that a pure tone with an RMS value
%   of 1 corresponds to 100 dB SPL.
%
%   SETDBSPL(lvl) returns a scaling constant that will scale a signal
%   with RMS value of 1 to the correct level. 
%
%   If the input is a matrix, it is assumed that each column is a signal.
%
%   SETDBSPL(insig,lvl,'ac') does the same, but considers only the AC
%   component of the signal (i.e. the mean is removed).
%
%   References:
%     B. C. J. Moore. An Introduction to the Psychology of Hearing. Academic
%     Press, 5th edition, 2003.
%     
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/setdbspl.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
%   See also: dbspl
  
%   Author: Peter L. Søndergaard, 2009

% ------ Checking of input parameters ---------

definput.flags.mean={'noac','ac'};
definput.keyvals.dboffset=100;
[flags,kv]=ltfatarghelper({'dboffset'},definput,varargin);

narginchk(1,5);

if ~isnumeric(inoutsig)
  error('%s: insig must be numeric.',upper(mfilename));
end;

% In the code below, "setdbspl" obtains the reference level from "dbspl"
% by calling "dbspl(1)", which will return only the offset measured in dB.

if nargin==1
  % Special mode, only the level has been given
  lvl=inoutsig;
  
  if ~isscalar(lvl) 
    error('%s: lvl must be a scalar.',upper(mfilename));
  end;

  inoutsig=gaindb(1,lvl-dbspl(1));
  return;
end;


if ~isnumeric(lvl) || ~isscalar(lvl) 
  error('%s: lvl must be a scalar.',upper(mfilename));
end;

% if (nargin<3) || (~ischar(options))
%   options='';
% end;


% ------ Computation --------------------------

if isvector(inoutsig)
  inoutsig = gaindb(inoutsig/rms(inoutsig,flags.mean),...
    lvl-dbspl(1,'dboffset',kv.dboffset));
else
  % If we have a matrix, set the level for every column.
  for ii=1:size(inoutsig,2);
    inoutsig(:,ii) = gaindb(inoutsig(:,ii)/rms(inoutsig(:,ii),flags.mean),...
      lvl-dbspl(1,'dboffset',kv.dboffset));
  end;
end;

