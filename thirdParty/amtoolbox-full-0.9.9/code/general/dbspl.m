function y = dbspl(insig,varargin)
%DBSPL RMS value of signal (in dB)
%   Usage: y = dbspl(insig);
%          y = dbspl(insig,'ac');
%
%   DBSPL(insig) computes the average SPL (sound pressure level) of the
%   input signal measured in dB, using the convention that a pure tone at
%   100 dB SPL has an RMS value of 1.
%
%   DBSPL(insig,'dboffset',dboffset) specifies a reference level to
%   convert between RMS and dB SPL. Default value is 100. Some commonly
%   used values are:
%
%    dboffset = 0. This convention was used for the development of the
%     lopezpoveda2001 and Breebaart models.
%
%    dboffset = -20*log10(20e-6) = 93.98. This corresponds to the common
%     convention of the reference being 20 micro Pa. Using this
%     convention, the numerical values of signals is the sound pressure
%     measured in Pascal.
%
%    dboffset = 100. This convention was used for the development of the
%     Dau and Jepsen models.
%
%   Globally changing the reference level for all functions in the
%   toolbox can be done by the following code:
%
%     ltfatsetdefaults('dbspl','dboffset',desired_dboffset);
%
%   and the currently used reference level is obtained by:
%
%     current_dboffset = dbspl(1);
%  
%   DBSPL takes the following flags at the end of the line of input
%   parameters:
%
%     'ac'      Consider only the AC component of the signal (i.e. the mean is
%               removed).
%
%     'dim',d   Work along specified dimension. The default value of []
%               means to work along the first non-singleton one.
%
%     'dboffset',dboffset
%               Specify offset in dB. Default value is 100.
%
%   See also: setdbspl
%
%   References:
%     B. C. J. Moore. An Introduction to the Psychology of Hearing. Academic
%     Press, 5th edition, 2003.
%     
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/dbspl.php

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

%   AUTHOR : Hagen Wierstorf

definput.keyvals.dim=[];
definput.flags.mean={'noac','ac'};
definput.keyvals.dboffset=100;
[flags,kv]=ltfatarghelper({'dim','dboffset'},definput,varargin);

  
% ------ Computation --------------------------

% The level of a signal in dB SPL is given by the following formula:
% level = 20*log10(p/p_0)
% To get to the standard used in the toolbox.
y = 20*log10( rms(insig,flags.mean) )+kv.dboffset;


