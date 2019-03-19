function s = adaptloop_init(nsigs,fs,varargin);
%ADAPTLOOP   Adaptation loops.
%   Usage: s = adaptloop(nsigs,fs,limit,minlvl,tau);
%          s = adaptloop(nsigs,fs,limit,minlvl);
%          s = adaptloop(nsigs,fs,limit);
%          s = adaptloop(nsigs,fs);
%
%   ADAPTLOOP(nsigs,fs,limit,minlvl,tau) applies non-linear adaptation to an
%   input signal insig sampled at a sampling frequency of fs Hz. limit is used
%   to limit the overshoot of the output, minlvl determines the lowest
%   audible threshhold of the signal and tau are time constants involved
%   in the adaptation loops. The number of adaptation loops is determined
%   by the length of tau.
%
%   ADAPTLOOP(nsigs,fs,limit,minlvl) does as above, but uses the values for
%   tau determined in Dau. 1996.
%
%   ADAPTLOOP(nsigs,fs,limit) does as above with a minimum threshhold minlvl
%   equal to 1e-5.
%
%   ADAPTLOOP(nsigs,fs) does as above with an overshoot limit of limit=10.
%
%   References:
%     T. Dau, D. Pueschel, and A. Kohlrausch. A quantitative model of the
%     effective signal processing in the auditory system. I. Model structure.
%     J. Acoust. Soc. Am., 99(6):3615-3622, 1996a.
%     
%     D. Pueschel. Prinzipien der zeitlichen Analyse beim Hoeren. PhD thesis,
%     Universitaet Goettingen, 1988.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/adaptloop_init.php

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

%   AUTHOR: Peter L. Søndergaard

% ------ Checking of input parameters and default parameters ---------

error(nargchk(2,5,nargin));

definput.keyvals.dim=[];
definput.import = {'adaptloop'};
[flags,keyvals,limit,minlvl_db,tau]  = ltfatarghelper({'limit','minlvl','tau'},definput,varargin);

% Convert minlvl
minlvl=setdbspl(minlvl_db);

% -------- Computation ------------------

s.loops=length(tau);

% Calculate filter coefficients for the loops.

% b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
% a1 coefficient of the upper IIR-filter
a1=exp(-1./(tau*fs));

% To get a range from 0 to 100 model units
s.corr = minlvl^(1/(2^s.loops));
s.mult = 100/(1-s.corr); 

% Determine steady-state levels. The values are repeated to fit the
% number of input signals.
state=repmat(minlvl.^(1./(2.^((1:s.loops).'))),nsigs);
  

% Overshoot Limitation.
  
% Max. possible output value
maxvalue = (1 - state.^2) * limit - 1;

% Factor in formula to speed it up 
factor = maxvalue * 2; 			

% Exponential factor in output limiting function
expfac = -2./maxvalue;
offset = maxvalue - 1;

% Load it all up

s.nsigs  = nsigs;
s.state  = state;
s.a1     = a1;
s.factor = factor;
s.expfac = expfac;
s.offset = offset;
s.minlvl = minlvl;
s.limit  = limit;
%OLDFORMAT

