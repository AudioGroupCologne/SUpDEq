function inoutsig = adaptloop(inoutsig,fs,varargin);
%ADAPTLOOP   Adaptation loops
%   Usage: outsig = adaptloop(insig,fs,limit,minlvl,tau);
%          outsig = adaptloop(insig,fs,limit,minlvl);
%          outsig = adaptloop(insig,fs,limit);
%          outsig = adaptloop(insig,fs);
%
%   ADAPTLOOP(insig,fs,limit,minlvl,tau) applies non-linear adaptation to an
%   input signal insig sampled at a sampling frequency of fs Hz. limit is
%   used to limit the overshoot of the output, minlvl determines the lowest
%   audible threshhold of the signal (in dB SPL) and tau are time constants
%   involved in the adaptation loops. The number of adaptation loops is
%   determined by the length of tau.
%
%   ADAPTLOOP(insig,fs,limit,minlvl) does as above, but uses the values for
%   tau determined in Dau. et al (1996a).
%
%   ADAPTLOOP(insig,fs,limit) does as above with a minimum threshhold minlvl*
%   equal to 0 dB SPL.
%
%   ADAPTLOOP(insig,fs) does as above with an overshoot limit of limit=10.
%
%   ADAPTLOOP takes the following flags at the end of the line of input
%   arguments:
%
%     'adt_dau'        Choose the parameters as in the Dau 1996 and 1997
%                      models. This consists of 5 adaptation loops with
%                      an overshoot limit of 10 and a minimum level of
%                      1e-5. This is a correction in regard to the model
%                      described in Dau et al. (1996a), which did not use 
%                      overshoot limiting. The adaptation loops have an 
%                      exponential spacing. This flag is the default.
%
%     'adt_puschel'    Choose the parameters as in the original Puschel 1988
%                      model. This consists of 5 adaptation loops without
%                      overshoot limiting. The adapation loops have a linear spacing.
%
%     'adt_breebaart'  As 'puschel'
%
%     'dim',d          Do the computation along dimension d of the input. 
%
%   See also: auditoryfilterbank, lopezpoveda2001
%
%   Demos: demo_adaptloop
%
%   References:
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. I. Model structure. J. Acoust.
%     Soc. Am., 110:1074-1088, August 2001.
%     
%     T. Dau, D. Pueschel, and A. Kohlrausch. A quantitative model of the
%     effective signal processing in the auditory system. I. Model structure.
%     J. Acoust. Soc. Am., 99(6):3615-3622, 1996a.
%     
%     D. Pueschel. Prinzipien der zeitlichen Analyse beim Hoeren. PhD thesis,
%     Universitaet Goettingen, 1988.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/adaptloop.php

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

% Copyright (c) 1999 - 2004 Stephan Ewert. All rights reserved.

%   AUTHOR : Stephan Ewert, Morten L. Jepsen, Peter L. Søndergaard

% ------ Checking of input parameters and default parameters ---------

if nargin<2
  error('Too few input parameters.');
end;

definput.import = {'adaptloop'};
definput.keyvals.dim=[];
[flags,keyvals,limit,minlvl_db,tau]  = ltfatarghelper({'limit','minlvl','tau'},definput,varargin);

% Convert minlvl from dB SPL to numerical value
minlvl=setdbspl(minlvl_db);

if ~isnumeric(tau) || ~isvector(tau) || any(tau<=0)
  error('%s: tau must be a vector with positive values.',upper(mfilename));
end;

if ~isnumeric(minlvl) || ~isscalar(minlvl) || minlvl<=0
  error('%s: minlvl must be a positive scalar.',upper(mfilename));
end;

if ~isnumeric(limit) || ~isscalar(limit) 
  error('%s: "limit" must be a scalar.',upper(mfilename));
end;  

% -------- Computation ------------------

[inoutsig,siglen,dummy,nsigs,dim,permutedsize,order]=assert_sigreshape_pre(inoutsig,[],keyvals.dim, ...
                                                  upper(mfilename));

% The current implementation of the adaptation loops works with
% dboffset=100, so we must change to this setting.
% The output is always the same, so there is no need for changing back.
 
% Obtain the dboffset currently used.
dboffset=dbspl(1);

% Switch the minimum level and signal to the correct scaling.
minlvl=gaindb(minlvl,dboffset-100);
inoutsig=gaindb(inoutsig,dboffset-100);

inoutsig=comp_adaptloop(inoutsig,fs,limit,minlvl,tau);

inoutsig=assert_sigreshape_post(inoutsig,dim,permutedsize,order);


