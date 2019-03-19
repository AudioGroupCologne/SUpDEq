function fb=filterbank_init(b,a,varargin);
%FILTERBANK_INIT  Wrapper around filter to multiple filters
%   Usage: outsig=filterbank_init(b,a);
%          outsig=filterbank_init(b,a,nsigs,hopsize);
%
%   fb=FILTERBANK_INIT(b,a) creates a filterbank structure fb for use with
%   FILTERBANK_BLOCK. The filterbank will filter the input signals with the
%   filters described in b and a.
%
%   fb=FILTERBANK_INIT(b,a,nsigs) does the same assuming that the input
%   to FILTERBANK_BLOCK will consist of nsigs signal at once.
%
%   See also: ufilterbankz, filterbank_block
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/filterbank_init.php

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

%   AUTHOR : Peter L. Søndergaard

% ------ Checking of input parameters ---------  

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

definput.keyvals.nsigs=1;
definput.keyvals.hopsize=1;

[flags,keyvals,fb.nsigs,fb.hopsize]  = ltfatarghelper({'nsigs','hopsize'},definput,varargin);

zilen=max(size(a,2),size(b,2))-1;

fb.b=b;
fb.a=a;
fb.nchannels=size(b,1);

fb.outstart = 0;
fb.outlen   = 0;
fb.outend   = 0;

% Initialize the initial conditions to zero.
fb.zi=zeros(zilen,fb.nsigs,fb.nchannels);

