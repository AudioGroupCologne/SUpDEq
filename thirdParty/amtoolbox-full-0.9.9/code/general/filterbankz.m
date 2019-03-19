function outsig=filterbankz(b,a,insig,hopsize)
%FILTERBANKZ  Filter bank with zero boundary condition
%   Usage: outsig=filterbankz(b,a,insig);
%          outsig=filterbankz(b,a,insig,hopsize);
%
%   FILTERBANKZ(b,a,insig,hopsize) filters the input signal with the
%   filters described in a and b. hopsize is a vector with a length
%   equal to the number of filters. Each channel is sub-sampled by the
%   corresponding hopsize.
%
%   If a and b are matrices then each row corresponds to a subband
%   channel.
%
%   If insig is a matrix then filtering is applied along the columns.
%
%   The output coefficients are stored a cell array. More precisely, the
%   n'th cell of c, c{m}, is a 2D matrix of size M(n) xW and
%   containing the output from the m'th channel subsampled at a rate of
%   a(m).  c{m}(n,l) is thus the value of the coefficient for time index
%   n, frequency index m and signal channel l.
%
%   See also: gammatone, ufilterbankz, auditoryfilterbank
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/filterbankz.php

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

%% ------ Checking of input parameters ---------  

error(nargchk(4,4,nargin));


%% ------ Computation --------------------------

[insig,siglen,dummy,nsigs,dim,permutedsize,order]=assert_sigreshape_pre(insig,[],[], ...
                                                  upper(mfilename));
nchannels=size(b,1);


outsig=cell(nchannels,1);

for ii=1:nchannels
  % Calculate the new length in the time domain of this channel
  outlen=ceil(siglen/hopsize);

  % Do the actual filtering.
  res = filter(b(ii,:),a(ii,:),insig);

  % Subsample the output, reshape a multidimensional array to the correct size and store.
  permutedsize(1)=outlen;
  outsig{ii} = assert_sigreshape_post(res(1:hopsize:siglen,:),dim,permutedsize,order);  
end;


