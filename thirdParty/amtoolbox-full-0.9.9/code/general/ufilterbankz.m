function outsig=ufilterbankz(b,a,insig,hopsize)
%UFILTERBANKZ  Uniform Filter bank with zero boundary condition
%   Usage: outsig=ufilterbankz(b,a,insig);
%          outsig=ufilterbankz(b,a,insig,hopsize);
%
%   UFILTERBANKZ(b,a,insig) filters the input signal with the filters
%   described in a and b.
%
%   UFILTERBANKZ(b,a,insig,hopsize) does the same, but only outputs every
%   hopsize sample in the time domain.
%
%   If a and b are matrices then each row corresponds to a subband
%   channel.
%
%   If insig is a matrix then filtering is applied along the columns.
%
%   If f is a single vector, then the output will be a matrix, where each
%   column in f is filtered by the corresponding filter in g. If f is
%   a matrix, the output will be 3-dimensional, and the third dimension will
%   correspond to the columns of the input signal
%   See also: gammatone, filterbankz, auditoryfilterbank
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/ufilterbankz.php

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

narginchk(3,4);

if nargin==3
  hopsize=1;
end;


%% ------ Computation --------------------------

[insig,siglen,dummy,nsigs,dim,permutedsize,order]=assert_sigreshape_pre(insig,[],[], ...
                                                  upper(mfilename));
nchannels=size(b,1);

outlen=ceil(siglen/hopsize);

outsig=zeros(outlen,nchannels,nsigs);

for ii=1:nchannels
  res = filter(b(ii,:),a(ii,:),insig);
  res = res(1:hopsize:siglen,:);  
  outsig(:,ii,:) = reshape(res,outlen,1,nsigs);
end;

% Modify permutedsize and order to reflect that first dimension has split
% in two
permutedsize=[outlen,nchannels,permutedsize(2:end)];
order=assert_groworder(order);

outsig=assert_sigreshape_post(outsig,dim,permutedsize,order);


