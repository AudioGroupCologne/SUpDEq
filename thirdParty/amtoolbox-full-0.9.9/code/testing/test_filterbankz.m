nsigs=3;

blocksize=4096;

nblocks=5;

fs=16000;

hopsize=2;

siglen=nblocks*blocksize+10*hopsize;



[b,a] = gammatone(erbspacebw(0,fs/2),fs,'complex');

fb=filterbank_init(b,a,nsigs,hopsize);

insig  = randn(siglen,nsigs);
outsig = zeros(siglen/hopsize,fb.nchannels,nsigs);

for ii=0:nblocks
  [outsig_block,fb]=filterbank_block(insig(1+ii*blocksize:min((ii+1)* ...
                                           blocksize,siglen),:),fb);
  
  outsig(fb.outstart:fb.outend,:,:)=outsig_block;
    

end;

outsig_ref=ufilterbankz(b,a,insig,hopsize);

norm(outsig(:)-outsig_ref(:))





%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/test_filterbankz.php

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

