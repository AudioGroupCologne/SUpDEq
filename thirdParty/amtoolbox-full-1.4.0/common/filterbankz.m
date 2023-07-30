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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/filterbankz.php


%   #Author : Peter L. Søndergaard (2012)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% ------ Checking of input parameters ---------  

narginchk(4, 4);

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



