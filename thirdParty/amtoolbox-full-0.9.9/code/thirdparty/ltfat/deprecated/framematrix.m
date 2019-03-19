function G=framematrix(F,L);
%FRAMEMATRIX  Frame synthesis operator matrix
%   Usage: G=framematrix(F,L);
%
%   G=frsynmatrix(F,L) returns the matrix representation G of the frame
%   synthesis operator for a frame F of length L. The frame object F*
%   must have been created using FRAME.
%
%   The frame synthesis operator matrix contains all the frame atoms as
%   column vectors. It has dimensions L xNcoef, where Ncoef is the
%   number of coefficients. The number of coefficients can be found as
%   Ncoef=framered(F)*L. This means that the frame matrix is usually
%   *very* large, and this routine should only be used for small values of
%   L.
%
%   The action of the frame analysis operator FRANA is equal to
%   multiplication with the Hermitean transpose of the frame
%   matrix. Consider the following simple example:
%
%     L=200;
%     F=frame('dgt','gauss',10,20);
%     G=frsynmatrix(F,L);
%     testsig = randn(L,1);
%     res = frana(F,testsig)-G'*testsig;
%     norm(res)
%
%   See also: frame, frana, frsyn
%
%   Url: http://ltfat.github.io/doc/deprecated/framematrix.html

% Copyright (C) 2005-2016 Peter L. Soendergaard <peter@sonderport.dk>.
% This file is part of LTFAT version 2.2.0
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

warning(['LTFAT: FRAMEMATRIX has been deprecated and will be removed',...
         ' in the future releases, please use FRSYNMATRIX instead.']);   

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

Lcheck=framelength(F,L);
if Lcheck~=L
    error('%s: Incompatible frame length.',upper(mfilename));
end;

if F.realinput
    
    %switch(F.type)
    %  case 'dgtreal'
        
    %  This code correctly reproduces the matrix represenation of the
    %  analysis operator, but not of the synthesis.
    %
    %    F2=frame('dgt',F.g,F.a,F.M);
    %    G2=frsynmatrix(F2,L);
    %    M2=floor(F.M/2)+1;
    %    N=L/F.a;
    %    G=zeros(L,M2*N);
    %    for n=0:N-1
    %        G(:,1+n*M2:(n+1)*M2)=G2(:,1+n*F.M:M2+n*F.M);
    %    end;
        
    %  otherwise
        error(['%s: The synthesis operator of real-valued-input frames does is ' ...
               'non-linear and does not have a matrix represenation.']);
        %end;
else
    
  % Generic code handles all frames where there are no extra coefficients
  % in the representation
  Ncoef = framered(F)*L;
  % sprintf for Octave compatibility
  assert(abs(Ncoef-round(Ncoef))<1e-3,sprintf('%s: There is a bug. Ncoef=%d should be an integer.',upper(mfilename),Ncoef));
  Ncoef=round(Ncoef);
  coef=eye(Ncoef);
  G = frsyn(F,coef);  
end;


