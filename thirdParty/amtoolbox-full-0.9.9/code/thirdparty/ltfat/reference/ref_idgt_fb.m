function [f]=ref_idgt_fb(coef,g,L,a,M)
%COMP_IDGT_FB  Filter bank IDGT.
%   Usage:  f=comp_idgt_fb(c,g,L,a,M);
%       
%   This is a computational routine. Do not call it directly.
%
%   See also: idgt
%
%   Url: http://ltfat.github.io/doc/reference/ref_idgt_fb.html

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

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

% Calculate the parameters that was not specified.
N=L/a;
b=L/M;

R=size(g,2);

W=prod(size(coef))/(M*N*R);

N=L/a;
b=L/M;

gl=length(g);
glh=floor(gl/2);  % gl-half

% Apply ifft to the coefficients.
coef=ifft(reshape(coef,M,N*W))*sqrt(M);

coef=reshape(coef,M,N,W);

% The fftshift actually makes some things easier.
g=fftshift(g);

g
f=zeros(L,W);

% Make multicolumn g by replication.
gw=repmat(g,1,W);

ff=zeros(gl,1);

% Rotate the coefficients, duplicate them until they have same
% length as g, and multiply by g.

for w=1:W
  
    
  % ----- Handle the first boundary using periodic boundary conditions. ---
  for n=0:ceil(glh/a)-1
    delay=mod(-n*a+glh,M);
    for ii=0:gl/M-1
      for m=0:delay-1
        ff(m+ii*M+1)=coef(M-delay+m+1,n+1,w)*g(m+ii*M+1);
      end;
      for m=0:M-delay-1
        ff(m+ii*M+delay+1)=coef(m+1,n+1,w)*g(m+delay+ii*M+1);
      end;
    end;
    
    sp=mod(n*a-glh,L);
    ep=mod(n*a-glh+gl-1,L);
    
    % Add the ff vector to f at position sp.
    for ii=0:L-sp-1
      f(sp+ii+1,w)+=ff(1+ii);
    end;
    for ii=0:ep
      f(1+ii,w)+=ff(L-sp+1+ii); 
    end;
  end;

  % ----- Handle the middle case. ---------------------
  for n=ceil(glh/a):floor((L-ceil(gl/2))/a)
    delay=mod(-n*a+glh,M);    
    for ii=0:gl/M-1
      for m=0:delay-1
        ff(m+ii*M+1)=coef(M-delay+m+1,n+1,w)*g(m+ii*M+1);
      end;
      for m=0:M-delay-1
        ff(m+ii*M+delay+1)=coef(m+1,n+1,w)*g(m+delay+ii*M+1);
      end;
    end;

    sp=mod(n*a-glh,L);
    ep=mod(n*a-glh+gl-1,L);
  
    % Add the ff vector to f at position sp.
    for ii=0:ep-sp
      f(ii+sp+1,w)+=ff(ii+1);
    end;
  end;

  % ----- Handle the last boundary using periodic boundary conditions. ---
  % This n is one-indexed, to avoid to many +1
  for n=floor((L-ceil(gl/2))/a)+1:N-1
    delay=mod(-n*a+glh,M);
    for ii=0:gl/M-1
      for m=0:delay-1
        ff(m+ii*M+1)=coef(M-delay+m+1,n+1,w)*g(m+ii*M+1);
      end;
      for m=0:M-delay-1
        ff(m+ii*M+delay+1)=coef(m+1,n+1,w)*g(m+delay+ii*M+1);
      end;
    end;
    
    sp=mod(n*a-glh,L);
    ep=mod(n*a-glh+gl-1,L);
    
    % Add the ff vector to f at position sp.
    for ii=0:L-sp-1
      f(sp+ii+1,w)+=ff(1+ii);
    end;
    for ii=0:ep
      f(1+ii,w)+=ff(L-sp+1+ii); 
    end;
  end;
end;

% Scale correctly.
f=sqrt(M)*f;







