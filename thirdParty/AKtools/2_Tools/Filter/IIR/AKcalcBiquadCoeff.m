% For documentation see AKfilter.m, and AKfilterDemo.m
% For python version from Frank Schulz (SONIBLE) see
% http://nbviewer.jupyter.org/github/sonible/nb/blob/master/iir_filter/index.ipynb
%
% Frank Schultz, FG Audiokommunikation, TU Berlin
% frank.schultz@tu-berlin.de, +49 175 15 49 763, Skype: j0shiiv
% 0.00 06.09.2010 taken from calcBiquadCoeff.c

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function [b,a] = AKcalcBiquadCoeff(B,A,fs)

	%2nd order analog filter:
	%B(1)=B0; B(2)=B1; B(3)=B2; A(1)=A0; A(2)=A1; A(3)=A2;
	%		  Y(s)		B0*s^2+B1*s+B2		B(1)*s^2+B(2)*s+B(3)
	%H(s)	= ---- =	--------------	=	--------------------
	%		  X(s)		A0*s^2+A1*s+A2		A(1)*s^2+A(2)*s+A(3)
	
	%bilinear transform H(s)->H(z) with s=2*fs*(z-1)/(z+1)
	
	%2nd order digital filter:
	%b(1)=b0; b(2)=b1; b(3)=b2; a(1)=1.; a(2)=a1; a(3)=a2;
	%		  Y(z)		b2*z^-2+b1*z^-1+b0		b(3)*z^-2+b(2)*z^-1+b(1)	
	%H(z)	= ---- =	------------------	=	------------------------
	%		  X(z)		a2*z^-2+a1*z^-1+1		a(3)*z^-2+a(2)*z^-1+a(1)

	fs2 = fs*fs;
	%following code could be more improved in future
	tmp =	 A(3) + 2*A(2)*fs	+ 4*A(1)*fs2;
	a(1) = 1.;
	b(1) =   B(3) + 2*B(2)*fs	+ 4*B(1)*fs2;
	b(2) = 2*B(3) - 8*B(1)*fs2;	
	a(2) = 2*A(3) - 8*A(1)*fs2;	
	b(3) =   B(3) - 2*B(2)*fs	+ 4*B(1)*fs2;	
	a(3) =   A(3) - 2*A(2)*fs	+ 4*A(1)*fs2;			
	b(1) = b(1)/tmp;
	b(2) = b(2)/tmp;
	b(3) = b(3)/tmp;
	a(2) = a(2)/tmp;
	a(3) = a(3)/tmp;
end