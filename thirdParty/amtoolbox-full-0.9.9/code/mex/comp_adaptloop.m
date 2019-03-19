function inoutsig = comp_adaptloop(inoutsig,fs,limit,minlvl,tau);
%COMP_ADAPTLOOP   Computation of adaptation loops.
%   Usage: outsig = comp_adaptloop(insig,fs,limit,minlvl,tau);
%
%   This is a computational routine. DO NOT call it directly. Call
%   ADAPTLOOP instead.
%
%   See also: adaptloop
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/mex/comp_adaptloop.php

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

% -------- Computation ------------------

% Determine the signal length and width
siglen=size(inoutsig,1);
nsigs=size(inoutsig,2);

% Determine the number of adaptation loops
nloops=length(tau);

% Calculate filter coefficients for the loops.

% b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
% a1 coefficient of the upper IIR-filter
b0=1./(tau*fs);
a1=exp(-b0);
b0=1-a1;

% To get a range from 0 to 100 model units
corr = minlvl^(1/(2^nloops));
mult = 100/(1-corr); 

% Apply minimum level to the input
inoutsig = max(inoutsig,minlvl);

% Determine steady-state levels. The values are repeated to fit the
% number of input signals.
state=minlvl.^(1./(2.^((1:nloops))));    

% Back up the value, because state is overwritten
stateinit=state;

if limit <=1 
  % No overshoot limitation

  for w=1:nsigs
    state=stateinit;
    
    for ii=1:siglen
      tmp1=inoutsig(ii,w);
      
      % Compute the adaptation loops.
      for jj=1:nloops
        tmp1=tmp1/state(jj);
        state(jj) = a1(jj)*state(jj) + b0(jj)*tmp1;         
      end;    
      
      % store the result.
      inoutsig(ii,w)=tmp1;
    end;
  
  end;
  
else 
  
  % Overshoot Limitation.
  
  % Max. possible output value
  maxvalue = (1 - state.^2) * limit - 1;
  
  % Factor in formula to speed it up 
  factor = maxvalue * 2; 			
  
  % Exponential factor in output limiting function
  expfac = -2./maxvalue;
  offset = maxvalue - 1;

  for w=1:nsigs
    
    % Determine steady-state levels. The values are repeated to fit the
    % number of input signals.
    state=stateinit;
  
    for ii=1:siglen

      tmp1=inoutsig(ii,w);
      
      for jj=1:nloops
        
        tmp1=tmp1/state(jj);
        
        if ( tmp1 > 1 )
          tmp1 = factor(jj)/(1+exp(expfac(jj)*(tmp1-1)))-offset(jj);
        end
        
        state(jj) = a1(jj)*state(jj) + b0(jj)*tmp1;
      
      end;
      
      % store the result.
      inoutsig(ii,w)=tmp1;    
      
    end;
  end
end;

% Scale to model units
inoutsig = (inoutsig-corr)*mult;


