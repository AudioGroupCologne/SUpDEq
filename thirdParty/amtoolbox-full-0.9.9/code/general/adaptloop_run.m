function [inoutsig,s] = adaptloop_run(s,inoutsig);
%COMP_ADAPTLOOP   Computation of adaptation loops.
%   Usage: [outsig,s] = adaptloop_run(s,insig);
%
%   See also: adaptloop
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/general/adaptloop_run.php

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

% -------- Computation ------------------

% Determine the signal length and width
siglen=size(inoutsig,1);

a1=s.a1;
b0=1-a1;

% Apply minimum level to the input
inoutsig = max(inoutsig,s.minlvl);

if s.limit <=1 
  % No overshoot limitation

  for ii=1:siglen
    tmp1=inoutsig(ii,:);
    
    % Compute the adaptation loops.
    for jj=1:s.loops
      tmp1=tmp1./s.state(jj,:);
      s.state(jj,:) = a1(jj)*s.state(jj,:) + b0(jj)*tmp1;         
    end;    
    
    % store the result.
    inoutsig(ii,:)=tmp1;
  end;  
else 
  
  % Overshoot Limitation. 
  
  for ii=1:siglen
      tmp1=inoutsig(ii,:);
      
      for jj=1:s.loops
        
        tmp1=tmp1./s.state(jj,:);

        for w=1:s.nsigs
          if ( tmp1(w) > 1 )
            tmp1(w) = s.factor(jj,w)/(1+exp(s.expfac(jj,w)*(tmp1(w)-1)))-s.offset(jj,w);
          end
          
        end;
        s.state(jj,:) = a1(jj)*s.state(jj,:) + b0(jj)*tmp1;
        
      end;
      
      % store the result.
      inoutsig(ii,:)=tmp1;    
      
  end;
end


% Scale to model units
inoutsig = (inoutsig-s.corr)*s.mult;


%OLDFORMAT

