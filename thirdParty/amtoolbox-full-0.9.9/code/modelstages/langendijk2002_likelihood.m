function [ la,le,ci ] = langendijk2002_likelihood( p,rang,tang,target,response )
%langendijk2002_likelihood   Likelihood estimation for evaluating model performance
%   Usage: [la,le,ci] = langendijk2002_likelihood(p,rang,tang,target,response)
%
%   Input parameters:
%     p         : pdf matrix
%     rang      : polar angles of possible response angles
%     tang      : polar angles of possible target angles
%     target    : target polar angles of localization test
%     response  : response polar angles of localization test
%
%   Output parameters:
%     la        :  actual likelihood
%     le        :  expected likelihood
%     ci        :  99% confidence interval for expected likelihood
%
%   XXX Describe the function.
%
%   See also: plot_langendijk2002_likelihood, langendijk2002
%
%   References:
%     E. Langendijk and A. Bronkhorst. Contribution of spectral cues to human
%     sound localization. J. Acoust. Soc. Am., 112:1583-1596, 2002.
%     
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/langendijk2002_likelihood.php

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

% AUTHOR : Robert Baumgartner
  
nt=length(target);

% pa represents pdf values of actual responses
pa=interp2([-90;rang(:);270],[-90;tang(:);270], ...
           [zeros(1,size(p,2)+2);[zeros(size(p,1),1),p,zeros(size(p,1),1)];...
            zeros(1,size(p,2)+2)] ,target,response);
la=-2*sum(log(pa))*55/nt; % actual likelihood

% random generator
lex=zeros(100,1);
for ind=1:100
  pe=zeros(size(target));
  for ind1=1:nt 
    post=find(tang>=target(ind1),1); % target position
                                     %         post=randi(size(p,2),1);
    posr = discreteinvrnd(p(:,post),1,1);
    pe(ind1)=p(posr,post);
  end
  lex(ind)=-2*sum(log(pe))*55/nt; 
end

le=mean(lex); % expected likelihood
err=2.58*std(lex);
ci=[le-err le+err]; % confidence interval

function [ X ] = discreteinvrnd(p,m,n)
% DISCRETEINVRND implements an inversion method for a discrete distribution
% with probability mass vector p and dimensions m,n
% Usage:    [ X ] = discreteinvrnd(p,m,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR : Robert Baumgartner, OEAW
% latest update: 2010-07-21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(m,n); 
for i = 1:m*n
    c = cumsum(p);
    u = max(c)*rand;
    X(i) = find(u < c ,1);
end


