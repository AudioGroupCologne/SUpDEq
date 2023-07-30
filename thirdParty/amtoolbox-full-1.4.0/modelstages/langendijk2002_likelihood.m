function [ la,le,ci ] = langendijk2002_likelihood( p,rang,tang,target,response )
%LANGENDIJK2002_LIKELIHOOD   Likelihood estimation
%
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
%   estimates the likelihood for evaluating model performance
%
%   See also: plot_langendijk2002_likelihood, langendijk2002
%
%   References:
%     E. Langendijk and A. Bronkhorst. Contribution of spectral cues to human
%     sound localization. J. Acoust. Soc. Am., 112:1583--1596, 2002.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/langendijk2002_likelihood.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: none
%   #Author : Robert Baumgartner (2013), OEAW Acoustical Research Institute

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
  
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



