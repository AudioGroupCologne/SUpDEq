function [ la,le,ci,lr,pvalue ] = baumgartner2014_likelistat( p,tang,rang,target,response,varargin )
%baumgartner2014_likelistat - Likelihood statistics for evaluation of model performance
%   Usage: [la,le,ci,lr,pvalue] = baumgartner2014_likelistat(p,tang,rang,target,response,varargin)
%   
%   Input arguments:
%     'p'              pdf matrix
%
%     'tang'           polar angles of possible target angles
%
%     'rang'           polar angles of possible response angles
%
%     'target'         target polar angles of localization test
%
%     'response'       response polar angles of localization test
%
%     'varargin'       Use 'normalize' for normalization of likelihoods.
%                      1 corresponds to unitary likelihood. This is the default.
%                      Use 'original' according to Langendijk et al. (2002).
% 
%   Output arguments:
%     'la'             actual likelihood
%
%     'le'             expected likelihood
%
%     'ci'             99% confidence interval for expected likelihood
%
%     'lr'             reference likelihoods
%                      1st dim: unimodal (1 gaussian dist.: std=17 deg, mu=0)
%                      2nd dim: bimodal  (2 gaussians: mu1=0, mu2=180)
%                      3rd dim: trimodal (mu1=0, mu2=90, mu3=180)
%                      4th dim: unitary
%   
%   See also: baumgartner2014
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/baumgartner2014_likelistat.php

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

definput.flags.normlikeli = {'normalize','original'};
definput.flags.plot = {'noplot','plot'};
definput.flags.jackknife = {'','jackknife'};

definput.keyvals.runs = 100;   

[flags,kv]=ltfatarghelper({'runs'},definput,varargin);

runs = kv.runs;

if iscell(p)
  
  pa = cell(length(p),1); % predicted probabilities at actual positions
  pe = cell(length(p),1); % predicted probabilities at expected positions
  for ll = 1:length(p)
    
    % Exclude target-response pairs not defined in prediction matrix
    tvalid = target{ll} < max(tang{ll}) & target{ll} > min(tang{ll});
    rvalid = response{ll} < max(rang{ll}) & response{ll} > min(rang{ll});
    idvalid = tvalid & rvalid;
    targets = target{ll}(idvalid);
    responses = response{ll}(idvalid);
    nt(ll) = sum(idvalid);
    pa{ll} = interp2(tang{ll},rang{ll},p{ll},targets,responses);
    pa{ll} = pa{ll}';
    
    % Define target angle boundaries necessary for computation of L_expect
    tangbound = tang{ll}(:)+0.5*diff([tang{ll}(1)-diff(tang{ll}(1:2));tang{ll}(:)]);
    
    % Indices of target positions
    post=zeros(size(targets));
    for ii = 1:length(post)
    	post(ii) = find(tangbound>=targets(ii),1);
    end
    
    % Generate response angles 
    pe{ll} = zeros(runs,length(post));
    for ind=1:runs
      posr=zeros(size(post));
      for jj = 1:length(post) 
        posr(jj) = discreteinvrnd(p{ll}(:,post(jj)),1);
        pe{ll}(ind,jj)=p{ll}(posr(jj),post(jj));
      end
    end
    
  end
  
  Nt = sum(nt);
  
  % Actual likelihood
  la=-2*sum(log([pa{:}]+eps));
  
  % Expected likelihood
  lex = -2*sum(log([pe{:}]+eps),2);
  le=mean(lex);
  err=2.58*std(lex); % standard deviation
%   err=2.626*std(lex); % 99%-quantile of Studen-t dist.
  ci=[le-err; le+err]; % 99% confidence interval
%   ci = [min(lex) max(lex)];
  
  % Probability that L_actual occurs as L_expected
  [muhat,sigmahat] = normfit(lex);
  if la > muhat
    pvalue = 2*(1-normcdf(la,muhat,sigmahat));
  else
    pvalue = 2*normcdf(la,muhat,sigmahat);
  end
  
  
% Reference Likelihoods ------------------------------------------------ 
  
  s = 17;   % std in degrees 
  idmed = ceil(length(p)/2);
  % unimodal
  pdfuni = exp(-0.5 * (rang{idmed}./s).^2); % normpdf
  pdfuni = pdfuni /sum(pdfuni);
  iduni = discreteinvrnd(pdfuni,Nt,runs);
  lr(1,1) = mean(-2*sum(log(pdfuni(iduni))));
  
  % unitary
  lr(4,1) = -2*Nt*log(1/length(rang{1}));
  
  
  % normalization of likelihoods, 1 corresponds to unitary
  
  if flags.do_normalize
      la = la/lr(4);
      le = le/lr(4);
      ci = ci/lr(4);
      lr = lr/lr(4);
  end

else

    
%   target = round(target);   % !!!!!!!!!!!!!!!!
  tangbound = tang(:)+0.5*diff([tang(1)-diff(tang(1:2));tang(:)]); % !!!!!!!!!

  
  % Actual Likelihood ----------------------------------------------------
  
  % Exclude undefined target angles
  tvalid = target < max(tang) & target > min(tang);
  rvalid = response < max(rang) & response > min(rang);
  idvalid = tvalid & rvalid;
  target = target(idvalid);
  response = response(idvalid);
  nt = sum(idvalid);
  pa = interp2(tang,rang,p,target,response);
         
  % Evaluate actual likelihood
  la=-2*sum(log(pa+eps)); % actual likelihood
  
  % Jackknife actual likelihood
  if flags.do_jackknife
    jackdiv = 2;
    njack = nt/jackdiv; % resample 
    lax = zeros(runs,1);
    idjack = zeros(njack,runs);
    for ind=1:runs
      idjack(:,ind) = randi(nt,njack,1);
      pax = interp2(tang,rang,p,target(idjack(:,ind)),response(idjack(:,ind)));
      lax(ind)=-2*sum(log(pax+eps));
    end
  end
  
  % Expected Likelihood --------------------------------------------------
  
  post=zeros(size(target)); % indices of target positions
  for ii = 1:nt
      post(ii) = find(tangbound>=target(ii),1);
  end
  
  lex=zeros(runs,1);
  for ind=1:runs
    
    if flags.do_jackknife
      postsel = post(idjack(:,ind));
    else
      postsel = post;
    end
    
    posr=zeros(size(postsel));
    pe=posr;
    for jj = 1:length(postsel) 
      posr(jj) = discreteinvrnd(p(:,postsel(jj)),1);
      pe(jj)=p(posr(jj),postsel(jj));
    end
    lex(ind)=-2*sum(log(pe+eps)); 
  end
  
  le=mean(lex); % expected likelihood

  err=2.58*std(lex); % 99%-quantile of standard normal dist.; std(lex) not standard error

%   err=2.626*std(lex);% 99%-quantile of Studen-t dist.; std(lex) as empirical standard deviation

  
  ci=[le-err le+err]; % 99% confidence interval
  
  % Bhattacharyya distance
  if flags.do_jackknife
    x = min([lax;lex]) : max([lax;lex]);
    plax = hist(lax,x);
    plax = plax/sum(plax);
    plex = hist(lex,x);
    plex = plax/sum(plex);
    Bdist = -log(sum(sqrt(plax .* plex)));
  else
    Bdist = nan;
  end
  
  
  % Reference Likelihoods ------------------------------------------------ 
  
  lr = zeros(4,1);
  s = 17;   % std in degrees 
  % unimodal
  pdfuni = normpdf(rang,0,s);
  pdfuni = pdfuni /sum(pdfuni);
  iduni = discreteinvrnd(pdfuni,nt,runs);
  lr(1) = mean(-2*sum(log(pdfuni(iduni))));
  
  % bimodal
  pdfbi = 0.5*normpdf(rang,0,s) + 0.5*normpdf(rang,180,s);
  pdfbi = pdfbi /sum(pdfbi);
  idbi = discreteinvrnd(pdfbi,nt,runs);
  lr(2) = mean(-2*sum(log(pdfbi(idbi))));
  
  % trimodal
  pdftri = 1/3*normpdf(rang,0,s) + 1/3*normpdf(rang,90,s) + 1/3*normpdf(rang,180,s);
  pdftri = pdftri /sum(pdftri);
  idtri = discreteinvrnd(pdftri,nt,runs);
  lr(3) = mean(-2*sum(log(pdftri(idtri))));
  
  % unitary
  lr(4) = -2*nt*log(1/length(rang));
  
  
  % normalization of likelihoods, 1 corresponds to unitary
  
  if flags.do_normalize
      la = la/lr(4);
      le = le/lr(4);
      ci = ci/lr(4);
      lr = lr/lr(4);
  end
  
end
  
if flags.do_plot
  plot_baumgartner2014_likelistat(la,le,ci,lr)
end
  
end
  
function [ X ] = discreteinvrnd(p,n,m)
% DISCRETEINVRND implements an inversion method for a discrete distribution
% with probability mass vector p for n trials
% Usage:    [ X ] = discreteinvrnd(p,n)
%
% AUTHOR : Robert Baumgartner

if ~exist('m','var')
    m=1;
end

p = p/sum(p);   % ensure probability mass vector
c = cumsum(p);
t = max(c)*rand(n,m); % rand returns ]0,1]
X = zeros(n,m);
for jj = 1:m
    for ii = 1:n
        X(ii,jj) = find(c >= t(ii,jj) ,1);
    end
end

end
