function plot_baumgartner2014_likelistat(la,le,ci,lr)
%plot_baumgartner2014_likelistat plots likelihood statistics according to Langendijk et al. (2002)
%   Usage:           plot_baumgartner2014_likelistat(la,le,ci)
%                    plot_baumgartner2014_likelistat(la,le,ci,lr)
%                    plot_baumgartner2014_likelistat(la,le,ci,lr)
%
%   Input arguments:
%       la:          actual likelihood
%
%       le:          expected likelihood
%
%       ci:          confidence interval for expected likelihood
%
%       lr:          reference likelihoods
%
%   See also: baumgartner2014
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/plot/plot_baumgartner2014_likelistat.php

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

% AUTHOR : Robert Baumgartner, OEAW Acoustical Research Institute


% if size(lr,1)==4
%     lr = mean(lr,2);
% else
%     lr = mean(lr,1);
% end


h=bar(la);
set(gca,'XTick',1:length(la))
set(gca,'YLim',[min(lr(:,1)) max(lr(:,4))-0.01],'Layer','top')
set(gca,'Box','on')
set(h,'FaceColor','white','BarWidth',0.6)

hold on
errorbar(le,(ci(:,1)-ci(:,2))/2,'k.');

% leg = legend('Actual likelihood','Expected likelihood');

% for ii = 1:size(lr,2)
%   plot([0,length(la)+1],[lr(:,ii),lr(:,ii)],'k:')
% end

ylabel({'Normalized';'Log-Likelihood'})
xlabel('Condition')
hold off



% hold on
% for ii = 1:length(lr)
%     plot(0.5:length(la)+0.5,lr(ii)*ones(length(la)+1,1),'k:')
% end
% 
% if length(la) == 1
%   h=bar(la);
%   set(gca,'Box','on')
%   set(h,'FaceColor','white','BarWidth',0.5)
% else % Jackknife
%   mla = mean(la);
%   stdla = std(la);
%   errorbar(mla,stdla,'k.');
% end
%   
% set(gca,'XTick',1:length(la))
% if lr(4)==1    % normalized likelihoods
%     set(gca,'YLim',[0.5 1.1],'Layer','top')
%     ylabel('Normalized Log-Likelihood')
% else
%     set(gca,'YLim',[200 550],'Layer','top')
%     ylabel('Log-Likelihood')
%     xlabel('Condition')
% end
% errorbar(le,(ci(:,1)-ci(:,2))/2,'k.');
% hold off

end
