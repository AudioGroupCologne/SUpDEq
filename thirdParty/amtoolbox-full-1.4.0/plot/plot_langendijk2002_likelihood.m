function plot_langendijk2002_likelihood(la,le,ci,tit,horlin)
%PLOT_LANGENDIJK2002_likelihood likelihood statistics according to Langendijk et al. (2002)
%   Usage:  plot_langendijk2002_likelihood(la,le,ci)
%           plot_langendijk2002_likelihood(la,le,ci,tit)
%           plot_langendijk2002_likelihood(la,le,ci,dynlin)
%
%   Input parameters:
%     la     : actual likelihood
%     le     : expected likelihood
%     ci     : confidence interval for expected likelihood
%     tit    : set to 'spatstrat' optionally
%     horlin : vector with dynamic values for reference lines
%
%   PLOT_LANGENDIJK2002_LIKELIHOOD(la,le,ci) plots likelihood statistics according to
%   Langendijk et al. (2002)
%
%   See also: langendijk2002_likelihood, langendijk2002
%
%   References:
%     E. Langendijk and A. Bronkhorst. Contribution of spectral cues to human
%     sound localization. J. Acoust. Soc. Am., 112:1583--1596, 2002.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_langendijk2002_likelihood.php


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


if ~exist('tit','var')
    tit='';
end

% figure('Name','Likelihood statistic','NumberTitle','off');
% clf
hold on
if ~exist('horlin','var')
    plot(0.5:length(la)+0.5,275*ones(length(la)+1,1),'k:') % unimodal
    plot(0.5:length(la)+0.5,350*ones(length(la)+1,1),'k:') % bimodal
    plot(0.5:length(la)+0.5,400*ones(length(la)+1,1),'k:') % trimodal
    plot(0.5:length(la)+0.5,437*ones(length(la)+1,1),'k:') % flat
else
    plot(0.5:length(la)+0.5,horlin(1)*ones(length(la)+1,1),'k:') % unimodal
    plot(0.5:length(la)+0.5,horlin(2)*ones(length(la)+1,1),'k:') % bimodal
    plot(0.5:length(la)+0.5,horlin(3)*ones(length(la)+1,1),'k:') % trimodal
    plot(0.5:length(la)+0.5,horlin(4)*ones(length(la)+1,1),'k:') % flat
end
h=bar(la);
set(gca,'XTick',1:length(la))
if strcmp('spatstrat',tit)==1
    set(gca,'XTickLabel',{'baseline';'dummy';'warped'})
end
set(gca,'YLim',[200 550],'Layer','top')
set(gca,'Box','on')
set(h,'FaceColor','white','BarWidth',0.5)
ylabel('Likelihood')
xlabel('Condition')
errorbar(le,(ci(:,1)-ci(:,2))/2,'k.');
hold off

end


