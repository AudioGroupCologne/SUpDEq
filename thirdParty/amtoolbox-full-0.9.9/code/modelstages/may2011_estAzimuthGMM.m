function [azEst,gridAz] = may2011_estAzimuthGMM(AZ_MODULE,intLL,nSources,bPlot)

if nargin < 2 || isempty(intLL);    intLL    = 'HIST'; end
if nargin < 3 || isempty(nSources); nSources = inf;    end
if nargin < 4 || isempty(bPlot);    bPlot    = false;  end

if isfinite(nSources)
    azEst = zeros(1,nSources);
end

gridAz = AZ_MODULE.rangeAZ;     
     
deltaAz = abs(diff(gridAz(1:2)));

% Select method for integrating localization information across time
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/may2011_estAzimuthGMM.php

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
switch(upper(intLL))
    
    case 'AVG'        
        % Integrate log-likelihood across channels
        prob_AF = exp(squeeze(sum(AZ_MODULE.loglik,1)));
        
        % Normalize such that probabilities sum up to one
        % for each frame
        xcorr_Freq = transpose(prob_AF ./ repmat(sum(prob_AF,2),[1 numel(gridAz)]));

        % Integrate pattern all frames 
        xcorr_ALL = sum(xcorr_Freq,2);
       
    case 'HIST'
        % Integrate log-likelihood across channels
        prob_AF = exp(squeeze(sum(AZ_MODULE.loglik,1)));
        
        % Normalize such that probabilities sum up to one
        % for each frame
        xcorr_Freq = transpose(prob_AF ./ repmat(sum(prob_AF,2),[1 numel(gridAz)]));
        
        % Find maximum lag per frame
        maxInt = argmax(xcorr_Freq,1);
        
        % Histogram
        xcorr_ALL = hist(gridAz(maxInt),gridAz);
end

% Find peaks, also consider endpoints as peak candidates
pIdx = may2011_findLocalPeaks([0; xcorr_ALL(:); 0]);
pIdx = pIdx - 1;

% Rank peaks
[temp,idx] = sort(xcorr_ALL(pIdx),'descend'); %#ok

% Number of azimuth estimates
nEst = min(numel(idx),nSources);

% Apply exponential interpolation to refine peak position
delta = may2011_interpolateParabolic(xcorr_ALL,pIdx(idx(1:nEst)));
        
% Azimuth estimate: Take most significant peaks
azEst(1:nEst) = gridAz(pIdx(idx(1:nEst))) + deltaAz * delta;


if bPlot || nargout == 0
    timeSec = 1:size(AZ_MODULE.loglik,2);
    
    mIdxPlot = argmax(xcorr_Freq,1);
    mDelta   = may2011_interpolateParabolic(xcorr_Freq,mIdxPlot); 
    
    figure;
    hSP1 = subplot(7,1,1:5);
    hpos = [0.1300    0.3690    0.7750    0.5560];
    set(hSP1,'position',hpos);
   
    imagesc(gridAz,timeSec,xcorr_Freq.');hold on;
    if isequal(upper(intLL),'HIST')
        h = plot(gridAz(mIdxPlot)+(5*mDelta),timeSec,'k.','MarkerSize',18,'linewidth',2);
        set(h,'color',[0.8 0.8 0.8])
    end
    hold off;
    ylabel('Time (s)')
    set(gca,'XTickLabel',[])
    hcb = colorbar;
    hpos = [0.9060    0.3690    0.0250    0.5560];

    set(hcb,'position',hpos);
    set(gca,'YTick',[0 0.5 1 1.5 2]);    
    set(gca,'YTickLabel',num2str([0 0.5 1 1.5 2].'))

    xlim([-93.5 93.5])
    xlim([-120 120])
    
    ht = title('GMM pattern');
    htpos = get(ht,'position');
    htpos(2) = htpos(2) * 0.25;
    set(ht,'position',htpos);
    
    xtickaz = {'' '' '-90' '' '-60'  '' '-30' '' '0' '' '30' '' '60' '' '90' '' ''};   
    set(gca,'xtick',-120:15:120,'xticklabel',xtickaz)
    
    hSP2 = subplot(7,1,6:7);
    hposS = get(hSP2,'position');
    hposS(2) = hposS(2) * 0.85;
    set(hSP2,'position',hposS);

    if isequal(upper(intLL),'HIST')
    h = bar(gridAz,xcorr_ALL,1);hold on;
    set(h,'FaceColor',[0.45 0.45 0.45])
    else
    h = plot([-95; gridAz; 95],[0; xcorr_ALL; 0],'-');hold on;
    set(h,'color',[0.25 0.25 0.25],'linewidth',1.5)
    end         
    plot(azEst,xcorr_ALL(pIdx(idx(1:nEst))),'kx','MarkerSize',18,'LineWidth',2.5)
    xlim([min(gridAz) max(gridAz)])
    ylim([0 1.35*max(xcorr_ALL)])
    xlabel('Azimuth (deg)')
    ylabel('Activity')
    set(gca,'YTick',[],'YTickLabel',[])
    
    xtickaz = {'' '' '-90' '' '-60'  '' '-30' '' '0' '' '30' '' '60' '' '90' '' ''};   
    set(gca,'xtick',-120:15:120,'xticklabel',xtickaz)
    
    box on;
    xlim([-93.5 93.5])
    xlim([-120 120])
        
    colormap(1-may2011_fireprint);
end


function maxidx = argmax(input, dim)

if nargin < 2 || isempty(dim)
    [temp, maxidx] = max(input);
else
    [temp, maxidx] = max(input,[],dim);
end

