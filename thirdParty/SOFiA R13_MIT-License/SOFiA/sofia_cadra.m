% /// ASAR/MARA Research Group
%  
% Technology Arts Sciences TH Köln
% Technical University of Berlin
% Deutsche Telekom Laboratories
% University of Rostock
% WDR Westdeutscher Rundfunk
% IOSONO GmbH Erfurt
% 
% SOFiA sound field analysis
% 
% CAD/R/A CAD Reflection Analysis R13-01-24
% 
% Copyright 2011-2017 Maximilian Rühl and Benjamin Bernschütz, rockzentrale 'AT' me.com  
%               
% This file is part of the SOFiA toolbox under MIT-License
%
%
% [] = sofia_cadra(detPeaks, cadPath, offset)
% ------------------------------------------------------------------------
% detPeaks          From sofia_pd()
%                   Listing of local maxima where each row represents a
%                   local maximum containing the following properties:
%                   1. Respective time slice number
%                   2. Time in ms (Refers to the center of the FFT block)
%                   3. Azimut in DEG    (DEG: Orientation Purpose)
%                   4. Elevation in DEG (DEG: Orientation Purpose)
%                   5. Intensity        (Ref. to the absolute maximum)
%                   6. Intensity in dB  (Ref. to the absolute maximum)
%                   7. Azimut in RAD
%                   8. Elevation in RAD
%
% cadPath           Path of the CAD Venue file used by realCADviewer
%
% offset            Free space around the origin. (Distance Compression)
%                   offset [0...1] referring to the maximum distance. 
%                   [Default 0.1] 
%
% This function uses a peak list generated with sofia_pd() and a simple CAD
% model of the respective room to visualize reflections within the CAD
% model.
%
% Dependencies: realCADviewer()
%

function [] = sofia_cadra(detPeaks, cadFile, offset)

if nargin <3, offset        =  0.1; end

if nargin <2
    [File,Path] = uigetfile('','Choose CAD data');
    cadPath=[Path,File];
end

if nargin <1
    error('detPeaks needed');
end

disp('SOFiA CAD/R/A CAD Reflection Analysis R13-01-24');

sphereSizeFactor   = 0.1;
bubblesize         = 50;
sizeOfFont         = 10;

load(cadFile,'venue')

if ~exist('venue','var')
    error('Invalid CAD data')    
end

maxDist = max(max(venue.obj3d(:,3,:)));
maxTime = max(detPeaks(:,1));

[Xm,Ym,Zm] = sph2cart(detPeaks(:,7), pi/2-detPeaks(:,8), detPeaks(:,1)*maxDist*(1-offset)/maxTime+maxDist*offset);

delete(gcf)

[x,y,z] = sphere;
surf(x*sphereSizeFactor, y*sphereSizeFactor, z*sphereSizeFactor) 


set(gcf,'Color', 'w');
axis equal
lighting phong;
camzoom(1.5);
col = colorbar();
xLabel = get(col,'xlabel');
set(xLabel,'String','dB'); 

[Xmb,Ymb,Zmb]=sph2cart(detPeaks(:,7), pi/2-detPeaks(:,8), maxDist);
colormapcopy = colormap;
colorspace = linspace(min(detPeaks(:,6)), max(detPeaks(:,6)), size(colormapcopy,1));
hold on

for i=1:size(detPeaks,1)    
    findMin = abs(colorspace - detPeaks(i,6));
    matchedColor = find(findMin == min(findMin)) ;    
    colorValue = [colormapcopy(matchedColor,1) colormapcopy(matchedColor,2) colormapcopy(matchedColor,3)];    
    line([0 Xmb(i)], [0 Ymb(i)], [0 Zmb(i)],'lineStyle','--','color',colorValue)
    scatter3(Xm(i),Ym(i),Zm(i), bubblesize, detPeaks(i,6), 'marker', 'o', 'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerFaceColor', colorValue)   
end

hold off

for i=1:size(detPeaks,1)    
    text(Xmb(i),Ymb(i),Zmb(i)+0.2,num2str(i), 'fontsize', sizeOfFont)
end

%3D-Model
realCADviewer(cadFile,'k','k')


