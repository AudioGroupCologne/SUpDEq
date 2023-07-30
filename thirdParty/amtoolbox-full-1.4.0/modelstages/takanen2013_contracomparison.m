function [activMap gainMtrx]= takanen2013_contracomparison(thetaL,thetaR,levels,energyL,energyR,fs,printFigs)
%TAKANEN2013_CONTRACOMPARISON Enhanance the contrast between the hemispheres
%   Usage: [activMap gainMtrx]= takanen2013_contracomparison(thetaL,thetaR,levels,energyL,energyR,fs,printFigs);
%
%   Input parameters:
%        theraL    : "where" cue of the left hemisphere
%        thetaR    : "where" cue of the right hemisphere
%        levels    : vector specifying the left/right location
%        energyL   : "what" cue of the left hemisphere
%        energyR   : "what" cue of the right hemisphere
%        fs        : sampling rate
%        printFigs : boolean value describing whether the computations are 
%                    illustrated or not
%
%   Output parameters:
%        activMap  : Matrix that describes in which of the six frequency
%                    ranges there is activation on a given location on
%                    map at a specific time instant
%        gainMtrx  : Matrix that describes the signal level dependent
%                    gains for the different activation values on the
%                    activity map
%        
%   The contralateral comparison is used to control the strengths of the
%   images on the binaural activity map. More specifically, the where cues
%   define positions of the images and an inhibitory coefficient is used to
%   reduce the strength of the image more closer to the center.
%   Additionally, inherent delays are introduced to the inhibition to
%   emulate the binaural adaptation found in the IC (McAlpine et. al.
%   2000).    
%
%   See also: takanen2013, takanen2013_formbinauralactivitymap, 
%             takanen2013_weightedaveragefilter
%
%   References:
%     D. McAlpine, D. Jiang, T. M. Shackleton, and A. R. Palmer. Responses of
%     neurons in the inferior colliculus to dynamic interaural phase cues:
%     Evidence for a mechanism of binaural adaptation. J. Neurophysiol.,
%     83(3):1356 -- 1365, Mar. 2000.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Visualization of functional
%     count-comparison-based binaural auditory model output. Hearing
%     research, 309:147--163, 2014. PMID: 24513586.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/takanen2013_contracomparison.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB M-Signal
%   #Author: Marko Takanen (2013)
%   #Author: Olli Santala 
%   #Author: Ville Pulkki

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% ------ Computations ----------------------------------------------------

thetaL = -1*thetaL;
dims= size(thetaL);
%difference between the activations on left and right side
difference = abs(thetaR-thetaL);

%compute the inhibitory coefficient based on the difference
gains = -1 + 2./(1+exp(-(difference/5)));
energyL(energyL>1) =1;
energyR(energyR>1)=1;

%store the original energies for plotting purposes
eL_orig = energyL;
eR_orig = energyR;

%delay and filter the inhibitory coefficient with a weighted average filter
gL = gains;
gR = gains;
gL((thetaL+thetaR)<0)=0;
gR((thetaL+thetaR)>0)=0;
multipL = takanen2013_weightedaveragefilter(gL,energyR,fs,0.01);
multipR = takanen2013_weightedaveragefilter(gR,energyL,fs,0.01);
multipL = 1-multipL;
multipR=1-multipR;
delay = round(0.0005*fs);
multipR = [ones(delay,dims(2));multipR(1:end-delay,:)];
multipL = [ones(delay,dims(2));multipL(1:end-delay,:)];

%multiply the energies with the inhibitory coefficient
energyL = energyL.*multipL;
energyR = energyR.*multipR;
%half-wave rectification
energyL(energyL<0) =0;
energyR(energyR<0) =0;

%optional plot
if(printFigs)
    band = 2;
    figure(98);
    g(1) = subplot(5,1,1);plot(thetaL(:,band));hold on;plot(thetaR(:,band),'r');hold off;
    legend('Left','Right');ylabel('Spatial cues');ylim([-95 95]);%ylim([-1.1 1.1]);
    g(2) = subplot(5,1,2);plot(gL(:,band),'b');hold on;plot(gR(:,band),'r');hold off;
    legend('Left','Right');ylabel('Inhibitory');ylim([-.3 1.1]);
    g(3) = subplot(5,1,3);plot(eL_orig(:,band),'b');hold on;plot(eR_orig(:,band),'r');hold off;
    legend('Left','Right');ylabel('Energy');ylim([-.3 1.1]);
    g(4) = subplot(5,1,4);plot(multipL(:,band),'b');hold on;plot(multipR(:,band),'r');hold off;
    legend('Left','Right');ylabel('Multiplier');ylim([-.3 1.1]);
    g(5) = subplot(5,1,5);plot(energyL(:,band),'b');hold on;plot(energyR(:,band),'r');hold off;
    legend('Left','Right');ylabel('Mult. energ.');ylim([-.3 1.1]);
    linkaxes(g,'x');
end

%% provide the parameters required for plotting the activity map
levelmtx = ones(dims(1),1)*levels;

%initialize the outputted matrices with zeros
activMap = zeros(dims(1),length(levels));
gainMtrx = activMap;
for i=1:dims(2)
    %find the closest match amongst the level values
    [~,tempL] = min(abs(levelmtx-(thetaL(:,i)*ones(1,length(levels)))),[],2);
    [~,tempR] = min(abs(levelmtx-(thetaR(:,i)*ones(1,length(levels)))),[],2);
    %edit the indexes so that they correspond to the right element in the matrix
    tempL2 = (tempL-1)*dims(1)+(1:dims(1))';
    tempR2 =(tempR-1)*dims(1)+(1:dims(1))';
    activMap(tempL2) = 1;
    activMap(tempR2)=1;
    gainMtrx(tempL2) = energyL(:,i);
    gainMtrx(tempR2) = energyR(:,i);
end


