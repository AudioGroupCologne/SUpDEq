function [thetaOut energyOut]= takanen2013_onsetenhancement(thetaIn,energyIn,fs,cfs)
%TAKANEN2013_ONSETENHANCEMENT Emphasize the onsets on direction analysis
%   Usage: [thetaOut energyOut]=takanen2013_onsetenhancement(thetaIn,energyIn,fs);
%
%   Input parameters:
%        thetaIn  : the "where" cue of the given hemisphere
%        energyIn : the ventral "what" stream output of the periphery model
%        fs       : sampling rate
%        cfs      : characterstic frequencies of the model
%
%   Output parameters:
%        thetaOut  : the resulting "where" cue applied in forming of the
%                    binaural activity map
%        energyOut : the resulting "what" cue applied in forming of the
%                    binaural activity map
%
%   This method aims to account for the emphasized role of onsets in
%   localization by the auditory system, and for the variation found in the
%   time resolution of the localization. The steps involved in this method
%   are listed below and a more detailed describtion about the method can
%   be found in in Takanen, Santala, Pulkki 2013 (Sec. 3.2.7)
%   
%   1) analyze the cues with two mechanisms
%
%       50-ms long time window to evaluate the average information 
%        over a longer time frame
%       3-ms long time frame to analyze the gradients of the
%        cues in order to be sensitive to the onsets
%
%   2) combine the informations obtained with the two mechanisms to
%
%       obtain the overall "where" and "what" cues
%
%       emphasize the information obtained with the mechanism
%        employing shorter time frame to give more weight to the
%        onsets
%
%   See also: takanen2013, takanen2013_periphery, takanen2013_weightedaveragefilter,
%             takanen2013_formbinauralactivitymap
%
%   References:
%     M. Takanen, O. Santala, and V. Pulkki. Visualization of functional
%     count-comparison-based binaural auditory model output. Hearing
%     research, 309:147--163, 2014. PMID: 24513586.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/takanen2013_onsetenhancement.php


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


%% ------ Set parameters for the computations ------------------------------
dims = size(thetaIn);
%within the range of 1 kHz and 2 kHz, none of the directional cues
%are considered accurate, hence the what cues are set to zero at those
%frequencies
limit1 = find(cfs>1000,1,'first');
limit2 = find(cfs<2000,1,'last');
%the original energy is however stored as the what cue has content also in
%that frequency range
originalEnergy = energyIn;
energyIn(:,limit1:limit2) =0;

%compute the rms value of the input energy so that the output energy can be
%scaled to a similar level
rmsOrig = (sqrt((ones(1,size(energyIn,1))*(energyIn.^2))./size(energyIn,1)));
%scale the "what" cues based on average values obtained with a pink noise
%signal reproduced at 60 dB SPL
x=amt_load('takanen2013','periphenergyaverages.mat');
averageEnerg=x.averageEnerg;
scaledEnergy = energyIn./(ones(dims(1),1)*averageEnerg);

%two different kinds of window lengths are used, namely 3 ms and 50 ms
win2 = hann(floor(0.003*fs))./floor(0.003*fs);
win1 = hann(floor(0.05*fs))./floor(0.05*fs);

%a lowpass filter is used to compute the %envelope of the signal in order 
%to compute the gradient
F = [200 500]; %band limits
A = [1 0];% band type: 0='stop', 1='pass'
dev=[10^(0.1/20)-1 0.0001]; % ripple/attenuation spec
[M,Wn,beta,typ]= kaiserord(F,A,dev,fs);  % window parameters
b=fir1(M,Wn,typ,kaiser(M+1,beta),'noscale');

%the mechanism employing shorter time frame inspects the overall gradient 
%that is computed across different frequency bands by highlighting the
%information at the lowest frequencies. A frequency-dependent gain-factor
%*g* is applied in this process.
l = dims(2):-1:1;
g = -1+2./exp(-l);
g=100*g./g(2);g(g>100)=100;g(g<1)=1;

%% ------ 1) Analysis with the two mechanisms --------------------------------

% 1.1) compute the gradient information
%compute the envelope with lowpass filtering of the input
envelope = filter(b,1,scaledEnergy);
envelope = [envelope((length(b)/2):end,:);zeros(length(b)/2-1,dims(2))];
%derive the instantaneous derivative by
grad = envelope-[zeros(1,dims(2));envelope(1:end-1,:)];
clear envelope
%employ half-wave rectification so that the analysis is sensitive only to
%the onsets, not to the offsets
grad = grad.*(grad>0);

%the average energy over long time window and the short time derivative
%need to be adjusted to the same range. This is done using the rms values
%of those energies for a pink noise burst in a concert hall
x=amt_load('takanen2013','onsetmultp.mat');
coeff=x.coeff;
%compute the average cue across frequencies from the gradient by weighting
%the lowest frequencies more
tauOfShortFrame = conv((thetaIn.*(grad.*(ones(dims(1),1)*g)))*ones(dims(2),1),win2,'same'...
    )./(conv((grad.*(ones(dims(1),1)*g))*ones(dims(2),1),win2,'same')+1e-30);
%set the same average cue to all frequency bands
tauOfShortFrame = tauOfShortFrame*ones(1,dims(2));
%compute also the average energy across frequencies
energOfShortFrame = conv(grad*ones(dims(2),1),win2,'same')*ones(1,dims(2));
clear grad

%and scale the gradient energy using the coefficient
energOfShortFrame = energOfShortFrame.*(ones(dims(1),1)*coeff);

% 1.2) compute the information over a longer time period
%initialize the where and what cue obtained with this mechanism
tauOfLongFrame = zeros(dims);
energOfLongFrame=tauOfLongFrame;
for i=1:dims(2)
   tauOfLongFrame(:,i) = (conv((thetaIn(:,i).*scaledEnergy(:,i)),win1,'same')...
        )./(conv(scaledEnergy(:,i),win1,'same')+1e-30);
   energOfLongFrame(:,i) =  conv(scaledEnergy(:,i),win1,'same');
end
clear thetaIn scaledEnergy
%% ------ Combination of the informations obtained ------------------------
%the where cues are combined by computing a weighted average in which the
%energies are used as the weight. Additionally, the information obtained by
%inspecting the gradients is multiplied with a scalar five in order to
%empasize the onsets on the localization
envelopeWeight =5;
thetaOut = (tauOfShortFrame.*energOfShortFrame*envelopeWeight+tauOfLongFrame.*energOfLongFrame)...
    ./(energOfLongFrame+energOfShortFrame*envelopeWeight+1e-30);
clear energOfLongFrame tauOfLongFrame tauOfShortFrame % reduce memory requirements

%the rms value of the energy obtained with the mechanism employing shorter
%time frame is computed
rms2 = (sqrt((ones(1,size(energOfShortFrame,1))*(energOfShortFrame.^2))./size(energOfShortFrame,1)));
%and scaled to a similar level as the original input. However, for
%visualization purposes, the mentioned energy is multiplied by two so that
%the onsets would be more visible in the resulting binaural activity map
energOfShortFrame = energOfShortFrame.*(ones(dims(1),1)*(2*rmsOrig./(rms2+1e-30)));
%the overall "what cue is formed by combining the original "what" stream
%input to the energy information obtained with the gradient mechanism
energyOut = max(originalEnergy,energOfShortFrame);


