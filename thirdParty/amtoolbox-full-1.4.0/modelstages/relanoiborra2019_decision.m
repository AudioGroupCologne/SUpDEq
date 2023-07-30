function [out] = relanoiborra2019_decision(clean_mfb, noisy_mfb, fs, fc, fc_mod,varargin)
%RELANOIBORRA2019_DECISION Decision stage of Relano-Iborra et al. (201)
%   Usage: [out, varargout] = relanoiborra2019_decision(clean_mfb, noisy_mfb, fs, fc, fc_mod,varargin)
%
%   This script builds the internal representations of the template and target signals
%   according to the CASP model (see references).
%   The code is based on previous versions of authors: Torsten Dau, Morten
%   L�ve Jepsen, Boris Kowalesky and Peter L. Soendergaard
%
%   The model has been optimized to work with speech signals, and the
%   preprocesing and variable names follow this principle. The model is
%   also designed to work with broadband signals. In order to avoid undesired
%   onset enhancements in the adaptation loops, the model expects to recive a
%   prepaned signal to initialize them.
%
%   Input parameter:
%
%      clean      :  clean speech template signal
%      noisy      :  noisy speech target signal
%      fs         :  Sampling frequency
%      flow       :  lowest center frequency of auditory filterbank
%      fhigh      :  highest center frequency of auditory filterbank
%      BMtype     :  model design: 'GT' gammatone (PEMO) or 'drnl'(CASP)
%      N_org      :  length of original sentence
%      sbj        :  subject profile for drnl definition
%
%   Outputs:
%      out           : correlation metric structure inlcuding:
%          .dint      : correlation values for each modulation band
%          .dsegments : correlation values from each time window and mod. band.
%          .dfinal    : final (averaged) correlation
%
%   References:
%     H. Relaño-Iborra, J. Zaar, and T. Dau. A speech-based computational
%     auditory signal processing and perception model. J. Acoust. Soc. Am.,
%     146(5), 2019.
%     
%     M. Jepsen, S. Ewert, and T. Dau. A computational model of human
%     auditory signal processing and perception. J. Acoust. Soc. Am., 124(1),
%     2008.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/relanoiborra2019_decision.php


%   #Author: Helia Relano Iborra (March 2019): v4.0 provided to the AMT team
%   #Author: Clara Hollomey (2021): adapted to the AMT
%   #Author: Piotr Majdak (2021): adapted to the AMT 1.0
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: M-Stats M-Signal

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Auditory filtering:
if isoctave
   warning(['Currently this model is only fully functional under MATLAB.']); 
end

definput.import={'relanoiborra2019'};
[flags,kv]  = ltfatarghelper({},definput,varargin);

if size(clean_mfb)~=size(noisy_mfb)
    error('internal representations must have the same size');
end

[Nsamp Naud_ch Nmod_ch] = size(noisy_mfb);

WinDurs = 1./fc_mod; % The window duration is the inverse of the centerfrequency of the modulation channel
WinDurs(1) = 1/2.5;

WinLengths = floor(WinDurs * fs);
Nsegments = floor(Nsamp./WinLengths)+ ones(1,length(fc_mod)); % The total number of segments is Nframes plus any additional "leftover"
    
      
if find(WinLengths == Nsamp)% If the duration of the stimulus is exactly equal to the window duration
        segIdx = find(WinLengths == Nsamp);
        Nsegments(segIdx) =  Nsegments(segIdx)-1;
end

for m=1:length(fc_mod)
    
    rule4th=find(fc > 4*fc_mod(m)); % Apply rule of cf_mod < 1/4 cf_aud
   
    % Build the modulation 2D matrices:                                       %
    speech=squeeze(clean_mfb(:, :, m));
    mix=squeeze(noisy_mfb(:, :, m));                                  
                                      
    % Delete discarded bands
     speech=speech(:, rule4th);
     mix=mix(:, rule4th);    
    
    tmp_ssnn = zeros(WinLengths(m), size(mix, 2), Nsegments(m)); % Allocate memory for multi-resolution
    tmp_ss = tmp_ssnn;
        
    segLengths = zeros(1,Nsegments(m)) ;
    
    % Find starting and ending points of the segments:
    
    for i = 1:Nsegments(m) % For each temporal segment of the signal
                               % find the start and end index of the frame
      if i > (Nsegments(m)-1)
          startIdx = 1 + (i-1)*WinLengths(m);
          endIdx = Nsamp;
      else
          startIdx = 1 + (i-1)*WinLengths(m);
          endIdx = startIdx + WinLengths(m)-1;
      end

      segment = startIdx:endIdx;
      segLengths(i) = length(segment);

      % internal representation of the temporal segments (samplesPerSegment x all bands x  number of segments)

      tmp_ss(1:segLengths(i), :, i) = speech(segment,:);
      tmp_ssnn(1:segLengths(i), :, i) = mix(segment,:);

      dint_temp(i, m) = corr2(tmp_ss(1:segLengths(i), :, i), tmp_ssnn(1:segLengths(i),:, i));
      dint_temp(dint_temp<0)=0; % Remove negative correlations
  
    end
dmod(m)=nansum(dint_temp(1:Nsegments(m), m))/Nsegments(m);

end

out.dint=dmod;
out.dsegments=dint_temp;
out.dfinal=mean(dmod);


