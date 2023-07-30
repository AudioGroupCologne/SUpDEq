function output = takanen2013_lso(ipsilateral, contralateral, fs, fc)
%TAKANEN2013_LSO Model of the lateral superior olive
%   Usage: output = takanen2013_lso(ipsilateral, contralateral, fs, fc);
%
%   Input parameters:
%        ipsilateral   : The ipsilateral "where" stream output from the
%                        model of the periphery
%        contralateral : The contralateral "where" stream output from the
%                        model of the periphery
%        fs            : sampling rate
%        fc            : characteristic frequencies
%
%   Output parameters:
%        output : Spatial cues for separate narrow bandwidths
%
%   This function models the lateral superior olive (LSO) by processing the 
%   output of the periphery model with the following steps:
%
%   1) First-order lowpass filter is applied to both the ipsilateral and
%      contralateral sides, and the contralateral side is delayed
%
%   2) The output is saturated at a certain dB level and limited
%      signals.
%
%   3) Weighted moving average filters with a short and longer time
%      constant are applied, the latter only below 1 kHz in order to slow
%      the response
%
%   See also: takanen2013, takanen2013_periphery, takanen2013_weightedaveragefilter
%
%   References:
%     V. Pulkki and T. Hirvonen. Functional count-comparison model for
%     binaural decoding. Acta Acustica united with Acustica, 95(5):883 --
%     900, Sept./Oct. 2009.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Visualization of functional
%     count-comparison-based binaural auditory model output. Hearing
%     research, 309:147--163, 2014. PMID: 24513586.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/takanen2013_lso.php


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

%% ------ The contralateral ear input is delayed by 0.2 ms ----------------
contradelay = round(0.0002*fs);
contralateral = [zeros(contradelay,size(contralateral,2));...
    contralateral(1:size(contralateral,1)-contradelay,:)];

%% ------ First-order lowpass filter with a time constant of 0.1 ms -------
tau = 0.0001; % seconds
B = 1-exp(-1/(fs*tau));
A = [1 -exp(-1/(fs*tau))];
ipsilateral = filter(B,A,ipsilateral);
contralateral = filter(B,A,contralateral);

%% ------ Level difference computation and limitation ---------------------
% The output is saturated at ILD of 18 dB
output = (10^(-18/20))*((ipsilateral)./(contralateral+1e-20));
% Limitation
output= output.*(output>0);
output(output>1) = 1;

%% ------ Spreading of short peaks ----------------------------------------
windowfunct = hann(round(0.001*fs));
windowfunct = windowfunct./sum(windowfunct);
for i=1:size(contralateral,2)
    output(:,i) = conv(output(:,i),windowfunct,'same');
end

%% ------ Weighted-moving average with a time constant of 4 ms ----------------
output = takanen2013_weightedaveragefilter(output,ipsilateral,fs,0.004);
% tau = 0.004; 
% B = 1-exp(-1/(tau*fs));A = [1 -exp(-1/(tau*fs))];
% output = (filter(B,A,(ipsilateral.^2).*output))./(filter(B,A,ipsilateral.^2)+1e-20);

%% ------ Self-weighted moving average ----------------
% Slow down the response above 1 kHz with a time constant of 50 ms
limit = find(fc<1000,1,'last');
output(:,1:limit) = takanen2013_weightedaveragefilter(output(:,1:limit),output(:,1:limit),fs,0.05);


