function [] = llado2022_weightsanalysis(NN_pretrained)
%LLADO2022_WEIGHTSANALYSIS analyses the neural networks' weights
%   Usage: llado2022_weightsanalysis(NN_pretrained);
%
%   Input parameters:
%     NN_pretrained : Pretrained network
%
%   llado2022_weightsAnalysis analyses the weights learnt by the NN and plots
%   them to understand the importance of the training features.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/llado2022_weightsanalysis.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB M - Communication Systems
%   #Author: Pedro Llado (2022)
%   #Author: Petteri Hyvärinen (2022)
%   #Author: Ville Pulkki (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% NN weights analysis: perceived direction
clear A B T TOTAL TOTALavg
for j = 1:8
    for i = 1:10
        A(:,:) = abs(NN_pretrained.preTrained_dir(1,j,i).net.IW{1}(:,:))';
        B(:) = abs(NN_pretrained.preTrained_dir(1,j,i).net.LW{2}(:,:))';
        T(i,:) = mean((A.*B)');
    end



    TOTAL(j,:) = mean(T);
end
TOTALavg = mean(TOTAL);
TOTALavg = TOTALavg./(sum(TOTALavg));
figure(3);

plot(TOTALavg(1:18),'b');
hold on;
figure(4);

plot(TOTALavg(19:end),'b');
hold on;


% NN weights analysis: position uncertainty
clear A B T TOTAL TOTALavg

for j = 1:8
    for i = 1:10
        A(:,:) = abs(NN_pretrained.preTrained_uncertainty(1,j,i).net.IW{1}(:,:))';
        B(:) = abs(NN_pretrained.preTrained_uncertainty(1,j,i).net.LW{2}(:,:))';
        T(i,:) = mean((A.*B)');

    end

    TOTAL(j,:) = mean(T);
end


TOTALavg = mean(TOTAL);
TOTALavg = TOTALavg./(sum(TOTALavg));
figure(3);
plot(TOTALavg(1:18),'r');
ylim([0.02 0.04])
xlim([0 19])
title("ITD weights")
legend("Perceived direction estimation", "Position uncertainty estimation",'Location','Southeast')
figure(4);
plot(TOTALavg(19:end),'r');
ylim([0.02 0.04])
xlim([0 19])
title("ILD weights")
legend("Perceived direction estimation", "Position uncertainty estimation",'Location','Southeast')
end


