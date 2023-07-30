function [net,tr] = llado2022_trainnn(x,y,hiddenLayerSize,augmentation_ratio,SNR)
%LLADO2022_TRAINNN trains the neural network
%   Usage: [net,tr] = llado2022_trainnn(x,y,hiddenLayerSize,augmentation_ratio,SNR);
%
%   Input parameters:
%     x                  : Features of the train subset
%     y                  : Labels for training the network
%     hiddenLayerSize    : Size of the hidden layer
%     augmentation_ratio : Ratio for data augmentation stage
%     SNR                : SNR of the augmented data
%
%   Output parameters:
%     net                : trained network
%     tr                 : training history
%
%   LLADO2022_TRAINNN trains the neural network
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/llado2022_trainnn.php


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
clear net;
net = fitnet(hiddenLayerSize);

net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 30/100;
net.divideParam.testRatio = 0/100;

%% Training data augmentation

Y_output_aug = y;
for aug_iter = 1:augmentation_ratio
    for id_col = 1:length(y(1,:))
        aux= y(:,id_col);
        Y_output_aug((aug_iter)*length(y(:,1))+1:(1+aug_iter)*length(y(:,1)),id_col) = aux;
    end
end

clear X_input_aug;
X_input_aug(:,:) = x(:,:);
for aug_iter = 1:augmentation_ratio
    for id_col = 1:length(x(1,:))
        %aux= awgn(x(:,id_col),SNR,'measured');
        auxnoise = randn(size(x(:,id_col)));
        aux = x(:,id_col) + scaletodbspl(auxnoise, dbspl(x(:,id_col)) - SNR);
        X_input_aug((aug_iter)*length(x(:,1))+1:(1+aug_iter)*length(x(:,1)),id_col) = aux;
    end
end

[net, tr] = train(net,X_input_aug',Y_output_aug');

end



