function [y_est] = llado2022_evaluatenn(x_test,NN_pretrained)
%LLADO2022_EVALUATENN evaluate the neural network
%   Usage: [y_est] = llado2022_evaluatenn(x_test,NN_pretrained)
%
%   Input parameters:
%     x_test             : Features of the test subset
%     NN_pretrained      : Pretrained network
%     hiddenLayerSize    : Size of the hidden layer
%
%   Output parameters:
%     y_est              : Estimated data
%
%   LLADO2022_EVALUATENN gives the estimation uncertainty of the neural
%   network
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/llado2022_evaluatenn.php


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

for iter = 1:NN_pretrained.nIter
    net_dir = NN_pretrained.preTrained_dir(1,end,iter).net;
    net_uncertainty = NN_pretrained.preTrained_uncertainty(1,end,iter).net;

    y_hat_dir(iter,:) = net_dir(x_test);
    y_hat_uncertainty(iter,:) = net_uncertainty(x_test);

end
clear clipPos;
clipPos = find(y_hat_dir < -90);
y_hat_dir(clipPos) = -90;
clear clipPos;
clipPos = find(y_hat_dir > 90);
y_hat_dir(clipPos) = 90;

y_est(:,1) = mean(y_hat_dir);
y_est(:,2) = mean(y_hat_uncertainty);
end


