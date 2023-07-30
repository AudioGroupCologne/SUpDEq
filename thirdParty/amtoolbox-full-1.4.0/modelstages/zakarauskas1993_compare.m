function [ out ] = zakarauskas1993_compare( in1,in2,do,bal )
% ZAKARAUSKAS1993_COMPARE Comparison process according to Zakarauskas et al. (1993)
%
%   Usage:      
%     [ out ] = zakarauskas1993_compare( in1,in2 )
%     [ out ] = zakarauskas1993_compare( in1,in2,do )
%     [ out ] = zakarauskas1993_compare( in1,in2,do,bal )
%
%   Input parameters:
%     in1:      modified DFT for one target position and both ears
%     in2:      stored DFT-templates for all positions and both ears
%     do:       differential order; default: 2
%     bal:      balance of left to right channel; default: 1
%
%   Output parameters:
%     out:      summed differences
%
%     This function compares two directional transfer functions with each other.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/zakarauskas1993_compare.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author: Robert Baumgartner (2010), OEAW Acoustical Research Institute

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


% default settings
if ~exist('do','var')
    do=2;
end
if ~exist('bal','var')
    bal=1;
end

ptemp=zeros(size(in2,2),2); % initialisation
for ch=1:2
    for ind=1:size(in2,2)
        z=diff(in1(:,ch),do)-diff(in2(:,ind,ch),do);
        ptemp(ind,ch)=sum(abs(z));
    end
end
p=(bal*ptemp(:,1)+1/bal*ptemp(:,2))/2; % balance
out=p;
end


