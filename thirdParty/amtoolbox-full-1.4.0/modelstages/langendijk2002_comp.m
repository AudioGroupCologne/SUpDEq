function [ si ] = langendijk2002_comp( target,template,varargin )
%LANGENDIJK2002_COMP Comparison process according to Langendijk et al. (2002)
%   Usage:        [ si ] = langendijk2002_comp( target,template )
%                 [ si ] = langendijk2002_comp( target,template,s,do,flags )
%
%   Input parameters:
%     target   :  modified DFT for one target position and both ears
%     template :  internal DTF-templates for all positions and both ears
%
%   Output parameters:
%     si       :  monaural similarity indices (one for each template entry)
%
%   LANGENDIJK2002_COMP compares the spectral features of an internal spectral
%   representation of a target sound with an internal template. As spectral
%   similarity criterion either the crosscorrelation coefficient ('xcorr') 
%   or a mapped version of the standard deviation of inter-spectral 
%   differences ('std') can be used.
%
%   The function accepts the following optional parameters:
%
%     's',s     Standard deviation of transforming Gaussian function 
%                (only for comparison process: std).
%                The default value is 2.
%
%     'do',do    Differential order. The default value is 0.
%
%     'std'      Use the 'std' comparison process. This is the default.
%
%     'xcorr'    Use the 'xcorr' comparison process. 
%
%   See also: langendijk2002
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/langendijk2002_comp.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: none
%   #Author : Robert Baumgartner (2013), OEAW Acoustical Research Institute

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


definput.import={'langendijk2002_comp'};
[flags,kv]=ltfatarghelper({'s','do'},definput,varargin);


%% Optional differentiation

if kv.do ~= 0
	target = diff(target,kv.do);
    template = diff(template,kv.do);
end


%% Comparison process

if flags.do_std
    
    idiff = repmat(target,[1,size(template,2),1]) - template;
    sigma = squeeze(std(idiff,1));
    si = exp(-0.5 * (sigma./kv.s).^2) ./ (sqrt(2*pi) .* kv.s); % normpdf
    
elseif flags.do_xcorr
    
    r = zeros(size(template,2),size(template,3));
    for ch = 1:size(template,3)
        for ii = 1:size(template,2)
            tmp = corrcoef(target(:,ch),template(:,ii,ch));
            r(ii,ch) = tmp(2);
        end
    end
    si = (r+1)/2;
    
end

end


