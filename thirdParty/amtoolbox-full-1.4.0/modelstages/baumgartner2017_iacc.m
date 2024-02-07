function [iacc,iacc_fc,fc] = baumgartner2017_iacc(binsig,varargin)
%BAUMGARTNER2017_IACC calculates the IACC
%   Usage:     [iacc,iacc_fc,fc] = baumgartner2017_iacc(binsig)
%
%   Input parameters:
%     binsig  : binaural time-domain signal (dimensions: time x channel)
%
%   Output parameters:
%     iacc    : normalized interaural cross-correlation coefficient (IACC)
%     iacc_fc : frequency-specific IACC
%     fc      : center frequencies
%
%   BAUMGARTNER2017_IACC calculates the interaural cross correlation
%   coefficient from a binaural signal.
%
%   See also: baumgartner2017
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/baumgartner2017_iacc.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: SOFA M-Stats M-Signal O-Statistics
%   #Author: Robert Baumgartner (2017), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.import={'baumgartner2017'};
definput.keyvals.maxlag = 0.001; % in sec

[flags,kv]=ltfatarghelper({},definput,varargin);

[mp,fc] = auditoryfilterbank(binsig(:,:),kv.fs,'flow',kv.flow,'fhigh',kv.fhigh);
if flags.do_ihc
  mp = ihcenvelope(mp,kv.fs,'ihc_dau1996');
end

maxlagn = round(kv.maxlag*kv.fs);
iacc_fc = nan(length(fc),1);
intens_fc = nan(length(fc),1);
for ifreq = 1:length(fc)
  intens_fc(ifreq) = sqrt(prod(sum(mp(:,ifreq,:).^2,1),3));
  iacnorm = xcorr(mp(:,ifreq,1),mp(:,ifreq,2),maxlagn,'coeff');
  iacc_fc(ifreq) = max(abs(iacnorm));
end
intWeight = intens_fc/sum(intens_fc); % intensity weighting
iacc = iacc_fc'*intWeight; % weighted average

end


