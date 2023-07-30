function outsig = takanen2013_weightedaveragefilter(insig,weight,fs,timeconst)
%TAKANEN2013_WEIGHTEDAVERAGEFILTER Compute the weighted or self-weighted average
%                 
%   Usage: outsig = takanen2013_weightedaveragefilter(insig,weight,fs,timeconst)
%
%   Input parameters:
%        insig     : signal from which the average is computed. Optionally,
%                    can be the same as the weight resulting in
%                    self-weighted average
%        weight    : signal to be used as weight in the computation
%        fs        : sampling rate
%        timeconst : time constant specifying the first-order IIR filter
%
%   Output parameters:
%        outsig    : output signal
%
%   This function computes the either the weighted or the self-weighted
%   average of the input signal using a first-order IIR filter whose time
%   constant is specified by the timeconst argument. More details about
%   the conputation can be found in Pulkki, Hirvonen 2009 (Sec. 3.2.3)
%
%   See also: takanen2013_mso, takanen2013_lso, takanen2013_wbmso,
%             takanen2013_onsetenhancement, takanen2013_contracomparison
%
%   References:
%     V. Pulkki and T. Hirvonen. Functional count-comparison model for
%     binaural decoding. Acta Acustica united with Acustica, 95(5):883 --
%     900, Sept./Oct. 2009.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Perceptually encoded signals and
%     their assessment. In J. Blauert, editor, The technology of binaural
%     listening. Springer, 2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/takanen2013_weightedaveragefilter.php


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


%% ------ Computation ----------------------------------------------------

%setting the parameters for a first-order IIR filter
B = 1-exp(-1/(fs*timeconst));
A = [1 -exp(-1/(fs*timeconst))];

%filter the weighted and self-weighted signals 
weighted = filter(B,A,(insig.*(weight.^2)));
selfWeighted = filter(B,A,(weight.^2));

%derive the output
outsig = weighted./(selfWeighted+1e-30);


