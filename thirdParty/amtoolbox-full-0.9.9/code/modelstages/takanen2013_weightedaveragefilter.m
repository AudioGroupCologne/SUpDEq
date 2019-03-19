function outsig = takanen2013_weightedaveragefilter(insig,weight,fs,timeconst)
%takanen2013_weightedaveragefilter Compute the weighted or self-weighted average
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
%     binaural decoding. Acta Acustica united with Acustica, 95(5):883 - 900,
%     Sept./Oct. 2009.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Perceptually encoded signals and
%     their assessment. In J. Blauert, editor, The technology of binaural
%     listening. Springer, 2013.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/takanen2013_weightedaveragefilter.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland


%% ------ Computation ----------------------------------------------------

%setting the parameters for a first-order IIR filter
B = 1-exp(-1/(fs*timeconst));
A = [1 -exp(-1/(fs*timeconst))];

%filter the weighted and self-weighted signals 
weighted = filter(B,A,(insig.*(weight.^2)));
selfWeighted = filter(B,A,(weight.^2));

%derive the output
outsig = weighted./(selfWeighted+1e-30);
