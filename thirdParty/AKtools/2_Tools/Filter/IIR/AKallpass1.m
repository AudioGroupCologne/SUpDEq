% For documentation see AKfilter.m, and AKfilterDemo.m
% A python implementation of most filters can be found here:
% http://nbviewer.jupyter.org/github/sonible/nb/blob/master/iir_filter/index.ipynb
%
% Frank Schultz, FG Audiokommunikation, TU Berlin
% frank.schultz@tu-berlin.de, +49 175 15 49 763, Skype: j0shiiv
% 0.00 06.09.2010 init dev
% see AKcalcBiquadCoeff.m for details
% to be checked!!!

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function [b,a]=AKallpass1(fg,fs,ai)
    wg=2*fs*tan(pi*fg/fs); %pre-warping
    B(1) = 0;			
    B(2) = -ai/wg;		
    B(3) = 1;
    A(1) = 0;	
    A(2) = +ai/wg;	
    A(3) = 1;
    [b,a] = AKcalcBiquadCoeff(B,A,fs);
end