% correctedPhase = AKcorrectPhase(phase,fs);
% applies phase correction as suggested by [1, p.38]
%
% this function is called from AKsweepFD.m
%
% WRAPPED phase expected!!!
% length = N/2+1 expected!!! dc..nqy part from even order FFT spectrum
%
% [1] Mueller, S. & Massarani, P. (2001): "Transfer-Function Measurement
%     with Sweeps. Directors Cut Including Previously Unreleased Material
%     And Some Corrections." J. Audio Eng. Soc., 49(6), 443-471.
%
% Andre Giese, Audio Communication Group, TU Berlin

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
function correctedPhase = AKcorrectPhase(phase,fs)
len=length(phase);

%get Parameters needed for calculation
phi_end=phase(end);
N=(len-1)*2;
df=fs/N;

%build a vector of length(phase) that keeps the correction
temp=ones(len,1);
temp=cumsum(temp)-1;
offset=df*phi_end/(fs/2);
correc=temp*offset;
%CORRECTING:
correctedPhase=phase-correc;


%DISPLAYING last 3 values
if 0
    fprintf('\tPhase @ NYQ-2: %10.4f  \n',correctedPhase(len-2));
    fprintf('\tPhase @ NYQ-1: %10.4f  \n',correctedPhase(len-1));
    fprintf('\tPhase @ NYQ: %10.4f  \n',correctedPhase(len));
end
