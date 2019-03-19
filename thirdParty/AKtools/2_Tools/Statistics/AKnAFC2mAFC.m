% hM = AKnAFC2mAFC(hN, N, M)
% transfroms detection rates from N-AFC tests to detection rates from M-AFC
% tests. With N-AFC tests, detection rates of 1/N represent chance level
% (i.e. guessing) and detection rates of 1/N+1/(2N) the threshold of
% perception (per definition).
%
% to obtain detection rates adjusted for the guessing probability call
% hM = AKnAFC2mAFC(hN, N, 1)
%
% I N P U T:
% hN - detection rates between 0 and 1. Can be scalar, vector or matrix of
%      any size, and will be transformed to a vector of size [1 numel(hN)]
% N  - Specifies the forced choice test underlying hN, e.g. N=2 would mean
%      that hN are detection rates from a 2-AFC test.
% M  - Specifies the target force choice test, e.g. M=3 will transform hN
%      to detection rates obtained from a 3-AFC test; M=1 will correct hN
%      for the guessing probability. M can be a scalar or vector.
%
% O U T P U T:
% hM - transformed probabilities of size [numel(M) numel(hN)]
%
% 11/2016  -  fabian.brinkmann@tu-berlin.de

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function hM = AKnAFC2mAFC(hN, N, M)

% reshape the input
hN = reshape(hN, [1 numel(hN)]);

% allocate output
hM = zeros(numel(M), numel(hN));

% transform to pure detection rate corrected for guessing probability
if N~=1
    h = (hN -1/N) * (1/(1-1/N));
else
    h = hN;
end

% transfomr to M-AFC detection rates
for mm = 1:numel(M)
    if M(mm) == 1
        hM(mm,:) = h;
    else
        hM(mm,:) = h*(1-1/M(mm)) + 1/M(mm);
    end
end