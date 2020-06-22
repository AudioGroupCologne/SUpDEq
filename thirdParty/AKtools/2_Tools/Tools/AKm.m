% [c, aMatched, bMatched] = AKm(a, b, operation, domain)
% matches the size of a and b, and applies element wise mathematical
% operations: e.g. if a is of size [N x M], and b is of size [1 N], b is
% repeatet to match a, i.e. b=repmat(b, [N 1]). c is than computet as
% c=a+B, c=a-b, c=a.*b, or c=a./b, according to operation:
%
% + AKm(a, b, '+')
% - AKm(a, b, '-')
% * AKm(a, b, '*')
% / AKm(a, b, '/')
% ^ AKm(a, b, '^')
%
% the operations can also be applied in different domains. In this case
% an FFT is performed before applying the operator and an IFFT is performed
% afterwards. This is done by appending 'cs' for operations on the complex
% spectra, or 'ms' for operations on the magnitude spectra, e.g.:
%
% AKm(a, b, '+', 'ms')
%
%
% I N P U T
% a, b       - scalar or matrix of any size (must be numeric)
% opearation - '+', '-', '*', '/', or '^'
% domain     - 'none' (default) the operatin is applied to a, and b without
%                  any transformation
%              'cs' the operation is applied on the complex spectra. In
%                  this case an FFT is performed before applying the
%                  operator and an IFFT is performed.
%              'ms' the operation is applied on the magnitude spectra
%
% O U T P U T
% c          - result 
% aMatched   - a after dimensionality matching
% bMatched   - b after dimensionality matching
%
% v1 7/2016 fabian.brinkmann@tu-berlin.de, Audio Communicatin Group,
%           TU Berlin

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
function [c, aMatched, bMatched] = AKm(a, b, operation, domain)

if nargin < 4
    domain = 'none';
end

% get size and dimensionality of input data
aS = size(a);
bS = size(b);
aD = numel(aS);
bD = numel(bS);

% check if any input is scalar
if aD>1 || bD>1
    
    % match dimension for dimension
    for nn = 1:max(aD, bD)
        
        % all input has minimally dimension nn
        if nn<=aD && nn<=bD
            if aS(nn)>bS(nn)
                if bS(nn)==1
                    
                    % match b to a
                    rep = ones(1,bD);
                    rep(nn) = aS(nn);
                    b = repmat(b, rep);
                    bS = size(b);
                    bD = numel(bS);
                    
                else
                    error('AKm:input', ['size(b,' num2str(nn) ') must be 1'])
                end
            elseif bS(nn)>aS(nn)
                if aS(nn)==1
                    
                    % match a to b
                    rep = ones(1,aD);
                    rep(nn) = bS(nn);
                    a = repmat(a, rep);
                    aS = size(a);
                    aD = numel(aS);
                    
                else
                    error('AKm:input', ['size(a,' num2str(nn) ') must be 1'])
                end
            end
            
        % b is of higher dimensionality
        elseif nn<=aD
            
            % match b to a
            rep = ones(1,nn);
            rep(nn) = aS(nn);
            b = repmat(b, rep);
            bS = size(b);
            bD = numel(bS);
            
        % a is of higher dimensionality
        elseif nn<=bD
            
            % match a to b
            rep = ones(1,nn);
            rep(nn) = bS(nn);
            a = repmat(a, rep);
            aS = size(a);
            aD = numel(aS);
            
        end
    end
    
end

aMatched = a;
bMatched = b;

% apply the domain
if strcmpi(domain, 'cs')
    a = fft(a);
    b = fft(b);
elseif strcmpi(domain, 'ms')
    a = abs(fft(a));
    b = abs(fft(b));
end

% do the math
switch  operation
    case '+'
        c = a+b;
    case '-'
        c = a-b;
    case '*'
        c = a.*b;
    case '/'
        c = a./b;
    case '^'
        c = a.^b;
end

% re-apply the domain
if strcmpi(domain, 'cs')
    c = ifft(c);
elseif strcmpi(domain, 'ms')
    c = circshift(ifft(c), round(size(c,1)/2));
end