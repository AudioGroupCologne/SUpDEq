% [y, h] = AKfractionalDelayCyclic(x, d, N, kaiserA)
%
% fractional delay of input vector x using a sinc filter [1].
% Delay is split into integer and fractional part, thus, the filter order N
% can be choosen independently from the delay d.
%
% I N P U T (default):
% x                - impulse responses of size [N M C], where N ist the
%                    number of samples, M the number of measurements, and C
%                    the number of channels.
% d                - fractional delays of size [C M]
% N (30)           - filter order (length = N+1)
% kaiser_A (60)    - side lobe rejection of kaiser window applied to sinc
%                    filter
%
% O U T P U T:
% y                - delayed output data
% h                - fractional delay filters. Will be zero if d is an
%                    integer value
%
% [1] Laakso et al. (1996): "Splitting the unit delay." In: IEEE Signal
%     Processing Magazine. 13(1):30-60.
%
% 04/2012 - fabian.brinkmann@tu-berlin.de
% 05/2018 - fabian.brinkmann@tu-berlin.de (added cyclic fractional delaying)

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
function [y, h] = AKfractionalDelayCyclic(x, d, N, kaiserA)

% -------------------------------------------- 0. define default parameters
if ~exist('N', 'var')
    N = 30;
end
if ~exist('kaiserA', 'var')
    kaiserA = 60;
end

% ------------------------------------------------------ 1. check the input
if size(d,1)==1 || size(d,2) == 1
    d = reshape(d, [1 numel(d)]);
end

if size(d,2)~=size(x,2) || size(d,1)~=size(x,3)
    error('AKfractionalDelay:Input', ['d must be of size [' num2str(size(x,3)) ' ' num2str(size(x,2)) ']'])
end

% -------------------------------------------- 2. handle multi channel data
if size(d,2)>1 || size(x,3)>1
    % allocate space
    h = zeros(N+1, size(x,2), size(x,3));
    y = zeros(size(x));
    
    % loop across measurements and channel
    for mm = 1:size(x,2)
        for cc = 1:size(x,3)
            [tmp, h(:,mm,cc)] = AKfractionalDelayCyclic(x(:,mm,cc), d(cc,mm), N, kaiserA);
            y(1:numel(tmp),mm,cc) = tmp;
        end
    end
    
    return
end


% ------------------------------------- 3. realise fractional part of delay
[L, C] = size(x);

    
if rem(abs(d), 1)~=0 % fractional delay
    if d > 0
        d_frac = d-floor(d);
    else
        d_frac = 1-rem(abs(d), 1);
    end
    
    % get discrete time vector
    % [1], eq. 21
    if ~mod(N, 2)
        % even N
        M_opt = round(d_frac)-N/2;
    else
        % odd N
        M_opt = floor(d_frac)-(N-1)/2;
    end
    
    n = M_opt-d_frac:M_opt-d_frac+N;
    
    % calculate filter
    h = sinc(n)';
    % get kaiser window
    kaiserWin = AKkaiser(N+1, kaiserA, d_frac);
    % apply kaiser window shifted by d_frac
    h  = h .* kaiserWin;
    
    % repeat channels to match input signal
    h = repmat(h, [1, C]);
    
    % filter (introduces additional delay)
    y = AKm(x, [h; zeros(L-N-1, C)], '*', 'cs');
    
    % remove delay from filtering, determined by filter order, d and d_frac
    if ~mod(N, 2)
        % even N
        if d > 0
            if d_frac < 0.5
                y = circshift(y, [-N/2 0]);
            else
                y = circshift(y, [-N/2+1 0]);
            end
        else
            if d_frac < 0.5
                y = circshift(y, [-N/2-1 0]);
            else
                y = circshift(y, [-N/2 0]);
            end
        end
    else
        % odd N
        if d > 0
            y = circshift(y, [-ceil(N/2)+1 0]);
        else
            y = circshift(y, [-ceil(N/2) 0]);
        end
    end
else
    y = x;
end

% ---------------------------------------- 4. realise integer part of delay
if d > 0
    y = circshift(y, [floor(d) 0]);
elseif d < 0
    y = circshift(y, [ceil(d) 0]);
else
    y = x;
end

% --- 5. generate zero filter output, if we did not apply fractional delays
if ~exist('h', 'var')
    h = zeros(N+1,1);
end