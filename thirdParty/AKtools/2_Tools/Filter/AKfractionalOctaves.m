% f_oct = AKfractionalOctaves(frac, f, fs)
% returns center frequencies, bandwidth and upper and lower cut-off
% frequencies for fractional octave filters according to DIN EN 61260.
%
%
% for example
% AKfractionalOctaves(3)
% returns frequencies of third octave filters
%
%
% I N P U T
% frac  - fraction of octave, e.g. 1 for octave frequencies, 3 for third
%         octave frequencies, etc. (default = 3)
% f     - two element vector specifiying the lower and upper limit of the
%         center frequencies in Hz (default = [32 16000]
% fs    - sampling rate in Hz. This is used to test if the highest upper
%         cut-off frequency of a filter is above fs/2. In case it is, it is
%         set to fs/2, i.e. the last filter band might have a smaller band-
%         width (default = 44100)
%
%
% O U T P U T
% f_oct - [N x 4] matrix that holds the N center frequencies in the first
%         column, the lower and upper cut-off frequencies in the second and
%         third column, and the bandwidth in the fourth column.
%
%
% 2018/03 - fabian.brinkmann@tu-berlin.de

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
function f_oct = AKfractionalOctaves(frac, f, fs)

% default input parameters
if nargin <= 2
    fs = 44100;
end
if nargin <= 1
    f = [32 16e3];
end
if nargin == 0
    frac = 3;
end

% get fractional octave center frequencies
if mod(frac,2)
    % DIN EN 61 260 Eq. (3)
    f_oct = 10^(3/10).^((-1000:500)'/frac)*1000;
else
    % DIN EN 61 260 Eq. (4)
    f_oct = 10^(3/10).^((2*(-1000:500)'+1)/(2*frac))*1000;
end

% restrict range of center frequencies according to user input
f_min = 10^(3/10).^(-1/(2*frac)).*min(f);
f_max = 10^(3/10).^( 1/(2*frac)).*max(f);
f_oct = f_oct(f_oct>=f_min & f_oct<=f_max);

% upper and lower cut-off frequency: DIN EN 61 260 Eq. (5-6)
f_oct(:,2) = 10^(3/10).^(-1/(2*frac)).*f_oct(:,1);
f_oct(:,3) = 10^(3/10).^( 1/(2*frac)).*f_oct(:,1);

% check if lowest cut-off freq. is > 0 and highest is < fs/2
id = find(f_oct(:,2)>0, 1, 'first');
f_oct = f_oct(id:end,:);
id = find(f_oct(:,3)<fs/2, 1, 'last');
if id == size(f_oct,1)
    f_oct = f_oct(1:id,:);
else
    f_oct = f_oct(1:id+1,:);
    f_oct(end, 3) = fs/2;
end

% band width
f_oct(:,4) = f_oct(:,3) - f_oct(:,2);