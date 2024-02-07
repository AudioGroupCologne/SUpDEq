% str = AKlistIntegers(vector,do_sort)
% takes an integer vector and builds a well readable string, by combining
% successive numbers via a "..".
%
% e.g. use
% AKlistIntegers([1:3,55:56,101,102:103])
%       = 1..3,55,56,101..103
%
% I N P U T
% vector     - integer vector (either [1,n] or [n,1])
% do_sort    - sort values before generating the list, default = true
%
% O U T P U T
% str        - generated character vector 
%
% v1 01/2018 helmholz@campus.tu-berlin.de, Audio Communicatin Group,
%            TU Berlin

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
function str = AKlistIntegers(vector,do_sort)

if nargin < 2
    do_sort = true;
end
if iscell(vector)
    vector = [vector{:}];
end
if ~isvector(vector)
    error('AKlistIntegers:input',...
        'Only row or column vectors allowed as input.');
end
if any(mod(vector,1))
    error('AKlistIntegers:input',...
        'Only Integers allowed as input.');
end

if do_sort
    vector = sort(vector);
end

str = []; %#ok<*AGROW>
i = 0;
while i < length(vector)
    i = i + 1;
    
    if ~isempty(str)
        str = [str,','];
    end
    
    str = [str,num2str(vector(i))];
    
    i_end = i;
    while i_end < length(vector)
        i_end = i_end + 1;
        if vector(i_end) ~= vector(i) + i_end - i
            i_end = i_end - 1;
            break;
        end
    end
    
    if i_end > i
        if i_end > i + 1
            str = [str,'..'];
        else
            str = [str,','];
        end
        str = [str,num2str(vector(i_end))];
        i = i_end;
    end
end

end
