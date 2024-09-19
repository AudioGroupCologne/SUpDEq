% clMsg = AKdisp(Msg, lastMsg)
% Prints the string Msg to the command window and erases the last message
% passed by argument lastMsg. Output argument clMsg is identical to Msg.
%
% Can be used for logging process in the command window
% Call AKdisp() for a demonstration
%
% 3/2015 - fabian.brinkmann@tu-berlin.de

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
function clMsg = AKdisp(Msg, lastMsg)

if nargin == 0
    clMsg = '';
    for n = 1:100
        Msg = ['AKdisp is counting ' num2str(n) '/100'];
        clMsg = repmat(sprintf('\b'), 1, numel(clMsg));
        fprintf([clMsg Msg])
        clMsg = Msg;
        pause(.05)
    end
elseif nargin == 1
    fprintf(Msg)
    clMsg = Msg;
    
    if nargout == 0
        clear clMsg
    end
elseif nargin == 2
    clMsg = repmat(sprintf('\b'), 1, numel(lastMsg));
    fprintf([clMsg Msg])
    clMsg = Msg;
    
    if nargout == 0
        clear clMsg
    end
else
    error('clDisp:nargin', 'Wrong number of input arguments')
end