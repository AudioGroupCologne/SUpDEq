% AKdependencies(type)
% checks the availablilty of third party code and throws an error if it
% does not exist in the Matlab search path.
%
% INPUT
% type - string or cell array containing multiple strings containing
%        {'auditory' 'FABIAN'}
%
% 10/2016 - fabian.brinkmann@tu-berlin.de

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
function AKdependencies(type)

    if iscellstr(type)
        for nn = 1:numel(type)
            AKcheckDep(type{nn})
        end
    elseif ischar(type)
        AKcheckDep(type)
    else
        error('''type'' must be a string or cell string.')
    end

end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
function AKcheckDep(typeString)

    switch lower(typeString)
        case 'auditory'
            if ~exist('ERBFilterBank.m', 'file')
                error('AKtools:dependencies', 'This part of AKtools needs the Auditory Toolbox which is available from: https://engineering.purdue.edu/~malcolm/interval/1998-010/');
            end
        case 'fabian'
            if ~exist('FABIAN_HRIR_measured_HATO_0.sofa', 'file')
                error('AKtools:dependencies', 'This part of AKtools needs the FABIAN HRIR data set which is available from: https://dx.doi.org/10.14279/depositonce-5718.2 inside the Matlab search path.');
            end
        otherwise
            error('AKtools:dependencies', ['''' typeString ''' is not a valid input to AKdependencies'])
    end

end
