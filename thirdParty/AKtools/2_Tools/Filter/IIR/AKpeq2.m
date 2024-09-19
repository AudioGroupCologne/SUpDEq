% For documentation see AKfilter.m, and AKfilterDemo.m
% For python version from Frank Schulz (SONIBLE) see
% http://nbviewer.jupyter.org/github/sonible/nb/blob/master/iir_filter/index.ipynb
%
% Frank Schultz, FG Audiokommunikation, TU Berlin
% frank.schultz@tu-berlin.de, +49 175 15 49 763, Skype: j0shiiv
% 0.00 07.09.2010 init dev
% see AKcalcBiquadCoeff.m for details
% do the prewarping for Q outside the function!
%
% Hannes Helmholz, FG Audiokommunikation, TU Berlin
% helmholz@campus.tu-berlin.de
% 21.02.2017 integration of former AKpeq2ABC.m and switch to standard types
% I, II and III

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
function [b,a]=AKpeq2(fm,fs,G,Q,type)

    if abs(G) < 1e-05 % flat        
        b = [1, 0, 0];
        a = [1, 0, 0];
        return;
    end
    
    % Q pre-warping should be done already
    w = 2*fs*tan(pi*fm/fs); % frequency pre-warping
    g = 10^(G/20);
    
    if strcmp(type,'I')
        gamma = g;
        delta = g;
    elseif strcmp(type,'II')
        gamma = 1;
        delta = g;
    elseif strcmp(type,'III')
        gamma = g^.5;
        delta = g^.5;
    % ELSE is prevented in AKfilter.m
    end

    if G > 0
        B = [1/w^2, delta/(Q*w),     1];
        A = [1/w^2, (delta/g)/(Q*w), 1];
    else
        B = [1/w^2, gamma/(Q*w),     1];
        A = [1/w^2, (gamma/g)/(Q*w), 1];
    end
    
    [b,a] = AKcalcBiquadCoeff(B,A,fs);
    
end
