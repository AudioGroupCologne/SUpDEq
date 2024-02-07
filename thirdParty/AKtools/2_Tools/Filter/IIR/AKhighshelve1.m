% For documentation see AKfilter.m, and AKfilterDemo.m
% For python version from Frank Schulz (SONIBLE) see
% http://nbviewer.jupyter.org/github/sonible/nb/blob/master/iir_filter/index.ipynb
%
% Frank Schultz, FG Audiokommunikation, TU Berlin
% frank.schultz@tu-berlin.de, +49 175 15 49 763, Skype: j0shiiv
% 0.00 06.09.2010 init dev
% see AKcalcBiquadCoeff.m for details
%
% Hannes Helmholz, FG Audiokommunikation, TU Berlin
% helmholz@campus.tu-berlin.de
% 21.02.2017 integration of former AKhighshelve1ABC.m and switch to
% standard types I, II and III

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
function [b,a]=AKhighshelve1(fg,fs,G,type)

    if abs(G) < 1e-05 % flat        
        b = [1, 0, 0];
        a = [1, 0, 0];
        return;
    end
    
    w = 2*fs*tan(pi*fg/fs); % frequency pre-warping
    g = 10^(G/20);
    
    if strcmp(type,'I')
        alpha = 1;
    elseif strcmp(type,'II')
        alpha = g^.5;
    elseif strcmp(type,'III')
        alpha = g^.25;
    % ELSE is prevented in AKfilter.m
    end
    
    if G > 0
        B = [0, g*alpha^-2/w, 1];
        A = [0, alpha^-2/w,   1];
    else
        B = [0, alpha^2/w,      1];
        A = [0, g^-1*alpha^2/w, 1];
    end
    
    [b,a] = AKcalcBiquadCoeff(B,A,fs);
    
end
