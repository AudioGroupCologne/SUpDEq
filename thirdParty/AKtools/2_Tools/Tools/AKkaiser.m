% kaiserWin = AKkaiser(L, A, frac)
% calculates kaiser window [1] of length L and side lobe attenuation A.
% Kaiser window is shifted by frac samples (-1 < frac < 1).
%
% For example usage see AKfractionalDelay.m
%
% [1] Oppenheim, A, R W Schafer, J R Buck (2004): Zeitdiskrete
% Signalverarbeitung. 2. ueberarbeitete Auflage. Pearson Education:Muenchen
% et al.
%
% 04/2012 - fabian.brinkmann@tu-berlin.de

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
function kaiserWin = AKkaiser(L, A, frac)

% calculate beta after [1], eq. 7.62
A = abs(A);

if A > 50
    beta = .1102*(A-8.7);
elseif A >= 21
    beta = .5842*(A-21)^.4 + .07886*(A-21);
else
    beta = 0;
end

% calculate Kaiser Window for diskrete points after [1], eq. 7.59
if mod(L,2)
    alpha = floor(L/2);
else
    alpha = L/2-.5;
end

I_0_b = besseli(0,beta);
L = (0+frac:L-1+frac)';
Z = beta*(1-((L-alpha)/alpha).^2).^.5;

% supress very small imaginary part
kaiserWin = real(besseli(0,Z))./I_0_b;
