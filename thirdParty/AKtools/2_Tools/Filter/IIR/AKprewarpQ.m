% For documentation see AKfilter.m, and AKfilterDemo.m
% For python version from Frank Schulz (SONIBLE) see
% http://nbviewer.jupyter.org/github/sonible/nb/blob/master/iir_filter/index.ipynb
%
% Frank Schultz, FG Audiokommunikation, TU Berlin
% frank.schultz@tu-berlin.de, +49 175 15 49 763, Skype: j0shiiv
% 0.00 07.09.2010 init dev
% Q prewarping for bandpass, bandstop and peq
%
% Hannes Helmholz, FG Audiokommunikation, TU Berlin
% helmholz@campus.tu-berlin.de
% 20.02.2017 extension of sources and sinus prewarping

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
function Qpre=AKprewarpQ(fm,fs,Q,type)
    if strcmp(type,'sin')
        % Robert Bristow-Johnson (1994): "The equivalence of various
        % methods of computing biquad coefficients for audio parametric
        % equalizers."  In: Proc. of 97th AES Convention, San Fransisco eq.
        % (14)
        w0 = 2*pi*fm/fs;
        BW = AKgetBWfromQ(Q);
        BW = BW*w0/sin(w0);
        Qpre = AKgetQfromBW(BW);
    elseif strcmp(type,'cos')
        % Thaden, R. (1997): "Entwicklung und Erprobung einer digitalen
        % parametrischen Filterbank." Diplomarbeit, RWTH Aaachen
        Qpre = Q*cos(pi*fm/fs);
    elseif strcmp(type,'tan')
        % Clark, R.J.; Ifeachor, E.C.; Rogers, G.M.; et al. (2000):
        % Techniques for Generating Digital Equalizer Coefficients.
        % In J. Aud. Eng. Soc. 48(4):281-298.
        Qpre = Q*(pi*fm/fs)/(tan(pi*fm/fs));
    else
        Qpre=Q; 
    end
end

function BW=AKgetBWfromQ(Q)
    BW = 2/log(2)*asinh(1/(2*Q));
end

function Q=AKgetQfromBW(BW)
    Q = 1/(2*sinh(log(2)/2*BW));
end
