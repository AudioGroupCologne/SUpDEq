% b = AKmoore1989(fc, g, q_bw, fs, N, M)
% calculated filters according to [1].
% For documentation see AKfilter.m, and AKfilterDemo.m
%
% [1] Brian C J Moore, Simon R Oldfield, Gary J Dooley: "Detection and
%     discrimination of spectral peaks and notches at 1 and 8 kHz," J.
%     Acoust. Soc. Am., vol. 85(2):820-836, (1989).
%
% 11/2015 - fabian.brinkmann@tu-berlin.de

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
function b = AKmoore1989(fc, g, q_bw, fs, N, M)

b = zeros(N, M);

% calculate filters
% NOTE: the filter order used is 16. In the original paper, a
% filter order of 8 was used, followed by a low-pass-to-band-pass
% transformation resulting in a filter order of 16. This is why the
% third noise in the original paper was highpassed with a 96 dB per
% oct. slope that equals order 96/6 = 16.
for m = 1:M
    
    % bandwidth
    f1 = fc(m)*q_bw(m)/2;
    
    % band pass
    d_bp = fdesign.bandpass('N,F3dB1,F3dB2', 16, (fc(m)-f1)/fs(m)*2, (fc(m)+f1)/fs(m)*2);
    H_bp = design(d_bp, 'butter');
    if isstable(H_bp)
        h_bp = impz(H_bp, N) * 10^(g(m)/20);
        H_bp = abs(fft(h_bp));
    else
        error('Filter not stabel. Change ''fc'' or ''bw''.')
    end
    
    % low pass
    d_lp = fdesign.lowpass('N,F3dB', 16, (fc(m)-f1)/fs(m)*2);
    H_lp = design(d_lp, 'butter');
    if isstable(H_lp)
        h_lp = impz(H_lp, N);
        H_lp = abs(fft(h_lp));
    else
        error('Filter not stabel. Change ''fc'' or ''bw''.')
    end
    
    % high pass
    d_hp = fdesign.highpass('N,F3dB', 16, (fc(m)+f1)/fs(m)*2);
    H_hp = design(d_hp, 'butter');
    if isstable(H_hp)
        h_hp = impz(H_hp, N);
        H_hp = abs(fft(h_hp));
    else
        error('Filter not stabel. Change ''fc'' or ''bw''.')
    end
    
    % power response of cascaded filters (-3 dB points add to 0 dB)
    B = sqrt(H_lp.^2 + H_bp.^2 + H_hp.^2);
    b(:,m) = ifft(B, 'symmetric');
    b(:,m) = AKphaseManipulation(b(:,m), fs(m), 'min', 2, 0);
end
