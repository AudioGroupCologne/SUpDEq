%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function cr = supdeq_convFFT(sig,filtr)
%
% Function to perform (2-channel) convolution in frequency domain.
%
% Output:
% cr                     - Convolution result
%
% Input:        
% sig                   - Audio signal to be filtered
% filtr                 - Single- or double-channel filter which should be
%                         applied to signal s
%
% Dependencies: -
%   
% (C) 2020 by JMA, Johannes M. Arend  
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function cr = supdeq_convFFT(sig,filtr)

if size(sig,2) > size(sig,1)
    sig = sig.';
    warning('Flipped dimensions of signal');
end

if size(filtr,2) > size(filtr,1)
    filtr = filtr.';
    warning('Flipped dimensions of filter');
end

L = length(sig)+length(filtr)-1;
NFFT = 2^nextpow2(L);

%For two channel filter
if size(filtr,2) > 1
    sig = [sig,sig];
end

cr = ifft(fft(sig,NFFT) .* fft(filtr,NFFT));
cr = cr(1:L,:);
