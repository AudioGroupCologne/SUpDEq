% f = n2freq(n,N, fs)
% 
% n2freq gets an index n and returns frequency in Hz. N is the length of
% the underlying data in samples, fs the sampling frequency.
% fs is optional, default is 44100. n can be a scalar or vector


function f = n2freq(n,N, fs)

if ~exist('fs', 'var')
    fs = 44100;
end

f = (n-1)*fs/N;