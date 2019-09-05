% n = freq2n(f,N, fs)
% 
% freq2n gets a frequency f and returns the index n. N is the length of the
% underlying data in samples, fs the sampling frequency. fs is optional,
% default is 44100. f can be a scalar or vector


function n = freq2n(f,N, fs)

if ~exist('fs', 'var')
    fs = 44100;
end

n = round(f /  (fs/N)) + 1;