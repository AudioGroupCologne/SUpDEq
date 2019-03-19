function out = breebaart2001_outmiddlefilter(in,fs)
%breebaart2001_outmiddlefilter simulates the outer- and middle ear transfer function from Breebaart et. al. 2001
%   
%   Usage: out = breebaart2001_outmiddlefilter(in,fs)
%
%   Input parameters:
%        in  : input acoustic signal.
%        fs     : sampling rate.
%
%   BREEBAART2001_OUTMIDDLEFILTER(in,fs) filters the input signal in with the
%   transfer function of the outer and middle ear  sampled with a frequency
%   of fs Hz as described in Breebaart (2001).
%
%   References:
%     J. Breebaart, S. van de Par, and A. Kohlrausch. Binaural processing
%     model based on contralateral inhibition. I. Model structure. J. Acoust.
%     Soc. Am., 110:1074-1088, August 2001.
%     
%
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/breebaart2001_outmiddlefilter.php

% Copyright (C) 2009-2015 Piotr Majdak and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 0.9.9
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   AUTHOR: Martina Kreuzbichler


% coefficients
q = 2 - cos(2*pi*4000/fs) - sqrt((cos(2*pi*4000/fs)-2).^2-1);
r = 2 - cos(2*pi*1000/fs) - sqrt((cos(2*pi*1000/fs)-2).^2-1);

% initialize vectors
len=length(in)+2; 
% get the length + overhead
y=zeros(len,1);    
% define y and x
x=zeros(len,1);     
x(3:len)=in;     

% main loop
for n=3:len
    y(n)=(1-q)*r*x(n) - (1-q)*r*x(n-1) + (q+r)*y(n-1) - q*r*y(n-2);
end

% set the output
out=y(3:end);

% to see bode
% imp=[1; zeros(1024,1)];
% out= outmiddleartransfunct(imp,48000);
% spect(fft(out),48000);




