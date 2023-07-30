function out = breebaart2001_outmiddlefilter(in,fs)
%BREEBAART2001_OUTMIDDLEFILTER simulates the outer- and middle ear transfer function from Breebaart et. al. 2001
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
%     Soc. Am., 110:1074--1088, August 2001.
%     
%
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/breebaart2001_outmiddlefilter.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: M-Signal
%   #Author: Peter L. Soendergaard (2011)
%   #Author: Martina Kreuzbichler (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 



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





