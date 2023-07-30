function outsig = viemeister1979(insig,fs)
%VIEMEISTER1979  The Viemeister (1979) leaky-integrator model
%   Usage: outsig=viemeister1979(insig,fs); 
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/viemeister1979.php


%
%   VIEMEISTER79(insig,fs) is a simpel amplitude modulation detection model
%   from the paper Viemeister 79. The input insig is a matrix containing
%   the intervals as column vectors. The first column is assumed to
%   always contain the target. The output is 1 if the model found that
%   the most amplitude modulated interval was the target, and 0 otherwise.
%
%   This model is included mostly as a test, as it is so simple.
%
%   References: viemeister1979tmt

%   #Author: C. Hollomey (2020): added four-pole Butterworth Bandpass
%   #StatusDoc: Submitted
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: M-Signal

%   AMT 1.0, PM (24.4.2021): Viemeister (1979) is about TMTFs, but the implementation is just a a leaky integrator. Is this the full model?

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.
  

narginchk(2, 2);
  
% 4-6 kHz four-pole Butterworth Bandpass (tentatively added)
[b, a] = butter (4, [4000/fs, 6000/fs]);
insig = filter(b, a, insig);
  
% halfwave rectification
insig = max(insig,0);

% first-order lowpass filter @ 65 Hz
[lp_b,lp_a] = butter(1,65/(fs/2));
insig = filter(lp_b, lp_a, insig);

% ac-coupled rms = std
%<<<<<<< HEAD:models/viemeister1979.m
%outsig = std(insig,1);

%=======
stddev = std(insig,1);

% Choose the interval with the highest standard deviation
[dummy,interval] = max(stddev);

% Answer is correct if we choose the first intervals
answer=(interval==1);


