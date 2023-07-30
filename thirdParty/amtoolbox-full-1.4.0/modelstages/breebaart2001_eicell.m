function z = breebaart2001_eicell(insig,fs,tau,ild,varargin)
%BREEBAART2001_EICELL  Excitation-inhibition cell computation for the Breebaart model
%   Usage: y = breebaart2001_eicell(insig,fs,tau,ild)
%
%   Input parameters:
%        insig	   : input signal, must be an [n by 2] matrix
%        fs        : sampling rate of input signal
%        tau       : characteristic delay in seconds (positive: left is leading)
%        ild       : characteristic ILD in dB (positive: left is louder)
%
%   Output parameters:
%        y	   : EI-type cell output as a function of time
%
%   BREEBAART2001_EICELL(insig,fs,tau,ild) compute the excitation-inhibition model on
%   the input signal insig.  The cell to be modelled responds to a delay
%   tau (measured in seconds) and interaural-level difference ild*
%   measured in dB.
%
%   BREEBAART2001_EICELL takes the following optional parameters:
%  
%     'tc',tc      Temporal smoothing constant. Default value is 30e-3.
%
%     'rc_a',rc_a  Parameter a for dynamic range compression.
%                  Default value is a=.1.
%
%     'rc_b',rc_b  Parameter b for dynamic range compression.
%                  Default value is b=0.00002.
%
%     'ptau',ptau  Time constant for p(tau) function. Default value is 2.2e-3.
%
%   See also: breebaart2001
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/breebaart2001_eicell.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: M-Signal
%   #Author: Jeroen Breebaart (2011)
%   #Author: Peter L. Soendergaard (2011)
%   #Author: Martina Kreuzbichler (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


if nargin<4
  error('%s: Too few input arguments.',upper(mfilename));
end;

definput.import={'breebaart2001_eicell'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% apply characteristic delay:
n = round( abs(tau) * fs );

l=insig(:,1);
r=insig(:,2);

if tau > 0,
    l = [zeros(n,1) ; l(1:end-n)];
else
    r = [zeros(n,1) ; r(1:end-n)];
end

% apply characteristic ILD:
l=gaindb(l, ild/2);
r=gaindb(r,-ild/2);

% compute instanteneous EI output:
x = (l - r).^2;

% temporal smoothing:
A=[1 -exp(-1/(fs*kv.tc))];
B=[1-exp(-1/(fs*kv.tc)) ];
y= filtfilt(B,A,x);% / ( (1-exp(-1/(fs*tc)))/2 );

% compressive I/O: Scale signal by 200. This approximately
% results in JNDs of 1 in the output
z = exp(-abs(tau)/kv.ptau) * kv.rc_a * log(kv.rc_b * y + 1);
% exp(-abs(tau)/0.0022) as in Larsen 2010
% 10^(-abs(tau)/0.005) as in Breebaart 2001a
% log10(kv.rc_b * y + 1) as in Davidson 2009


