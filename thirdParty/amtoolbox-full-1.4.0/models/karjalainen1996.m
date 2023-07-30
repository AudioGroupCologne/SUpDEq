function [slow,fast]=karjalainen1996(insig,fs)
%KARJALAINEN1996  Non-linear adapation network
%   Usage: [slow,fast]=karjalainen1996(insig,fs)
%
%   [slow,fast]=KARJALAINEN1996(insig,fs) computes the non-linear
%   adaptation network from Karjalainen et al. (1996) on the input signal
%   insig sampled with a sampling frequency of fs Hz.
% 
%    XXX What are the assumptions on the input? The example below
%     generates NaN's
%
%    XXX What is slow and fast*? In which units are they defined?
%
%    XXX Which level convention is used ?
%
%    XXX Bibtex entry for the correct paper to cite.
%
%   Examples:
%   ---------
%
%   The following show the adapation to a simple test signal:
%
%     [insig,fs] = greasy;
%     [slow,fast]=karjalainen1996(insig,fs);
%
%     subplot(1,2,1);
%     plot(slow);
%
%     subplot(1,2,2);
%     plot(fast);
%
%   This file (and the corresponding mex file) where originally part of
%   HUTear- Matlab toolbox for auditory modeling. The toolbox is available at 
%   http://www.acoustics.hut.fi/software/HUTear/>`
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/karjalainen1996.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Author: Aki Härmä, Helsinki University of Technology
%   #Author: Peter L. Sondergaard (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
 

[slow,fast]=comp_karjalainen1996(insig,fs);

