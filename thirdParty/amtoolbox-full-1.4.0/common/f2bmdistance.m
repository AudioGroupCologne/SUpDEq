function [f, varargout] = f2bmdistance( segnum, varargin)
%F2BMDISTANCE calculates the resonance frequency of each cochlear segment
%    
%   Usage: f = f2bmdistance( segnum );
%          f = f2bmdistance( flow, fhigh, segnum );
%
%   F2BMDISTANCE uses the Greenwood position-frequency map function (1990)
%   for calculating the characteristic frequencies (CFs) and the cochlear 
%   positions (Pos) corresponding to the distribution of segnum channels. The
%   the low-frequency (flow) and the high-frequency (fhigh) range can be
%   optionally specified.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/common/f2bmdistance.php


%   #Author: The AMT Team (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if nargin == 1
    %assume we have the segment index given
    n=0:1/segnum:1;
    f=zeros(1,length(n)-1);
    seg=zeros(1,segnum);
    for i=1:segnum
        seg(i)=(n(i)+n(i+1))/2;
    end

    %% Greenwoood parameters
    A = 165.44; % for human auditory span it is A=160.1377
    alfa = 2.1;
    k = 1;
    %% Calculate the Greenwood function
    fe=zeros(size(seg));
    fe=A.*(10.^(alfa.*seg)-k);
    f=fliplr(fe);
    %% Display if needed
    % figure(1)
    % stem(seg,f,':r');
elseif nargin == 3
    %assume we have flow, fhigh, and number of channels given
    Flow = segnum;
    Fhigh = varargin{1};
    N = varargin{2};
    %% Greenwoood parameters
    A = 165.44; % for human auditory span it is A=160.1377
    alfa = 2.1;
    k = 1;
    Plow=(log10((Flow/A)+k)/alfa);
    Phigh=(log10((Fhigh/A)+k)/alfa);
    buff=zeros(1,N+1);Pos=zeros(1,N);
    for i=1:N
        Pos(i)=((Plow-Phigh)/(N-1))*(i-1)+Phigh;
    end

    %%
    CF=zeros(size(Pos));
    CF=A.*(10.^(alfa.*Pos)-k);
    Pos=1-Pos;
    %%
    f = CF;
    varargout{1} = Pos;
else
    error('f2bmdistance either requires the segment index, or flow, fhigh, and channelnumber as input arguments.');
    
end


