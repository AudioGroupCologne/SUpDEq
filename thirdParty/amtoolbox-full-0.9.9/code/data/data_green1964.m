function [Pc_est SNRenv] = data_green1964(SNRenv_lin,material)
%data_green1964  Converts the overall SNRenv to percent correct
%   Usage: [Pc_est SNRenv ] = data_green1964(SNRenv_lin,material);
%
%   Input parameters:
%     SNRenv_lin :  matrix with the SNRenv values (not in dB),
%                   one for each SNR and processing combination.
% 
%     material   :  Specify the speech material used.
%                   Current: 'CLUE', 'DANTALE' or 'DANTALE II'
%
%   XXX Description is missing
%  
%   References:
%     D. M. Green and T. G. Birdsall. The effect of vocabulary size. In J. A.
%     Swets, editor, Signal Detection and Recognition by Human Observers,
%     pages 609-619. John Wiley & Sons, 1964.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_green1964.php

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

%  AUTHOR: Søren Jørgensen august 2010


if nargin < 2
    amt_disp('you have to specify the spech material used: CLUE, DANTALE or DANTALEII')
end

% ---------- Determine the model parameters based on the speech material used
switch material
    case 'CLUE'
        m = 8000; 
        sigma_s = .6;
        
    case 'DANTALE'
        m = 8000 ; 
        sigma_s = .9;
     
    case 'DANTALEII'
        m = 50; 
        sigma_s = .9; 
    
end
% -------- The general conversion constants. They are the same for alle materials
    k = .82; 
    q = .5;

%--------  Only used for output:
   SNRenv = 10*log10(SNRenv_lin);

% ---------- Converting from SNRenv to d_prime  --------------
    d_prime = k*(SNRenv_lin).^q; %SNRm_lin; %

%----------- Converting from d_prime to Percent correct, Green and Birdsall (1964)----------
    Un = 1*norminv(1-(1/m)); 
    mn = Un + (.577 /Un);% F^(-1)[1/n] Basically gives the value that would be drawn from a normal destribution with probability p = 1/n.
    sig_n=  1.28255/Un; 
    Pc_est = normcdf(d_prime,mn,sqrt(sigma_s.^2+sig_n.^2))*100;  

        

