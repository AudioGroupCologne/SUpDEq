function [Pc_est, SNRenv] = data_green1964(SNRenv_lin,material)
%DATA_GREEN1964  Converts the overall SNRenv to percent correct
%   Usage: [Pc_est SNRenv ] = data_green1964(SNRenv_lin,material);
%
%   Input parameters:
%     SNRenv_lin :  matrix with the SNRenv values (not in dB),
%                   one for each SNR and processing combination.
% 
%     material   :  Specify the speech material used.
%                   Current: 'CLUE', 'DANTALE' or 'DANTALE II'
%
%
%   for 3 speech corpi, this function converts their SNRenv to dprime values
%   and to percent correct
%  
%   References:
%     D. M. Green and T. G. Birdsall. The effect of vocabulary size. In J. A.
%     Swets, editor, Signal Detection and Recognition by Human Observers,
%     pages 609--619. John Wiley & Sons, 1964.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_green1964.php


%   #Author: Søren Jørgensen (2010)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


if nargin < 2
    amt_disp('Specify the spech material: CLUE, DANTALE or DANTALEII');
end

% ---------- Determine the model parameters based on the speech material used
switch material
    case 'CLUE'XXX Description is missing
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

        


