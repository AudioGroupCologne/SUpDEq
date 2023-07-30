function output = joergensen2011_combineinformation(input,SNRs,conditions,Nsp)
%JOERGENSEN2011_COMBINEINFORMATION  Combine information
%   Usage: output = joergensen2011_combineinformation(input,SNRs,conditions,Nsp);
%
%   Input parameters:
%     input      : Cell array with the SNRenv results for each
%                  processing condition (n),SNR (k), and speech sample (q)
%     SNRs       : Vector with the SNRs used
%     conditions : Vector with the procssing conditions used.
%     Nsp        : Number of speech samples
%
%   Output parameters:
%     output     : Structure containing the SNRenv
%
%   JOERGENSEN2011_COMBINEINFORMATION(input,SNRs,conditions,Nsp) combines
%   the SNRenv across modulation and audio filters. It is also possible to
%   extracts other information such as the excitation patterns or long-term
%   spectra.
%
%   The output struct contains the following fields:
%
%     .combined_aud      [n,k] Matrix with overall SNRenv for each processing condition and SNR
%
%     .SNRenvs           Cell array {n,k,q} with entries for each condition, SNR, and speech sample. Each entry
%                        is Matrix with an SNRenv value for each corpus
%
%
%   See also: joergensen2011
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/joergensen2011_combineinformation.php


%   #StatusDoc: Submitted
%   #StatusCode: Submitted
%   #Verification: Untrusted
%   #Requirements: M-Signal M-Stats
%   #Author: Søren Jørgensen (2010)
%   #Author: Peter L. Sondergaard (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if nargin <6
    TimeSpectra = 'no';
end

if nargin <5
    ExPatterns = 'no';
end
res = input;

for q = 1:Nsp %for each of the speech samples
    for n = conditions %for each of the processing conditions
        for k = 1:length(SNRs) %for each SNR
            
            Output_tmp = res{n,k,q};
            % extracting the SNRenvs from the internal representation:
            %           
            SNRenvs{n,k,q} = Output_tmp.outSNRenvs(:,:);
% ------------------ Combining information --------------------------------
            
            tmp = SNRenvs{n,k,q};
            
            %                Converting to linear values:
            linear = 10.^(tmp*.1);
            linear(find(tmp == 0)) = 0;
            
            combined_mod_tmp = sqrt(sum(linear.^2,1)); % Combining across modulation filters using integration model (Green and Swets 1988)
            combined_aud(n,k,q) = (sqrt(sum(combined_mod_tmp.^2)));% Combining across auditory filters
            clear combined_mod_tmp;
            
        end
    end
    
end
%
% saving output:
output.combined_aud = combined_aud;
output.SNRenvs = SNRenvs;




