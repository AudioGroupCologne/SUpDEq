function output = joergensen2011_combineinformation(input,SNRs,conditions,Nsp)
%joergensen2011_combineinformation  Combine information
%   Usage: output = joergensen2011_combineinformation(input,SNRs,conditions,Nsp);
%
%   Input parameters:
%     input      : Cell array with the SNRenv results for each
%                  processing condition (n),SNR (k), and speech sample (q)
%     SNRs       : Vector with the SNRs used
%     conditions : Vector with the procssing conditions used.
%     Nsp        : Number the number of speech samples in the XXX
%
%   Output parameters:
%     output     : Structure containing the following fields:
%
%                  output.combined_aud 
%                     [n,k] Matrix with overall SNRenv
%                     for each processing condition and SNR
% 
%                  output.SNRenvs
%                     Cell array {n,k,q} with entries for each
%                     condition, SNR, and speech sample. Each entry
%                     is Matrix with an SNRenv value for each XXX
%
%   JOERGENSEN2011_COMBINEINFORMATION(input,SNRs,conditions,Nsp) combines
%   the SNRenv across modulation and audio filters. It is also possible to
%   extracts other information such as the excitation patterns or long-term
%   spectra.
%
%   XXX Better description
%
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/modelstages/joergensen2011_combineinformation.php

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

% Author: Søren Jørgensen august 2010

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



