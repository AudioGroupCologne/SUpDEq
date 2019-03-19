%DEMO_TAKANEN2013 Demo of the binaural model by Takanen, Santala and Pulkki
%
%   This script generates a figure showing the result of the binaural
%   auditory model by Takanen, Santala and Pulkki (2013) for sound source
%   distributions consisting of different number of sound sources simulated
%   with HRTFs to emit incoherent samples of pink noise. The resulting
%   activity map shows that the activation spreads as the width of the
%   distribution increases, which is in accordance with the results of the
%   psychoacoustical experiment by Santala and Pulkki (2011).
%
%   Optionally, pre-computed cochlear model outputs can be applied to 
%   significantly reduce the required computation time. The pre-computed 
%   cochlear model outputs can be obtained from the authors.
%
%   Figure 1: Output of the audiory model
%
%      The activity map.
%
%   See also: takanen2013
%
%   References:
%     O. Santala and V. Pulkki. Directional perception of distributed sound
%     sources. J. Acoust. Soc. Am., 129:1522 - 1530, 2011.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Visualization of functional
%     count-comparison-based binaural auditory model output. Hearing
%     research, 309:147-163, 2014. PMID: 24513586.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Perceptually encoded signals and
%     their assessment. In J. Blauert, editor, The technology of binaural
%     listening. Springer, 2013.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/demos/demo_takanen2013.php

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

%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland

%% Starting of the script
% Use pre-computed cochlear model outputs, otherwise set preComp=0;
preComp = 1;

compType = 1;
printFigs = 0;
printMap = 1;

%if the user wishes to use pre-computed cochlea model outputs to reduce the
%required computation time
if preComp ==1
    filename='demo_cochlea.mat';
    data=amt_load('takanen2013',filename);
    output= takanen2013(data.tests.cochlea,data.tests.fs,compType,printFigs,printMap);
    title(data.tests.scenario);
    set(gca,'Ytick',data.tests.ytickPos);set(gca,'YtickLabel',data.tests.ytickLab(end:-1:1));
    ylabel(data.tests.ylab);
%otherwise, binaural input signals are used
else
    filename='demo_binsig.mat';
    data=amt_load('takanen2013',filename);
    output= takanen2013(data.tests.insig,data.tests.fs,compType,printFigs,printMap);
    title(data.tests.scenario);
    set(gca,'Ytick',data.tests.ytickPos);set(gca,'YtickLabel',data.tests.ytickLab(end:-1:1));
    ylabel(data.tests.ylab);
end
