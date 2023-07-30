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
%   Normall, pre-computed cochlear model outputs can be applied to 
%   significantly reduce the required computation time. 
%
%   Start the AMT in the 'redo' mode to re-calculate the cochlear model
%   (be patient, this might take a few hours)
%
%   Figure 1: Output of the audiory model: the activity map.
%
%   See also: takanen2013
%
%   References:
%     O. Santala and V. Pulkki. Directional perception of distributed sound
%     sources. J. Acoust. Soc. Am., 129:1522 -- 1530, 2011.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Visualization of functional
%     count-comparison-based binaural auditory model output. Hearing
%     research, 309:147--163, 2014. PMID: 24513586.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Perceptually encoded signals and
%     their assessment. In J. Blauert, editor, The technology of binaural
%     listening. Springer, 2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/demos/demo_takanen2013.php


%   #Author: Marko Takanen (2013)
%   #Author: Olli Santala (2013)
%   #Author: Ville Pulkki (2013)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% Starting of the script
% Use pre-computed cochlear model outputs, otherwise set preComp=0;
flags=amt_configuration;

compType = 1;
printFigs = 0;
printMap = 1;

if flags.do_redo, %use binaural input signals and compute cochlear model
    filename='demo_binsig.mat';
    data=amt_load('takanen2013',filename);
    output= takanen2013(data.tests.insig,data.tests.fs,compType,printFigs,printMap);
    title(data.tests.scenario);
    set(gca,'Ytick',data.tests.ytickPos);set(gca,'YtickLabel',data.tests.ytickLab(end:-1:1));
    ylabel(data.tests.ylab);
else % use pre-computed cochlea model outputs to reduce the computation time
    filename='demo_cochlea.mat';
    data=amt_load('takanen2013',filename);
    output= takanen2013(data.tests.cochlea,data.tests.fs,compType,printFigs,printMap);
    title(data.tests.scenario);
    set(gca,'Ytick',data.tests.ytickPos);set(gca,'YtickLabel',data.tests.ytickLab(end:-1:1));
    ylabel(data.tests.ylab);
end

