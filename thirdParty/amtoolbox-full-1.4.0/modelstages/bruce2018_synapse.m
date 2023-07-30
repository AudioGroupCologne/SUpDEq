function [psth, meanrate, varrate, synout, trd_vector, trel_vector] = bruce2018_synapse(vihc, fc, nrep, dt, noisetype, implnt, spont, tabs, trel)
%BRUCE2018_SYNAPSE synapse model proposed by Bruce et al. (2018)
%
%   Usage:
%     [psth, meanrate, varrate, synout, trd_vector, trel_vector] = bruce2018_synapse(vihc, fc, nrep, dt, noisetype, implnt, spont, tabs, trel)
%
%   Input parameters:
%     vihc      : the inner hair cell (IHC) relative transmembrane potential [V]
%     fc        : vector of center frequencies [Hz]
%     nrep      : number of stimulus repetitions (about 10 - 200)
%     dt        : discrete time distance, 1/sampling frequency, needs to be either 100, 200, or 500 kHz
%     noiseType : "variable" or "fixed (frozen)" fGn: 1 for variable fGn and 0 for fixed (frozen) fGn
%     implnt    : "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
%     spont     : nerve fibres
%     tabs      : absolute timing info
%     trel      : relative timing info
%
%   Output parameters:
%     psth        : the peri-stimulus time histogram (PSTH) (or a spike train if nrep = 1)
%     meanrate    : mean spiking rate
%     varrate     : variable spiking rate
%     synout      : synapse output
%     trd_vector  : refractory period
%     trel_vector : relative refractory period
%
%   BRUCE2018_SYNAPSE calculates the output of the synapse model by Bruce et al. (2018) and several associated parameters
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/bruce2018_synapse.php


%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: MATLAB MEX M-Signal
%   #Author: Ian Bruce: basic code of the model
%   #Author: Alejandro Osses (2020): original implementation
%   #Author: Clara Hollomey (2021): adapted to the AMT 1.0
%   #Author: Piotr Majdak (2021): adaptations to exp_osses2022; specificSRautoTiming added

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 



 [psth, meanrate, varrate, synout, trd_vector, trel_vector] = ...
        comp_bruce2018_Synapse(vihc(:)',fc,nrep,dt,noisetype,implnt,spont,tabs,trel);
      
  meanrate=meanrate';
  varrate=varrate';
  psth=psth';
  synout=synout';
end


