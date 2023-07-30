function [ANresp,var_rate,psth] = zilany2014_synapse(vihc,fc,nrep,tdres,fiberType,noiseType,implnt)
%ZILANY2014_SYNAPSE Auditory nerve (AN) synapse model
%   Usage: [ANresp,var_rate, psth] = zilany2014_synapse(vihc,fc,nrep,tdres,fiberType,noiseType,implnt);
%
%   Input parameters:
%     vihc       : Output from inner hair cells (IHCs) in Volts
%     fc         : Center frequencies (Hz)
%     nrep       : Number of repetitions for the mean rate, rate variance 
%                  & psth calculation. Default is 1.
%     tdres      : simulation time resolution, fs_mod^(-1)
%     fiberType  : Type of the fiber based on spontaneous rate (SR) in spikes/s
%                  1: Low SR; 2: Medium SR (default); 3: High SR.
%     noiseType  : Fractional Gaussian noise will be different in every
%                  simulation (1), or will be always the same (0, default)
%     implnt     : 0...Use approxiate implementation of the power-law (default). 
%                  1...Use actual implementation of the power-law functions.
%
%   Output parameters:
%     ANresp     : AN response in terms of the estimated instantaneous mean 
%                  spiking rate (incl. refractoriness) in nf different AN 
%                  fibers spaced equally on the BM
%     var_rate   : var rate
%     psth       : Spike histogram
%
%   ZILANY2014_SYNAPSE returns modeled responses of one AN fibers to a specific inner haircell potential.
%
%   Please cite the references below if you use this model.
%
%
%   Demos: demo_zilany2014
%
%   References:
%     M. S. A. Zilany, I. C. Bruce, and L. H. Carney. Updated parameters and
%     expanded simulation options for a model of the auditory periphery. The
%     Journal of the Acoustical Society of America, 135(1):283--286, Jan.
%     2014.
%     
%     M. Zilany, I. Bruce, P. Nelson, and L. Carney. A phenomenological model
%     of the synapse between the inner hair cell and auditory nerve:
%     Long-term adaptation with power-law dynamics. J. Acoust. Soc. Am.,
%     126(5):2390 -- 2412, 2009.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/zilany2014_synapse.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB MEX M-Signal
%   #Author: Muhammad Zilany 
%   #Author: Robert Baumgartner: adapted to the AMT
%   #Author: Clara Hollomey (2020): adapted to AMT 1.0
%   #Author: Piotr Majdak (2021): C1 and C2 outputs

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


      [ANresp,var_rate,psth] = comp_zilany2014_synapse(vihc(:)',fc,nrep,tdres,fiberType,noiseType,implnt);
      ANresp=ANresp';
      var_rate=var_rate';
      psth=psth';
    
end


