function [vihc, varargout] = bruce2018_innerhaircells(insig, fc, nrep, dt, duration, cohc, cihc, species)
%BRUCE2018_INNERHAIRCELLS inner hair cell potential
%
%   Usage:
%     vihc = bruce2018_innerhaircells(insig, fc, nrep, dt, duration, cohc, cihc, species)
%
%   Input parameters:
%     insig    : input stimulus [time x 1]
%     fc       : vector of center frequencies [Hz]
%     dt       : discrete time distance, 1/sampling frequency, needs to be either 100, 200, or 500 kHz
%     duration : stimulus pause duration [s]
%     nrep     : number of stimulus repetitions (about 10 - 200)
%     cohc     : outer hair cell coefficient (1.0...NH)
%     cihc     : inner hair cell coefficient (1.0...NH)
%     species  : 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
%
%   Output parameters:
%     vihc     : the inner hair cell (IHC) relative transmembrane potential [V] [time x 1]
%     C1       : chirp filter C1 output  [time x 1]
%     C2       : Wideband filter C2 output  [time x 1]
%
%   BRUCE2018_INNERHAIRCELLS calculates the inner hair cells' relative transmembrane potential [V]
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/bruce2018_innerhaircells.php


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


[vihc, C1, C2] = comp_bruce2018_IHC(insig(:)',fc,nrep,dt,duration,cohc,cihc,species);

vihc=vihc'; % AMT 1.0: time is first dimension
if nargout >=1
    varargout{1} = C1';
    varargout{2} = C2';
end

end


% Author notes:
% model_IHC_BEZ2018 - Bruce, Erfani & Zilany (2018) Auditory Nerve Model
%
%     vihc = model_IHC_BEZ2018(pin,CF,nrep,dt,reptime,cohc,cihc,species);
%
% vihc is the inner hair cell (IHC) relative transmembrane potential (in volts)
%
% pin is the input sound wave in Pa sampled at the appropriate sampling rate (see instructions below)
% CF is the characteristic frequency of the fiber in Hz
% nrep is the number of repetitions for the psth
% dt is the binsize in seconds, i.e., the reciprocal of the sampling rate (see instructions below)
% reptime is the time between stimulus repetitions in seconds - NOTE should be equal to or longer than the duration of pin
% cohc is the OHC scaling factor: 1 is normal OHC function; 0 is complete OHC dysfunction
% cihc is the IHC scaling factor: 1 is normal IHC function; 0 is complete IHC dysfunction
% species is the model species: "1" for cat, "2" for human with BM tuning from Shera et al. (PNAS 2002),
%    or "3" for human BM tuning from Glasberg & Moore (Hear. Res. 1990)
%
% For example,
%
%    vihc = model_IHC_BEZ2018(pin,1e3,10,1/100e3,0.2,1,1,2); **requires 8 input arguments
%
% models a normal human fiber of high spontaneous rate (normal OHC & IHC function) with a CF of 1 kHz, 
% for 10 repetitions and a sampling rate of 100 kHz, for a repetition duration of 200 ms, and
% with approximate implementation of the power-law functions in the synapse model.
%
%
% NOTE ON SAMPLING RATE:-
% Since version 2 of the code, it is possible to run the model at a range
% of sampling rates between 100 kHz and 500 kHz.
% It is recommended to run the model at 100 kHz for CFs up to 20 kHz, and
% at 200 kHz for CFs> 20 kHz to 40 kHz.


