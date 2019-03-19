function [ANresp,fc,varargout] = zilany2014(spl,stim,fsstim,varargin)
%ZILANY2014 Auditory nerve (AN) model
%   Usage: [ANresp,fc] = zilany2014(lvl,stim,fsstim);
%          [ANresp,fc,vihc,psth] = zilany2014(lvl,stim,fsstim);
%
%   Input parameters:
%     spl         : Sound pressure level (SPL; re 20e-6 Pa) of stimulus in dB
%     stim        : Pressure waveform of stimulus (timeseries)
%     fsstim      : Sampling frequency of stimulus  
%
%   Output parameters:
%     ANresp     : AN response in terms of the estimated instantaneous mean 
%                  spiking rate (incl. refractoriness) in nf different AN 
%                  fibers spaced equally on the BM
%     fc         : Frequency vector containing the nf center frequencies
%     vihc       : Output from inner hair cells (IHCs) in Volts
%     psth       : Spike histogram
%
%   ZILANY2014(...) returns modeled responses of multiple AN fibers tuned to 
%   various characteristic frequencies characterstic frequencies evenly spaced 
%   along ERB scale (see audspace).
%
%   Please cite the references below if you use this model.
%
%   This function takes the following optional key/value pairs:
%
%     'flow',flow     Lowest centre frequency. Default value is 100.
%
%     'fhigh',fhigh   Highest centre frequency. Default value is 16000.
%
%     'nfibers',nf    Number of fibers between lowest and highest
%                     frequency. The fibers will be equidistantly spaced
%                     on the basilar membrane. Default value is 500.
%
%     'fsmod',fsmod   Model sampling rate. It is possible to run the model 
%                     at a range of fsmod between 100 kHz and 500 kHz.
%                     Default value is 200kHz for cats and 100kHz for humans.
%
%     'fiberType',fT  Type of the fiber based on spontaneous rate (SR) in spikes/s
%                     fT=1 for Low SR; fT=2 for Medium SR (default); 
%                     fT=3 for High SR.
%
%     'cohc',cohc     OHC scaling factor: 1 denotes normal OHC function (default); 
%                     0 denotes complete OHC dysfunction.
%
%     'cihc',cihc     IHC scaling factor: 1 denotes normal IHC function (default); 
%                     0 denotes complete IHC dysfunction.
%
%     'nrep',nrep     Number of repetitions for the mean rate, 
%                     rate variance & psth calculation. Default is 1.
%
%   ZILANY2014 accepts the following flag:
%
%     'human'         Use model parameters for humans. This is the default.
%
%     'cat'           Use model parameters for cats.
%
%     'fixedFGn'      Fractional Gaussian noise will be the same in every 
%                     simulation. This is the default.
%
%     'varFGn'        Fractional Gaussian noise will be different in every 
%                     simulation.
%
%     'approxPL'      Use approxiate implementation of the power-law
%                     functions. This is the default.
%
%     'actualPL'      Use actual implementation of the power-law functions.
% 
%
%   Demos: demo_zilany2014
%
%   References:
%     M. S. A. Zilany, I. C. Bruce, and L. H. Carney. Updated parameters and
%     expanded simulation options for a model of the auditory periphery. The
%     Journal of the Acoustical Society of America, 135(1):283-286, Jan.
%     2014.
%     
%     M. Zilany, I. Bruce, P. Nelson, and L. Carney. A phenomenological model
%     of the synapse between the inner hair cell and auditory nerve:
%     Long-term adaptation with power-law dynamics. J. Acoust. Soc. Am.,
%     126(5):2390 - 2412, 2009.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/models/zilany2014.php

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

% AUTHOR: code provided by Muhammad Zilany, AMT compatibility adapted by Robert Baumgartner

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

% Define input flags and values
definput.flags.species = {'human','cat'};
definput.flags.noiseType = {'fixedFGn','varFGn'};
definput.flags.powerLawImp = {'approxPL','actualPL'};
definput.keyvals.flow    = 100;
definput.keyvals.fhigh   = 16000;
definput.keyvals.nfibers = 500;
definput.keyvals.fsmod   = 100e3;
definput.keyvals.fiberType = 2;
definput.keyvals.cohc    = 1;
definput.keyvals.cihc    = 1;
definput.keyvals.nrep    = 1;
[flags,kv]  = ltfatarghelper({'flow','fhigh','nfibers','fsmod',...
  'fiberType','cohc','cihc','nrep'},definput,varargin);

% Species settings
if flags.do_cat
  kv.fsmod = max(kv.fsmod,200e3);
end

% Stimulus settings
stim	= resample(stim,kv.fsmod,fsstim);	% stim fs = mod fs
idnz	= stim ~= 0;                    % ignore pauses
lvlref = 20*log10(1/20e-6);           % Reference level: 20 micro Pa
stim(idnz) = setdbspl(stim(idnz),spl,'dboffset',lvlref); % Calibrate level
stim = stim(:)'; % stim must be a row vector

% characteristic frequencies evenly spaced along ERB scale
fc = audspace(kv.flow,kv.fhigh,kv.nfibers,'erb');

% tdres is the binsize in seconds, i.e., the reciprocal of the sampling rate
tdres   = 1/kv.fsmod;   

% reptime is the time between stimulus repetitions in seconds; -> set 
% twice the duration of stim
reptime = (2*length(stim))/kv.fsmod;

% species is either "cat" (1) or "human" (2): "1" for cat and "2" for human
species = flags.do_human+1;

% noiseType is for fixed fGn (noise will be same in every simulation) or variable fGn: "0" for fixed fGn and "1" for variable fGn
noiseType = flags.do_varFGn + 0;

% implnt is for "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
implnt = flags.do_actualPL + 0;

% Call AN model and loop for all fibers tuned to different CFs
vihc = zeros(kv.nfibers,round(reptime*kv.fsmod));
ANresp = vihc;
psth = vihc;
for jj = 1:kv.nfibers
  
  % Call IHC model (mex'ed C model)
  vihc(jj,:) = comp_zilany2014IHC(stim,fc(jj),kv.nrep,tdres,reptime,kv.cohc,kv.cihc,species);
  
  % Call Synapse model
  [ANresp(jj,:),varrate,psth(jj,:)] = comp_zilany2014Synapse(...
    vihc(jj,:),fc(jj),kv.nrep,tdres,kv.fiberType,noiseType,implnt);
end



if nargout >= 3
  varargout{1} = vihc;
  if nargout == 4
    varargout{2} = psth;
  end
end

% AUTHOR notes:

% vihc = model_IHC(pin,CF,nrep,tdres,reptime,kv.cohc,kv.cihc,flags.do_human+1);
% [meanrate,varrate,psth] = model_Synapse(vihc,CF,nrep,tdres,kv.fiberType,flags.do_varFGn,flags.do_actualPL);
%
% Output parameters:
% vihc is the inner hair cell (IHC) potential (in volts)
% meanrate is the estimated instantaneous mean rate (incl. refractoriness)
% varrate is the estimated instantaneous variance in the discharge rate (incl. refractoriness)
% psth is the peri-stimulus time histogram 
%
% Input parameters:
% pin is the input sound wave in Pa sampled at the appropriate sampling rate (see instructions below)
% CF is the characteristic frequency of the fiber in Hz
% nrep is the number of repetitions for the mean rate, rate variance & psth calculation
% tdres is the binsize in seconds, i.e., the reciprocal of the sampling rate (see instructions below)
% reptime is the time between stimulus repetitions in seconds - NOTE should be equal to or longer than the duration of pin
% cohc is the OHC scaling factor: 1 is normal OHC function; 0 is complete OHC dysfunction
% cihc is the IHC scaling factor: 1 is normal IHC function; 0 is complete IHC dysfunction
% species is either "cat" (1) or "human" (2): "1" for cat and "2" for human
% fiberType is the type of the fiber based on spontaneous rate (SR) in spikes/s - "1" for Low SR; "2" for Medium SR; "3" for High SR
% noiseType is for fixed fGn (noise will be same in every simulation) or variable fGn: "0" for fixed fGn and "1" for variable fGn
% implnt is for "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
%
% Example:
%
%   This example shows how to model a normal fiber of high spontaneous rate 
%   (normal OHC & IHC function) with a CF of 1 kHz, for 10 repititions 
%   and a sampling rate of 100 kHz, for a repetition duration of 200 ms, 
%   with a fixed fractional Gaussian noise (fGn) and also with approximate 
%   implementation of the power-law functions in the synapse model. :::
%
%    vihc = model_IHC(pin,1e3,10,1/100e3,0.200,1,1,1); **requires 8 input arguments
%    [meanrate,varrate,psth] = model_Synapse(vihc,1e3,10,1/100e3,3,0,0); **requires 7 input arguments
%
