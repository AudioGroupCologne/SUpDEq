function [r_mean,psth,ihc,c1,c2,r_var,output] = zilany2014(stim,fsstim,fc,varargin)
%ZILANY2014 Auditory-nerve filterbank (improved)
%
%   Usage: [ANresp,fc] = zilany2014(stim,fsstim, fc);
%          [ANresp,fc,psth,ihc,t_var] = zilany2014(stim,fsstim, fc);
%
%   Input parameters:
%     stim        : Pressure waveform of stimulus (timeseries)
%     fsstim      : Sampling frequency of stimulus  
%     fc         : Frequency vector containing the CFs. 
%                  Use fc=audspace(lo,hi,numCF,'erb'); to space equally on the 
%                  ERB frequency scale. 
%
%   Output parameters:
%     r_mean     : Instantaneous mean spiking rate (incl. refractoriness) 
%                  of different AN fibers at corresponding CFs. Size: [time CFs]
%     psth       : Spike histogram
%     ihc        : Output from inner hair cells (IHCs) in Volts
%     c1         : Output from the chirping filter C1
%     c2         : Output from the chirping filter C2
%     r_var      : Instananeous variance in the discharge rate of the ANs
%
%
%   This function takes the following optional key/value pairs:
%
%     'fsmod',fsmod   Model sampling rate. It is possible to run the model 
%                     at a range of fsmod between 100 kHz and 500 kHz.
%                     Default value is 200kHz for cats and 100kHz for humans.
%
%     'fiberType',fT  Type of the fiber based on spontaneous rate (SR)
%                     1: Low SR, SR fixed to 0.1 spikes/s
%                     2: Medium SR, SR fixed to 4 spikes/s
%                     3: High SR, SR fixed to 100 spikes/s
%                     4: Custom, defined by the fibre numbers in numH, numM and numL
%
%     'numH'          Number of high SR fibres. Only if fiberType is 4. 
%     'numM'          Number of medium SR fibres. Only if fiberType is 4. 
%     'numL'          Number of low SR fibres. Only if fiberType is 4.
%
%     'cohc',cohc     OHC scaling factor: 1 denotes normal OHC function (default); 
%                     0 denotes complete OHC dysfunction.
%
%     'cihc',cihc     IHC scaling factor: 1 denotes normal IHC function (default); 
%                     0 denotes complete IHC dysfunction.
%
%     'nrep',nrep     Number of repetitions for the mean rate, 
%                     rate variance & psth calculation. 
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
%     'shera2002'     Selects the BM tuning from Shera et al. (2002) (default).
% 
%     'glasberg1990'  Selects the BM tuning from Glasberg & Moore (1990)
%
%   ZILANY2014(...) returns modeled responses of multiple AN fibers tuned to 
%   various characteristic frequencies (CFs). Middle-ear filtering is 
%   included and corresponds to middleearfilter(...,'zilany2009');
%   Please cite the references below if you use this model.
%
%   See also: demo_zilany2014 zilany2014_synapse zilany2014_innerhaircells
%             zilany2014_ffgn zilany2014 zilany2007
%             data_baumgartner2016 plot_roenne2012 plot_roenne2012_chirp
%             plot_roenne2012_tonebursts demo_carney2015 baumgartner2016_spectralanalysis
%             roenne2012_click roenne2012_chirp carney2015_generateneurogram
%             roenne2012_tonebursts exp_osses2022
%             bruce2018 roenne2012 baumgartner2013
%             baumgartner2016
%
%   Demos: demo_zilany2014
%
%   References:
%     C. A. Shera, J. J. J. Guinan, and O. A. J. Revised estimates of human
%     cochlear tuning from otoacoustic and behavioral measurements.
%     Proceedings of the National Academy of Sciences of the United States of
%     America, 99(5):3318--3323, 2002.
%     
%     B. R. Glasberg and B. Moore. Derivation of auditory filter shapes from
%     notched-noise data. Hearing Research, 47(1-2):103--138, 1990.
%     
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
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/zilany2014.php


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

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading default values, reading input arguments:

if nargout >= 7
    output = [];
end

% Define input flags and values
definput.import = {'zilany2014'}; % load defaults from arg_zilany2014

[flags,kv]  = ltfatarghelper({'fsmod','fiberType','cohc','cihc','nrep'},definput,varargin);

% Species settings
if flags.do_cat, kv.fsmod = max(kv.fsmod,200e3); end

% Stimulus settings
stim	= resample(stim,kv.fsmod,fsstim); % stim fs = mod fs
fs = kv.fsmod;
tdres   = 1/fs; % time interval of processing (in s), i.e., reciprocal of the sampling rate

stim = stim(:)'; % stim will be a row vector

numCF=length(fc);

% reptime is the time between stimulus repetitions in seconds; -> set 
% twice the duration of stim
reptime = kv.reptime*length(stim)/fs;

% species must be 1 for cat, 2 for human (shera2002) or 3 for human (glasberg1990)
if flags.do_cat, species=1; end
if flags.do_human
  if flags.do_shera2002, species=2; end
  if flags.do_glasberg1990, species=3; end
end

% noiseType is for fixed fGn (noise will be same in every simulation) or variable fGn: "0" for fixed fGn and "1" for variable fGn
noiseType = flags.do_varFGn + 0;

% implnt is for "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
implnt = flags.do_actualPL + 0;

% Call AN model and loop for all fibers tuned to different CFs
ihc = zeros(round(reptime*fs)*kv.nrep,numCF);
r_mean = zeros(round(reptime*fs),numCF);
r_var = zeros(round(reptime*fs),numCF);
psth = zeros(round(reptime*fs),numCF);

if nargout >= 7
    % Memory allocation:
    meanrate_LSR = zeros(size(r_mean));
    meanrate_MSR = zeros(size(r_mean));
    meanrate_HSR = zeros(size(r_mean));
    
    psth_LSR = zeros(size(psth));
    psth_MSR = zeros(size(psth));
    psth_HSR = zeros(size(psth));
end

if nargout>=5,
  c1 = zeros(round(reptime*fs*kv.nrep),numCF);
  c2 = zeros(round(reptime*fs*kv.nrep),numCF);
end

if kv.fiberType==4,
  fiberTypes=[ones(kv.numL,1);ones(kv.numM,1)*2;ones(kv.numH,1)*3]; % custom fiber types
else
  fiberTypes=kv.fiberType;
end
for jj = 1:numCF
  if nargout>=5,
    [ihc(:,jj),c1(:,jj),c2(:,jj)] = zilany2014_innerhaircells(stim,fc(jj),kv.nrep,tdres,reptime,kv.cohc,kv.cihc,species);
  else
    ihc(:,jj) = zilany2014_innerhaircells(stim,fc(jj),kv.nrep,tdres,reptime,kv.cohc,kv.cihc,species);    
  end

  for ii=1:length(fiberTypes)
	  amt_disp(['CF = ' int2str(jj) '/' int2str(numCF) '; spont = ' int2str(ii) '/' int2str(length(fiberTypes))],'volatile');
    [r_mean1,r_var1,psth1] = zilany2014_synapse(ihc(:,jj)',fc(jj),kv.nrep,tdres,fiberTypes(ii),noiseType,implnt);
    r_mean(:,jj)=r_mean(:,jj)+r_mean1;
    r_var(:,jj)=r_var(:,jj)+r_var1;
    psth(:,jj)=psth(:,jj)+psth1;
    
    if nargout >= 7
        switch fiberTypes(ii)
            case 1
                meanrate_LSR(:,jj) = meanrate_LSR(:,jj)+r_mean1;
                psth_LSR(:,jj) = psth_LSR(:,jj)+psth1;
                
            case 2
                meanrate_MSR(:,jj) = meanrate_MSR(:,jj)+r_mean1;
                psth_MSR(:,jj) = psth_MSR(:,jj)+psth1;
                
            case 3
                meanrate_HSR(:,jj) = meanrate_HSR(:,jj)+r_mean1;
                psth_HSR(:,jj) = psth_HSR(:,jj)+psth1;
        end
    end
    
  end
  if ~isempty(fiberTypes), amt_disp(); end
  if length(fiberTypes)>1
    
    r_mean(:,jj)=r_mean(:,jj)/length(fiberTypes);
    r_var=r_var/length(fiberTypes);
    psth(:,jj)=psth(:,jj)/length(fiberTypes);
    
    if nargout >= 7
        if kv.numL ~= 0
            meanrate_LSR = meanrate_LSR / kv.numL;
            psth_LSR(:,jj)=psth_LSR(:,jj)/kv.numL;
        end
        output.meanrate_LSR = meanrate_LSR;
        if kv.numM ~= 0
            meanrate_MSR = meanrate_MSR / kv.numM;
            psth_MSR(:,jj)=psth_MSR(:,jj)/kv.numM;
        end
        output.meanrate_MSR = meanrate_MSR;
        if kv.numH ~= 0
            meanrate_HSR = meanrate_HSR / kv.numH;
            psth_HSR(:,jj)=psth_HSR(:,jj)/kv.numH;
        end
        output.meanrate_HSR = meanrate_HSR;
    end
  end
  % [sum(psth) sum(psth_LSR) sum(psth_MSR) sum(psth_HSR) (4*sum(psth_LSR)+4*sum(psth_MSR)+12*sum(psth_HSR))/20]
end

if ~isempty(kv.psth_binwidth),
  
    bin = round(kv.psth_binwidth*fs);
    psth_bin = zeros(size(buffer(psth(:,1),bin),2),numCF);
    for ii=1:numCF
        psth_bin(:,ii)=sum(buffer(psth(:,ii),bin))/(kv.nrep*kv.psth_binwidth);
    end
    psth=psth_bin;
  
    if nargout >= 7
        %%%
        psth_bin = zeros(size(buffer(psth_LSR(:,1),bin),2),numCF);
        for ii=1:numCF
          psth_bin(:,ii)=sum(buffer(psth_LSR(:,ii),bin))/(kv.nrep*kv.psth_binwidth);
        end
        psth_LSR=psth_bin;
        output.psth_LSR = psth_LSR;

        %%%
        psth_bin = zeros(size(buffer(psth_MSR(:,1),bin),2),numCF);
        for ii=1:numCF
          psth_bin(:,ii)=sum(buffer(psth_MSR(:,ii),bin))/(kv.nrep*kv.psth_binwidth);
        end
        psth_MSR=psth_bin;
        output.psth_MSR = psth_MSR;

        %%%
        psth_bin = zeros(size(buffer(psth_HSR(:,1),bin),2),numCF);
        for ii=1:numCF
          psth_bin(:,ii)=sum(buffer(psth_HSR(:,ii),bin))/(kv.nrep*kv.psth_binwidth);
        end
        psth_HSR=psth_bin;
        output.psth_HSR = psth_HSR;
    end
  
end


