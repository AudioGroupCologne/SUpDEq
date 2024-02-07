function [output,info] = bruce2018(stim,fsstim, fc, varargin) 
%BRUCE2018 Auditory-nerve filterbank (improved synapse)
%   Usage: [output] = bruce2018(stim,fsstim, fc);
%          [output,info] = bruce2018(stim, fsstim, fc, varargin);
%
%
%   Input parameters:
%     stim        : Pressure waveform of stimulus (timeseries)
%
%     fsstim      : Sampling frequency of stimulus
%
%     fc         : Frequency vector containing the CFs. 
%                  Use logspace(log10(flow), log10(fhigh),numCF) to
%                  replicate the results from Bruce et al. (2018).
%
%     varargin   : various flags and key-value pairs.
%
%
%   Output parameters:
%     output       : a struct containing the various modelstage outputs
%                    The struct output contains:
%
%                    - cohcs: actually used COHCs in the simulations
%                    - cihcs : actually used CIHCs in the simulations
%                    - sponts : actually used spontanous rates in the simulations
%                    - tabss : actually used absolute timings in the simulations
%                    - trels : actually used relative timings in the simulations
%                    - neurogram_ft : (fine-timing) neurogram calculated from the synapse output [time CFs]
%                    - t_ft : time axis of neurogram_ft (s)
%                    - neurogram_mr : average of spike count in each PSTH bin [time CF]
%                    - t_mr : time axis of neurogram_mr (s)
%                    - neurogram_Sout : synapse output [time CF]
%                    - t_Sout : time axis of neurogram_Sout (s)
%                    - fc : actually used CFs (Hz)
%                    - psth_ft : fine-timing PSTH [time CF]
%                    - meanrate : average spiking rate (spikes/s) [time CF]
%                    - varrate : variance of the spiking rate (spikes/s) [time CF]
%
%
%   BRUCE2018(...) returns modeled responses of multiple AN fibers tuned to 
%   various characteristic frequencies characterstic. 
%   
%   Please cite the references below if you use this model.
%
%   This function takes the following optional key/value pairs:
%
%     'ag_fs',ag_fs    Frequencies at which the audiogram should be
%                      evaluated. Used only in fitaudiogram.
%
%     'ag_db',ag_db    Hearing loss (dB) at frequencies ag_fs. 
%                      Used only in fitaudiogram.
%
%     'cohcs',cohcs    OHC scaling factors: 1 denotes normal OHC function (default); 
%                      0 denotes complete OHC dysfunction. Can be a vector
%                      of the size of fc or a scalar. Not used in fitaudiogram.
%
%     'cihcs',cihcs    IHC scaling factors: 1 denotes normal IHC function (default); 
%                      0 denotes complete IHC dysfunction. Can be a vector
%                      of the size of fc or a scalar. Not used in fitaudiogram.
%
%     'numL',nl        number of nerve fibres with low SR. Used only in autoSR or specificSRautoTiming.
%
%     'numM',nm        number of nerve fibres with medium SR. Used only in autoSR or specificSRautoTiming.
%
%     'numH',nh        number of nerve fibres with high SR. Used only in autoSR or specificSRautoTiming.
%
%     'lossL',lls      loss of low-SR fibres ranging from 0 (no fibres)  
%                      to 1 (healthy, all fibres). Used only in autoSR or specificSRautoTiming.
%
%     'lossM',lms      loss of medium-SR fibres ranging from 0 (no fibres)
%                      to 1 (healthy, all fibres). Used only in autoSR or specificSRautoTiming.
%
%     'lossH',lhs      loss of high-SR fibres ranging from 0 (no fibres)
%                      to 1 (healthy, all fibres). Used only in autoSR or specificSRautoTiming.
%
%     'SRL',SRL        SR of the low-SR fibres (spikes/s). Used only in specificSRautoTiming.
%
%     'SRM',SRM        SR of the medium-SR fibres (spikes/s). Used only in specificSRautoTiming.
%
%     'SRH',SRH        SR of the high-SR fibres (spikes/s). Used only in specificSRautoTiming.
%
%     'numsponts',n    Overall numbers of fibers. Used only in specificSR.
%
%     'spont',spont    SR (spikes/s). Can be scalar or size of fc. Used only in specificSR.
%
%     'tabs',tabs      Absolute timings (s). Can be scalar or size of fc. Used only in specificSR.
%
%     'trel',trel      Relative timings (s). Can be scalar or size of fc. Used only in specificSR.
%
%     'psthbinwidth_mr',psthbw   mean-rate binwidth (s).
%
%     'windur_ft',winft          fine-timing neurogram window length.
%
%     'windur_mr',winmr          mean-rate neurogram window length.
%
%     'nrep',nrep                Number of repetitions for the mean rate, 
%                                rate variance & psth calculation. Default is 1.
%
%     'reptime',rt               length of one repetition of the stimuli with pause.
%                                Default is 1.2  stimulus duration.
%
%     'fsmod',fsmod              Model sampling rate. It is possible to run the model 
%                                at a range of fsmod between 100 kHz and 500 kHz.
%                                Default value is 200 kHz for cats and 100 kHz for humans.
%
%   BRUCE2018 accepts the following flags:
%
%     'fitaudiogram'     Calculate the hearing-loss factors cihcs and 
%                        cohcs from the frequencies ag_fs and threshold 
%                        shifts ag_dbloss (dB) by using BRUCE2018_fitaudiogram.
%                        The default parameters reflect a healthy cochlea.
%
%     'no_fitaudiogram'  Default. Use directly provided cihcs and cohcs
%                        either as scalars or size of fc.
%                        The default parameters reflect a healthy cochlea. 
%
%     'human'          Default. Use model parameters for humans.
%
%     'cat'            Use model parameters for cats.
%
%     'fixedFGn'       Default. Fractional Gaussian noise will be the same in every 
%                      simulation.
%
%     'varFGn'         Fractional Gaussian noise will be different in every 
%                      simulation.
%
%     'approxPL'       Default. Use approxiate implementation of the power-law
%                      functions. 
%
%     'actualPL'       Use actual implementation of the power-law functions.
%
%     'outputPerSynapse' Output the synapse output of each individual 
%                        nerve fibre.
%                        This can considerably slow down the calculations.
%     
%     'outputPerCF'      Default. Output the average results over all synapses. 
%                        This mode is faster than outputPerSynapse.
%
%     'autoSR'            Generate the parameters for the AN 
%                         population wil be generated by the function 
%                         BRUCE2018_generateanpopulation based on 
%                         the number of low, medium, and high SR fibres
%                         numL, numM, numH. Default. 
%
%     'specificSRautoTiming'  Generate the timing parameters for the AN 
%                             population by the function 
%                             BRUCE2018_generateanpopulation based on 
%                             the number of low, medium, and high SR fibres
%                             numL, numM, numH. The spontanous rates are
%                             provided in SRL, SRM, and SRH. 
%
%     'specificSR'   Do not generate AN parameters. The overall number 
%                    of fibres numsponts, spontanouse rate spont, 
%                    absolute timing tabs, and relative timing trel
%                    must be provided. 
%
%
%   If 'outputPerSynapse' is specified, psth_ft, meanrate, and
%   varrate have the dimensions [time CF Syn] (with Syn as the number 
%   of synapses) and output additionally contains:
%
%     'synout'   the output of a synapse [time CF Syn]
%
%     'info'     a struct containing the parameter settings applied
%
%
%
%   See also: plot_bruce2018 demo_bruce2018_auditorynervemodel
%             demo_bruce2018 bruce2018_synapse bruce2018_generateanpopulation
%             bruce2018_innerhaircells bruce2018_fitaudiogram bruce2018_ffgn
%             exp_bruce2018 demo_carney2015 carney2015_fitaudiogram carney2015_generateneurogram
%             exp_osses2022 zilany2014
%
%   References:
%     I. C. Bruce, Y. Erfani, and M. S. R. Zilany. A phenomenological model
%     of the synapse between the inner hair cell and auditory nerve:
%     Implications of limited neurotransmitter release sites. Hearing
%     Research, 360:40--54, 2018.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/bruce2018.php


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

%
%   The other two modelstages 'bruce2018_innerhaircells' and
%   'bruce2018_synapse' are always active. '_innerhaircells' is called for
%   each element in 'fcs' (each characteristic frequency), and '_synapse'
%   is called for each nerve fiber. Per default and for execution speed, 
%   only results from outside of the nerve fiber loop are written to the 
%   struct 'output'. To retrieve results per fiber, set the flag
%   'ouputPerSynapse'. The actual parameters used in bruce2018 are
%   output to the 'info' struct.

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

% Define input flags and values
definput.import = {'bruce2018'}; % load defaults from arg_bruce2018

[flags,kv]  = ltfatarghelper({},definput,varargin);

% Derive the number of CF fibres
numCF = length(fc);

fs=kv.fsmod;  % sampling rate (Hz) of the model
tdres   = 1/fs; % sampling interval (s), i.e., the reciprocal of fs

stim	= resample(stim,fs,fsstim);	% resample the stimulus to the model fs
T  = length(stim)/fs;  % actual stimulus duration in seconds
smw_ft = hamming(kv.windur_ft);
smw_mr = hamming(kv.windur_mr);
simdur = ceil(T*kv.reptime/kv.psthbinwidth_mr)*kv.psthbinwidth_mr; % time interval (s) between stimulus repetitions

% species is either "cat" (1) or "human" (2): "1" for cat and "2" for human
if flags.do_cat, species = 1; else species = 2; end    

% noiseType is for fixed fGn (noise will be same in every simulation) or 
% variable fGn: "0" for fixed fGn and "1" for variable fGn
if flags.do_varFGn, noiseType = 0; else noiseType = 1; end

% implnt is for "approxiate" or "actual" implementation of the power-law functions: 
%"0" for approx. and "1" for actual implementation
if flags.do_actualPL, implnt = 1; else implnt = 0; end

if flags.do_fitaudiogram
  % mixed loss
  dbloss = interp1(kv.ag_fs,kv.ag_dbloss,fcs,'linear','extrap');
  [cohcs,cihcs,OHC_loss, IHC_loss]=bruce2018_fitaudiogram(fc,dbloss,species);
  info.OHC_loss = OHC_loss;
  info.IHC_loss = IHC_loss;
else
  if length(kv.cihcs) ~= numCF && length(kv.cihcs) ~= 1
    error('cihc needs to either be a scalar or a vector of the same size as fc.')
  end
  if length(kv.cohcs) ~= numCF && length(kv.cohcs) ~= 1
    error('cohc needs to either be a scalar or a vector of the same size as fc.')
  end

  cohcs = kv.cohcs.* ones(1, numCF);
  cihcs = kv.cihcs.* ones(1, numCF);
  info.OHC_loss = [];
  info.IHC_loss = [];
end

% generate AN parameters? 
if flags.do_autoSR,
    % generate SR and timing for a population based on numL, numM, and numH
  [sponts,tabss,trels] = bruce2018_generateanpopulation(numCF,[kv.numL kv.numM kv.numH]);
  numsponts = round([kv.lossL kv.lossM kv.lossH].*[kv.numL kv.numM kv.numH]);
  cntsponts = sum(numsponts);
end
if flags.do_specificSRautoTiming,
  % use specific SRs but calculate the timing based on numL, numM, and numH
  [~,tabss,trels] = bruce2018_generateanpopulation(numCF,[kv.numL kv.numM kv.numH]);
  sponts.LS = kv.SRL*ones(1,kv.numL); 
  sponts.MS = kv.SRM*ones(1,kv.numM); 
  sponts.HS = kv.SRH*ones(1,kv.numH);
  numsponts = round([kv.lossL kv.lossM kv.lossH].*[kv.numL kv.numM kv.numH]);
  cntsponts = sum(numsponts);
end
if flags.do_specificSR,
  % use predefined SR parameters for all fibers, used in e.g., exp_bruce2018('fig8b');
  if ~isscalar(kv.numsponts), error('numsponts needs to a scalar when running bruce2018 in specificSR mode.'); end
  cntsponts = kv.numsponts;
  sponts = kv.spont;
  tabss = kv.tabs;
  trels = kv.trel;
end

% initializations
init_size = simdur*fs;
output_len=floor(simdur/tdres+0.5)*kv.nrep; % length of vihc, c1, c2, synout
psth_len=output_len/kv.nrep; % length of vr, mr, psth_ft
psthbins = round(kv.psthbinwidth_mr*fs);  % number of psth_ft bins per psth bin

neurogram_ft = zeros(round(init_size),numCF);
neurogram_Sout = zeros(round(init_size)*kv.nrep,numCF);
neurogram_mr = zeros(round(init_size)/round(kv.psthbinwidth_mr*fs),numCF);
output.vihc=zeros(output_len,numCF);
output.C1=zeros(output_len,numCF);
output.C2=zeros(output_len,numCF);
if flags.do_outputPerSynapse
  output.psth_ft=zeros(psth_len,numCF,cntsponts);
  output.meanrate=zeros(psth_len,numCF,cntsponts);
  output.varrate=zeros(psth_len,numCF,cntsponts);
  output.psth=zeros(psth_len/psthbins,numCF,cntsponts);
end
if flags.do_outputPerCF
  % Memory allocation:
  output.psth_ft=zeros(psth_len,numCF);
  output.meanrate=zeros(psth_len,numCF);
  output.varrate=zeros(psth_len,numCF);
  output.psth=zeros(psth_len/psthbins,numCF);
  
  output.meanrate_LSR = zeros(size(output.meanrate));
  output.meanrate_MSR = zeros(size(output.meanrate));
  output.meanrate_HSR = zeros(size(output.meanrate));
  
  output.psth_LSR = zeros(size(output.psth));
  output.psth_MSR = zeros(size(output.psth));
  output.psth_HSR = zeros(size(output.psth));
end

for ii = 1:numCF
  %for each CF, collect the appropriate tabs, trels and sponts, and calculate
  %  the inner hair cell potential (for the whole stimulus)
    % FC = fc(ii);
    % cohc = cohcs(ii);
    % cihc = cihcs(ii);
    
    [vihc, C1, C2] = bruce2018_innerhaircells(stim,fc(ii),kv.nrep,tdres,simdur,cohcs(ii),cihcs(ii),species);  
    % [vihc, C1, C2] = zilany2014_innerhaircells(stim,FC,kv.nrep,tdres,simdur,cohc,cihc,species); equivalent to bruce2018
    output.vihc(:,ii) = vihc;
    output.C1(:,ii) = C1;
    output.C2(:,ii) = C2;

    if cntsponts>0

      if flags.do_autoSR || flags.do_specificSRautoTiming, % use generated SR parameters
          spont = [sponts.LS(ii,1:numsponts(1)) sponts.MS(ii,1:numsponts(2)) sponts.HS(ii,1:numsponts(3))];
          tabs = [tabss.LS(ii,1:numsponts(1)) tabss.MS(ii,1:numsponts(2)) tabss.HS(ii,1:numsponts(3))];
          trel = [trels.LS(ii,1:numsponts(1)) trels.MS(ii,1:numsponts(2)) trels.HS(ii,1:numsponts(3))];
      end   
      if flags.do_specificSR % use the same SR parameters for all fibers
          spont = sponts.*ones(1, cntsponts);
          tabs = tabss.*ones(1, cntsponts);
          trel = trels.*ones(1, cntsponts); 
      end    
      
      for jj = 1:cntsponts  % calculate the synapse output for each fibre type     

          amt_disp(['CF = ' int2str(ii) '/' int2str(numCF) '; spont = ' int2str(jj) '/' int2str(cntsponts) '; SR = ' num2str(spont(jj)) ' spikes/s'],'volatile');

          [psth_ft,mr,vr,synout] = bruce2018_synapse(vihc,fc(ii),kv.nrep,tdres,noiseType,implnt,spont(jj),tabs(jj),trel(jj));

          neurogram_Sout(:,ii) = neurogram_Sout(:,ii)+synout;
          cnt = sum(reshape(psth_ft,psthbins,length(psth_ft)/psthbins))';
          neurogram_ft(:,ii) = neurogram_ft(:,ii)+filter(smw_ft,1,psth_ft);
          neurogram_mr(:,ii) = neurogram_mr(:,ii)+filter(smw_mr,1,cnt);

          if flags.do_outputPerSynapse
              output.psth_ft(:, ii, jj) = psth_ft;
              output.meanrate(:, ii, jj) = mr;
              output.varrate(:, ii, jj) = vr;
              output.synout(:, ii, jj) = synout;
              output.psth(:,ii,jj)=cnt;
          end
          if flags.do_outputPerCF,
              output.psth_ft(:,ii)=output.psth_ft(:,ii)+psth_ft;
              output.meanrate(:,ii)=output.meanrate(:,ii)+mr;
              output.varrate(:,ii)=output.varrate(:,ii)+vr;
              output.psth(:,ii)=output.psth(:,ii)+cnt;
              
              fiberType = find(jj <= cumsum([kv.numL kv.numM kv.numH]),1,'first'); 
              
              switch fiberType
                  case 1
                      output.meanrate_LSR(:,ii) = output.meanrate_LSR(:,ii) + mr;
                      output.psth_LSR(:,ii) = output.psth_LSR(:,ii) + cnt;
                       
                  case 2
                      output.meanrate_MSR(:,ii) = output.meanrate_MSR(:,ii) + mr;
                      output.psth_MSR(:,ii) = output.psth_MSR(:,ii) + cnt;
                      
                  case 3
                      output.meanrate_HSR(:,ii) = output.meanrate_HSR(:,ii) + mr;
                      output.psth_HSR(:,ii) = output.psth_HSR(:,ii) + cnt;
              end
              
          end
      end
      amt_disp();
      if flags.do_outputPerCF
        output.psth_ft(:,ii)=output.psth_ft(:,ii)/cntsponts;
        output.meanrate(:,ii)=output.meanrate(:,ii)/cntsponts;
        output.varrate(:,ii)=output.varrate(:,ii)/cntsponts;
        
        if kv.numL ~= 0
            output.meanrate_LSR(:,ii)=output.meanrate_LSR(:,ii)/kv.numL;
        end
        if kv.numM ~= 0
            output.meanrate_MSR(:,ii)=output.meanrate_MSR(:,ii)/kv.numM;
        end
        if kv.numH ~= 0
            output.meanrate_HSR(:,ii)=output.meanrate_HSR(:,ii)/kv.numH;
        end
        
      end
    end
end
% amt_disp();
neurogram_ft = neurogram_ft(1:kv.windur_ft/2:end,:); % 50% overlap in Hamming window
t_ft = 0:kv.windur_ft/2/fs:(size(neurogram_ft,1)-1)*kv.windur_ft/2/fs; % time vector for the fine-timing neurogram
neurogram_mr = neurogram_mr(1:kv.windur_mr/2:end,:); % 50% overlap in Hamming window
t_mr = 0:kv.windur_mr/2*kv.psthbinwidth_mr:(size(neurogram_mr,1)-1)*kv.windur_mr/2*kv.psthbinwidth_mr; % time vector for the mean-rate neurogram
t_Sout = 0:1/fs:(size(neurogram_Sout,1)-1)/fs; % time vector for the synapse output neurogram    
           
%--------------------------------------------------------------------------
% write everything into an output struct
output.cohcs = cohcs;
output.cihcs = cihcs;
output.sponts = sponts;
output.tabss = tabss;
output.trels = trels;
output.t_ft = t_ft;
output.neurogram_ft = neurogram_ft;
output.neurogram_mr = neurogram_mr;
output.neurogram_Sout = neurogram_Sout;
output.t_mr = t_mr;
output.t_Sout = t_Sout;
output.fc = fc;
info.keyvals = kv;
info.flags = flags;


