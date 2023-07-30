function [output,cf] = verhulst2018(insig,fs,fc_flag,varargin)
%VERHULST2018 Cochlear transmission-line model (improved, incl. brainstem)
%
%   Usage: 
%     output = verhulst2018(insig,fs,fc_flag)
%     [output, cf] = verhulst2018(insig,fs,fc_flag, ...)
%
%   Input parameters:
%
%     insig          : the input signal to be processed. Each column is 
%                      processed in parallel, so it is possible to run several simulations in parallel
%     fs             : sampling rate (Hz)
%     fc_flag        : list of frequencies specifying the probe positions 
%                      along the basilar membrane, or all to probe all 
%                      1000 cochlear sections, or abr to probe 401 locations
%                      between 112 and 12000 Hz.
%     numH           : number of high-SR fibers. Must be larger than zero.
%     numM           : number of medium-SR fibers. Can be zero.
%     numL           : number of low-SR fibers. Can be zero.
%     ic             : calculate the IC responses
%
%   Output parameters:
%
%     cf        : Center frequencies (Hz) of the probed basiliar membrane sections.
%     output    : Structure with the following fields:
%     fs_an     : Sample rate (Hz) of the output.
%     fs_abr    : Sample rate (Hz) of the brainstem sections (IC, CN, W1, W3, and W5).
%     w1        : Wave 1, output of the AN model.
%     w3        : Wave 3, output of the CN model.
%     w5        : Wave 5, output of the IC model.
%     an_summed : Sum of HSR, MSR and LSR responses (per channel) and the input to the CN (modelled by verhulst2015_cn). Provided by default. Can be disabled by the flag no_an.
%     ihc       : IHC receptor potential. Provided by default. Can be disabled by the flag no_ihc.
%     cn        : Detailed output of the CN. Provided by default. Can be disabled by the flag no_cn.
%     ic        : Detailed output of the IC. Provided by default. Can be disabled by the flag no_ic.
%
%   This function computes the basilar membrane displacement.
%
%
%   The output can optionally provide the following information:
%
%     'anfH'   responses of the high-SR fibers. Optional, only when called with flag anfH.
%
%     'anfM'   responses of the medium-SR fibers. Optional, only when called with flag anfM.
%
%     'anfL'   responses of the low-SR fibers. Optional, only when called with flag anfL.
%
%     'v'      velocity of the basilar membrane sections [time section channel] and input to the AN (modelled by VERHULST2018_ihctransduction. Optional, provided only when called with flag v.
%
%     'y'     displacement of the basilar membrane sections [time section channel]. Can be disabled by the flag no_y.
%
%     'oae'   ottoacoustic emission as sound pressure at the middle ear. Can be disabled by the flag oae.
%
%   References:
%     S. Verhulst, A. Altoè, and V. Vasilkov. Functional modeling of the
%     human auditory brainstem response to broadband stimulation.
%     hearingresearch, 360:55--75, 2018.
%     
%
%   License
%   --------
%
%   This model is licensed under the UGent Academic License. For non-commercial academic research, 
%   you can use this file and/or modify it under the terms of that license. Further usage details 
%   are provided in the in the AMT directory "licences".
%
%   See also: verhulst2015 verhulst2018 demo_verhulst2018 demo_verhulst2012
%             verhulst2018_ihctransduction verhulst2015_cn
%             verhulst2015_ic verhulst2018_auditorynerve exp_verhulst2012
%             verhulst2012 verhulst2015
%             verhulst2018 middleearfilter data_takanen2013 takanen2013_periphery
%             exp_osses2022 exp_takanen2013 takanen2013
%
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/verhulst2018.php


%   #License: ugent
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal PYTHON C
%   #Author: Alejandro Osses (2020): primary implementation based on https://github.com/HearingTechnology/Verhulstetal2018Model
%   #Author: Piotr Majdak (2021): adaptations for the AMT 1.0

% This file is licenced under the terms of the UGent Academic License, which details can be found in the AMT directory "licences" and at <https://raw.githubusercontent.com/HearingTechnology/Verhulstetal2018Model/master/license.txt>.
% For non-commercial academic research, you can use this file and/or modify it under the terms of that license. This file is distributed without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. 

amt_info('once');
% ------ Checking of input parameters ------------

if nargin<3
	error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig)
	error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
    error('%s: fs must be a positive scalar.',upper(mfilename));
end;

definput.import={'verhulst2018'};  % load defaults from arg_verhulst2018
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

IrrPct = keyvals.IrrPct;
storeflag = ' ';
if flags.do_oae, storeflag = [storeflag 'e']; end
if flags.do_y  , storeflag = [storeflag 'y']; end
if flags.do_ihc, storeflag = [storeflag 'i']; end

[Nr_signals,idx]=min(size(insig));
irregularities=keyvals.irr_on.*ones(1,Nr_signals);

%%% Fix model parameters so far: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject = keyvals.subject; % seed
non_linear_type = keyvals.non_linear_type;

% Number of neurones: 13-3-3 => no synaptopathy, any loss of neurones => synaptopathy
numH = keyvals.numH;
numM = keyvals.numM;
numL = keyvals.numL;

if idx == 1
    error('If multiple signals are used as input make sure each column is one signal')
end

%%% Stage 1A. Outer ear: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_outerear
    % Same outer ear filter as in other peripheral models
    hp_fir = headphonefilter(fs);% Getting the filter coefficients at fs
    N = ceil(length(hp_fir)/2);  % group delay for a FIR filter of order length(hp_fir)
    M = size(insig,2);
    insig = [insig; zeros(N,M)]; % group delay compensation: step 1 of 2.
    insig = filter(hp_fir,1,insig); % filtering
    insig = insig(N+1:end,1:M); % group delay compensation: step 2 of 2
end

%%% Stage 1B. Middle ear: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[b,a] = middleearfilter(fs,'verhulst2018');
if flags.do_middleear
    insig = filter(b,a,insig);
end

if flags.do_no_middleear
    % If no middle-ear filter is applied, a gain/attenuation equivalent to
    % the level in the filter band-pass is applied to the input signal. This
    % is important to get an appropriate compression in the subsequent filter
    % bank.
    K = 8192; % arbitrary
    h = ones(K,1);
    h = h.*freqz(b,a,K);

    me_gain_TF = max( 20*log10(abs(h)) ); % max of the filter response
    insig = gaindb(insig,me_gain_TF);
end

%%% Stage 2. Cochlear filter bank - Transmission line model: %%%%%%%%%%%%%%
%         2.1 Signal preprocessing:
fs_in = fs;
if fs ~= keyvals.fs_up
    insig = resample(insig,keyvals.fs_up,fs);
    fs = keyvals.fs_up;
end
stim=transpose(insig); % transpose it (python C-style row major order)

if(ischar(fc_flag)) %if probing all sections 1001 output (1000 sections plus the middle ear)
    switch  fc_flag
        case 'all'
            Nr_sections = 1000;
        case 'half'
            Nr_sections = 500;
        case 'abr'
            Nr_sections = 401;
            % Nothing to do, it is correctly formatted
        otherwise
            error('fc flag not recognised it should be either ''all'',''abr'',''half'', or a numeric array');
    end
    fc_str = fc_flag;
else
    % then fc is a numeric array that will be passed as a column vector
    fc_flag=round(fc_flag(:));
    fc_str = 'custom CF(s)';
    Nr_sections = length(fc_flag);
end
probes=fc_flag;

  % Load poles for the selected hearing profile
poles=amt_load('verhulst2018','Poles.mat','Poles');
if ~isfield(poles.Poles,keyvals.hearing_profile)
  error(['Hearing profile ' keyvals.hearing_profile ' not available.']);
end
sheraPo=poles.Poles.(keyvals.hearing_profile);

dir_model = fullfile(amt_basepath,'environments','verhulst2018',filesep);
dir_data  = fullfile(dir_model,'out',filesep);
version_year = 2018; 

% if ~exist(fullfile(dir_model,'tridiag.so'),'file')
% 	error('/environments/verhulst2018/tridiag.so library is missing. Run amt_mex');
% end

channels = Nr_signals;
L_samples = length(stim(1,:));

amt_disp('VERHULST 2018: The following parameters will be passed to Python:',flags.disp);
amt_disp(['  Number of signals to be processed: ',num2str(channels),' channels'],flags.disp);
amt_disp(['  Cochlear hearing profile: ',keyvals.hearing_profile,' (poles between ',num2str(sheraPo(1)),' and ',num2str(sheraPo(end)),')'],flags.disp);
amt_disp(['  Number of auditory nerve fibres: (',num2str([numH,numM,numL]),') (HSR,MSR,LSR)'],flags.disp);
amt_disp(['  Number of cochlear sections to be stored: flag ',num2str(fc_str),' (all=1000, half=500, abr=401)'],flags.disp);
amt_disp(['  Seed number: ',num2str(subject),' (subject ''variable'')'],flags.disp);
amt_disp(['  Irregularities=',num2str(irregularities(1)),' (1=on,0=off), IrrPrct=',num2str(IrrPct)],flags.disp);
amt_disp(['  Non linear type=',num2str(non_linear_type),'(''vel'' or ''disp'')'],flags.disp);
amt_disp(['  Extra variables to be saved: ',num2str(storeflag),' (storeflag)'],flags.disp);

%%% Storing input, change dir, run the model, and come back
amt_disp('Cochlear processing...','volatile');

in.stim=stim; in.fs=fs; in.channels=channels; in.subject=subject;
in.sheraPo=sheraPo; in.irregularities=irregularities;
in.probes=probes; in.storeflag=storeflag; in.IrrPct=IrrPct;

in.non_linear_type=non_linear_type; in.version_year=version_year;
out.cf=[Nr_sections,1,Nr_signals];
out.v=[Nr_sections,L_samples,Nr_signals];
if flags.do_y,   out.y=[Nr_sections,L_samples,Nr_signals]; end
if flags.do_oae, out.e=[L_samples,1,Nr_signals]; end
out=amt_extern('Python','verhulst2018','run_cochlear_model.py',in,out); 

v = out.v;

if flags.do_no_v, out = rmfield(out,'v'); end % reduce memory usage if v not used.

fs_abr = keyvals.subfs; % default subfs = 20000 Hz
fs_an  = keyvals.subfs;
DECIMATION = fs/keyvals.subfs;

output(1:Nr_signals) = struct('fs_bm',fs_in);

for ii=1:Nr_signals
  
    amt_disp(['Neural processing ' num2str(ii) ' of ' num2str(Nr_signals)],'volatile');

    v_here = squeeze(v(:,:,ii))';

    output(ii).cf=squeeze(out.cf(:,1,ii));
    if flags.do_v,   output(ii).v=resample(v_here,fs_in,fs); end
    if flags.do_y  , output(ii).y=resample(squeeze(out.y(:,:,ii))',fs_in,fs); end
    if flags.do_oae, output(ii).oae=resample(squeeze(out.e(:,:,ii)),fs_in,fs); end
    
    %%% Stage 3. Inner hair cell model: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flags.do_ihc || flags.do_an || flags.do_cn || flags.do_ic
        outsig = verhulst2018_ihctransduction(v_here,fs,version_year,keyvals); % unscaled input
        output(ii).ihc = resample(outsig,fs_in,fs);
        output(ii).fs_ihc=fs_in;
    end

    if flags.do_an || flags.do_cn || flags.do_ic
        %%% Stage 4: Auditory nerve model and Wave I: %%%%%%%%%%%%%%%%%%%%%%%%%
        % Reducing the sampling frequency: Decimation.
      n = 30;
      for jj = 1:size(outsig,2)
        Vm_res(:,jj) = decimate(outsig(:,jj),DECIMATION,n,'fir'); % n-th order FIR filter before decimation
      end
      Vm_res(1:DECIMATION,:) = repmat(outsig(1,:),DECIMATION,1);
      
      if numH ~= 0 % HSR neurones
          anfH = verhulst2018_auditorynerve(Vm_res,fs_abr,keyvals.kSR_H,keyvals.kmax_H,version_year);
      else
          error('There should be at least one HSR neurone, set numH to a non-null value...')
      end
      
      if numM ~= 0 % MSR neurones
          anfM = verhulst2018_auditorynerve(Vm_res,fs_abr,keyvals.kSR_M,keyvals.kmax_M,version_year);
      else
          anfM = zeros(size(Vm_res));
      end
        
      if numL ~= 0 % LSR neurones
          anfL = verhulst2018_auditorynerve(Vm_res,fs_abr,keyvals.kSR_L,keyvals.kmax_L,version_year);
      else
          anfL = zeros(size(Vm_res)); % empty array if no anfL, saves some computation power
      end

      if flags.do_anfH, output(ii).anfH = anfH; end
      if flags.do_anfM, output(ii).anfM = anfM; end
      if flags.do_anfL, output(ii).anfL = anfL; end

      outsig = numL*anfL+numM*anfM+numH*anfH;
      if flags.do_an, output(ii).an_summed = outsig; end

      % Loading scaling constants for Wave I, III, and V:
      M1 = keyvals.M1;
      M3 = keyvals.M3;
      M5 = keyvals.M5;

      switch Nr_sections
          case {401,500}                
              cal_factor = 1; % same cochlear tuning in both configurations
          case 1000                
              cal_factor = 0.5; % more dense cochlear resolution (twice as many channels)
          otherwise
              cal_factor = 1;
      end

      output(ii).w1 = cal_factor*M1*sum(outsig,2);
    end
    
    if flags.do_cn || flags.do_ic
      %%% Stage 5: Cochlear nucleus: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Tex = keyvals.Tex_cn;
      Tin = keyvals.Tin_cn;
      dly = keyvals.dly_cn;
      A   = keyvals.Acn;
      S   = keyvals.Scn;
      outsig = verhulst2015_cn(outsig,fs_abr,Tex,Tin,dly,A,S);
      if flags.do_cn
        if ~iscell(outsig)
          output(ii).cn = outsig;
        else
          output(ii).cn_mfb = outsig;
          output(ii).cn = il_sum_cell(outsig);
        end
      end

      if flags.do_no_mfb
        % Only one CN filter: outsig is numeric
        output(ii).w3 = cal_factor*M3*sum(outsig,2);
      end
      if flags.do_mfb
        % CN from the modulation filter bank: outsig is a cell variable
        output(ii).w3 = cal_factor*M3*sum(il_sum_cell(outsig),2);
      end
    end

    if flags.do_ic
      %%% Stage 6: Inferior Colliculus: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Tex = keyvals.Tex_ic;
      Tin = keyvals.Tin_ic;
      dly = keyvals.dly_ic;
      A   = keyvals.Aic;
      S   = keyvals.Sic;
      outsig = verhulst2015_ic(outsig,fs_abr,Tex,Tin,dly,A,S);
      if flags.do_ic
        if ~iscell(outsig)
          output(ii).ic = outsig;
        else
          output(ii).ic_mfb = outsig;
          output(ii).ic = il_sum_cell(outsig);
        end

      end
      if flags.do_no_mfb
        % Only one IC filter: outsig is numeric
        output(ii).w5 = cal_factor*M5*sum(outsig,2);
      end
      if flags.do_mfb
        % IC from the modulation filter bank: outsig is a cell variable
        output(ii).w5 = cal_factor*M5*sum(il_sum_cell(outsig),2);
      end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    output(ii).fs_an  = fs_an;
    output(ii).fs_abr = fs_abr;
end

output(1).keyvals = keyvals;
amt_disp();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inoutsig = il_sum_cell(inoutsig)

for i = 2:length(inoutsig)

    inoutsig{1} = inoutsig{1}+inoutsig{i};

end
inoutsig = inoutsig{1}; % keeps only the first cell, and converts it into a double array


