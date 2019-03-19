function data = data_baumgartner2016(varargin)
%DATA_BAUMGARTNER2016  Data from Baumgartner et al. (2016)
%   Usage: data = data_baumgartner2016(flag)
%
%   Input parameters:
%
%   Output parameters:
%     data : data structure contains the following fields:
%
%            .id         listener ID
%
%            .S          listener-specific sensitivity parameter
%
%            .mrs        listener-specific task-induced response scatter (derived
%                        from central lateral response precision in baseline condition)
%
%            .Obj        DTF data in SOFA Format
%
%            .pe_exp     experimental local polar RMS error
%
%            .qe_exp     experimental quadrant error rate
%
%            .target     experimental target angles
%
%            .response   experimental response angles
%
%            .itemlist   experimental item list. Columns denote:
%                        1:4 ... azi_target,ele_target,azi_response,ele_response
%                        5:8 ... lat_target,pol_target,lat_response,pol_response
%                        9   ... F/B-Confusion resolved pol_response
%
%          If the 'model'-falg is set the output contains also the following fields
%
%            .S          listener-specific sensitivity parameter.
%
%            .Obj        DTF data in SOFA Format.
%
%            .pe_exp     experimental local polar RMS error in baseline condition.
%
%            .qe_exp     experimental quadrant error rate in baseline condition.
%
%            .target     experimental target angles.
%
%            .response   experimental response angles.
%
%            .stim       target stimulus.
%
%            .fsstim     sampling rate of target stimulus.
%
%   DATA_BAUMGARTNER2016(flag) returns data from Baumgartner et al. (2016)
%   describing a model for sound localization in sagittal planes (SPs)
%   on the basis of listener-specific directional transfer functions (DTFs).
%
%   DATA_BAUMGARTNER2016 accepts the following flags:
%
%     'baumgartner2014' data of the pool from Baumgartner et al. (2014). This is the default.
%     'Long'            300ms at 50+-5dB SL.
%     '10dB'            3ms at 10+-5dB SL.
%     '20dB'            3ms at 20+-5dB SL.
%     '30dB'            3ms at 30+-5dB SL.
%     '40dB'            3ms at 40+-5dB SL.
%     '50dB'            3ms at 50+-5dB SL.
%     '60dB'            3ms at 60+-5dB SL.
%     '70dB'            3ms at 70+-5dB SL.
%     'all'             All conditions stated above. Itemlists in cell array.
%     'ConditionNames'  To receive cell array with all condition names.
%
%     'model'     DTFs, sensitivities and test stimuli necessary for model 
%     	          predictions. Sensitivity paramters will be calibrated if
%     	          calibration data does not exist or does not match the
%     	          current setting of baumgartner2016.
%
%
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%   2) Data in hrtf/baumgartner2014 and hrtf/baumgartner2016
%
%   3) Data in auxdata/baumgartner2016
%
%   Examples:
%   ---------
%
%   To get all listener-specific data of the pool from Baumgartner et al. (2014), use:
%
%     data_baumgartner2016;
%
%   To get all listener-specific data of the LocaLevel study, use:
%
%     data_baumgartner2016('Long');
%
%   See also: baumgartner2016, exp_baumgartner2016
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_baumgartner2016.php

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

% AUTHOR : Robert Baumgartner

% TODO: explain Data in description;

%% ------ Check input options --------------------------------------------

definput.import={'baumgartner2016'};
definput.flags.condition = {'baumgartner2014';'Long';'10dB';'20dB';'30dB';'40dB';'50dB';'60dB';'70dB';'all'};
definput.flags.conditionNames = {'';'ConditionNames'};
definput.flags.modeldata = {'','model'};

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_recalib
  flags.cachemode = 'redo';
else
  flags.cachemode = 'normal';
end
  
%% Output Only Condition Names of LocaLevel
if flags.do_ConditionNames
  data = definput.flags.condition(2:end-1);
  return
end

%% Listener pool (listener-specific SP-DTFs) from Baumgartner et al. (2014)
if flags.do_baumgartner2014  
  listeners = {'NH12';'NH15';'NH21';'NH22';'NH33';'NH39';'NH41';'NH42';... % NumChan
               'NH43';'NH46';'NH53';'NH55';'NH58';'NH62';'NH64';'NH68';...
               'NH71';'NH72';'NH14';'NH16';'NH17';'NH18';'NH57';...
               };
  data=cell2struct(listeners,'id',2);

  for ii = 1:length(data)

    data(ii).S = 0.5; % default sensitivity

    filename = fullfile(SOFAdbPath,'baumgartner2014',...
      ['ARI_' data(ii).id '_hrtf_M_dtf 256.sofa']);

    data(ii).Obj = SOFAload(filename);
    data(ii).fs = data(ii).Obj.Data.SamplingRate;

  end

  data = loadBaselineData(data,kv.latseg,kv.dlat);

  %% prior distribution
  data = addpriordist(data);

  %% Calibration of S

  if flags.do_SPLtemAdapt
    kv.SPLtem = kv.SPL;
  end

  %%% Define cache name according to settings for auditory periphery model
  cachename = ['calibration_g' num2str(kv.gamma,'%u') ...
      '_mrs' num2str(kv.mrsmsp,'%u') ...
      '_do' num2str(kv.do,'%u') ...
      '_tar' num2str(kv.SPL,'%u') 'dB_tem' num2str(kv.SPLtem,'%u') 'dB_'...
      flags.fbank];
  if flags.do_gammatone
    cachename = [cachename '_'  num2str(1/kv.space,'%u') 'bpERB'];
    if flags.do_middleear; cachename = [cachename '_middleear']; end
    if flags.do_ihc; cachename = [cachename '_ihc']; end
  else % zilany
    cachename = [cachename '_' flags.fibertypeseparation];
  end
  if kv.prior > 0 
    cachename = [cachename '_prior' num2str(kv.prior,'%u')];
  end
  if kv.tiwin < 0.5
    cachename = [cachename '_tiwin' num2str(kv.tiwin*1e3) 'ms']; 
  end
  cachename = [cachename '_mgs' num2str(kv.mgs)]; 

  c = amt_cache('get',cachename,flags.cachemode);
  if isempty(c) %|| not(isequal(c.kv,kv))

    % reset listener-specific MRS to default
    for ii = 1:length(data)
      data(ii).mrs = kv.mrsmsp;
    end

    amt_disp('Calibration procedure started. Please wait!','progress')
    data = baumgartner2016_calibration(data,'argimport',flags,kv);

    c.data = rmfield(data,{'Obj','fs','itemlist','target','response'}); % reduce filesize
    c.kv = kv;
    amt_cache('set',cachename,c)

  else

    for ss = 1:length(data)
      for ii = 1:length(c.data)
        if strcmp(data(ss).id,c.data(ii).id)
          data(ss).S = c.data(ii).S;
          data(ss).mrs = c.data(ii).mrs;
          if isfield(c.data,'prior')
            data(ss).prior = c.data(ii).prior;
          else
            data(ss).prior = kv.prior;
          end
        end
      end
    end

  end 
  
else % Loca Level
  
  %% Extract localization data

  d = amt_load('baumgartner2016','data.mat');

  if flags.do_all

    for ll = 1:length(d.subject)

      data(ll).condition = d.condition;
      data(ll).id = d.subject(ll).id;

      data(ll).SL = [50,10:10:70];
      data(ll).SPL = data(ll).SL + d.subject(ll).SPLtoSLoffset;
      data(ll).SPL(2:end) = data(ll).SPL(2:end) + d.subject(ll).LongToShortOffset;

      for C = 1:length(d.condition)

        data(ll).itemlist{C} = real(d.subject(ll).expData{C}(:,1:8));
        data(ll).pe_exp(C) = localizationerror(data(ll).itemlist{C},'rmsPmedianlocal');
        data(ll).qe_exp(C) = localizationerror(data(ll).itemlist{C},'querrMiddlebrooks');

      end
    end

  else

    C = find(ismember(d.condition,flags.condition));

    for ll = 1:length(d.subject)

      data(ll).itemlist = real(d.subject(ll).expData{C}(:,1:8));
      data(ll).id = d.subject(ll).id;
      if flags.do_Long
        data(ll).SL = 50;
        data(ll).SPL = data(ll).SL + d.subject(ll).SPLtoSLoffset;
      else % Short
        data(ll).SL = str2num(flags.condition(1:2));
        data(ll).SPL = data(ll).SL + d.subject(ll).SPLtoSLoffset + d.subject(ll).LongToShortOffset;
      end

      data(ll).pe_exp = localizationerror(data(ll).itemlist,'rmsPmedianlocal');
      data(ll).qe_exp = localizationerror(data(ll).itemlist,'querrMiddlebrooks');

    end

  end


  %% Listener-specific SP-DTFs
  if flags.do_model

      for ii = 1:length(data)

        filename = fullfile(SOFAdbPath,'baumgartner2016',...
          ['ARI_' data(ii).id '_hrtf_M_dtf 256.sofa']);

        data(ii).Obj = SOFAload(filename);
        data(ii).fs = data(ii).Obj.Data.SamplingRate;

        if flags.do_Long
          data(ii).stim = d.subject(ii).stim.long;
        else % short
          data(ii).stim = d.subject(ii).stim.short;
        end
        data(ii).fsstim = d.subject(ii).stim.fs;

      end

      %% prior districution
      data = addpriordist(data);

      %% Calibration of S
 
      if flags.do_SPLtemAdapt
        kv.SPLtem = kv.SPL;
      end
      
      cachename = ['calibration_localevel_g' num2str(kv.gamma,'%u') ...
          '_mrs' num2str(kv.mrsmsp,'%u') ...
          '_do' num2str(kv.do,'%u') ...
          '_tar' num2str(kv.SPL,'%u') 'dB_tem' num2str(kv.SPLtem,'%u') 'dB_'...
          flags.fbank];
      if flags.do_gammatone
        cachename = [cachename '_'  num2str(1/kv.space,'%u') 'bpERB'];
        if flags.do_middleear; cachename = [cachename '_middleear']; end
        if flags.do_ihc; cachename = [cachename '_ihc']; end
      else % zilany
        cachename = [cachename '_' flags.fibertypeseparation];
      end
      if kv.prior > 0 
        cachename = [cachename '_prior' num2str(kv.prior,'%u')];
      end
      if kv.tiwin < 0.5
        cachename = [cachename '_tiwin' num2str(kv.tiwin*1e3) 'ms']; 
      end

      c = amt_cache('get',cachename,flags.cachemode);
      if isempty(c) %|| not(isequal(c.kv,kv))

        c.SL = 50; % dB SL of targets
        c.SPL = c.SL + [d.subject.SPLtoSLoffset];
        for ii = 1:length(data)
          c.stim{ii} = d.subject(ii).stim.long;
        end

        %% Baseline data for calibration
        baseline = data_baumgartner2016('Long');

        for ll = 1:length(data)  

          data(ll).pe_exp = localizationerror(baseline(ll).itemlist,'rmsPmedianlocal'); % s(ll).baseline.pe_exp
          data(ll).qe_exp = localizationerror(baseline(ll).itemlist,'querrMiddlebrooks'); % s(ll).baseline.qe_exp
          data(ll).mrs = localizationerror(data(ll).itemlist,'precLcentral');

          for ii = 1:length(kv.latseg)

            latresp = baseline(ll).itemlist(:,7);
            idlat = latresp <= kv.latseg(ii)+kv.dlat & latresp > kv.latseg(ii)-kv.dlat;
            mm2 = baseline(ll).itemlist(idlat,:);

            data(ll).target{ii} = mm2(:,6); % polar angle of target
            data(ll).response{ii} = mm2(:,8); % polar angle of response
            data(ll).Nt{ii} = length(data(ll).target{ii});

          end

        end
        %%

        % reset listener-specific MRS to default
        for ii = 1:length(data)
          data(ii).mrs = kv.mrsmsp;
        end

        amt_disp('Calibration procedure started. Please wait!','progress')
        data = baumgartner2016_calibration(data,'argimport',flags,kv,'c',c);

        c.data = rmfield(data,{'Obj','fs','itemlist','target','response'}); % reduce filesize
        c.kv = kv;
        amt_cache('set',cachename,c)

      else

        for ii = 1:length(data)
          idx = find(ismember({c.data.id},data(ii).id));
          data(ii).S = c.data(idx).S;
        end

      end

  end

end
end


function s = loadBaselineData(s,latseg,dlat)

  % latseg = 0;%[-20,0,20]; 
  % dlat = 30;%10;

  % Experimental baseline data
  numchan = data_goupell2010('BB');
  methods = data_majdak2010('Learn_M');
  spatstrat = data_majdak2013('BB');
  ctcL = data_majdak2013ctc('Learn');

  for ll = 1:length(s)

    s(ll).itemlist = [];

    s(ll).itemlist = [s(ll).itemlist ; numchan(ismember({numchan.id},s(ll).id)).mtx];
    s(ll).itemlist = [s(ll).itemlist ; methods(ismember({methods.id},s(ll).id)).mtx];
    s(ll).itemlist = [s(ll).itemlist ; spatstrat(ismember({spatstrat.id},s(ll).id)).mtx];
    s(ll).itemlist = [s(ll).itemlist ; ctcL(ismember({ctcL.id},s(ll).id)).mtx];

    s(ll).pe_exp = localizationerror(s(ll).itemlist,'rmsPmedianlocal');
    s(ll).qe_exp = localizationerror(s(ll).itemlist,'querrMiddlebrooks');
    s(ll).mrs = localizationerror(s(ll).itemlist,'precLcentral');

    for ii = 1:length(latseg)

      latresp = s(ll).itemlist(:,7);
      idlat = latresp <= latseg(ii)+dlat & latresp > latseg(ii)-dlat;
      mm2 = s(ll).itemlist(idlat,:);

      s(ll).pe_exp_lat(ii) = localizationerror(mm2,'rmsPmedianlocal');
      s(ll).qe_exp_lat(ii) = localizationerror(mm2,'querrMiddlebrooks');

      s(ll).target{ii} = mm2(:,6); % polar angle of target
      s(ll).response{ii} = mm2(:,8); % polar angle of response
      s(ll).Ntar{ii} = length(s(ll).target{ii});

    end
    
  end
end


function data = addpriordist(data)
  dang = 30; % angular width of segments
  Tmin = 5; % min. # of targets to estimate prior distribution
  edges = -90:dang:270;
  for ii = 1:length(data)
    try
      T = histcounts(data(ii).itemlist(:,6),edges);
      R = histcounts(data(ii).itemlist(:,8),edges);
    catch
      centers = edges(1:end-1)+diff(edges)/2;
      T = hist(data(ii).itemlist(:,6),centers);
      R = hist(data(ii).itemlist(:,8),centers);
    end
    T(T<Tmin) = nan;
    RvT = R./T;
    RvT(isnan(RvT)) = 1;
    data(ii).priordist.y = RvT;
    data(ii).priordist.x = edges(1:end-1)+dang/2;
  end
end
