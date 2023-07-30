function data = data_macpherson2000(varargin)
%DATA_MACPHERSON2000 results of Experiment I from macpherson2000level
%   Usage: data = data_macpherson2000
%
%   Output parameters:
%     data.id       : listener identification
%     data.SL       : across-listener average of probed sensation levels (Exp. I)
%     data.dur      : probed stimulus duration in msec
%     data.err      : localization error
%     data.itemlist : raw target-response data (only for flag 'original')
%     data.pgf      : polar angle gain for the front (only for flag 'retrieved') 
%     data.pgr      : polar angle gain for the back (only for flag 'retrieved')
%
%   DATA_MACPHERSON2000 returns listeners' polar angle gains obtained in 
%   free-field localization experiments. Data was retrieved from FIG.4 in 
%   macpherson2000level. Five listeners were tested.
%
%   DATA_MACPHERSON2000 accepts the following optional flags:
%
%     'fig4'      Polar gain and scatter of responses to noise-burst targets. 
%                 (a,b) Front-hemisphere targets. (c,d) Rear-hemisphere targets. 
%                 The dotted line at approximately 26 degrees shows the scatter 
%                 expected if responses were uniformly distributed throughout 
%                 the quasi-veridical region rather than concentrated near 
%                 the regression line. This is the default.
%
%     'fig7'      Effect of level on polar angle gain for 3- and 100-ms Gaussian 
%                 noise bursts. Each panel shows data for one listener. 
%                 Open circles, 3-ms stimuli; Diamonds, 100-ms stimuli. 
%                 Error bars show the standard error of the regression 
%                 coefficient (polar angle gain).
%
%     '3ms'       Short stimuli. This is the default. 
%
%     '100ms'     Long stimuli.
%
%     'low'       Low-level stimuli around 30 dB SL. This is the default. 
%
%     'high'      High-level stimuli around 55 dB SL.  
%
%     'original'  Data from Ewan Macpherson. This is the default. 
%
%     'retrieved' Data retrieved from Figure in the paper. 
%
%     'errflag'   One of the error flags defined in �localizationerror�. 
%
%   Examples:
%   ---------
%
%   To get the long-duration, high-level data shown in Fig. 4 use :
%
%     data = data_macpherson2000('fig4','100ms','high');
%
%   References:
%     E. A. Macpherson and J. C. Middlebrooks. Localization of brief sounds:
%     Effects of level and background noise. J. Acoust. Soc. Am., 108(1834),
%     2000.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_macpherson2000.php


%   #Author: Robert Baumgartner (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%% Check input options

definput.import = {'localizationerror'};

% Define input flags
definput.flags.fig = {'fig4','fig7'};
definput.flags.dur = {'3ms','100ms'};
definput.flags.SL = {'','low','high'};
definput.flags.source = {'original','retrieved'};

% Parse input options
[flags,kv]  = ltfatarghelper({},definput,varargin);


%% Extract data

if flags.do_original
  
  d = amtload('macpherson2000','data.mat');
  
  if flags.do_fig4
    
    data = d.fig4.data;
    
    for ii=1:length(data)
      
      idstim = data(ii).stimtype == 0; % Gaussian burst
      data(ii).stimtype = data(ii).stimtype(idstim,:);
      data(ii).itemlist = data(ii).itemlist(idstim,:);
      data(ii).SL = data(ii).SL(idstim,:);
      data(ii).dur = data(ii).dur(idstim,:);
      
      if flags.do_3ms
        iddur = data(ii).dur == 3;
      else
        iddur = data(ii).dur == 100;
      end
      data(ii).stimtype = data(ii).stimtype(iddur,:);
      data(ii).itemlist = data(ii).itemlist(iddur,:);
      data(ii).SL = data(ii).SL(iddur,:);
      data(ii).dur = data(ii).dur(iddur,:);
      
      if flags.do_low
        idsl = data(ii).SL <= 42.5;
      elseif flags.do_high
        idsl = data(ii).SL >= 42.5;
      else % all
        idsl = data(ii).SL > 0;
      end
      data(ii).stimtype = data(ii).stimtype(idsl,:);
      data(ii).itemlist = data(ii).itemlist(idsl,:);
      data(ii).SL = data(ii).SL(idsl,:);
      data(ii).dur = data(ii).dur(idsl,:);
      
      data(ii).err = localizationerror(data(ii).itemlist,flags.errorflag);
      
    end
    
    data = rmfield(data,'stimtype');
    
  elseif flags.do_fig7
    
    data = d.fig7.data;
    
    for ii=1:length(data)
      
      if flags.do_3ms
        iddur = data(ii).dur == 3;
      else
        iddur = data(ii).dur == 100;
      end
      data(ii).itemlist = data(ii).itemlist(iddur,:);
      data(ii).SL = data(ii).SL(iddur,:);
      data(ii).dur = data(ii).dur(iddur,:);
      
      if flags.do_low
        idsl = data(ii).SL <= 42.5;
      elseif flags.do_high
        idsl = data(ii).SL >= 42.5;
      else % all
        idsl = data(ii).SL > 0;
      end
      data(ii).itemlist = data(ii).itemlist(idsl,:);
      data(ii).SL = data(ii).SL(idsl,:);
      data(ii).dur = data(ii).dur(idsl,:);
      
    end
    
  end
  
end
  
if flags.do_retrieved % summary statistics retrieved from fig 3
      
  amtdisp('Summary statistics retrieved from Fig. 3.')
  
    Ns = 5; % # of listeners

    % SL ranges reported in Sec. II.A.2
    lowSL = mean([25 35]);
    highSL = mean([50 60]);

    short = 3;  % msec
    long = 100; % msec

    pgf(:,1) = [2.2 2.75 3.45 2.7 2.9]/4.65*1.5;
    pgf(:,2) = [2.4 2.8 2.15 2.65 2.9]/4.65*1.5;
    pgf(:,3) = [1.45 2.35 1.7 2.75 2.7]/4.65*1.5;
    pgf(:,4) = [2.4 1.55 0.75 2.0 1.55]/4.65*1.5;

    pgr(:,1) = [2.0 3.4 3.05 3.08 3.8]/4.65*1.5;
    pgr(:,2) = [2.15 2.9 2.95 2.65 3.5]/4.65*1.5;
    pgr(:,3) = [1.55 2.35 2.45 2.6 3.2]/4.65*1.5;
    pgr(:,4) = [1.15 2.05 2.1 2.4 2.7]/4.65*1.5;

    IDs = {'S30';'S18';'S53';'S52';'S51'};

    for ii = 1:Ns
      data(ii).SL = [lowSL,highSL,lowSL,highSL];
      data(ii).dur = [long,long,short,short];
      data(ii).id = IDs{ii};
      data(ii).pgf = pgf(ii,:)';
      data(ii).pgr = pgr(ii,:)';
    end

end

end

