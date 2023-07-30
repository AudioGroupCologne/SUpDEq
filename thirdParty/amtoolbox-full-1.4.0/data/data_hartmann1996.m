function data = data_hartmann1996(varargin)
%DATA_HARTMANN1996 data from Hartmann and Wittenberg (1996) 
%   Usage: data = data_hartmann1996(condition,f0)
%
%   Output parameters:
%     data    : structure with individual and average data
%
%   The condition flag may be one of:
%
%     'ILD' or 'fig7'   ILDs up to nprime set to zero. This is the default.
%
%     'ISLD' or 'fig8'  ISLDs (interaural spectral level differeces) maintained 
%                       while flattening right-ear HRTF up to nprime.
%
%   The f0 flag may be one of:
%
%     '125'   Fundamental frequency of 125 Hz and highest harmonic at 4750 Hz. 
%             This is the default.
%
%     '250'   Fundamental frequency of 250 Hz and highest harmonic at 9500 Hz. 
%             This option is only available for condition ILD/fig7.
%
%
%   Data from Hartmann and Wittenberg (1996) with interaural cue
%   alterations up to a certain harmonic nprime.
%
%   References:
%     W. M. Hartmann and A. Wittenberg. On the externalization of sound
%     images. J. Acoust. Soc. Am., 99(6):3678--88, June 1996.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/data/data_hartmann1996.php


%   #Author: Robert Baumgartner (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

definput.flags.f0 = {'125','250'};
definput.flags.expirement = {'ILD','fig7','ISLD','fig8'};

[flags]=ltfatarghelper({},definput,varargin);

%% Actual data
if flags.do_ILD || flags.do_fig7
  
  if flags.do_125
    % listeners C
    C = [1.0, 3.0; 8.0, 3.0;14, 2.0;19, 2.3;        25, 2.0; 38, 2.0];
    % listeners R
    R = [          8.0, 2.7;14, 2.0;19, 1.5;22, 1.0;25, 1.0; 38, 0.50];
    % overlap for avg
    idC = 2:6;
    idR = [1:3,5:6];
    
  elseif flags.do_250 
    % listeners C
    C = [8.0, 2.8;        19, 2.0; 22, 2.0;  25, 1.8; 38, 0.58];
    % listeners R
    R = [8.0, 2.5;14, 1.9;19, 1.0; 22, 0.88; 25, 0.48;38, 0.38];
    % overlap for avg
    idC = 1:5;
    idR = [1,3:6];
    
  end
  
  data.subj(1).ID = 'C';
  data.subj(1).nprime = C(:,1);
  data.subj(1).Escore = C(:,2);
  
  data.subj(2).ID = 'R';
  data.subj(2).nprime = R(:,1);
  data.subj(2).Escore = R(:,2);
  
  data.avg.nprime = data.subj(1).nprime(idC);
  data.avg.Escore = mean([data.subj(1).Escore(idC),data.subj(2).Escore(idR)],2);
  
elseif flags.do_ISLD || flags.do_fig8
  
  ID = {'A','C','R','W'};
  indData = {...
    [5.0, 2.9;8.0, mean([2.7,0.98]);14, mean([2.5,1.5]);19, 0.38;38, 0.46];... % averaged in case of "split images"
  	[8.0, 3.0;14, mean([2.9, 0.48]);19, 0.97;25, 1.8;38, 0.96];...
  	[1.0, mean([3.1, 1.0]);5.0, mean([3.1, 1.0]);8.0, 1.0;19, 1.0;25, 1.7;38, 1.0];...
  	[5.0, 1.6;8.0, 2.3;19, 1.5;25, 2.0;38, 1.6]};
  
  Ns = length(ID);
  for ii = 1:Ns
    data.subj(ii).ID = ID{ii};
    data.subj(ii).nprime = indData{ii}(:,1);
    data.subj(ii).Escore = indData{ii}(:,2);
  end
  
  data.avg.nprime = [5;8;19;38];
  Escore = nan(length(data.avg.nprime),Ns); % initialize with data form W
  for ii = 1:Ns
    for ee = 1:length(data.avg.nprime)
      id = find(data.subj(ii).nprime == data.avg.nprime(ee));
      if isscalar(id)
        Escore(ee,ii) = data.subj(ii).Escore(id);
      end
    end
  end
  data.avg.Escore = nanmean(Escore,2);
  
end
end


