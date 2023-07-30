function output = exp_langendijk2002(varargin)
%EXP_LANGENDIJK2002  Experiment from Langendijk & Bronkhorst (2002)
%   Usage: output = exp_langendijk2002(flag);
%
%   EXP_LANGENDIJK2002(flags) recreates figures from Langendijk &
%   Bronkhorst (2002)
%
%   The following flags can be specified;
%
%     'plot'    Plot the output of the experiment. This is the default.
%
%     'no_plot'  Don't plot, only return data.
%
%     'fig7'    Listener P6
%
%     'fig9'    Listener P3
%
%   You can choose between two of his listeners P3 and P6. The required
%   data (DTF data and response patterns) will be provided by precalculated
%   mat-files due to high computing time (optionally data can be 
%   recalculated by using DATA_LANGENDIJK2002. 
%
%   The following subfigures show the probability density functions (PDFs)
%   and actual response angles as functions of the target angle for the 
%   exemplary listener and different conditions. The shading of each cell 
%   codes the response probability (light/dark denotes high/low probability):
%
%    Subfigure 1: Baseline condition
%
%    Subfigure 2: 2-octave condition (4-16kHz)
%
%    Subfigure 3: 1-octave condition (low 4-8kHz)
%
%    Subfigure 4: 1-octave condition (middle 5.7-11.3kHz)
%
%    Subfigure 5: 1-octave condition (high 8-16kHz)
%   
%   Subfigure 6 shows the likelihood statistics for the actual responses
%   (bars) in comparison to the means (dots) and the 99% confidence  
%   intervals (whiskers) of the expected likelihood. See paper for further 
%   details about the likelihood statistics.
%
%   The output are the pdfs for the baseline condition.
%
%   Examples:
%   ---------
%
%   To display Figure 7 use :
%
%     exp_langendijk2002('fig7');
%
%   To display Figure 9 use :
%
%     exp_langendijk2002('fig9');
%
%   See also: langendijk2002, langendijk2002_likelihood, plot_langendijk2002, 
%             plot_langendijk2002_likelihood, data_langendijk2002
%
%   References:
%     E. Langendijk and A. Bronkhorst. Contribution of spectral cues to human
%     sound localization. J. Acoust. Soc. Am., 112:1583--1596, 2002.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_langendijk2002.php


%   #Author: Robert Baumgartner (2016)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% ------ Check input options --------------------------------------------

definput.flags.type = {'missingflag','fig7','fig9'};
definput.flags.plot = {'plot','no_plot'};

% Parse input options
[flags]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SETTINGS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig7
  listener='P6';  % ID of listener (P3 or P6)
elseif flags.do_fig9
  listener='P3';
end

fs = 48000;     % sampling frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dtfdata=amt_load('langendijk2002',[listener '_data.mat']);
% loads hM data for all conditions 
% data can be recalculated by calling data_langendijk2002('expdata')

% pdf calcualtion
%h = waitbar(0,'Please wait...');
pb  = langendijk2002( dtfdata.medir,dtfdata.medir,fs); % baseline
%waitbar(1/5)
p2o = langendijk2002( dtfdata.medir2o,dtfdata.medir,fs); % 2-oct (4-16kHz)
%waitbar(2/5)
p1ol= langendijk2002( dtfdata.medir1ol,dtfdata.medir,fs); % 1-oct (low:4-8kHz)
%waitbar(3/5)
p1om= langendijk2002( dtfdata.medir1om,dtfdata.medir,fs); % 1-oct (middle:5.7-11.3kHz)
%waitbar(4/5)
p1oh= langendijk2002( dtfdata.medir1oh,dtfdata.medir,fs); % 1-oct (high:8-16kHz)
%waitbar(5/5)

% likelihood estimations
la=zeros(5,1);le=zeros(5,1);ci=zeros(5,2);
idb=1:2:length(dtfdata.targetb); % in order to get comparable likelihoods
[la(1),le(1),ci(1,:)] = langendijk2002_likelihood( pb,dtfdata.pol,dtfdata.pol,dtfdata.targetb(idb),dtfdata.responseb(idb) );
[la(2),le(2),ci(2,:)] = langendijk2002_likelihood( p2o,dtfdata.pol,dtfdata.pol,dtfdata.targetc,dtfdata.response2o );
[la(3),le(3),ci(3,:)] = langendijk2002_likelihood( p1ol,dtfdata.pol,dtfdata.pol,dtfdata.targetc,dtfdata.response1ol );
[la(4),le(4),ci(4,:)] = langendijk2002_likelihood( p1om,dtfdata.pol,dtfdata.pol,dtfdata.targetc,dtfdata.response1om );
[la(5),le(5),ci(5,:)] = langendijk2002_likelihood( p1oh,dtfdata.pol,dtfdata.pol,dtfdata.targetc,dtfdata.response1oh );
%close(h)

output = pb;

if flags.do_plot
  figure('Name',listener)
  clf
  % pdf plots with actual responses
  subplot(2,3,1)
  hold all;    
  plot_langendijk2002(pb,dtfdata.pol,dtfdata.pol,'no_colorbar');
  title(['Baseline']);    
  h=plot( dtfdata.targetb, dtfdata.responseb, 'ko');
  set(h, 'MarkerFaceColor','w');

  subplot(2,3,2)
  hold all;
  plot_langendijk2002(p2o,dtfdata.pol,dtfdata.pol,'no_colorbar');
  title(['2-oct (4-16kHz)']);
  h=plot( dtfdata.targetc, dtfdata.response2o, 'ko');
  set(h, 'MarkerFaceColor','w');

  subplot(2,3,3)
  hold all;
  plot_langendijk2002(p1ol,dtfdata.pol,dtfdata.pol,'no_colorbar');
  title(['1-oct (low: 4-8kHz)']);
  h=plot( dtfdata.targetc, dtfdata.response1ol, 'ko');
  set(h, 'MarkerFaceColor','w');

  subplot(2,3,4)
  hold all;
  plot_langendijk2002(p1om,dtfdata.pol,dtfdata.pol,'no_colorbar');
  title(['1-oct (middle: 5.7-11.3kHz)']);
  h=plot( dtfdata.targetc, dtfdata.response1om, 'ko');
  set(h, 'MarkerFaceColor','w');

  subplot(2,3,5)
  hold all;
  plot_langendijk2002(p1oh,dtfdata.pol,dtfdata.pol,'no_colorbar');
  title(['1-oct (high: 8-16kHz)']);
  h=plot( dtfdata.targetc, dtfdata.response1oh, 'ko');
  set(h,'MarkerFaceColor','w')

  % likelihood statistic
  subplot(2,3,6)
  plot_langendijk2002_likelihood(la,le,ci);
  set(gca,'XLim',[0.5 5.5])
end

end

