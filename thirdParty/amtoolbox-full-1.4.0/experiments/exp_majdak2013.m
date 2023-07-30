function data = exp_majdak2013(varargin)
%EXP_MAJDAK2013 Reproduce figures from Majdak et al. (2013)
%   Usage: data = exp_majdak2013(flag) 
%
%   The following flags can be specified
%
%     'fig6'   Participants localization performance before, 
%              during, and after the training. Quadrant errors.
%              Polar errors. Lateral errors.
% 
%
%   Requirements: 
%   -------------
%
%   In Matlab, Statistics and Machine Learning Toolbox and 
%   Curve Fitting Toolbox are required.
%
% 
%   See also: data_majdak2013
%
%   References:
%     P. Majdak, T. Walder, and B. Laback. Effect of long-term training on
%     sound localization performance with spectrally warped and band-limited
%     head-related transfer functions. J. Acoust. Soc. Am., 134:2148--2159,
%     2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/experiments/exp_majdak2013.php


%   #Requirements: M-statistics M-curve
%   #Author: David Poirier-Quinot (2022): first implementation
%   #Author: Clara Hollomey (2023): integration in the AMT
%   #Author: Piotr Majdak (2023): various fixes

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose.

% parse inputs
%definput.flags.type = {'fig6_quadrant', 'fig6_polar', 'fig6_lateral'};
definput.flags.datatype = {'missingflag', 'fig6'};
[flags,kv]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

if flags.do_fig6
% define data to display based on parsed input
%tmp = strsplit(flags.type, '_'); errorTypeStr = tmp{2};
errorTypeStr = {'quadrant', 'polar', 'lateral'};

for ii = 1:numel(errorTypeStr)
    
switch errorTypeStr{ii}

    case 'quadrant'
        err='querrMiddlebrooks';
    
    case 'polar'
        err='rmsPmedianlocal';
    
    case 'lateral'
        err='precL';
end

% load data
data = data_majdak2013('fig6');

% init locals
errbar = 'std';  errbaridx = 2;
% errbar = 'ci';  errbaridx = 3;
N = 4200;
nonTrainingSessionNames = {'learn', 'post_learn', 'pre_test_dummy', 'post_test_dummy', 'pre_test_warped', 'post_test_warped'};
traininSessionNames = {'training_1', 'training_2', 'training_3', 'training_4', 'training_5', 'training_6', 'training_7', 'training_8', 'training_9', 'training_10', 'training_11', 'training_12', 'training_13', 'training_14', 'training_15', 'training_16', 'training_17', 'training_18', 'training_19', 'training_20', 'training_21'};
configs = struct('sessionNames', {{}}, 'stat', [], 'isControl', 0);

% loop over control conditions
for iControl = 0:1

    % loop over pre/post training sessions
    for iSession = 1:length(nonTrainingSessionNames)
        configs(end+1) = struct('sessionNames', {{nonTrainingSessionNames{iSession}}}, 'stat', [], 'isControl', iControl);
    end

    % add training sessions (grouped)
    configs(end+1) = struct('sessionNames', {traininSessionNames}, 'stat', [], 'isControl', iControl);

end

% remove first (dummy)
configs(1) = [];


%% compute stats

amt_disp('');

% loop over config (sessions)
for iConfig = 1:length(configs)
    
    % init local 
    config = configs(iConfig);
    meanallid = [];

    % loop over subsession (for training)
    for iSession = 1:length(config.sessionNames)
        
        % log
        amt_disp(sprintf('process stat \n  session:    %d/%d \n  subsession: %d/%d', iConfig, length(configs), iSession, length(config.sessionNames)), 'volatile');

        % init filter
        selVect = data.is_control == config.isControl;
        selVect = selVect & ismember(data.session_str, config.sessionNames{iSession});
        
        % discard first 50 trials for learning session 
        if( ismember(config.sessionNames, 'learn') ); selVect = selVect & data.trial > 50; end
        
        % stats
        meanallid(iSession, :) = bootstrp(sum(selVect), @(m) localizationerror(m,err), data.pos(selVect, :));

    end

    % compute stats
    configs(iConfig).stat = getStat(meanallid.');
    
end


%% model training performances

% control group
x=(200:200:N)';
yhat=[configs(14).stat(:,1); ];
yfit=fit(x, yhat-min(yhat),'exp1');
[yfitC, r, J, COVB]=nlinfit(x,yhat,@exponential,[yfit.a yfit.b min(yhat)]);
ci2=nlparci(real(yfitC),r,'covar',COVB);

% warped group
x=(200:200:N)';
yhatW=[configs(7).stat(:,1); ];
yfitW=fit(x, yhatW-min(yhatW),'exp1');
[yfitW, rW, J, COVBW]=nlinfit(x,yhatW,@exponential,[yfitW.a yfitW.b min(yhatW)]);
ciW=nlparci(real(yfitW),rW,'covar',COVBW);


%% plot data

% init 
figure; hold on;
set(gcf,'Position',[10   226   616   484],'PaperType','A4');
prepos=-200;
postpos=N+400;

statPP = reshape([configs(8:13).stat], 4, 6).';
statPPW = reshape([configs(1:6).stat], 4, 6).';
stat = configs(14).stat;
statW = configs(7).stat;

% model: control
newx=(100:50:N+100)';

% model: warped
[ypred, delta]=nlpredci(@exponential,newx,real(yfitW),real(rW),'covar',real(COVBW));
hp=patch([newx; flipud(newx)], [ypred+delta; flipud(ypred-delta)],[1 1 1]*0.7);
set(hp, 'Edgecolor', [1 1 1]*0.7);
leg_model=plot(newx,ypred,'k','LineWidth',1);
[ypred, delta]=nlpredci(@exponential,newx,real(yfitC),real(r),'covar',real(COVB));
hp=patch([newx; flipud(newx)], [ypred+delta; flipud(ypred-delta)],[1 1 1]*0.8);
set(hp, 'Edgecolor', [1 1 1]*0.8);
plot(newx,ypred,'k','LineWidth',1);

% control: pre- & post
h=errorbar([prepos; N+100], [statPP(1,1); NaN], [statPP(2,errbaridx); NaN], 'g^'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','g'); leg_BBcontrol=h;
h=errorbar([postpos; N+100], [statPP(2,1); NaN], [statPP(2,errbaridx); NaN], 'g^'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','g');
h=errorbar([prepos; N+100], [statPP(3,1); NaN], [statPP(3,errbaridx); NaN], 'ro'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','r'); leg_BLcontrol=h;
h=errorbar([postpos; N+100], [statPP(4,1); NaN], [statPP(4,errbaridx); NaN], 'ro'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','r');
h=errorbar([prepos; N+100], [statPP(5,1); NaN], [statPP(5,errbaridx); NaN], 'bs'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','b'); leg_WPcontrol=h;
h=errorbar([postpos; N+100], [statPP(6,1); NaN], [statPP(6,errbaridx); NaN], 'bs'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','b');

% warped: pre- & post
h=errorbar([prepos+50; N+150], [statPPW(1,1); NaN], [statPPW(2,errbaridx); NaN], 'g^'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','w'); leg_BBtarget=h;
h=errorbar([postpos+50; N+150], [statPPW(2,1); NaN], [statPPW(2,errbaridx); NaN], 'g^'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','w');
h=errorbar([prepos+50; N+100], [statPPW(3,1); NaN], [statPPW(3,errbaridx); NaN], 'ro'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','w'); leg_BLtarget=h;
h=errorbar([postpos+50; N+100], [statPPW(4,1); NaN], [statPPW(4,errbaridx); NaN], 'ro'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','w');
h=errorbar([prepos+50; N+100], [statPPW(5,1); NaN], [statPPW(5,errbaridx); NaN], 'bs'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','w'); leg_WPtarget=h;
h=errorbar([postpos+100; N+100], [statPPW(6,1); NaN], [statPPW(6,errbaridx); NaN], 'bs'); errorbar_tick(h,80,'UNITS');
set(h,'LineWidth', 1,  'LineStyle','none','MarkerFaceColor','w');

% training: control group
h1=errorbar([(200:200:N)'; N+100],[stat(:,1); NaN],[stat(:,errbaridx); NaN],'ro'); %errorbar_tick(h1,80,'UNITS');
set(h1,'LineWidth', 1, 'LineStyle','none','MarkerFaceColor','r','LineStyle','-');
% training: warped group
h1=errorbar([(200:200:N)'; N+100],[statW(:,1); NaN],[statW(:,errbaridx); NaN],'bs'); % errorbar_tick(h1,80,'UNITS');
set(h1,'LineWidth', 1, 'LineStyle','none','MarkerFaceColor','w','LineStyle','-');


xlabel('Training day number','FontName','Arial');
[mxx,meta]=localizationerror(zeros(1, 9), err);
ylabel(meta.ylabel,'FontName','Arial');
box on;
set(gca,'XLim',[prepos-300 postpos+300],'FontName','Arial');
set(gca,'XTick',[200:2*200:N ]);
celli=1:N/200';
c=num2cell(celli(1:2:end))';
set(gca,'XTickLabel',c);
set(gca, 'TickLength', [0.02 0.05]);
set(gca,'LineWidth',1);
yaxval=get(gca,'YLim');


switch err
    case 'querrMiddlebrooks'
        set(gca,'YLim', [2 34.9]);
    case 'rmsPmedianlocal'
        set(gca,'YLim', [26 49]);
    case 'precL'
        set(gca,'YLim', [10.2 20.2]);
        set(gca, 'YTick', 10:2:20);
        h=legend([leg_BBcontrol leg_BLcontrol leg_WPcontrol leg_BBtarget leg_BLtarget leg_WPtarget leg_model], ...
            'Control Broadband', 'Control Band-Limited', 'Control Warped', ...
            'Target Broadband', 'Target Band-Limited', 'Target Warped');
        set(h,'Fontsize',10,'LineWidth',1);
    case 'rmsL'
        set(gca,'YLim', [7.2 19.2]);
end

% save figure
%file = mfilename('fullpath');% current filename as results_path

amt_disp(real([yfitC(1) yfitW(1) -1/yfitC(2) -1/yfitW(2) yfitC(3) yfitW(3)]));

end
end
end

%% local functions

function [stat] = getStat(meanallid)
    
% init locals
cialpha = 0.05; % 95% confidence interval

% compute stats
[muprepost,sigprepost,muciprepost,sigciprepost] = normfit(meanallid, cialpha);
stat(:,1) = muprepost;
stat(:,2) = sigprepost;
stat(:,3) = diff(muciprepost)/2;
stat(:,4) = diff(sigciprepost)/2;

end


function y=exponential(a,x)
y=real(a(1).*exp(a(2)*x)+a(3));
end

function errorbar_tick(h,w,xtype)
%ERRORBAR_TICK Adjust the width of errorbars
%   ERRORBAR_TICK(H) adjust the width of error bars with handle H.
%      Error bars width is given as a ratio of X axis length (1/80).
%   ERRORBAR_TICK(H,W) adjust the width of error bars with handle H.
%      The input W is given as a ratio of X axis length (1/W). The result 
%      is independent of the x-axis units. A ratio between 20 and 80 is usually fine.
%   ERRORBAR_TICK(H,W,'UNITS') adjust the width of error bars with handle H.
%      The input W is given in the units of the current x-axis.
%
%   See also ERRORBAR
%

% Author: Arnaud Laurent
% Creation : Jan 29th 2009
% MATLAB version: R2007a
%
% Notes: This function was created from a post on the french forum :
% http://www.developpez.net/forums/f148/environnements-developpement/matlab/
% Author : Jerome Briot (Dut) 
%   http://www.mathworks.com/matlabcentral/newsreader/author/94805
%   http://www.developpez.net/forums/u125006/dut/
% It was further modified by Arnaud Laurent and Jerome Briot.

% Check numbers of arguments
error(nargchk(1,3,nargin))

% Check for the use of V6 flag ( even if it is depreciated ;) )
flagtype = get(h,'type');

% Check number of arguments and provide missing values
if nargin==1
	w = 80;
end

if nargin<3
   xtype = 'ratio';
end

% Calculate width of error bars
if ~strcmpi(xtype,'units')
    dx = diff(get(gca,'XLim'));	% Retrieve x limits from current axis
    w = dx/w;                   % Errorbar width
end

% Plot error bars
if strcmpi(flagtype,'hggroup') % ERRORBAR(...)
    
    hh=get(h,'children');		% Retrieve info from errorbar plot
    x = get(hh(2),'xdata');		% Get xdata from errorbar plot
    
    x(4:9:end) = x(1:9:end)-w/2;	% Change xdata with respect to ratio
    x(7:9:end) = x(1:9:end)-w/2;
    x(5:9:end) = x(1:9:end)+w/2;
    x(8:9:end) = x(1:9:end)+w/2;

    set(hh(2),'xdata',x(:))	% Change error bars on the figure

else  % ERRORBAR('V6',...)
    
    x = get(h(1),'xdata');		% Get xdata from errorbar plot
    
    x(4:9:end) = x(1:9:end)-w/2;	% Change xdata with respect to the chosen ratio
    x(7:9:end) = x(1:9:end)-w/2;
    x(5:9:end) = x(1:9:end)+w/2;
    x(8:9:end) = x(1:9:end)+w/2;

    set(h(1),'xdata',x(:))	% Change error bars on the figure
    
end

end

