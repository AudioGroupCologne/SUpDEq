function [h]=plot_ziegelwanger2013(Obj,data,type,color,ele,ch,time,style,width)
%PLOT_ZIEGELWANGER2013 plot time-of-arrival
%   Usage: plot_ziegelwanger2013(Obj,data,type,color,ele,ch,time,style,width)
%
%
%   Input parameters:
%       Obj   : SOFA object
%       data  : Input data
%       type  : Plot type
%       color : Line color
%       ele   : Elevation [deg]
%       ch    : Channel Index
%       time  : Flag (1...time vector [ms] 0...[samples])
%       style : Linestyle
%       width : Linewidth
% 
%   Output parameters:
%       h: figure handle
%
%   Estimates the Time-of-Arrival for each column in input data hM and corrects 
%   the results with a geometrical model of the head.
%
%   Examples:
%   ---------
% 
%   To plot the modelled TOA in the horizontal plane after using
%   ziegelwanger2013, use:
%
%       plot_ziegelwanger2013(Obj,Obj.Data.Delay,1,'b',0,1,1);
%
%   See also: ziegelwanger2013_onaxis, ziegelwanger2013_offaxis,
%             data_ziegelwanger2013, exp_ziegelwanger2013
%
%   References:
%     P. Majdak and H. Ziegelwanger. Continuous-direction model of the
%     broadband time-of-arrival in the head-related transfer functions. In
%     ICA 2013 Montreal, volume 19, page 050016, Montreal, Canada, 2013. ASA.
%     
%     H. Ziegelwanger and P. Majdak. Modeling the broadband time-of-arrival
%     of the head-related transfer functions for binaural audio. In
%     Proceedings of the 134th Convention of the Audio Engineering Society,
%     page 7, Rome, 2013.
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/plot/plot_ziegelwanger2013.php


%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Verified
%   #Requirements: SOFA M-Signal M-Optimization
%   #Author: Harald Ziegelwanger (2013), Acoustics Research Institute, Vienna, Austria

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 


%% ----------------------------check variables-----------------------------
if time
    data=data/Obj.Data.SamplingRate*1000;
end

if exist('style','var')
    if ischar(style)
        style={style};
    end
end

if ~ischar(color)
    rgbcolor=color;
    color='b'*ones(1,size(color,1));
end

%% --------------------------------plot TOA--------------------------------
[~,idx1]=sort(Obj.SourcePosition(:,1));
idx2=find(Obj.SourcePosition(idx1,2)==ele);

for ii=1:length(type)
    switch type(ii)
        case 1
            h=plot(Obj.SourcePosition(idx1(idx2),1),data(idx1(idx2),ii*2-1+(ch-1)),[color(ii) '-']);
        case 2
            h=plot(Obj.SourcePosition(idx1(idx2),1),data(idx1(idx2),ii*2-1+(ch-1)),[color(ii) '--']);
        case 3
            idx3=find(Obj.SourcePosition(idx1(idx2),2)==ele & data(idx1(idx2),ii*2-1+(ch-1))~=0);
            if exist('style','var')
                h=plot(Obj.SourcePosition(idx1(idx2(idx3)),1),data(idx1(idx2(idx3)),ii*2-1+(ch-1)),[color(ii) style{ii}],'MarkerSize',width(ii));%,'MarkerFaceColor',color(ii)
            else
                h=plot(Obj.SourcePosition(idx1(idx2(idx3)),1),data(idx1(idx2(idx3)),ii*2-1+(ch-1)),[color(ii) '.']);
            end
        case 4
            h=plot(Obj.SourcePosition(idx1(idx2),1),data(idx1(idx2),ii*2-1+(ch-1)),[color(ii) style{ii}],'LineWidth',width(ii));
    end
    if exist('rgbcolor','var')
        set(h,'color',rgbcolor(ii,:));
    end
    hold on
end

title(['Elevation: ' num2str(ele) 'deg'])
xlim([0 359])
xlabel('Azimuth (deg)')
if time
    ylabel('Time (ms)')
else
    ylabel('Time (samples)')
end
grid on
set(gca,'XTick',[0 90 180 270 360]);

end %of funciton


