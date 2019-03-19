function [h]=plot_ziegelwanger2013(Obj,data,type,color,ele,ch,time,style,width)
%plot_ziegelwanger2013 XXX Headline missing
%   Usage: plot_ziegelwanger2013(Obj,data,type,color,ele,ch,time,style,width)
%
%   PLOT_ZIEGELWANGER2013(Obj,data,type,color,ele,ch,time,style,width)
%   plots TOA-data in horizontal planes.
%
%   Input:
%       Obj   : SOFA object
%       data  : XXX Description missing
%       type  : XXX Description missing
%       color : XXX Description missing
%       ele   : XXX Description missing
%       ch    : XXX Description missing
%       time  : XXX Description missing
%       style : XXX Description missing
%       width : XXX Description missing
% 
%   Output:
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
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/plot/plot_ziegelwanger2013.php

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

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna,
% Austria


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
