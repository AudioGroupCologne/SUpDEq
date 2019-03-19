% AKtightenFigure(hfig,extention)
% resizes figure (either current or by given figure handle) to have no
% margin around any contained objects.
% Since the content detection might not work perfectly in every case an
% extention can be given (exact values need to be individually evaluated).
% The extention can also be given as the first and only parameter.
%
% e.g. use
% AKtightenFigure;
% AKtightenFigure(hfig);
% AKtightenFigure([0,.5,.5,0]);
% AKtightenFigure(hfig,[0,.5,.5,0]);
%
% I N P U T
% hfig       - figure handle, default = gcf
% extention  - 4 element vector targeting the directions
%              [ left , bottom , right , top ]
%              negative values can be given to shrink the figure even more,
%              default = [.01,.01,.01,.01]
%
% v1 03/2017 helmholz@campus.tu-berlin.de, Audio Communicatin Group,
%            TU Berlin
% v2 01/2018 helmholz@campus.tu-berlin.de, Audio Communicatin Group,
%            TU Berlin

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function AKtightenFigure(hfig,extention)

if nargin < 1
    hfig = gcf;
    extention = [.01,.01,.01,.01];
elseif nargin < 2
    if isa(hfig,'matlab.ui.Figure')
        extention = [.01,.01,.01,.01];
    else
        extention = hfig;
        hfig = gcf;
    end
end
if length(extention) ~= 4
    error('AKtightenFigure:extention','Extention vector with wrong size given.');
end

% find all the axes in the figure
hax = findall(hfig,'type','axes');

% compute the tighest box that includes all axes
tighest_box = [Inf,Inf,-Inf,-Inf]; % left bottom right top
for i=1:length(hax)
    set(hax(i),'units','centimeters');
    
    p = get(hax(i),'position');
    ti = get(hax(i),'tightinset');
    
    % get position as left, bottom, right, top
    p = [p(1),p(2),p(1)+p(3),p(2)+p(4)] + ti.*[-1,-1,1,1];
 
    tighest_box(1) = min(tighest_box(1),p(1));
    tighest_box(2) = min(tighest_box(2),p(2));
    tighest_box(3) = max(tighest_box(3),p(3));
    tighest_box(4) = max(tighest_box(4),p(4));
end

% add extension
tighest_box = tighest_box + extention.*[-1,-1,1,1];

% move all axes to left-bottom
for i=1:length(hax)
    if strcmp(get(hax(i),'tag'),'legend')
        continue;
    end
    p = get(hax(i),'position');
    set(hax(i),'position',...
        [p(1)-tighest_box(1),p(2)-tighest_box(2),p(3),p(4)]);
end

% resize figure to fit tightly
set(hfig,'units','centimeters');
p = get(hfig, 'position');

width = tighest_box(3)-tighest_box(1);
height = tighest_box(4)-tighest_box(2); 
set(hfig,'position',[p(1),p(2),width,height]);

% set papersize
set(hfig,'PaperUnits','centimeters');
set(hfig,'PaperSize',[width,height]);
set(hfig,'PaperPositionMode','manual');
set(hfig,'PaperPosition',[0,0,width,height]);

end

