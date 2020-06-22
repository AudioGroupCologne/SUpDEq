%MAXIMIZE  Maximize a figure window to fill the entire screen
%
% Examples:
%   maximize
%   maximize(hFig)
%
% Maximizes the current or input figure so that it fills the whole of the
% screen that the figure is currently on. This function is platform
% independent.
%
%IN:
%   hFig - Handle of figure to maximize. Default: gcf.

function maximize(hFig)
if nargin < 1
    hFig = gcf();
end
% Java-less method (since R2018a)
if exist('verLessThan', 'file') && ~verLessThan('matlab', '9.4')
    set(hFig, 'WindowState', 'maximized');
    return;
end
% Old Java version
drawnow(); % Required to avoid Java errors
set(hFig, 'WindowStyle', 'normal');
% Suppress warning about Java obsoletion
oldState = warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFig = get(handle(hFig), 'JavaFrame'); 
jFig.setMaximized(true);
warning(oldState);
end