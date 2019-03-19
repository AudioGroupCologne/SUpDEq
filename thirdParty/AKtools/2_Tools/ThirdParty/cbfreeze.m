function CBH = cbfreeze(varargin)
    
    if verLessThan('matlab', '8.4')
        CBH = cbfreeze2014a_orLess(varargin);
    else
        CBH = cbfreeze2014b_orHigher(varargin);
    end
    
    % Output?:
    if ~nargout
        clear CBH
    else
        CBH(~ishandle(CBH)) = [];
    end

end


function CBH = cbfreeze2014a_orLess(varargin)
%CBFREEZE   Freezes the colormap of a colorbar.
%
%   SYNTAX:
%           cbfreeze
%           cbfreeze('off')
%           cbfreeze(H,...)
%     CBH = cbfreeze(...);
%
%   INPUT:
%     H     - Handles of colorbars to be freezed, or from figures to search
%             for them or from peer axes (see COLORBAR).
%             DEFAULT: gcf (freezes all colorbars from the current figure)
%     'off' - Unfreezes the colorbars, other options are:
%               'on'    Freezes
%               'un'    same as 'off'
%               'del'   Deletes the colormap(s).
%             DEFAULT: 'on' (of course)
%
%   OUTPUT (all optional):
%     CBH - Color bar handle(s).
%
%   DESCRIPTION:
%     MATLAB works with a unique COLORMAP by figure which is a big
%     limitation. Function FREEZECOLORS by John Iversen allows to use
%     different COLORMAPs in a single figure, but it fails freezing the
%     COLORBAR. This program handles this problem.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * If no colorbar is found, one is created.
%     * The new frozen colorbar is an axes object and does not behaves
%       as normally colorbars when resizing the peer axes. Although, some
%       time the normal behavior is not that good.
%     * Besides, it does not have the 'Location' property anymore.
%     * But, it does acts normally: no ZOOM, no PAN, no ROTATE3D and no
%       mouse selectable.
%     * No need to say that CAXIS and COLORMAP must be defined before using
%       this function. Besides, the colorbar location. Anyway, 'off' or
%       'del' may help.
%     * The 'del' functionality may be used whether or not the colorbar(s)
%       is(are) froozen. The peer axes are resized back. Try: 
%        >> colorbar, cbfreeze del
%
%   EXAMPLE:
%     surf(peaks(30))
%     colormap jet
%     cbfreeze
%     colormap gray
%     title('What...?')
%
%   SEE ALSO:
%     COLORMAP, COLORBAR, CAXIS
%     and
%     FREEZECOLORS by John Iversen
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cbfreeze.m
%   VERSION: 1.1 (Sep 02, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jun 08, 2009)
%   1.1      Fixed BUG with image handle on MATLAB R2009a. Thanks to Sergio
%            Muniz. (Sep 02, 2009)

%   DISCLAIMER:
%   cbfreeze.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Parameters:
cbappname = 'Frozen';         % Colorbar application data with fields:
                              % 'Location' from colorbar
                              % 'Position' from peer axes befor colorbar
                              % 'pax'      handle from peer axes.
axappname = 'FrozenColorbar'; % Peer axes application data with frozen
                              % colorbar handle.
 
% Set defaults:
S = 'on';                   Sopt = {'on','un','off','del'};
H = get(0,'CurrentFig');

% Check inputs:
if nargin==2 && (~isempty(varargin{1}) && all(ishandle(varargin{1})) && ...
  isempty(varargin{2}))
 
 % Check for CallBacks functionalities:
 % ------------------------------------
 
 varargin{1} = double(varargin{1});
 
 if strcmp(get(varargin{1},'BeingDelete'),'on') 
  % Working as DeletFcn:

  if (ishandle(get(varargin{1},'Parent')) && ...
      ~strcmpi(get(get(varargin{1},'Parent'),'BeingDeleted'),'on'))
    % The handle input is being deleted so do the colorbar:
    S = 'del'; 
    
   if ~isempty(getappdata(varargin{1},cbappname))
    % The frozen colorbar is being deleted:
    H = varargin{1};
   else
    % The peer axes is being deleted:
    H = ancestor(varargin{1},{'figure','uipanel'}); 
   end
   
  else
   % The figure is getting close:
   return
  end
  
 elseif (gca==varargin{1} && ...
                     gcf==ancestor(varargin{1},{'figure','uipanel'}))
  % Working as ButtonDownFcn:
  
  cbfreezedata = getappdata(varargin{1},cbappname);
  if ~isempty(cbfreezedata) 
   if ishandle(cbfreezedata.ax)
    % Turns the peer axes as current (ignores mouse click-over):
    set(gcf,'CurrentAxes',cbfreezedata.ax);
    return
   end
  else
   % Clears application data:
   rmappdata(varargin{1},cbappname) 
  end
  H = varargin{1};
 end
 
else
 
 % Checks for normal calling:
 % --------------------------
 
 % Looks for H:
 if nargin && ~isempty(varargin{1}) && all(ishandle(varargin{1}))
  H = varargin{1};
  varargin(1) = [];
 end

 % Looks for S:
 if ~isempty(varargin) && (isempty(varargin{1}) || ischar(varargin{1}))
  S = varargin{1};
 end
end

% Checks S:
if isempty(S)
 S = 'on';
end
S = lower(S);
iS = strmatch(S,Sopt);
if isempty(iS)
 error('CVARGAS:cbfreeze:IncorrectStringOption',...
  ['Unrecognized ''' S ''' argument.' ])
else
 S = Sopt{iS};
end

% Looks for CBH:
CBH = cbhandle(H);

if ~strcmp(S,'del') && isempty(CBH)
 % Creates a colorbar and peer axes:
 pax = gca;
 CBH = colorbar('peer',pax);
else
 pax = [];
end


% -------------------------------------------------------------------------
% MAIN 
% -------------------------------------------------------------------------
% Note: only CBH and S are necesary, but I use pax to avoid the use of the
%       "hidden" 'Axes' COLORBAR's property. Why... ??

% Saves current position:
fig = get(  0,'CurrentFigure');
cax = get(fig,'CurrentAxes');

% Works on every colorbar:
for icb = 1:length(CBH)
 
 % Colorbar axes handle:
 h  = double(CBH(icb));
 
 % This application data:
 cbfreezedata = getappdata(h,cbappname);
 
 % Gets peer axes:
 if ~isempty(cbfreezedata)
  pax = cbfreezedata.pax;
  if ~ishandle(pax) % just in case
   rmappdata(h,cbappname)
   continue
  end
 elseif isempty(pax) % not generated
  try
   pax = double(get(h,'Axes'));  % NEW feature in COLORBARs
  catch
   continue
  end
 end
 
 % Choose functionality:
 switch S
  
  case 'del'
   % Deletes:
   if ~isempty(cbfreezedata)
    % Returns axes to previous size:
    oldunits = get(pax,'Units');
    set(pax,'Units','Normalized');
    set(pax,'Position',cbfreezedata.Position)
    set(pax,'Units',oldunits)
    set(pax,'DeleteFcn','')
    if isappdata(pax,axappname)
     rmappdata(pax,axappname)
    end
   end
   if strcmp(get(h,'BeingDelete'),'off') 
    delete(h)
   end
   
  case {'un','off'}
   % Unfrozes:
   if ~isempty(cbfreezedata)
    delete(h);
    set(pax,'DeleteFcn','')
    if isappdata(pax,axappname)
     rmappdata(pax,axappname)
    end
    oldunits = get(pax,'Units');
    set(pax,'Units','Normalized')
    set(pax,'Position',cbfreezedata.Position)
    set(pax,'Units',oldunits)
    CBH(icb) = colorbar(...
     'peer'    ,pax,...
     'Location',cbfreezedata.Location);
   end
 
  otherwise % 'on'
   % Freezes:
 
   % Gets colorbar axes properties:
   cb_prop  = get(h);
   
   % Gets colorbar image handle. Fixed BUG, Sep 2009
   hi = findobj(h,'Type','image');
    
   % Gets image data and transform it in a RGB:
   CData = get(hi,'CData'); 
   if size(CData,3)~=1
    % It's already frozen:
    continue
   end
  
   % Gets image tag:
   Tag = get(hi,'Tag');
  
   % Deletes previous colorbar preserving peer axes position:
   oldunits = get(pax,'Units');
              set(pax,'Units','Normalized')
   Position = get(pax,'Position');
   delete(h)
   cbfreezedata.Position = get(pax,'Position');
              set(pax,'Position',Position)
              set(pax,'Units',oldunits)
  
   % Generates new colorbar axes:
   % NOTE: this is needed because each time COLORMAP or CAXIS is used,
   %       MATLAB generates a new COLORBAR! This eliminates that behaviour
   %       and is the central point on this function.
   h = axes(...
    'Parent'  ,cb_prop.Parent,...
    'Units'   ,'Normalized',...
    'Position',cb_prop.Position...
   );
  
   % Save location for future call:
   cbfreezedata.Location = cb_prop.Location;
  
   % Move ticks because IMAGE draws centered pixels:
   XLim = cb_prop.XLim;
   YLim = cb_prop.YLim;
   if     isempty(cb_prop.XTick)
    % Vertical:
    X = XLim(1) + diff(XLim)/2;
    Y = YLim    + diff(YLim)/(2*length(CData))*[+1 -1];
   else % isempty(YTick)
    % Horizontal:
    Y = YLim(1) + diff(YLim)/2;
    X = XLim    + diff(XLim)/(2*length(CData))*[+1 -1];
   end
  
   % Draws a new RGB image:
   image(X,Y,ind2rgb(CData,colormap),...
    'Parent'            ,h,...
    'HitTest'           ,'off',...
    'Interruptible'     ,'off',...
    'SelectionHighlight','off',...
    'Tag'               ,Tag...
   )  

   % Removes all   '...Mode'   properties:
   cb_fields = fieldnames(cb_prop);
   indmode   = strfind(cb_fields,'Mode');
   for k=1:length(indmode)
    if ~isempty(indmode{k})
     cb_prop = rmfield(cb_prop,cb_fields{k});
    end
   end
   
   % Removes special COLORBARs properties:
   cb_prop = rmfield(cb_prop,{...
    'CurrentPoint','TightInset','BeingDeleted','Type',...       % read-only
    'Title','XLabel','YLabel','ZLabel','Parent','Children',...  % handles
    'UIContextMenu','Location',...                              % colorbars
    'ButtonDownFcn','DeleteFcn',...                             % callbacks
    'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle',...
    'PlotBoxAspectRatio','DataAspectRatio','Position',... 
    'XLim','YLim','ZLim'});
   
   % And now, set new axes properties almost equal to the unfrozen
   % colorbar:
   set(h,cb_prop)

   % CallBack features:
   set(h,...
    'ActivePositionProperty','position',...
    'ButtonDownFcn'         ,@cbfreeze,...  % mhh...
    'DeleteFcn'             ,@cbfreeze)     % again
   set(pax,'DeleteFcn'      ,@cbfreeze)     % and again!  
  
   % Do not zoom or pan or rotate:
   setAllowAxesZoom  (zoom    ,h,false)
   setAllowAxesPan   (pan     ,h,false)
   setAllowAxesRotate(rotate3d,h,false)
   
   % Updates data:
   CBH(icb) = h;   

   % Saves data for future undo:
   cbfreezedata.pax       = pax;
   setappdata(  h,cbappname,cbfreezedata);
   setappdata(pax,axappname,h);
   
 end % switch functionality   

end  % MAIN loop


% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

% Output?:
if ~nargout
 clear CBH
else
 CBH(~ishandle(CBH)) = [];
end

% Returns current axes:
if ishandle(cax) 
 set(fig,'CurrentAxes',cax)
end

end % [EOF]   cbfreeze.m

function CBH = cbfreeze2014b_orHigher(varargin)
%CBFREEZE   Freezes the colormap of a colorbar.
%
%   SYNTAX:
%           cbfreeze
%           cbfreeze off
%           cbfreeze(CMAP)          % Freezes it with this colormap!
%           cbfreeze(CMAP,'off')
%           cbfreeze(H,...)
%     CBH = cbfreeze(...);
%
%   INPUTS:
%     CMAP  - Colormap matrix or name to be used at the colorbar.
%             DEFAULT: (uses the figure one)
%     H     - Handles of colorbars to be frozen, or from figures to search
%             for them or from peer axes (see COLORBAR).
%             DEFAULT: gcf (freezes all colorbars from the current figure)
%     'off' - Unfreezes the colorbars, other options are:
%               'on'    Freezes
%               'un'    same as 'off'
%               'del'   Deletes the colorbars.
%             DEFAULT: 'on' (of course)
%
%   OUTPUTS (all optional):
%     CBH - Color bar handle(s).
%
%   DESCRIPTION:
%     MATLAB works with a unique COLORMAP by figure which is a big
%     limitation. Function FREEZECOLORS by John Iversen allows to use
%     different COLORMAPs in a single figure, but it fails freezing the
%     COLORBAR. This program handles this problem.
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * If no colorbar is found, one is created on each of the axes.
%     * But, if you need to creat it with LCOLORBAR instead of COLORBAR,
%       use:
%        >> lcolorbar(...,'Tag','Colorbar')
%     * The new frozen colorbar is an axes object and does not behaves
%       as normally colorbars when resizing the peer axes. Although, some
%       time the normal behavior is not that good.
%     * Besides, it does not have the 'Location' property anymore.
%     * But, it does acts normally: no ZOOM, no PAN, no ROTATE3D and no
%       mouse selectable.
%     * No need to say that CAXIS and COLORMAP must be defined before using
%       this function. Besides, the colorbar location. Anyway, 'off' or
%       'del' may help.
%     * The 'del' functionality may be used whether or not the colorbar(s)
%       is(are) froozen. The peer axes are resized back. Try:
%        >> colorbar, cbfreeze del
%
%   EXAMPLE:
%     surf(peaks(30))
%     colormap jet
%     cbfreeze
%     colormap gray
%     title('What...?')
%
%   SEE ALSO:
%     COLORMAP, COLORBAR, CAXIS
%     and
%     FREEZECOLORS by John Iversen
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   cbfreeze.m
%   VERSION: 2.1 (Jul 03, 2014) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>)
%   MATLAB:  8.2.0.701 (R2013b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   REVISIONS:
%   1.0      Released. (Jun 08, 2009)
%   1.1      Fixed BUG with image handle on MATLAB R2009a. Thanks to Sergio
%            Muniz. (Sep 02, 2009)
%   2.0      Fixed several BUGs about scientific notation, thanks to Rafa
%            and Jenny from the FileExchange. Changed application name to
%            'cbfreeze' on both, the colorbar and the peer axes. New
%            optional input CMAP. (Jun 05, 2014) 
%   2.1      Fixed BUGs about reading inputs, thanks to Maxime Desbiens. 
%            (Jul 03, 2014) 

%   DISCLAIMER:
%   cbfreeze.m is provided "as is" without warranty of any kind, under the
%   revised BSD license.

%   Copyright (c) 2009-2014 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Parameters:
appName = 'cbfreeze';

% Set defaults:
OPT  = 'on';
H    = get(get(0,'CurrentFigure'),'CurrentAxes');
CMAP = [];

% Checks inputs:
assert(nargin<=3,'CAVARGAS:cbfreeze:IncorrectInputsNumber',...
    'At most 3 inputs are allowed.')
assert(nargout<=1,'CAVARGAS:cbfreeze:IncorrectOutputsNumber',...
    'Only 1 output is allowed.')

% Checks from where CBFREEZE was called:
if (nargin~=2) || (isempty(varargin{1}) || ...
        ~all(reshape(ishandle(varargin{1}),[],1)) ...
        || ~isempty(varargin{2}))
    % CBFREEZE called from Command Window or M-file:

    % Reads H in the first input: Version 2.1
    if ~isempty(varargin) && ~isempty(varargin{1}) && ...
            all(reshape(ishandle(varargin{1}),[],1)) 
        H = varargin{1};
        varargin(1) = [];
    end
    
    % Reads CMAP in the first input: Version 2.1
    if ~isempty(varargin) && ~isempty(varargin{1})
        if isnumeric(varargin{1}) && (size(varargin{1},2)==3) && ...
                (size(varargin{1},1)==numel(varargin{1})/3)
            CMAP = varargin{1};
            varargin(1) = [];
        elseif ischar(varargin{1}) && ...
                (size(varargin{1},2)==numel(varargin{1}))
            temp = figure('Visible','off');
            try
                CMAP = colormap(temp,varargin{1});
            catch
                close temp
                error('CAVARGAS:cbfreeze:IncorrectInput',...
                    'Incorrrect colormap name ''%s''.',varargin{1})
            end
            close temp
            varargin(1) = [];
        end
    end
    
    % Reads options: Version 2.1
    while ~isempty(varargin)
        if isempty(varargin{1}) || ~ischar(varargin{1}) || ...
                (numel(varargin{1})~=size(varargin{1},2))
            varargin(1) = [];
            continue
        end
        switch lower(varargin{1})
            case {'off','of','unfreeze','unfreez','unfree','unfre', ...
                    'unfr','unf','un','u'}
                OPT = 'off';
            case {'delete','delet','dele','del','de','d'}
                OPT = 'delete';
            otherwise
                OPT = 'on';
        end
    end
    
    % Gets colorbar handles or creates them:
    CBH = cbhandle(H,'force');
    
else
    
    % Check for CallBacks functionalities:
    % ------------------------------------
    
    varargin{1} = double(varargin{1});
    
    if strcmp(get(varargin{1},'BeingDelete'),'on')
        % CBFREEZE called from DeletFcn:
        
        if (ishandle(get(varargin{1},'Parent')) && ...
                ~strcmpi(get(get(varargin{1},'Parent'),'BeingDeleted'),'on'))
            % The handle input is being deleted so do the colorbar:
            OPT = 'delete';
            
            if strcmp(get(varargin{1},'Tag'),'Colorbar')
                % The frozen colorbar is being deleted:
                H = varargin{1};
            else
                % The peer axes is being deleted:
                H = ancestor(varargin{1},{'figure','uipanel'});
            end
        else
            % The figure is getting close:
            return
        end
        
    elseif ((gca==varargin{1}) && ...
            (gcf==ancestor(varargin{1},{'figure','uipanel'})))
        % CBFREEZE called from ButtonDownFcn:
        
        cbdata = getappdata(varargin{1},appName);
        if ~isempty(cbdata)
            if ishandle(cbdata.peerHandle)
                % Sets the peer axes as current (ignores mouse click-over):
                set(gcf,'CurrentAxes',cbdata.peerHandle);
                return
            end
        else
            % Clears application data:
            rmappdata(varargin{1},appName)
        end
        H = varargin{1};
    end
    
    % Gets out:
    CBH = cbhandle(H);
    
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Keeps current figure:
cfh = get(0,'CurrentFigure');

% Works on every colorbar:
for icb = 1:length(CBH)
    
    % Colorbar handle:
    cbh = double(CBH(icb));
    
    % This application data:
    cbdata = getappdata(cbh,appName);
    
    % Gets peer axes handle:
    if ~isempty(cbdata)
        peer = cbdata.peerHandle;
        if ~ishandle(peer)
            rmappdata(cbh,appName)
            continue
        end
    else
        % No matters, get them below:
        peer = [];
    end
    
    % Choose functionality:
    switch OPT
        
        case 'delete'
            % Deletes:
            if ~isempty(peer)
                % Returns axes to previous size:
                oldunits = get(peer,'Units');
                set(peer,'Units','Normalized');
                set(peer,'Position',cbdata.peerPosition)
                set(peer,'Units',oldunits)
                set(peer,'DeleteFcn','')
                if isappdata(peer,appName)
                    rmappdata(peer,appName)
                end
            end
            if strcmp(get(cbh,'BeingDelete'),'off')
                delete(cbh)
            end
            
        case 'off'
            % Unfrozes:
            if ~isempty(peer)
                delete(cbh);
                set(peer,'DeleteFcn','')
                if isappdata(peer,appName)
                    rmappdata(peer,appName)
                end
                oldunits = get(peer,'Units');
                set(peer,'Units','Normalized')
                set(peer,'Position',cbdata.peerPosition)
                set(peer,'Units',oldunits)
                CBH(icb) = colorbar(...
                    'Peer'    ,peer,...
                    'Location',cbdata.cbLocation);
            end
            
        case 'on'
            % Freezes:

            % Gets colorbar axes properties:
            cbprops = get(cbh);
            
            % Current axes on colorbar figure:
            fig = ancestor(cbh,{'figure','uipanel'});
            cah = get(fig,'CurrentAxes');
            
            % Gets colorbar image handle. Fixed BUG, Sep 2009
            himage = findobj(cbh,'Type','image');
            
            % Gets image data and transforms them to RGB:
            CData = get(himage,'CData');
            if size(CData,3)~=1
                % It's already frozen:
                continue
            end
            
            % Gets image tag:
            imageTag = get(himage,'Tag');
            
            % Deletes previous colorbar preserving peer axes position:
            if isempty(peer)
                peer = cbhandle(cbh,'peer');
            end
            oldunits = get(peer,'Units');
            set(peer,'Units','Normalized')
            position = get(peer,'Position');
            delete(cbh)
            oldposition = get(peer,'Position');
            
            % Seves axes position
            cbdata.peerPosition = oldposition;
            set(peer,'Position',position)
            set(peer,'Units',oldunits)
            
            % Generates a new colorbar axes:
            % NOTE: this is needed because each time COLORMAP or CAXIS
            %       is used, MATLAB generates a new COLORBAR! This
            %       eliminates that behaviour and is the central point
            %       on this function.
            cbh = axes(...
                'Parent'  ,cbprops.Parent,...
                'Units'   ,'Normalized',...
                'Position',cbprops.Position...
                );
            
            % Saves location for future calls:
            cbdata.cbLocation = cbprops.Location;
            
            % Move ticks because IMAGE draws centered pixels:
            XLim = cbprops.XLim;
            YLim = cbprops.YLim;
            if     isempty(cbprops.XTick)
                % Vertical:
                X = XLim(1) + diff(XLim)/2;
                Y = YLim    + diff(YLim)/(2*length(CData))*[+1 -1];
            else % isempty(YTick)
                % Horizontal:
                Y = YLim(1) + diff(YLim)/2;
                X = XLim    + diff(XLim)/(2*length(CData))*[+1 -1];
            end
            
            % Gets colormap:
            if isempty(CMAP)
                cmap = colormap(fig);
            else
                cmap = CMAP;
            end
            
            % Draws a new RGB image:
            image(X,Y,ind2rgb(CData,cmap),...
                'Parent'            ,cbh,...
                'HitTest'           ,'off',...
                'Interruptible'     ,'off',...
                'SelectionHighlight','off',...
                'Tag'               ,imageTag)
            
            % Moves all '...Mode' properties at the end of the structure,
            % so they won't become 'manual':
            % Bug found by Rafa at the FEx. Thanks!, which also solves the
            % bug found by Jenny at the FEx too. Version 2.0
            cbfields = fieldnames(cbprops);
            indmode  = strfind(cbfields,'Mode');
            temp     = repmat({'' []},length(indmode),1);
            cont     = 0;
            for k = 1:length(indmode)
                % Removes the '...Mode' properties:
                if ~isempty(indmode{k})
                    cont = cont+1;
                    temp{cont,1} = cbfields{k};
                    temp{cont,2} = getfield(cbprops,cbfields{k});
                    cbprops = rmfield(cbprops,cbfields{k});
                end
            end
            for k=1:cont
                % Now adds them at the end:
                cbprops = setfield(cbprops,temp{k,1},temp{k,2});
            end
            
            % Removes special COLORBARs properties:
            cbprops = rmfield(cbprops,{...
                'CurrentPoint','TightInset','BeingDeleted','Type',...       % read-only
                'Title','XLabel','YLabel','ZLabel','Parent','Children',...  % handles
                'UIContextMenu','Location',...                              % colorbars
                'ButtonDownFcn','DeleteFcn',...                             % callbacks
                'CameraPosition','CameraTarget','CameraUpVector', ...
                'CameraViewAngle',...
                'PlotBoxAspectRatio','DataAspectRatio','Position',...
                'XLim','YLim','ZLim'});
            
            % And now, set new axes properties almost equal to the unfrozen
            % colorbar:
            set(cbh,cbprops)
            
            % CallBack features:
            set(cbh,...
                'ActivePositionProperty','position',...
                'ButtonDownFcn'         ,@cbfreeze,...  % mhh...
                'DeleteFcn'             ,@cbfreeze)     % again
            set(peer,'DeleteFcn'        ,@cbfreeze)     % and again!
            
            % Do not zoom or pan or rotate:
            %if isAllowAxesZoom(fig,cbh)
            setAllowAxesZoom  (    zoom(fig),cbh,false)
            %end
            %if isAllowAxesPan(fig,cbh)
            setAllowAxesPan   (     pan(fig),cbh,false)
            %end
            %if isAllowAxesRotate(fig,cbh)
            setAllowAxesRotate(rotate3d(fig),cbh,false)
            %end
            
            % Updates data:
            CBH(icb) = cbh;
            
            % Saves data for future undo:
            cbdata.peerHandle = peer;
            cbdata.cbHandle   = cbh;
            setappdata(cbh ,appName,cbdata);
            setappdata(peer,appName,cbdata);
            
            % Returns current axes:
            if ishandle(cah)
                set(fig,'CurrentAxes',cah)
            end
            
    end % switch functionality
    
end  % MAIN loop

% Resets the current figure
if ishandle(cfh)
    set(0,'CurrentFigure',cfh)
end

% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

% Output?:
if ~nargout
    clear CBH
else
    CBH(~ishandle(CBH)) = [];
end

end

% [EOF] CBFREEZE.M by Carlos A. Vargas A.
