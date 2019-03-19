function varargout=amtcache(cmd,name,varargin)
S=warning;
warning('on');
warning('Warning: AMTCACHE will be removed in a future release. Use AMT_CACHE instead. ');  
warning(S);

%   Based on amt_cache from 9.7.2017.
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/legacy/amtcache.php

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

persistent CacheURL CacheMode;
if isempty(CacheURL)
  CacheURL=['http://www.sofacoustics.org/data/amt-' amt_version('version') '/cache'];
end
if isempty(CacheMode)
  CacheMode='normal';
end

switch cmd
  case 'set'
    f=dbstack('-completenames');
    fn=f(2).file;
    token=urlencode(strrep(fn(length(amt_basepath)+1:end),'\','/'));
    tokenpath=fullfile(amt_basepath,'cache',token);
    tokenfn=fullfile(tokenpath,[name '.mat']);
    
    if ~exist(tokenpath,'dir'); mkdir(tokenpath); end
    for ii=3:nargin
      cache(ii-2).name=inputname(ii);
      cache(ii-2).value=varargin{ii-2};
    end
    save(tokenfn,'cache','-v6');
    varargout{1}=tokenfn;
    
  case 'get'
    if nargin<3, varargin{1}='global'; end  % if not provided: global
    if strcmp(varargin{1},'global'), varargin{1}=CacheMode; end  % if global: use the stored cache mode
      % now let's parse the cache mode
    switch varargin{1}
      case 'redo' % force recalculations in any case
        for ii=1:nargout, varargout{ii}=[]; end

      case 'cached' % use local cache. If not available download from the internet. If not available throw an error.
        f=dbstack('-completenames');
        fn=f(2).file;
        token=urlencode(strrep(fn(length(amt_basepath)+1:end),'\','/'));
        tokenpath=fullfile(amt_basepath,'cache',token);
        tokenfn=fullfile(tokenpath,[name '.mat']);
        if ~exist(tokenfn,'file'),
          webfn=[CacheURL '/' urlencode(token) '/' name '.mat'];
          amt_disp(['Cache: Downloading ' name '.mat for ' token],'progress');
          if ~exist(tokenpath,'dir'); mkdir(tokenpath); end
          [~,stat]=urlwrite(webfn,tokenfn);
          if ~stat
            error(['Unable to download file from remote cache: ' webfn]);
          end          
        end
        load(tokenfn);
        for ii=1:nargout
          varargout{ii}=cache(ii).value;
        end
        
      case 'localonly' % use local cache only. If not available, enforce recalculation
        f=dbstack('-completenames');
        fn=f(2).file;
        token=urlencode(strrep(fn(length(amt_basepath)+1:end),'\','/'));
        tokenpath=fullfile(amt_basepath,'cache',token);
        tokenfn=fullfile(tokenpath,[name '.mat']);
        if ~exist(tokenfn,'file'),
            amt_disp(['Cached data not found: ' tokenfn],'progress');
            amt_disp('Enforce recalculation...','progress');
            for ii=1:nargout, varargout{ii}=[]; end % enforce recalculation
        else
          load(tokenfn);
          for ii=1:nargout
            varargout{ii}=cache(ii).value;
          end
        end

      case 'normal' % use local cache. If not available download from the internet. If not available recalculate.
        f=dbstack('-completenames');
        fn=f(2).file;
        token=urlencode(strrep(fn(length(amt_basepath)+1:end),'\','/'));
        tokenpath=fullfile(amt_basepath,'cache',token);
        tokenfn=fullfile(tokenpath,[name '.mat']);
        if ~exist(tokenfn,'file'),
          webfn=[CacheURL '/' urlencode(token) '/' name '.mat'];
          amt_disp(['Cache: Downloading ' name '.mat for ' token],'progress');
          if ~exist(tokenpath,'dir'); mkdir(tokenpath); end
          [~,stat]=urlwrite(webfn,tokenfn);
          if ~stat
            amt_disp(['Cached data not found: ' webfn],'progress');
            amt_disp('Enforce recalculation...','progress');
            for ii=1:nargout, varargout{ii}=[]; end % enforce recalculation
          else         
            load(tokenfn);  % downloaded to local cache. Load...
            for ii=1:nargout
              varargout{ii}=cache(ii).value;
            end
          end 
        else
          load(tokenfn);  % Locally available, load...
          for ii=1:nargout
            varargout{ii}=cache(ii).value;
          end
        end
        
    end
  case 'setURL'
    CacheURL=name;
  case 'getURL'
    varargout{1}=CacheURL;
  case 'clearAll'
    cachepath=fullfile(amt_basepath,'cache');
    if strcmp(input(['clearAll clears ' strrep(cachepath,'\','\\') '. Type YES for confirmation: '],'s'),'YES'), 
      amt_disp(['Clearing ' cachepath ' ...'],'progress');
      rmdir(cachepath, 's');       
    end
  case 'setMode'
    CacheMode = name;
  otherwise
    error('Unsupported command');
end
