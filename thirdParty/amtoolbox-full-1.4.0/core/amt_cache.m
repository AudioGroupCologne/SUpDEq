function varargout=amt_cache(cmd,name,varargin)
%AMT_CACHE  Caches variables for later or retrieves variables from cache
%   Usage: var = amt_cache('get', package, flags);
%          amt_cache('set', package, variables);
%   
%   Input parameters:
%
%     get      : gets the content of a package from the cache.
%                variables = AMT_CACHE('get', package) reads a package from the cache
%                and outputs its content in var. package must a be a string
%                identifying the package of variables. If the package contains multiple
%                variables, variables can be a list of variables like 
%                [var1, var2, ... , varN] = .... The order of returned variables is the
%                same as that used for saving in cache.
%                ... = AMT_CACHE('get', package, flags) allows to control the 
%                behaviour of accessing the cache. 
%
%     set      : stores variables as a package in the cache. 
%                AMT_CACHE('set', package, variables) saves variables in the cache using
%                the name package. variables can be a list of variables separated by
%                comma.
%                
%     getURL   : outputs the URL of the cache in the internet. 
%
%     setURL   : sets the URL of the internet cache to a new URL. 
%
%     clearAll : clears the cache directory. An interactive confirmation is required.
%
%
%   AMT_CACHE is a utility function for avoiding long computation times by retrieving
%   pre-calculated intermediate results from the AMT's online repository.
%
%   The following flags can be specified
%
%     'normal'    Use cached package.
%                 If the cached package is locally not available, it will be downloaded from the internet. 
%                 If it is remotely not available, enforce recalculation of the package. 
%                 Note that this method may by-pass the actual processing and thus does not always test the
%                 actual functionality of a model.
%
%     'cached'    Enforce to used cached package.
%                 If the cached package is locally not available, it will be 
%                 downloaded from the internet. 
%                 If it is remotely not available, enforce recalculation of the package. 
%                 Note that this method may by-pass the actual processing and thus 
%                 does not always test the actual functionality of a model. 
%                 It is, however, very convenient for fast access of results 
%                 like plotting figures. On the internet, the cached packages 
%                 are available for the release versions only.
%
%     'redo'    Enforce the recalculation of the package.
%               [..] = AMT_CACHE('get', [..]) outputs empty variables always. 
%
%     'localonly'    Package will be recalculated when locally not available.
%                    Do not connect to the internet.
% 
%
%
%   This is an example of using the cache in a function. 
%   In this example, we store the variables x, y, and z in the package
%   fileABC*:
%
%     definput.import={'amt_cache'}; % load default flags from arg_amt_cache
%     [flags,~]  = ltfatarghelper({},definput,varargin); % overwrite flags if user-provided
%     [x,y,z] = amt_cache('get', 'fileABC', flags);
%     if isempty(x)
%           % calculate your variables x,y,z here
%         amt_cache('set', 'fileABC', x, y, z); % save calculated variables in fileABC
%     end
%     %  use your variables x, y, z here
%
%   Note that in this example, the flags indicating the mode of caching are
%   stored in flags.cachemode which can be achieved by:
% 
%     definput.import={'amt_cache'};
%     [flags, keyvals] = ltfatarghelper({}, definput, varargin); 
%
%   at the begin of the function. This way, the cache mode can be provided by the 
%   user if required. 
%
%
%   See also: data_ziegelwanger2013
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/core/amt_cache.php


%   #Author: Piotr Majdak (2015)
%   #Author: Clara Hollomey (2021)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

%%   `amt_cache` supports the following commands:

[~, kv] = amt_configuration;
mlock;
persistent CacheURL CacheMode;
caller = dbstack;
last = 2;

if isempty(CacheURL)
  CacheURL = kv.cacheURL;
  if ( numel(caller)==1 ) || ( numel(caller) > 1 && ~strcmp('amt_configuration', caller(last).name) )
    amt_configuration('cacheURL', CacheURL);
  end  
end

if isempty(CacheMode)
  CacheMode='normal'; 
  if ( numel(caller)==1 ) || ( numel(caller) > 1 && ~strcmp('amt_configuration', caller(last).name) )
    amt_configuration('cachemode', CacheMode);
  end  
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
        if ~exist(tokenfn,'file')
          webfn=[CacheURL '/' urlencode(token) '/' name '.mat'];
          amt_disp(['Cache: Downloading ' name '.mat for ' token]);
          if ~exist(tokenpath,'dir'); mkdir(tokenpath); end
          if isoctave
            [~, stat]=urlwrite(webfn, tokenfn);
          else
            try
              outfilename = websave(tokenfn,webfn);
              stat = 1;
            catch
              stat = 0;
            end
          end

          if ~stat     
          for ii = 2:numel(kv.version)
            webfn=[kv.downloadURL kv.version{ii} '/cache/' urlencode(token) '/' name '.mat'];
            if isoctave
              [~,stat]=urlwrite(webfn,tokenfn);
            else
              try
                outfilename = websave(tokenfn,webfn);
                stat = 1;
              catch
                stat = 0;
              end
            end  

            if stat
              disp(['Found data in version: ' kv.version{ii}]);
              break;
            end
          
            if ii == numel(kv.version)
              error(['Unable to download file: ' webfn]);
            end
          end
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
            amt_disp(['Cached data not found: ' tokenfn]);
            amt_disp('Enforce recalculation...');
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
          amt_disp(['Cache: Downloading ' name '.mat for ' token]);
          if ~exist(tokenpath,'dir'); mkdir(tokenpath); end
          if isoctave
            [~, stat]=urlwrite(webfn, tokenfn);
          else
            try
              outfilename = websave(tokenfn,webfn);
              stat = 1;
            catch
              stat = 0;
            end
          end

          if ~stat     
            amt_disp('Data not found. Searching in previous AMT versions...');
          for ii = 2:numel(kv.version)
            webfn=[kv.downloadURL kv.version{ii} '/cache/' urlencode(token) '/' name '.mat'];
            if isoctave
              [~,stat]=urlwrite(webfn,tokenfn);
            else
              try
                outfilename = websave(tokenfn,webfn);
                stat = 1;
              catch
                stat = 0;
              end
            end  

            if stat
              amt_disp(['Found data in version: ' kv.version{ii}]);
              %load(tokenfn);  % downloaded to local cache. Load...
              %for ii=1:nargout
              %varargout{ii}=cache(ii).value;
              %end
              break;
            end
          
            if ii == numel(kv.version)
                amt_disp(['Cached data not found: ' webfn]);
                amt_disp('Enforce recalculation...');
            %else         
            %    load(tokenfn);  % downloaded to local cache. Load...
            %    for ii=1:nargout
            %    varargout{ii}=cache(ii).value;
            %    end
            end 
          end
          end
          %end
          if stat
            load(tokenfn);  % downloaded to local cache. Load...
            for ii=1:nargout
              varargout{ii}=cache(ii).value;
            end
          else
            for ii=1:nargout
              varargout{ii}=[];
            end    
          end
          %if ~stat
          %  amt_disp(['Cached data not found: ' webfn]);
          %  amt_disp('Enforce recalculation...');
          %  for ii=1:nargout, varargout{ii}=[]; end % enforce recalculation
          %else         
          %  load(tokenfn);  % downloaded to local cache. Load...
          %  for ii=1:nargout
          %    varargout{ii}=cache(ii).value;
          %  end
          %end 
        %else
        %  load(tokenfn);  % Locally available, load...
        %  for ii=1:nargout
        %    varargout{ii}=cache(ii).value;
        %  end
       % end
        
        %end
    else
        load(tokenfn);  % Locally available, load...
        for ii=1:nargout
           varargout{ii}=cache(ii).value;
        end
    end    
    end
  case 'setURL'
    CacheURL=name;
    
    if ( numel(caller)==1 ) || ( numel(caller) > 1 && ~strcmp('amt_configuration', caller(last).name) )
      amt_configuration('cacheURL', name);
    end
    
  case 'getURL'
    varargout{1}=CacheURL;
  case 'clearAll'
    cachepath=fullfile(amt_basepath,'cache');
    if strcmp(input(['clearAll clears ' strrep(cachepath,'\','\\') '. Type YES for confirmation: '],'s'),'YES'), 
      amt_disp(['Clearing ' cachepath ' ...']);
      rmdir(cachepath, 's');       
    end
  case 'setMode'
    CacheMode = name;
    %if exist(flags.cachemode) && ~strcmp(flags.cachemode, name)
    if ( numel(caller)==1 ) || ( numel(caller) > 1 && ~strcmp('amt_configuration', caller(last).name) )
      amt_configuration('cachemode', name);
    end

  otherwise
    error('Unsupported command');
    end
end




