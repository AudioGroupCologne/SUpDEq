function [k,v]=may2011_findlocalpeaks(x,m,w)
%MAY2011_FINDLOCALPEAKS finds peaks with optional quadratic interpolation
%
%   Input parameters:  
%     X       : is the input signal
%     M       : mode ('q' performs quadratic interpolation, 'v' finds
%               valleys instead of peaks)
%     W       : is the width tolerance; a peak will be eliminated if there is
%               a higher peak within +-W samples
%
%   Output parameters:  
%     K       : are the peak locations in X (fractional if M='q')
%     V       : are the peak amplitudes: if M='q' the amplitudes will be interpolated
%               whereas if M~='q' then V=X(K). 
%
%   Outputs are column vectors regardless of whether X is row or column.
%   If there is a plateau rather than a sharp peak, the routine will place the
%   peak in the centre of the plateau. When the W input argument is specified,
%   the routine will eliminate the lower of any pair of peaks whose separation
%   is <=W; if the peaks have exactly the same height, the second one will be eliminated.
%   All peak locations satisfy 1<K<length(X).
%
%   If no output arguments are specified, the results will be plotted.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/may2011_findlocalpeaks.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: MATLAB M-Signal
%   #Author: Mike Brookes (2005)
%   #Author: Tobias May (2014)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

if nargin<2
    m=' ';
end
nx=length(x);
if any(m=='v')
    x=-x(:);        % invert x if searching for valleys
else
    x=x(:);        % force to be a column vector
end
dx=x(2:end)-x(1:end-1);
r=find(dx>0);
f=find(dx<0);

if length(r)>0 & length(f)>0    % we must have at least one rise and one fall
    dr=r;
    dr(2:end)=r(2:end)-r(1:end-1);
    rc=repmat(1,nx,1);
    rc(r+1)=1-dr;
    rc(1)=0;
    rs=cumsum(rc); % = time since the last rise
    
    df=f;
    df(2:end)=f(2:end)-f(1:end-1);
    fc=repmat(1,nx,1);
    fc(f+1)=1-df;
    fc(1)=0;
    fs=cumsum(fc); % = time since the last fall
    
    rp=repmat(-1,nx,1);
    rp([1; r+1])=[dr-1; nx-r(end)-1];
    rq=cumsum(rp);  % = time to the next rise
    
    fp=repmat(-1,nx,1);
    fp([1; f+1])=[df-1; nx-f(end)-1];
    fq=cumsum(fp); % = time to the next fall
    
    k=find((rs<fs) & (fq<rq) & (floor((fq-rs)/2)==0));   % the final term centres peaks within a plateau
    v=x(k);
    
    if any(m=='q')         % do quadratic interpolation
        b=0.5*(x(k+1)-x(k-1));
        a=x(k)-b-x(k-1);
        j=(a>0);            % j=0 on a plateau
        v(j)=x(k(j))+0.25*b(j).^2./a(j);
        k(j)=k(j)+0.5*b(j)./a(j);
        k(~j)=k(~j)+(fq(k(~j))-rs(k(~j)))/2;    % add 0.5 to k if plateau has an even width
    end
    
    % now purge nearby peaks
    
    if nargin>2
        j=find(k(2:end)-k(1:end-1)<=w);
        while any(j)
            j=j+(v(j)>=v(j+1));
            k(j)=[];
            v(j)=[];
            j=find(k(2:end)-k(1:end-1)<=w);
        end
    end
else
    k=[];
    v=[];
end
if any(m=='v')
    v=-v;    % invert peaks if searching for valleys
end
if ~nargout
    if any(m=='v')
        x=-x;    % re-invert x if searching for valleys
        ch='v';
    else
        ch='^';
    end
    plot(1:nx,x,'-',k,v,ch);
end


