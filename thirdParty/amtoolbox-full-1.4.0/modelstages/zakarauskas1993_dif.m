function [ out ] = zakarauskas1993_dif( in,do_it )
%ZAKARAUSKAS1993_DIF differentiates input signal
%
%   Usage:
%     [ out ] = zakarauskas1993_dif( in,do_it )
%
%   Input parameters:
%     in:       magnitude of spectrum
%     do_it:       differential order; default: 0
%
%   Output parameters:
%     out:      derivative of input signal
%
%   This function differentiates the input signal.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/zakarauskas1993_dif.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author: Robert Baumgartner (2010), OEAW Acoustical Research Institute

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

% default settings
if ~exist('do_it','var')
    do_it=0;
end

if length(size(in))==1
    switch do_it
        case {1}
            x1=zeros(length(in)-1,1);
            for ind=1:size(x1,1)
                x1(ind) = in(ind+1) - in(ind);
            end
            out=x1;
        case {2}
            x1=zeros(length(in)-1,1);
            for ind=1:size(x1,1)
                x1(ind) = in(ind+1) - in(ind);
            end
            x2=zeros(length(x1)-1,1);
            for ind=1:size(x2,1)
                x2(ind) = x1(ind+1) - x1(ind);
            end
            out=x2;
        otherwise
            out=in;
    end
elseif length(size(in))==2
    switch do_it
        case {1}
            x1=zeros(size(in,1)-1,size(in,2));
            for ii=1:size(in,2)
                for ind=1:size(x1,1)
                    x1(ind,ii) = in(ind+1,ii) - in(ind,ii);
                end
            end
            out=x1;
        case {2}
            x1=zeros(size(in,1)-1,size(in,2));
            for ii=1:size(in,2)
                for ind=1:size(x1,1)
                    x1(ind,ii) = in(ind+1,ii) - in(ind,ii);
                end
            end
            x2=zeros(size(x1,1)-1,size(in,2));
            for ii=1:size(in,2)
                for ind=1:size(x2,1)
                    x2(ind,ii) = x1(ind+1,ii) - x1(ind,ii);
                end
            end
            out=x2;
        otherwise
            out=in;
    end
else
    disp('Please insert each channel seperately!');
end
end



