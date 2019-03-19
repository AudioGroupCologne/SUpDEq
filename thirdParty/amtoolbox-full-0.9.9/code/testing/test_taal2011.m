f=greasy;
f=[greasy;greasy];

L=numel(f);

x=f;
y=setdbspl(randn(L,1),dbspl(f)-20);

d=taal2011(x,y,16000);
d_ref=ref_stoi(x,y,16000);
d_ref_1=ref_stoi_1(x,y,16000);

% Exact reference
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/testing/test_taal2011.php

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
(d_ref-d_ref_1)/d_ref

% Modified with respect to the reference
(d-d_ref)/d_ref

% Test for multi-signal
x=randn(20000,2);
y=x+0.1*randn(20000,2);

d_ref(1)= ref_stoi_1(x(:,1),y(:,1),10000);
d_ref(2)= ref_stoi_1(x(:,2),y(:,2),10000);
d       = taal2011(x,y,10000);

(d-d_ref)./d_ref
