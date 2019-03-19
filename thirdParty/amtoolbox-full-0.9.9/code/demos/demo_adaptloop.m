%DEMO_ADAPTLOOP  Show the effect of adaptation
%
%   This script demonstrates the effect of adaptation applied to a test
%   signal with and without noise.
%
%   The test signal is made of a sinosoidal ramp up and down between 0
%   and 1.
%
%   Figure 1: Clean test signal
%
%      This figure shows the effect of adaptation on the clean test signal with and
%      without overshoot limiting.
%
%   Figure 2: Noisy test signal
%
%      This figure shows the effect of adaptation on the noisy test signal
%      with and without overshoot limiting. Notice that in the second plot,
%      the initial spike at the beginning of the signal caused from the sharp
%      transition from complete silence to noise is magnitudes larger than
%      the values in the rest of the output.
%
%   See also: adaptloop
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/demos/demo_adaptloop.php

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


siglen=10000;
fs=10000;

% This is the default minimum level (0 dB) of the adaptation loops. The
% loops assume that a signal is never silent, and sets all values below
% minlvl equal to minlvl. For plotting purposes, we do the same explicitly.
minlvl=setdbspl(0);

part=siglen/10;

insig=[zeros(2*part,1);
       rampup(part);
       ones(2*part,1);
       rampdown(part);
       zeros(4*part,1)];

insig=max(insig,minlvl);

figure;

x=(0:siglen-1)/fs;
subplot(3,1,1);
plot(x,20*log10(insig));
title('Input signal');
xlabel('time / s');
ylabel('level / Db');

subplot(3,1,2);
plot(x,adaptloop(insig,fs,0));
title('Adaptation.');
xlabel('time / s');
ylabel('level / model units');

subplot(3,1,3);
plot(x,adaptloop(insig,fs));
title('Adaptation w. limiting.');
ylabel('level / model units');
xlabel('time / s');

% Add a low level of noise
insig=abs(insig+0.001*randn(siglen,1));
insig=max(insig,minlvl);

figure;

subplot(3,1,1);
plot(x,20*log10(insig));
title('Input signal with added Gaussian noise.');
ylabel('level / Db');
xlabel('time / s');

subplot(3,1,2);
plot(x,adaptloop(insig,fs,0));
title('Adaptation.');
ylabel('level / model units');
xlabel('time / s');

subplot(3,1,3);
plot(x,adaptloop(insig,fs));
title('Adaptation w. limiting.');
ylabel('level / model units');
xlabel('time / s');
