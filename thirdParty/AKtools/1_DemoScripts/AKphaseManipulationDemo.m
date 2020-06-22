% Demo script for manipulating the phase response of impulse responses.
% This might be usefull from time to time...
%
% 09/2013 - fabian.brinkmann@tu-berlin.de

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under  the License. 
close all; clear; clc

% load HRIR for manipulating it's phase
x = AKhrirInterpolation(0, 0, 0, 'measured_sh');

% try with odd and even number of samples and see how linear phase IR
% changes...
Ntrunc = 256;    % 256, or 255
x = x(1:Ntrunc);

% set the group delay for linear phase generation
% (The value of (Ntrunc-1)/2 means that the peak of the linear phase IR
%  will be in it's middle - play with this)
Ngroup = (Ntrunc-1)/2;

% get IRs with manipulated phase behaviour
y_min  = AKphaseManipulation(x, 44100, 'min', 0);      % the last argument zero-pads x before generating the minimum phase
                                                       % to achieve better results. You should increase this until differences
                                                       % between minimum-phase and original magnitude spectrum are small enough.
                                                       % There is output on the command window that helps you with this.
                                                       % In this case the differences are 0.07 dB in the important
                                                       % range - which looks ok to me
                                                            
y_lin  = AKphaseManipulation(x, 44100, 'lin', Ngroup); % the last parameter sets the group delay in samples

y_zero = AKphaseManipulation(x, 44100, 'zero');

% plot
AKf(30,20)
subplot(2,3,1)
    AKp(x, 't2d', 'c', [.7 .7 .7])
    AKp(y_min, 'tc2d', 'dr', [-1.5 2.5])
    title('Minimum phase')
    legend('orignial', 'minimum phase', 'location', 'NorthEast')
subplot(2,3,2)
    AKp(x, 't2d', 'c', [.7 .7 .7])
    AKp(y_lin, 'tc2d', 'dr', [-1.5 2.5])
    title('Linear phase')
    legend('orignial', 'linear phase', 'location', 'NorthEast')
subplot(2,3,3)
    AKp(x, 't2d', 'c', [.7 .7 .7])
    AKp(y_zero, 'tc2d', 'dr', [-1.5 2.5])
    title('Zero phase')
    legend('orignial', 'zero phase', 'location', 'NorthEast')
subplot(2,1,2)
    AKp(x, 'm2d', 'c', [.7 .7 .7])
    AKp(y_min, 'm2d', 'dr', [-20 20])
    AKp(y_lin, 'm2d', 'dr', [-20 20])
    AKp(y_zero, 'm2d', 'dr', [-20 20])
    title('Magnitude responses')
