function [spOut, vOut] = ashida2016_LSOmodelCOC(spEx, spIn, DT, varargin)
%ASHIDA2016_LSOMODELCOC Coincidence counting model of lateral superior olive 
%
%   Input parameters:
%     spEx : excitatory input spike vector (number of inputs at each time step)
%     spIn : inhibitory input spike vector (number of inputs at each time step)
%     DT   : size of time step [ms] 
%
%   Output parameters:
%     spOut : number of output spikes at each time step 
%     vOut  : internal state (virtual membrane potential) at each time step 
%
% Notes on the input/output arguments 
%  - spEx and spIn should have the same length. Otherwise, the shorter 
%    of them determines the total simulation time length. 
%
% Notes on the model 
%  The model was originally created by Ashida et al. (2016). This model 
%  counts the number of synchronized excitatory inputs in the coincidence 
%  window and generates a spike if the count reaches the threshold. 
%  Effects of inhibition is modeled as a (weighted) subtraction of 
%  excitatory inputs in the inhibition time window. For a comparison 
%  with other LSO models, see Ashida et al. (2017). 
%
% References 
% Ashida G, Kretzberg J, Tollin DJ (2016) PLoS Comput Biol 12: e1004997
%  "Roles for coincidence detection in coding amplitude-modulated sounds" 
% Ashida G, Tollin DJ, Kretzberg J (2017) PLoS Comput Biol 13(12): e1005903. 
% "Physiological models of the lateral superior olive" 
% https://doi.org/10.1371/journal.pcbi.1005903
%
% Revisions
% #Author: Go Ashida (ver 0.9, May 16, 2017)
% #Author: Go Ashida (ver 1.0, Dec 28, 2017) : updated references
% #Author: Jonas Klug (June 2, 2019) : Edit to call from a framework
% #Author: Alejandro Osses (2023) : Integration in AMT 1.4
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/ashida2016_LSOmodelCOC.php


%%%%%% Copyright 2017 Go Ashida (go.ashida@uni-oldenburg.de) %%%%%%%%%%%%%
% Permission is hereby granted under the Apache License, Version 2.0; 
% Users of this file must be in compliance with this license, a copy of 
% which may be obtained at http://www.apache.org/licenses/LICENSE-2.0
% This file is provided on an "AS IS" basis, WITHOUT WARRANTIES OR 
% CONDITIONS OF ANY KIND, either express or implied.

definput.import = {'ashida2016_LSOmodelCOC'};
[flags,kv]  = ltfatarghelper({},definput,varargin);
lso = kv.lso;

%% pre-defined parameters
Tref = lso.Tref;    %1.6; % [ms] 
ThEx = lso.ThEx;    %8;   % threshold 
WiEx = lso.WiEx;    %0.8; % [ms] coincidence window
gIn = lso.gIn;    %2;   % threshold increase by inhibition 
WiIn = lso.WiIn;    %1.6; % [ms] inhibition window 

Nref = round(Tref/DT); % steps
NwEx = round(WiEx/DT); % steps
NwIn = round(WiIn/DT); % steps

%% data vectors
Nsteps = min([length(spEx), length(spIn)]);
spOut = zeros(1,Nsteps,'int8');
vOut = zeros(1,Nsteps,'single');
spSum = 0; % input spike counter
refCount = 0; % refractory counter

%% main loop 
for t = 1:Nsteps

  % add excitatory inputs at time t 
  spSum = spSum + spEx(t); 
  % subtract inhibitory inputs at time t  
  spSum = spSum - spIn(t) * gIn; 
  % remove just-expired excitatory inputs 
  if(t-NwEx>0); spSum = spSum - spEx(t-NwEx); end
  % remove just-expired inhibitory inputs 
  if(t-NwIn>0); spSum = spSum + spIn(t-NwIn) * gIn ; end

  % check for threshold crossing 
  % (1) if in refractory period, then decrement counter 
  if(refCount>0) 
   spOut(t)=0; refCount = refCount-1;  
  % (2) if threshold is reached, generate a spike and set refractory counter
  elseif(spSum>=ThEx) 
   spOut(t) = 1; refCount = Nref; 
  % (3) if no threshold crossing happened, then no spike output 
  else 
   spOut(t)=0; 
  end 

  % save internal state at time t 
  vOut(t) = spSum;

end
