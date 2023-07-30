function definput=arg_laback2023(definput)
% ARG_LABACK2023
%
% Gets the default values for the experiment by B. Laback (2023)
%       
% See also: exp_laback2023
%
%   #Author: Bernhard Laback (2023)
%   #Author: Clara Hollomey (2023)
%   #Author: Alejandro Osses (2023): Integration into AMT 1.4
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/defaults/arg_laback2023.php


%Abbreviations
%   T Target stimulus
%   P Precursor stimulus

definput.flags.MOC = {'no_MOC','MOC'}; % new flag by AO

definput.keyvals.stimdb = 60;        %Stimulus intensity (in dB SPL): [default: 60] (MidLev= 60, LowLev=45, HighLev=75)
definput.keyvals.model = 1;          %if 1: MOC model; if 2: Zilany2014 model
definput.keyvals.Experiment = 2;     %Index of Experiment [1, 2, or 3], as in Laback, B. (2020). Peripheral and Central Factors in Contextual Lateralization based on Interaural Level Difference, Proceedings of Forum Acusticum, 2020, Lyon.
definput.keyvals.F0cue = 0;          %Only relevant for Exp. 2: if 0: T and P have same envelope rate, if 1: T has double envelope rate as P
definput.keyvals.ANmode = 3;         %Only used for Zilany2014 model; if 1: LSO based on Zilany2014; if 2; Zilany2014-AN-spike based; 3: Zilany2014-AN-synapse based;
definput.keyvals.SubExp = 1;         %Only relevant for Exp. 2: if 0: Exp 2A (Target varied); if 1: Exp2B (Precursor varied)
definput.keyvals.AdaptMode = 0;      %if 0: standard mode; if 1: a Pre-Precursor with a duration equal to the precrusor (600 ms) and gap of 500 ms is added before the precursor

% definput.keyvals.modelreps = 5;      %number of repetitions of each random noise stimulus
definput.keyvals.PrecILD = 10;       %ILD (in dB) of P (has to be positive; a condition with negative ILD is automatically added) 
definput.keyvals.PrecPrecILD = 0;    %ILD (in dB) of Pre-P (a positive value produces an Ipsi ILD)
definput.keyvals.PrecITD = 0;        %ITD (in degrees) of P (has to be positive; a condition with negative ILD is automatically added) [default: 2000]
% definput.keyvals.step=2;             %Stepsize of ILDs of T [default: 2]
% definput.keyvals.NrILD=6;            %Number of ILDs of T, starting from zero) [default: 6]
definput.keyvals.frozennoise = 0;    %if 1: use frozen noise; if 0: newly generated noises for different P condition
definput.keyvals.Onset = 0;

definput.keyvals.length_stimTarg = [];
definput.keyvals.noiseILDvec = []; % noise ILD vector from sig_laback2023.m

definput.keyvals.DurT  = 301e-3;     % T duration (in sec) [default: 301e-3]
definput.keyvals.DurN  = 600e-3;     % P duration (in sec) [default: 600e-3]
definput.keyvals.DurG  = 10e-3;      % Duration of gap between P and T (in sec) [default: 10e-3]
definput.keyvals.RampDur = 50;       % Duration of raised-cosine ramps for stimuli (in ms) [default: 50]
definput.keyvals.factor = 1;         % Amplification factor for transposed stimuli => not needed any more

definput.keyvals.precmode = 0;       % Precursor mode
definput.keyvals.AdaptMode = 0;
definput.keyvals.Modeldelay = 240;   % Delay between stimulus and model response (in samples) [default: 240]
definput.keyvals.Onset = [];         % Defined in the groups 'exp1' and 'exp2', below (integrate better)
definput.keyvals.ANmode = 3;

%AN model parameters
definput.keyvals.fs = 100e3;             %AN sampling rate in Hz (must be 100, 200 or 500 kHz)
definput.keyvals.nrep = 20;              %Number of stimulus repetitions of AN model [not relevant with synapse implementation]
definput.keyvals.shocks = [0 0];         %Set to 1 to enable shocking MOCR ([ipsi contra])
definput.keyvals.psthbinwidth = 0.5e-3;  % binwidth in seconds; %default: 0.5e-3
definput.keyvals.gapvec=zeros(definput.keyvals.fs*definput.keyvals.DurG,1);
definput.keyvals.CF = 4e3;                %Center frequency of P and T stimuli (in Hz)
definput.keyvals.CF_Model = definput.keyvals.CF;                 %Charecteristic frequency (in Hz) of AN neurons to be predicted [default: 4e3]  
definput.keyvals.spontvec = [0.8 10 130];% spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects stimulus parameters (default: [0.8 10 130])
definput.keyvals.cohc  = 1.0;            %OHC gain factor; if 1: normal OHC function; if 0: completely lost OHC function 
definput.keyvals.cihc  = 1.0;            %IHC gain factor; if 1: normal IHC function; if 0: completely lost IHC function 
definput.keyvals.BinMOC = 1;             %if 1: binaural MOC, with the amount determined by mocr_max_def; if 0, contralateral MOC, but no ipsilateral MOC
definput.keyvals.mocr_max_def = 1.0;     %Factor for MOC-induced OHC gain reduction; if 1: full reduction; if 0: no gain reduction (MOC activation)

% Experiment smalt=========================================================
%General parameters
stimdb = 60;        %Stimulus intensity (in dB SPL): [default: 60] (MidLev= 60, LowLev=45, HighLev=75)
model = 1;          %if 1: MOC model; if 2: Zilany2014 model
Experiment = 2;     %Index of Experiment [1, 2, or 3], as in Laback, B. (2020). Peripheral and Central Factors in Contextual Lateralization based on Interaural Level Difference, Proceedings of Forum Acusticum, 2020, Lyon.
F0cue = 0;          %Only relevant for Exp. 2: if 0: T and P have same envelope rate, if 1: T has double envelope rate as P

% modelreps = 5;      %number of repetitions of each random noise stimulus
PrecILD = 10;       %ILD (in dB) of P (has to be positive; a condition with negative ILD is automatically added) 
PrecITD = 0;        %ITD (in degrees) of P (has to be positive; a condition with negative ILD is automatically added) [default: 2000]
NrILD=6;            %Number of ILDs of T, starting from zero) [default: 6]
frozennoise = 0;    %if 1: use frozen noise; if 0: newly generated noises for different P condition

DurT  = 301e-3;     %T duration (in sec) [default: 301e-3]
DurN  = 600e-3;     %P duration (in sec) [default: 600e-3]
DurG  = 10e-3;      %Duration of gap between P and T (in sec) [default: 10e-3]
RampDur = 50;       %Duration of raised-cosine ramps for stimuli (in ms) [default: 50]
factor = 1;         %Amplification factor for transposed stimuli => not needed any more


%AN model parameters
fs = 100e3;             %AN sampling rate in Hz (must be 100, 200 or 500 kHz)
nrep = 20;              %Number of stimulus repetitions of AN model [not relevant with synapse implementation]
shocks = [0 0];         %Set to 1 to enable shocking MOCR ([ipsi contra])
psthbinwidth = 0.5e-3;  % binwidth in seconds; %default: 0.5e-3
gapvec=zeros(fs*DurG,1);
CF = 4e3;                %Center frequency of P and T stimuli (in Hz)
CF_Model = CF;                 %Charecteristic frequency (in Hz) of AN neurons to be predicted [default: 4e3]  
spontvec = [0.8 10 130];% spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects stimulus parameters (default: [0.8 10 130])
cohc  = 1.0;            %OHC gain factor; if 1: normal OHC function; if 0: completely lost OHC function 
cihc  = 1.0;            %IHC gain factor; if 1: normal IHC function; if 0: completely lost IHC function 
BinMOC = 1;             %if 1: binaural MOC, with the amount determined by mocr_max_def; if 0, contralateral MOC, but no ipsilateral MOC
mocr_max_def = 1.0;     %Factor for MOC-induced OHC gain reduction; if 1: full reduction; if 0: no gain reduction (MOC activation)

definput.keyvals.fiberType = [];

definput.groups.exp1  = {'stimdb', stimdb, 'model', model, 'Experiment', Experiment,...
    'F0cue', F0cue, 'PrecILD', PrecILD, 'PrecITD', PrecITD,...
    'frozennoise', frozennoise, 'DurT', DurT, 'DurN', DurN,...
    'DurG', DurG, 'RampDur', RampDur, 'factor', factor,'fs', fs, 'nrep', nrep,...
    'shocks', shocks, 'psthbinwidth', psthbinwidth, 'gapvec', gapvec, 'CF', CF, ...
    'CF_Model', CF_Model, 'spontvec', spontvec, 'cohc', cohc, 'cihc', cihc, 'BinMOC',...
    BinMOC, 'mocr_max_def', mocr_max_def };

                                                               
%Experiment zilany=========================================================
model = 2;          %if 1: MOC model by Smalt (2014) => should no be used; if 2: AN model by Zilany2014 
                    %in case of model 2, the two MOC options correspond to OHC gain = 0 and 1, respectively
ANmode = 3;         %Only used for Zilany2014 model; if 1: LSO based on Zilany2014; if 2; Zilany2014-AN-spike based; 3: Zilany2014-AN-synapse based;
Experiment = 1;     %Index of Experiment [1, 2, or 3], as in Laback, B. (2020). Peripheral and Central Factors in Contextual Lateralization based on Interaural Level Difference, Proceedings of Forum Acusticum, 2020, Lyon.
SubExp = 1;         %Only relevant for Exp. 2: if 0: Exp 2A (Target varied); if 1: Exp2B (Precursor varied)
F0cue = 1;          %Only relevant for Exp. 2: if 0: 0T and P have same envelope rate, if 1: T has double envelope rate as P

% modelreps = 5;      %number of repetitions of each random noise stimlus (Synapse: 5; LSO: 30)
% NrILD = 6;          %Number of ILDs of T, starting from zero) [default: 6]
% step=2;             %Stepsize of ILDs of T [default: 2]


PCond = 1;          %four different stimulus conditions (1: no T, 2:diotic P, 3: ipsi P, 4: contra P)  
ILDStartInd = 1;    %start index for ILD increment; if 2: ILD corresponds to "step"
MOCStartInd = 1;    %if 2: MOC on, if 2: startin with MOC off
                    %t stimulus conditions (1: no T, 2:diotic P, 3: ipsi P, 4: contra P)   
cohc = 1;           %OHC gain factor; if 1: normal OHC function; if 0: completely lost OHC function 
FibInd = 1;         %Starting index for running different fiber types
nrepExc =20;        %Number of stimulus repetitions at ipsilateral (excitatory) ear [Default: 20] 
nrepInh = 8;      	%Number of stimulus repetitions at contralateral (inhibitory) ear [Default: 8]
AdaptMode = 0;      %if 0: standard mode; if 1: a Pre-Precursor with a duration equal to the precrusor (600 ms) and gap of 500 ms is added before the precursor
Onset = 0;          %Onset of response window re start of target
frozennoise = 0;    %if 0: newly generated noises for different P condition; if 1: use frozen noise 
PrecILD = 10;       %ILD (in dB) of lateral P (has to be positive; a condition with negative ILD is automatically added) 
PrecPrecILD = 0;    %ILD (in dB) of Pre-P (a positive value produces an Ipsi ILD)
PrecITD = 0;        %ITD (in us) of P (has to be positive; a condition with negative ILD is automatically added) [default: 2000]
stimdb = 60;        %Stimulus intensity (in dB SPL): [default: 60] (MidLev= 60, LowLev=45, HighLev=75)
DurT  = 301e-3;     %T duration (in sec) [default: 301e-3]
DurN  = 600e-3;     %P duration (in sec) [default: 600e-3] and Prec-precursor
DurG  = 10e-3;      %Duration of gap between P and T (in sec) [default: 10e-3]
DurGLong  = 500e-3; %Duration (in sec) of gap between Pre-precursor and precursor [default: 500e-3]
RampDur = 50;       %Duration of raised-cosine ramps for stimuli (in ms) [default: 50]
factor = 1;         %amplification factor for transposed stimuli => not needed any more
%==========================================================================

definput.groups.exp2  = {   'model', model, ...
                            'ANmode', ANmode, ...
                            'Experiment', Experiment,...
                            'SubExp', SubExp, ...
                            'F0cue', F0cue, ...
                            'cohc', cohc, ...
                            'AdaptMode', AdaptMode, ...
                            'frozennoise', frozennoise,...
                            'PrecILD', PrecILD, ...
                            'PrecPrecILD', PrecPrecILD, ...
                            'PrecITD', PrecITD, ...
                            'stimdb', stimdb, ...
                            'DurT', DurT, 'DurN', DurN, 'DurG', DurG, ...
                            'RampDur', RampDur, 'factor',factor, ...
                            'fs', fs, ...
                            'nrep', nrep, ...
                            'shocks', shocks, ...
                            'psthbinwidth', psthbinwidth,...
                            'gapvec', gapvec, ...
                            'CF', CF, ...
                            'CF_Model', CF_Model, ...
                            'spontvec', spontvec, ...
                            'cihc', cihc, ...
                            'BinMOC', BinMOC, ...
                            'mocr_max_def', mocr_max_def, ...
                            'Onset', Onset};
