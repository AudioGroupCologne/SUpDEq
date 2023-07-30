function [timeout,meout,mocr,c1filterout,c2filterout,c1vihc,c2vihc,vihc,synout,psth]  =  smalt2014(pin,CF,nrep,binwidth,T,cohc,cihc,spont,mocr_max,mocr_threshold,mocr_slope,mocr_binauralratio,shocks)
%SMALT2014 Medial olivocochlear reflex in auditory nerve responses
%   Usage: [timeout,meout,mocr,c1filterout,c2filterout,c1vihc,c2vihc,vihc,synout,psth] = smalt2014(pin,CF,nrep,binwidth,T,cohc,cihc,spont,mocr_max,mocr_threshold,mocr_slope,mocr_binauralratio,shocks)
%
%
%   Input parameters:
%     pin                : the input sound wave in Pa sampled at the appropriate sampling rate (see instructions below)
%     cf                 : the characteristic frequency of the fiber in Hz
%     nrep               : the number of repetitions for the psth
%     binwidth           : the binsize in seconds, i.e., the reciprocal of the sampling rate (see instructions below)
%     reptime            : the time between stimulus repetitions in seconds - NOTE should be equal to or longer than the duration of pin
%     cohc               : the ohc scaling factor: 1 is normal OHC function; 0 is complete OHC dysfunction
%     cihc               : the ihc scaling factor: 1 is normal IHC function; 0 is complete IHC dysfunction
%     spont              : the spontaneous rate of the fiber in spikes/s - NOTE a value of 50 was used in Zilany and Bruce (2006)
%     mocr_max           : the maximum ohc gain reduction possible: 1 ohc gain can be completely reduced; 0 means ohc gain cannot be changed by efferent path
%     mocr_threshold     : the threshold of the LD block before gain reduction starts to occur
%     mocr_slope         : the slope of the LD block (how quickly gain reduction starts to occur as level increases)
%     mocr_binauralratio : the ipsi / contra efferent gain reduction ratio
%     shocks             : a length 2 vector with value 0 or 1 for ipsi/contra.  0 indicates do not shock the system, 1 indicates shock the system (effectively reduces gain by mocr_max
%
%
%   Output parameters:
%     timeout     : an array of times in seconds
%     meout       : the output of the middle-ear filter
%     mocr        : the output of the efferent block of the model (MOCR timecourse)
%     c1filterout : the output of the C1 (signal path) BM filter
%     c2filterout : the output of the C2 (parallel path) BM filter
%     c1vihc      : the output of the C1 IHC transduction function
%     c2vihc      : the output of the C2 IHC transduction function
%     vihc        : the IHC potential
%     synout      : the synapse output in spikes/s
%     psth        : the peri-stimulus time histogram
%
%
%   [timeout,meout,mocr,c1filterout,c2filterout,c1vihc,c2vihc,vihc,synout,psth] = smalt2014(pin,1e3,10,1/500e3,0.200,1,1,50,1,-130.5,0.01638,[.4654 .2332],[0 0]);
%   
%   models a normal fiber of spontaneous rate 50 spikes/sec (normal OHC & IHC function)
%   with a CF of 1 kHz, for 10 repititions and a sampling rate of 500kHz, 
%   for a repetition duration of 200 ms and a binaural stimulus pin = [2xN]
%   for a monaural stimulus, use pin [1xN] (row vector)
%   
%   Note on the sampling rate: In this version of the code, only 100kHz sampling rate is supported
%   
%   References:
%     C. J. Smalt, M. G. Heinz, and E. A. Strickland. Modeling the
%     Time-Varying and Level-Dependent Effects of the Medial Olivocochlear
%     Reflex in Auditory Nerve Responses. Journal of the Association for
%     Research in Otolaryngology, 15(2):159--173, Apr. 2014. [1]http ]
%     
%     References
%     
%     1. https://doi.org/10.1007/s10162-013-0430-z
%     
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/models/smalt2014.php


%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Verified
%   #Author: Christopher Smalt (2013)
%   #Author: Michael Heinz (2013)
%   #Author: Elisabeth Strickland (2013)
%   #Author: Clara Hollomey (2023): Integration in the AMT
%   #Author: Piotr Majdak (2023): Documentation fixes

S = size(pin);
if S(1)==2 && S(2) > 1
    % Binaural
    [timeout,meout,mocr,c1filterout,c2filterout,c1vihc,c2vihc,vihc,synout,psth] ...
        =  comp_smalt2014_binaural(pin,CF,nrep,binwidth,T*2,cohc,cihc,spont,mocr_max,mocr_threshold,mocr_slope,mocr_binauralratio,shocks);
        
else
    % Monaural
    [timeout,meout,mocr,c1filterout,c2filterout,c1vihc,c2vihc,vihc,synout,psth] ...
        = comp_smalt2014_monaural(pin,CF,nrep,binwidth,T,cohc,cihc,spont,mocr_max(1),mocr_threshold,mocr_slope,mocr_binauralratio,shocks(1));
end
