function output = data_takanen2013(varargin)
%DATA_TAKANEN2013 Data applied in the model by Takanen, Santala and Pulkki
%   Usage: output = data_takanen2013(flag)
%
%   DATA_TAKANEN2013(flag) returns the data specified by the flag. The
%   optional datasets are either matrices or structs containing data
%   required in the different processing steps of the binaural auditory
%   model by Takanen, Santala and Pulkki (2013). 
%
%   The following flags can be specified:
%
%     'cochleardelays'        Frequency-dependent delays of the cochlea 
%                             model by Verhulst et al (2012) are 
%                             compensated for in the TAKANEN2013_PERIPHERY
%                             function. The delays were determined by 
%                             analyzing the the cross-correlation- 
%                             functions between different frequency bands.
%
%     'msolimits'             In TAKANEN2013_MSO the 
%                             ipsilateral input signals are divided by 
%                             scaling values and thereafter limited between
%                             0 and 1. The scaling values were obtained by
%                             computing the average level of the
%                             TAKANEN2013_PERIPHERY output for a pink noise
%                             signal reproduced at 30 dB SPL
%
%     'wbmsomultp'            In TAKANEN2013_WBMSO the energy output is multiplied
%                             so that the energy of the wide-band MSO model is
%                             in a similar level as compared to the energies
%                             of the narrowband MSO and LSO models. The
%                             values with which the energy output is
%                             multiplied were obtained through an iterative
%                             process.
%
%     'lookuptable'           In TAKANEN2013_DIRECTIONMAPPING, the 
%                             outputs of the different MSO and LSO models
%                             are mapped into directions ranging from -20
%                             to 90. The mapping is implemented following
%                             the idea of self-calibration, where the outputs
%                             of the models are compared separately to sets
%                             of reference values computed with HRTF-processed
%                             samples. The sample employed for the MSO and the
%                             LSO models was a 80-ms-long pink noise burst,
%                             whereas the sample for the wide-band MSO model
%                             was a impulse response of a first-order 
%                             Butterworth lowpass filter with a cut-off 
%                             frequency of 500 Hz
%
%     'onsetmultp'            In TAKANEN2013_ONSETENHANCEMENT
%                             the energies of the short-term and of the 
%                             long-term directional cues are scaled to a 
%                             similar level with the help of a set of pre-
%                             computed values, the values which were
%                             obtained by the computing the average levels
%                             of the two energies for a pink noise burst 
%                             convolved with binaural room impulse response
%                             of a concert hall
%
%     'periphenergyaverages'  In TAKANEN2013_FORMBINAURALACTIVITYMAP
%                             the levels of the what cues are 
%                             scaled in order to visualize the evoked 
%                             activations on the binaural activity map.
%                             The values were obtained by computing the 
%                             average levels of the what cue for a pink 
%                             noise burst reproduced at 60 dB SPL.
%
%     'noplot'                Don't plot, only return data. This is the default.
%
%     'plot'                  Plot the data.
%
%   If no flag is given, the function will print the list of valid flags.
%
%   Examples:
%   ---------
%
%   To load and plot the pre-computed frequency-dependent cochlear model delays use:
%
%     data_takanen2013('cochleardelays','plot');
%
%   To load and plot the pre-computed ipsilateral limits for MSO model use:
%
%     data_takanen2013('msolimits','plot');
%
%   To load and plot the pre-computed wide-band MSO energy multiplier use:
%
%     data_takanen2013('wbmsomultp','plot');
%
%   To load and plot the reference values for the MSO and LSO models at different
%   frequencies for directions use:
%
%     data_takanen2013('lookuptable','plot');
%
%   References:
%     M. Takanen, O. Santala, and V. Pulkki. Visualization of functional
%     count-comparison-based binaural auditory model output. Hearing
%     research, 309:147-163, 2014. PMID: 24513586.
%     
%     M. Takanen, O. Santala, and V. Pulkki. Perceptually encoded signals and
%     their assessment. In J. Blauert, editor, The technology of binaural
%     listening. Springer, 2013.
%     
%     S. Verhulst, T. Dau, and C. A. Shera. Nonlinear time-domain cochlear
%     model for transient stimulation and human otoacoustic emission. J.
%     Acoust. Soc. Am., 132(6):3842 - 3848, 2012.
%     
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/data/data_takanen2013.php

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

%   AUTHOR: Marko Takanen, Olli Santala, Ville Pulkki
%
%   COPYRIGHT (C) 2013 Aalto University
%                      School of Electrical Engineering
%                      Department of Signal Processing and Acoustics
%                      Espoo, Finland


definput.import={'amtredofile'};
definput.flags.type={'missingflag','cochleardelays','msolimits','wbmsomultp','lookuptable','onsetmultp','periphenergyaverages'};
definput.flags.plot = {'noplot','plot'};
[flags,keyvals]  = ltfatarghelper({},definput,varargin);


if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;
%% Data employed in the model of periphery
if flags.do_cochleardelays
    %frequency-dependent delays of the cochlea model by Verhulst et al. (2012)
    %are compensated for in the *takanen2013_periphery.m*-function. The delays
    %were determined by analyzing the the cross-correlation-functions between
    %different frequency bands.
    x=amt_load('takanen2013','cochleardelays.mat');
    output = x.velocitydelays;
    if flags.do_plot
        figure;plot(output);xlabel('Frequency band');ylabel('Delay in samples');
        title('Frequency-dependent delays of the cochlea model');
    end
end
%% Data employed in the MSO model
if flags.do_msolimits
    %in *takanen2013_mso.m*-function the ipsilateral input signals are
    %divided by scaling values and thereafter limited between 0 and 1. The
    %scaling values were obtained by computing the average level of the
    %*periphOutput.left* for a pink noise signal reproduced at 30 dB SPL
    x=amt_load('takanen2013','msolimits.mat');
    output = x.limits;
    if flags.do_plot
        figure;plot(output);xlabel('Frequency band');ylabel('Limit value');
        title('Frequency-dependent limits for the ipsilateral input');
    end
end
%% Data employed in the Wide-band MSO model
if flags.do_wbmsomultp
    %in the *takanen2013widebandmso.m*-function the energy output is multiplied
    %so that the energy of the wide-band MSO model is in a similar level as
    %compared to the energies of the narrowband MSO and LSO models. The values
    %with which the energy output is multiplied were obtained through an
    %iterative process.
    x=amt_load('takanen2013','wbmsomultp.mat');
    output = x.multp;
    if flags.do_plot
        figure;plot(output);xlabel('Frequency band');ylabel('Multiplying coefficient');
        title('Frequency-dependent multiplier for the wb-MSO energy');
    end
end 
%% Data employed in the direction mapping
if flags.do_lookuptable
    %in the *takanen2013_directionmapping.m*, the outputs of the different MSO
    %and LSO models are mapped into directions ranging from -20 to 90. The
    %mapping is implemented following the idea of self-calibration, where 
    %theoutputs of the models are compared separately to sets of reference
    %values computed with HRTF-processed samples. The sample employed for 
    %the narrowband MSO and the LSO models was a 80-ms-long pink noise burst,
    %whereas the sample for the wide-band MSO model was a impulse response 
    %of a first-order Butterworth lowpass filter with a cut-off frequency 
    %of 500 Hz
    x=amt_load('takanen2013','lookuptable.mat');
    output = x.referencevalues;
    if flags.do_plot
        figure;plot(output.angles,output.mso(:,8),output.angles,output.lso(:,8),'r',...
            output.angles,output.wbmso(:,8),'k');xlabel('Azimuth angle [degrees]');
        ylabel('Output values');legend('MSO','LSO','Wide-band MSO','location','NorthWest');
        title('Rerefence values for the different models');
    end
end 
%% Data employed in the onset contrast enhancement
if flags.do_onsetmultp
    %in the *takanen2013_onsetenhancement.m*-function, the energies of the
    %short-term and of the long-term directional cues are scaled to a similar
    %level with the help of a set of pre-computed values, the values which were
    %obtained by the computing the average levels of the two energies for a
    %pink noise burst convolved with binaural room impulse response of a
    %concert hall.
    x=amt_load('takanen2013','onsetmultp.mat');
    output = x.coeff;
    if flags.do_plot
        figure;plot(output);xlabel('Frequency band');ylabel('Multiplying coefficient');
        title('Frequency-dependent multiplier for short-term cue energy');
    end
end
%% Data employed both in onset contrast enhancement and in the forming of the binaural activity map
if flags.do_periphenergyaverages
    %in the *takanen2013_formbinauralactivitymap.m*-function the levels of the
    %what cues are scaled in order to visualize the evoked activations on the
    %binaural activity map. The values were obtained by computing the average
    %levels of the what cue for a pink noise burst reproduced at 60 dB SPL.
    x=amt_load('takanen2013','periphenergyaverages.mat');
    output = x.averageEnerg;
    if flags.do_plot
        figure;plot(output);xlabel('Frequency band');ylabel('Average output');
        title('Average outputs of the periphery model for a pink noise at 60 dB SPL');
    end
end

