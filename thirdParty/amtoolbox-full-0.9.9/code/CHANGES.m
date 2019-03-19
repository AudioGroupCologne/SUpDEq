%Changes throughout the release history of the AMT
%   
%   Version 0.9.9 (September 2017)
%   ==============================
% 
%   - data_baumgartner2017 provided
%   - data_baumgartner2017looming provided)
%   - kelavas2015 improved (appears in the code and documentation now)
%   - bug #78 fixed: hohmann2002_process did not work for multiband input signals.
%
%   Version 0.9.8 (12.7.2017)
%   =========================
% 
%   Compatibility breaks:
%     - amt --> amt_*
%         - Use amt_start to start the AMT
%         - Legacy files provided but will be removed in the future. 
%         - amthelp --> amt_version
%     - gfb --> hohmann2002. Legacy files provided. See note #11 for more details.
%     - zilany2007humanized --> zilany2007. Legacy file provided.
%     - drnl --> lopezpoveda2001. Legacy file provided.
%     - modfilterbankepsm --> ewert2000. No legacy file. 
%     - HRTFs: only SOFA files allowed now.
%     - Functions creating/modifying signals have the prefix sig_ now:
%         - Legacy files provided but will be removed in the future.
%        	- irns --> sig_yost1996
%        	- whitenoiseburst --> sig_whitenoiseburst
%        	- transposedtone --> sig_transposedtone
%       	- perfectsweep --> sig_linsweep
%        	- notchednoise --> sig_notchednoise
%        	- bmsin --> sig_lindemann1986
%        	- simulatedimpulseresponse --> sig_joergensen2011
%       	- itdsin --> sig_itdsin
%       	- ildsin --> sig_ildsin
%       	- itdildsin --> sig_itdildsin
%       	- competingtalkers -> sig_competingtalkers
%       	- breebaart2001siggen --> sig_breebaart2001
%       	- bincorrnoise --> sig_bincorrnoise
%       	- bandpassnoisefreq --> sig_bandpassnoise
%     - Model stages in modelstages have the format modelXX_stageYY now:
%         - modelXXstageYY --> modelXX_stageYY
%         - ffGn --> zilany2014_ffGn
%         - breebaart2001preproc --> breebaart2001_preproc
%         - No legacy files provided. Adapt your code if directly calling model stages from your code. 
% 
%   New
%     - baumgartner2017: sound externalization model
%     - baumgartner2016: level-dependent sagittal-plane sound localization model for NH and HI listeners
%     - hohmann2002: the gfb_ framework of Gammatone filterbank integrated as hohmann2002 framework
%     - kelvasa2015: sound localization in cochlear-implant listeners
%     - emuexp: emulation of experiments using interative runs like 3-AFC
%     - breebaart2001_centralproc: decision stage from Breebaart et al. (2001). 
%     - exp_breebaart2001: reproduces results from Breebaart et al. (2001) based on emuexp
%     - model initiative: interface to the model initiative (Dietz et al. 2016)
% 
% 
%   Structural changes:
%     - baumgartner2014 decomposed into model stages
%     - exp_spille2013 merged into exp_dietz2011
%     - directories re-structured:
%         - main scripts of a model go to: model. They are called nameyear
%         - scripts with model stages (model scripts other than the main one) go to: modelstages. They are called nameyear_postfix
%         - scripts not being part of a specific model go to: general. They do not (!) start with model name.
%         - scripts generating audio signals go to: signals. They are called sig_nameyear or sig_functionality
%         - scripts returning measured data go to: data. They are called data_nameyear
%         - scripts with default parameters of other scripts go to: defaults. They are called arg_callingfunction
%         - files being compiled to mex files go to: mex. In order to be compiled, they must have the prefix comp_
%         - files being compiled to oct files go to: oct
%         - other files requiring compilation/installation go to: bin. Their compilation must be considered in make.bat (Windows) and Makefile (Linux/MacOS).
%         - HRTFs go to hrtfs. They are in SOFA and will be downloaded on the fly. 
%         - Lagecy files providing backwards compatibility fo to: legacy
%         - binaural, monaural, speech, filters removed
% 
%   Other updates:
%     - dietz2011: minor bug fixes
%     - data_joergensen2011: completion of data
%     - exp_baumgartner2014: new figures, compatibility improved
%     - demo_hohmann2002: new figures
%     - installation simplified (amt_mex does all the installation now)
% 
%   Version 0.9.7 (10.6.2015)
%   =========================
% 
%   New:
%     - Caching of data: see amt_cache and the cache directory.
%     - Automatic download of auxiliary data. See amt_load and the auxdata directory.
%     - Control for messages output in the command line. See amt_disp.
%     - SOFA files for HRTFs: requires SOFA API, see SOFAload.
%     - Models: zilany2014, joergensen2011, joergensen2013, georganti2013
%     - Signals: sig_joergensen2011
%     - Clean documentation (no errors, no warnings).
%     - makefile for Linux, compiling of cpp files
% 
%   Structure changes:
%     - arg_ functions moved to arg directory, comp_ functions moved to mex, directory comp removed
%     - plot functions renamed to plot_ and moved to plot directory
%     - amt_start and amt_mex improved
%     - readme file for sourceforge added
%     - reference directory removed (it was a directory with original contributions to the AMT)
% 
%   Other changes:
%     - interpolation for various polar-angle samplings  
%     - added new experiment in exp_baumgartner2014: fig5_baumgartner2015aro
%     - stability improvements in baumgartner2014
%     - minor bugfix in demo_baumgartner2013 and doc update in baumgartner2014
%     - exp_lindemann1986: fig 14b disabled: it takes ages and is wrong anyway...
%     - added reference for verhulst2012
%     - wierstorf2013: additional files for HRTF handling removed, load the itd-to-angle look-up table with data_wierstorf2013.m now. 
%     - hrtf/enzner2008 removed (enzner2008 data are in auxdata now)
%     - langendijk2002: data and HRTFs removed from repository (are required data)
%     - changed the order of announcements on amt_start
%     - jelfs2011: removed dependency on read_hrir
%     - plotjelfs2011 moved to demo_jelfs2011.m (plotjelfs2011 was actually a demo).
%     - progress output supressed in the documentation
%     - amt_disp introduced for displaying information depending on the start-up condition of the AMT.
%     - enzner2008 and exp_enzner2008 split in the model and experiment part.
%     - exp_georganti2013 works for me. Documentation is missing yet.
%     - may2011 documentation integrated
%     - 2014 version of may2011 added. demo_may2011 works but documentation invalid yet.
%     - extractsp: stability improvement
%     - minor documentation and stability updates, new function baumgartner2014parametrization and functionalities in localizationerror added.
%     - Fixed imag ILD in dietz2011
%     - Added function to load some simulated monaural room impulses responses.
%     - documentation updates and use of SOFA's remote load functionality in data_baumgartner2014.
%     - major style overhaul of the Joergsen 2011 and 2013 models. Experiments included etc. Does not yet pass mat2doc, and sound files are missing.
%     - localizationerror: new performance measures added
%     - data_majdak2010 and ...2013ctc: Angles forced to be real valued.
% 
% 
%   Version 0.9.6
%   =============
% 
%   New: 
%     - Gammatone validation provided, including exp_gammatone, demo_gammatone, exp_hohmann2002, and demo_hohmann2002
%     - Wierstorf et al. (2013) provided, including wierstorf2013, exp_wierstorf2013, and integration with the SFS toolbox
% 
%   Fixes:
%     - data_goode1994: more details provided
%     - jelfs2011 works now with SOFA HRTFs stored in hrtf/jelfs2011/, e.g., kemar.sofa
%     - hohmann2007 naming resolved. hohmann2007 renamed to herzke2007, the primary model is called hohmann2002 now
%     - localizationerror: missing error types added
%     - exp_spille2013: uses lowpass f_inst
%     - dietz2011 improved to better reflect the corresponding publications
%
%   Url: http://amtoolbox.sourceforge.net/amt-0.9.9/doc/CHANGES.php

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

