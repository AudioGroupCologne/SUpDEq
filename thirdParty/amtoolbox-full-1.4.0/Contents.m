% AMT - Online documentation
%
%
%  General
%  -------
%
%  This is the most complete, and up-to-date description of the AMT. It consists of auto-generated documentation of all models and other functions included in the AMT. This documentation is directly included in the M-files. Note that because of the automatic generation, the appearance on the website may suffer some details. 
%  
%  The current status of the models can be found in the section "Models" at this website.
%
%  Changes from previous AMT versions ca be found in CHANGES. 
%
%  New to the AMT? 
%  ---------------
%
%    AMT_START            - Start the AMT. Only the third-party toolboxes available on your system will be loaded. Note that the LTFAT is an essential toolbox and still will be downloaded and installed.
%
%    CHANGES              - Changes to previous versions
%
%  Note that the AMT full package provides all third-party toolboxes and pre-compiled binaries. For more details, read the documentation of AMT_START and AMT_MEX.
% 
%   Requirements to run the AMT 1.x
%   -------------------------------
%
%   Matlab 2018b (or more recent) or Octave 8.2 (or more recent). Further, the large time-frequency analysis toolbox (*LTFAT*) is essential. It will be downloaded and installed by amt_start on the fly if not available on your system.
%
%   Further requirements are optional and depend on the model. Type amt_info('X'); with X being the model name to display the model-specific requirements. See also AMT_INFO for more details. 
%
%   The model-specific requirements are as follows:
%
%   1) Many models require compiled binaries. We provide pre-compiled binaries within the AMT full package, however, compilation on your system may still be required. To this end, the AMT needs to access the GNU Compiler Collection (GCC). The availability can be checked by executing system('gcc --version'); in Matlab or Octave. To compile the AMT binaries, execute AMT_MEX. 
%
%   2) The models VERHULST2012, VERHULST2015, and VERHULST2018 require Python (version >= 2.6) with the packages numpy and scipy. On Linux, install the packages in the shell: sudo apt-get install python-scipy python-numpy. On Windows, install Python from <https://www.python.org/> and add python.exe to the Windows search path. Then install the packages by executing python -m pip install numpy scipy in the Windows Command Window. In Matlab/Octave check the Python version with system('python -V').
%
%   3) System-dependent toolboxes may be required for some models. On Matlab, mostly the Signal Processing and the Statistics Toolbox may be required. On Octave, the signal, statistics, netcdf, and optim packages may be required. 
%
%   4) Third-party toolboxes may be required, depending on the model. All third-party toolboxes are provided in the AMT full package. In addition to that, amt_start('install'); downloads and installs all third-party toolboxes. 
%
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/Contents.php





