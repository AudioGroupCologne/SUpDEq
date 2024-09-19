V1.2.1: February 26, 2024
-------------------------
- Ported to GIT (without history)
- Fixed bug in AKp when doing balloon plots
- Added support for playrec and new Mac Silicon processors

V1.2.0: March 8, 2018, revision 115
-----------------------------------
- Added
   + AKroomSimulationDemo.m, and 2_Tools/RoomSimulation/
   + AKairAttenuation.m
   + AKsphericalHeadDemo.m, AKsphericalHead.m, and AKsphericalHead.pdf
   + AKsphericalHeadDiffuse.m
   + AKshEqualization.m
   + AKshRadial.m
   + AKshEnergy.m
   + AKsoftLimit.m
   + AKairAttenuation
   + AKfractionalOctaves.m
   + AKfrequencyTicks.m
   + AKcollectFiles.m
   + AKtightenFigure.m
   + AKlistIntegers.m
   + More demo data in 0_DemoData/
- AKmeasureIR.m can now measure absolute sound levels after calibration
- AKmeasureDemo.m has a structure that is easier to follow.
- AKmeasureIR.m can now measure channels in any order and non-successive channles, i.r. [3 1 5:8]
- AKfilter.m is more efficient and the structure was improved
- AKp.m can now use spherical spline interpolation for displaying spherical data
- AKspeedOfSound.m can now account for humidity, pressure and CO_2 concentration

V1.1.0: February 9, 2017, revision 29
-------------------------------------
- The ‘FABIAN head-related transfer function data base’ is now available from: https://dx.doi.org/10.14279/depositonce-5718.2
- Added
   + AKio, and AKioDemo
   + playrec external for Windows (64 bit only)
   + AKsofa2ak
- Minor bug fixes and improvements in the documentation and robustness

V1.0.1: December 19, 2016, revision 15
--------------------------------------
- Added
   + AKboxplot
   + AKboxplotDemo
   + AKaverage
   + AKtoa
   + AKtoaDemo
   + AKcrossCorr
   + AKaverage
   + AKfractionalDelay
   + AKkaiser
   + AKdisp
   + getRatings.m
   + some small statistics tools
- Minor bug fixes and improvements in the documentation and robustness

----------------------------------
V1.0: November 8, 2016, revision 0
- Initial public release
