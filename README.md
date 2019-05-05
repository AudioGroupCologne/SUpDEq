## SUpDEq - Spatial Upsampling by Directional Equalization
This toolbox is a MATLAB implementation of the SUpDEq method, as presented by Pörschmann et al. (2019) [1]. In general, the SUpDEq method is an approach to generate dense HRTF datasets based on sparsely measured HRTF datasets. For example, applying the SUpDEq method, a (technically) appropriate full-spherical dense HRTF dataset with 2702 directions can be derived from only 38 actually measured HRTFs. 

Basically, the method attempts to remove direction-dependent temporal and spectral components of the sparse HRTF dataset by a spectral division (equalization) of the HRTFs with corresponding rigid sphere transfer functions (STFs). This equalized HRTF dataset is then transformed to the spherical harmonics (SH) domain by means of a spherical Fourier transform [2][3]. Next, a de-equalization is performed by extracting the equalized HRTFs at the spatial sampling points of the dense grid (SH-interpolation / spatial upsampling) and multiplying these HRTFs with the respective STFs in frequency domain.

## Requirements
The toolbox is implemented in MATLAB R2015b and R2018a and requires the Signal Processing Toolbox. Older versions of MATLAB might also work. The following third party toolboxes are required for full functionality. All listed toolboxes are part of the SUpDEq toolbox and are stored in the folder "thirdParty".

- SOFiA (Sound Field Analysis Toolbox) [4]  
 https://github.com/AudioGroupCologne/SOFiA
- AKtools [5]   
https://www.ak.tu-berlin.de/menue/digitale_ressourcen/research_tools/aktools/
- SOFA API (Spatially Oriented Format for Acoustics) [6]  
https://sourceforge.net/projects/sofacoustics/
- AMT (Auditory Modeling Toolbox) [7]  
http://amtoolbox.sourceforge.net  

SOFiA and AMT partly work with MEX files. The required MEX files for Mac and Windows (64 Bit) are pre-compiled and stored in the respective folders. 


## Installation
Simply clone or zip download this repository and run the "supdeq_start" script. This will add all required paths and install SOFA, if not previously installed on your system. 

## First steps
After installation, we recommend to go through the "supdeq_demo" script step by step. This will lead you through the basic processing. The resulting dense HRTF dataset can then be perceptually and technically compared to the (dense) reference HRTF dataset with the provided evaluation functions.

## Documentation
A HTML documentation of the MATLAB code can be found here: http://audiogroup.web.th-koeln.de/SUpDEq_doc/index.html  
The third-party toolboxes are excluded from the documentation.

## About
Christoph Pörschmann and Johannes M. Arend^a  
TH Köln - University of Applied Sciences  
Institute of Communications Engineering  
Department of Acoustics and Audio Signal Processing  
Betzdorfer Str. 2, D-50679 Cologne, Germany  
https://www.th-koeln.de/akustik  


^a Also at: Technical University of Berlin  
Audio Communication Group  
Einsteinufer 17c, D-10587 Berlin, Germany  

Thanks to Fabian Brinkmann (Audio Communication Group, TU Berlin) for useful discussions and contributions!  

The SUpDeq method is patent protected.



## References
[1] C. Pörschmann*, J.M. Arend*, and F. Brinkmann, “Directional Equalization of Sparse Head-Related Transfer Function Sets for Spatial Upsampling,” IEEE/ACM Trans. Audio, Speech, Lang. Process., vol. 27, no. 6, pp. 1060–1071, 2019, *These authors contributed equally to this work. 
[2] B. Bernschütz, “Microphone Arrays and Sound Field Decomposition for Dynamic Binaural Recording,” TU Berlin, 2016.  
[3] B. Rafaely, Fundamentals of Spherical Array Processing. Berlin Heidelberg: Springer-Verlag, 2015.  
[4] B. Bernschütz, C. Pörschmann, S. Spors, and S. Weinzierl, “SOFiA Sound Field Analysis Toolbox,” in Proceedings of the International Conference on Spatial Audio - ICSA 2011, 2011, pp. 8–16.  
[5] F. Brinkmann and S. Weinzierl, “AKtools - an open software toolbox for signal acquisition, processing, and inspection in acoustics,” in Proceedings of the 142nd AES Convention, Berlin, Germany, 2017, pp. 1–6.  
[6] P. Majdak, Y. Iwaya, T. Carpentier, R. Nicol, M. Parmentier, A. Roginska, Y. Suzuki, K. Watanabe, H. Wierstorf, H. Ziegelwanger, and M. Noisternig, “Spatially Oriented Format for Acoustics: A Data Exchange Format Representing Head-Related Transfer Functions,” in Proceedings of the 134th AES Convention, Rome, Italy, 2013, pp. 1–11.  
[7] P. Søndergaard and P. Majdak, "The Auditory Modeling Toolbox," in The Technology of Binaural Listening, edited by J. Blauert, Berlin Heidelberg: Springer-Verlag, pp. 33-56, 2013.
