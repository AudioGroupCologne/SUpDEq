## SUpDEq - Spatial Upsampling by Directional Equalization
This toolbox is a MATLAB implementation of the SUpDEq method, as presented in [1-6]. In general, the SUpDEq method is an approach to generate dense HRTF datasets based on sparsely measured HRTF datasets. For example, applying the SUpDEq method, a (technically) appropriate full-spherical dense HRTF dataset with 2702 directions can be derived from only 38 actually measured HRTFs. 

Basically, the method attempts to remove direction-dependent temporal and spectral components of the sparse HRTF dataset by a spectral division (equalization) of the HRTFs with corresponding rigid sphere transfer functions (STFs). This equalized HRTF dataset is then transformed to the spherical harmonics (SH) domain by means of a spherical Fourier transform [7][8]. Next, a de-equalization is performed by extracting the equalized HRTFs at the spatial sampling points of the dense grid (SH-interpolation / spatial upsampling) and multiplying these HRTFs with the respective STFs in frequency domain.

Various extentions of the toolbox, as presented in [2-6], allow e.g. the synthesis of near-field HRTFs based on far-field datasets (distance variation), distance error compensation of measured HRTF datasets, or low frequency extention of (equalized) measured HRTFs.

## Requirements
The toolbox is implemented in MATLAB R2015b and R2018a and requires the Signal Processing Toolbox. Older versions of MATLAB might also work. The following third party toolboxes are required for full functionality. All listed toolboxes are part of the SUpDEq toolbox and are stored in the folder "thirdParty".

- SOFiA (Sound Field Analysis Toolbox) [9]  
 https://github.com/AudioGroupCologne/SOFiA
- AKtools [10]   
https://www.ak.tu-berlin.de/menue/digitale_ressourcen/research_tools/aktools/
- SOFA API (Spatially Oriented Format for Acoustics) [11]  
https://sourceforge.net/projects/sofacoustics/
- AMT (Auditory Modeling Toolbox) [12]  
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
[2] C. Pörschmann and J. M. Arend, “A Method for Spatial Upsampling of Directivity Patterns of Human Speakers by Directional Equalization,” in Proceedings of the 45th DAGA, 2019, pp. 1458–1461.
[3] J. M. Arend and C. Pörschmann, “Synthesis of Near-Field HRTFs by Directional Equalization of Far-Field Datasets,” in Proceedings of the 45th DAGA, 2019, pp. 1454–1457.
[4] C. Pörschmann and J. M. Arend, “Obtaining Dense HRTF Sets from Sparse Measurements in Reverberant Environments,” in Proceedings of the AES International Conference on Immersive and Interactive Audio, York, UK, 2019, pp. 1–10.
[5] C. Pörschmann, J. M. Arend, and F. Brinkmann, “Spatial upsampling of individual sparse head-related transfer function sets by directional equalization,” in Proceedings of the 23rd International Congress on Acoustics, 2019, pp. 4870–4877.
[6] J. M. Arend and C. Pörschmann, “Spatial upsampling of sparse head-related transfer function sets by directional equalization - Influence of the spherical sampling scheme,” in Proceedings of the 23rd International Congress on Acoustics, 2019, pp. 2643–2650.
[7] B. Bernschütz, “Microphone Arrays and Sound Field Decomposition for Dynamic Binaural Recording,” TU Berlin, 2016.  
[8] B. Rafaely, Fundamentals of Spherical Array Processing. Berlin Heidelberg: Springer-Verlag, 2015.  
[9] B. Bernschütz, C. Pörschmann, S. Spors, and S. Weinzierl, “SOFiA Sound Field Analysis Toolbox,” in Proceedings of the International Conference on Spatial Audio - ICSA 2011, 2011, pp. 8–16.  
[10] F. Brinkmann and S. Weinzierl, “AKtools - an open software toolbox for signal acquisition, processing, and inspection in acoustics,” in Proceedings of the 142nd AES Convention, Berlin, Germany, 2017, pp. 1–6.  
[11] P. Majdak, Y. Iwaya, T. Carpentier, R. Nicol, M. Parmentier, A. Roginska, Y. Suzuki, K. Watanabe, H. Wierstorf, H. Ziegelwanger, and M. Noisternig, “Spatially Oriented Format for Acoustics: A Data Exchange Format Representing Head-Related Transfer Functions,” in Proceedings of the 134th AES Convention, Rome, Italy, 2013, pp. 1–11.  
[12] P. Søndergaard and P. Majdak, "The Auditory Modeling Toolbox," in The Technology of Binaural Listening, edited by J. Blauert, Berlin Heidelberg: Springer-Verlag, pp. 33-56, 2013.
