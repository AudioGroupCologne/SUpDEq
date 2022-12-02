## SUpDEq - Spatial Upsampling by Directional Equalization
This toolbox is a MATLAB implementation of the SUpDEq method, as presented in [1-7]. In general, the SUpDEq method is an approach to generate dense HRTF datasets based on sparsely measured HRTF datasets. For example, applying the SUpDEq method, a (technically) appropriate full-spherical dense HRTF dataset with 2702 directions can be derived from only 38 actually measured HRTFs. 

Basically, the method attempts to remove direction-dependent temporal and spectral components of the sparse HRTF dataset by a spectral division (equalization) of the HRTFs with corresponding rigid sphere transfer functions (STFs). This equalized HRTF dataset is then transformed to the spherical harmonics (SH) domain by means of a spherical Fourier transform [9-10]. Next, a de-equalization is performed by extracting the equalized HRTFs at the spatial sampling points of the dense grid (SH-interpolation / spatial upsampling) and multiplying these HRTFs with the respective STFs in frequency domain.

Various extentions of the toolbox, as presented in [2-6], allow e.g., the synthesis of near-field HRTFs based on far-field datasets (distance variation), distance error compensation of measured HRTF datasets, or low frequency extention of (equalized) measured HRTFs.

A comparison of various methods for spherical harmonics interpolation of time-aligned (in terms of SUpDEq "equalized") HRTFs is presented in [7]. The study also contains a perceptual evaluation of the SUpDEq method. 

The function `supdeq_interpHRTF` furthermore contains various pre-processing and interpolation methods discussed in [7-8], which can be freely combined. The function also contains the post-interpolation magnitude correction `MCA` (Magnitude-Corrected and Time-Aligned Interpolation) presented in [11] to further improve interpolation results obtained with time-aligned interpolation. 

## Requirements
The toolbox is implemented in MATLAB R2015b, R2018a, and R2020a and requires the Signal Processing Toolbox. Older versions of MATLAB might also work. The following third party toolboxes are required for full functionality. All listed toolboxes are part of the SUpDEq toolbox and are stored in the folder "thirdParty".

- [SOFiA](https://github.com/AudioGroupCologne/SOFiA) (Sound Field Analysis Toolbox) [12]  
- [AKtools](https://www.ak.tu-berlin.de/menue/publications/open_research_tools/aktools/) [13]   
- [SOFA API](https://sourceforge.net/projects/sofacoustics/) (Spatially Oriented Format for Acoustics) [14]  
- [AMT](http://amtoolbox.sourceforge.net) (Auditory Modeling Toolbox) [15]   
- [SFS](https://github.com/sfstoolbox/sfs-matlab) (Sound Field Synthesis Toolbox) [16]   
- [Triangle/Ray Intersection](https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection) [17]   

SOFiA and AMT partly work with MEX files. The required MEX files for Mac and Windows (64 Bit) are pre-compiled and stored in the respective folders. 


## Installation
Simply clone or zip download this repository and run the `supdeq_start` script. This will add all required paths and install SOFA, if not previously installed on your system. 

## First Steps
After installation, we recommend to go through the `supdeq_demo` script step by step. This will lead you through the basic processing. The resulting dense HRTF dataset can then be perceptually and technically compared to the (dense) reference HRTF dataset with the provided evaluation functions.

## HRTF Interpolation

For an extended selection of pre-processing, interpolation, and post-processing methods, we recommend checking the `supdeq_interpHRTF` function. The function provides various SUpDEq pre-processing methods as well as other time-alignment methods as described in [7]. The function allows SH, Natural Neighbor, or Barycentric interpolation of the (pre-processed) HRTFs [8]. Pre-processing and interpolation can be freely combined to obtain the best results depending on the specific conditions. Moreover, the function allows to perform post-interpolation magnitude-correction using the MCA method (MCA - Magnitude-Corrected and Time-Aligned Interpolation) presented in [11], further improving the interpolation results. The script `supdeq_demo_MCA` offers an easy introduction to MCA interpolation.

## Near-Field HRTF Synthesis

For our studies on near-field HRTFs and nearby sound sources (e.g., [18][19]), we implemented various methods for near-field HRTF synthesis based on far-field data, which we integrated in the SUpDEq toolbox. The function `supdeq_dvf` allows synthesizing near-field HRTFs using distance variation functions (DVFs) in combination with acoustic parallax effects, as used for the listening experiments in [19]. With the function `supdeq_rangeExt`, HRTFs in SH domain can be shifted in distance using radial functions. 

## Documentation
A HTML documentation of the MATLAB code can be found here:  
[http://audiogroup.web.th-koeln.de/SUpDEq_doc/index.html](http://audiogroup.web.th-koeln.de/SUpDEq_doc/index.html)  
The third-party toolboxes are excluded from the documentation.

## About
Johannes M. Arend and Christoph Pörschmann  
TH Köln - University of Applied Sciences  
Institute of Communications Engineering  
Department of Acoustics and Audio Signal Processing  
Betzdorfer Str. 2, D-50679 Cologne, Germany  
[https://www.th-koeln.de/akustik](https://www.th-koeln.de/akustik) 

J. M. Arend is now at Technical University of Berlin  
Audio Communication Group  
Einsteinufer 17c, D-10587 Berlin, Germany  
[https://www.ak.tu-berlin.de/menue/fachgebiet_audiokommunikation](https://www.ak.tu-berlin.de/menue/fachgebiet_audiokommunikation/)

Thanks to Fabian Brinkmann (Audio Communication Group, TU Berlin) for useful discussions and contributions!  

## References
[1] C. Pörschmann\*, J.M. Arend\*, and F. Brinkmann, “Directional Equalization of Sparse Head-Related Transfer Function Sets for Spatial Upsampling,” IEEE/ACM Trans. Audio, Speech, Lang. Process., vol. 27, no. 6, pp. 1060–1071, 2019, \*These authors contributed equally to this work.  
[2] C. Pörschmann and J. M. Arend, “A Method for Spatial Upsampling of Directivity Patterns of Human Speakers by Directional Equalization,” in Proceedings of the 45th DAGA, 2019, pp. 1458–1461.  
[3] J. M. Arend and C. Pörschmann, “Synthesis of Near-Field HRTFs by Directional Equalization of Far-Field Datasets,” in Proceedings of the 45th DAGA, 2019, pp. 1454–1457.  
[4] C. Pörschmann and J. M. Arend, “Obtaining Dense HRTF Sets from Sparse Measurements in Reverberant Environments,” in Proceedings of the AES International Conference on Immersive and Interactive Audio, York, UK, 2019, pp. 1–10.  
[5] C. Pörschmann, J. M. Arend, and F. Brinkmann, “Spatial upsampling of individual sparse head-related transfer function sets by directional equalization,” in Proceedings of the 23rd International Congress on Acoustics, 2019, pp. 4870–4877.  
[6] J. M. Arend and C. Pörschmann, “Spatial upsampling of sparse head-related transfer function sets by directional equalization - Influence of the spherical sampling scheme,” in Proceedings of the 23rd International Congress on Acoustics, 2019, pp. 2643–2650.   
[7] J. M. Arend, F. Brinkmann, and C. Pörschmann, “Assessing Spherical Harmonics Interpolation of Time-Aligned Head-Related Transfer Functions,” J. Audio Eng. Soc., vol. 69, no. 1/2, pp. 104–117, 2021.  
[8] C. Pörschmann, J. M. Arend, D. Bau, and T. Lübeck, “Comparison of Spherical Harmonics and Nearest-Neighbor based Interpolation of Head-Related Transfer Functions,” in Proceedings of the AES International Conference on Audio for Virtual and Augmented Reality (AVAR), Redmond, WA, USA, 2020, pp. 1–10.  
[9] B. Bernschütz, “Microphone Arrays and Sound Field Decomposition for Dynamic Binaural Recording,” TU Berlin, 2016.  
[10] B. Rafaely, Fundamentals of Spherical Array Processing. Berlin Heidelberg: Springer-Verlag, 2015.   
[11] J. M. Arend, C. Pörschmann, S. Weinzierl, F. Brinkmann, "Magnitude-Corrected and Time-Aligned Interpolation of Head-Related Transfer Functions," (Manuscript submitted for publication).   
[12] B. Bernschütz, C. Pörschmann, S. Spors, and S. Weinzierl, “SOFiA Sound Field Analysis Toolbox,” in Proceedings of the International Conference on Spatial Audio - ICSA 2011, 2011, pp. 8–16.  
[13] F. Brinkmann and S. Weinzierl, “AKtools - an open software toolbox for signal acquisition, processing, and inspection in acoustics,” in Proceedings of the 142nd AES Convention, Berlin, Germany, 2017, pp. 1–6.  
[14] P. Majdak, Y. Iwaya, T. Carpentier, R. Nicol, M. Parmentier, A. Roginska, Y. Suzuki, K. Watanabe, H. Wierstorf, H. Ziegelwanger, and M. Noisternig, “Spatially Oriented Format for Acoustics: A Data Exchange Format Representing Head-Related Transfer Functions,” in Proceedings of the 134th AES Convention, Rome, Italy, 2013, pp. 1–11.  
[15] P. Søndergaard and P. Majdak, "The Auditory Modeling Toolbox," in The Technology of Binaural Listening, edited by J. Blauert, Berlin Heidelberg: Springer-Verlag, pp. 33-56, 2013.  
[16] H. Wierstorf and S. Spors, "Sound Field Synthesis Toolbox," in Proceedings of the 132th AES Convention, Budapest, Hungary, 2012, pp. 1–4.  
[17] Jaroslaw Tuszynski (2021). Triangle/Ray Intersection (https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection), MATLAB Central File Exchange. Retrieved October 29, 2021.   
[18] J. M. Arend, H. R. Liesefeld, and C. Pörschmann, “On the influence of non-individual binaural cues and the impact of level normalization on auditory distance estimation of nearby sound sources,” Acta Acust., vol. 5, no. 10, pp. 1–21, 2021.  
[19] J. M. Arend, M. Ramírez, H. R. Liesefeld, and C. Pörschmann, “Do near-field cues enhance the plausibility of non-individual binaural rendering in a dynamic multimodal virtual acoustic scene?,” Acta Acust., vol. 5, no. 55, pp. 1–14, 2021.  

