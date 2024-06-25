## SUpDEq - Spatial Upsampling by Directional Equalization (Description of the Original Branch)
This toolbox is a MATLAB implementation of the SUpDEq method, as presented in [1-7]. In general, the SUpDEq method is an approach to generate dense HRTF datasets based on sparsely measured HRTF datasets. For example, applying the SUpDEq method, a (technically) appropriate full-spherical dense HRTF dataset with 2702 directions can be derived from only 38 actually measured HRTFs. 

Basically, the method attempts to remove direction-dependent temporal and spectral components of the sparse HRTF dataset by a spectral division (equalization) of the HRTFs with corresponding rigid sphere transfer functions (STFs). This equalized HRTF dataset is then transformed to the spherical harmonics (SH) domain by means of a spherical Fourier transform [9-10]. Next, a de-equalization is performed by extracting the equalized HRTFs at the spatial sampling points of the dense grid (SH-interpolation / spatial upsampling) and multiplying these HRTFs with the respective STFs in frequency domain.

A comparison of various methods for spherical harmonics interpolation of time-aligned (in terms of SUpDEq "equalized") HRTFs is presented in [7]. The study also contains a perceptual evaluation of the SUpDEq method. 

The function `supdeq_interpHRTF` furthermore contains various pre-processing and interpolation methods discussed in [7,8,11,18], which (in most cases) can be freely combined. The function also contains the post-interpolation magnitude correction `MCA` (Magnitude-Corrected and Time-Aligned Interpolation) presented in [11] to further improve interpolation results obtained with time-aligned interpolation. 

Various extensions of the toolbox, as presented in [2-8,20], allow e.g., the synthesis of near-field HRTFs based on far-field datasets (distance variation), distance error compensation of measured HRTF datasets, or low frequency extension of (equalized) measured HRTFs.

## SUpDEq - Spatial Upsampling by Directional Equalization (MCA Upsampling Optimisation, Actual Branch)
Complementing the work of Arend et al. (2023) [11] in the original branch, our focus in this project has been on enhancing interpolation methods for HRTFs. We have implemented and optimized the SUpDEq method along with the Magnitude-Corrected and Time-Aligned Interpolation (MCA) technique to achieve this goal. The objective is to significantly improve the accuracy and perceived quality of spatial audio reproduction


## Requirements
The toolbox is implemented in MATLAB R2015b, R2018a, R2020a, and R2022b and requires the Signal Processing Toolbox. Older versions of MATLAB might also work. The following third party toolboxes are required for full functionality. All listed toolboxes are part of the SUpDEq toolbox and are stored in the folder "thirdParty".

- [SOFiA](https://github.com/AudioGroupCologne/SOFiA) (Sound Field Analysis Toolbox) [12]  
- [AKtools](https://www.ak.tu-berlin.de/menue/publications/open_research_tools/aktools/) [13]   
- [SOFA API](https://sourceforge.net/projects/sofacoustics/) (Spatially Oriented Format for Acoustics) [14]  
- [AMT](http://amtoolbox.sourceforge.net) (Auditory Modeling Toolbox) [15]   
- [SFS](https://github.com/sfstoolbox/sfs-matlab) (Sound Field Synthesis Toolbox) [16]   
- [Triangle/Ray Intersection](https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection) [17]   
- [SARITA](https://github.com/AudioGroupCologne/SARITA) (Spherical Array Interpolation by Time Alignment) [18]   

SOFiA and AMT partly work with MEX files. The required MEX files for Mac and Windows (64 Bit) are pre-compiled and stored in the respective folders. 


## Installation
For a simple utilisation simply clone or zip download this repository and add it to your Matlab path.

However, if you want to execute the MCA demo with the HUTUBS database as HRIRS (script `supdeq_demo_MCA_HUTUBS.`) you will have to download this database directly from [the official page](https://depositonce.tu-berlin.de/items/dc2a3076-a291-417e-97f0-7697e332c960) (The HUTUBS head-related transfer function (HRTF) database) [21]  and copy the unzipped content of the HIRIRs folder directly into the folder `materials/HUTUBSHRIRs`.

## First Steps
After installation, we recommend to go through the `supdeq_demo` script step by step for the original SUpDEq implementation as presented in [1]. This will lead you through the basic processing. The resulting dense HRTF dataset can then be perceptually and technically compared to the (dense) reference HRTF dataset with the provided evaluation functions.

## MCA Upsampling Optimisation Usage
To test the MCA with ILD compensation use the demo files `supdeq_demo_MCA` or `supdeq_demo_MCA_HUTUBS`. Please, mind that the latter needs the HUTUBS database to be manually downloaded as described previously. Further analysis of the results of the `supdeq_demo_MCA` and `supdeq_demo_MCA_HUTUBS` script is provided by the `MAG_error_dILD/dILD_error_ku100.m` and `MAG_error_dILD/magnitude_error_ku100.m` scripts. The Interim results are automaticly stored into the `mat_data` folder. The computed errors can eventually be used by the plotting scripts in the `plots` folder.

Notes: For the HUTUBS-subject91 dataset, magnitude and ILD errors cannot be computed due to the reference having more sampling grid points than the interpolated data.

## HRTF Interpolation
For an extended selection of pre-processing, interpolation, and post-processing methods, we recommend checking the `supdeq_interpHRTF` function. The function provides various SUpDEq pre-processing methods as well as other time-alignment methods as described in [7]. The function allows SH, Natural Neighbor, Barycentric, or SARITA interpolation of the (pre-processed) HRTFs [7,8,18]. Pre-processing and interpolation can (in most cases) be freely combined to obtain the best results depending on the specific conditions. Moreover, the function allows to perform post-interpolation magnitude-correction using the MCA method (MCA - Magnitude-Corrected and Time-Aligned Interpolation) presented in [11], further improving the interpolation results. The script `supdeq_demo_MCA` offers an easy introduction to MCA interpolation.

## Near-Field HRTF Synthesis

For our studies on near-field HRTFs and nearby sound sources (e.g., [19][20]), we implemented various methods for near-field HRTF synthesis based on far-field data, which we integrated in the SUpDEq toolbox. The function `supdeq_dvf` allows synthesizing near-field HRTFs using distance variation functions (DVFs) in combination with acoustic parallax effects, as used for the listening experiments in [20]. With the function `supdeq_rangeExt`, HRTFs in SH domain can be shifted in distance using radial functions. 

## Original Documentation
A HTML documentation of the MATLAB code can be found here:  
[http://audiogroup.web.th-koeln.de/SUpDEq_doc/index.html](http://audiogroup.web.th-koeln.de/SUpDEq_doc/index.html)  
The third-party toolboxes are excluded from the documentation.

## MCA Upsampling Optimisation Main Files

    1. supdeq_interpHRTF.m
        - This script has been extended to incorporate ILD compensation. It implements the algorithm to generate ILD filters based on interpolated HRTF datasets.

    2. supdeq_demo_MCA.m
        - A script to generate .mat files containing interpolated data for the ku100 dataset.

    3. supdeq_demo_MCA_HUTUBS.m
        - Similar to `supdeq_demo_MCA.m`, but specifically for the HUTUBS-subject91 dataset.


## Directory Structure and Files in `MAG_error_dILD` folder

- `magnitude_error_ku100.`
Computes the magnitude error for the ku100 dataset and saves the result as a .mat file in the `mat_data` folder.

- `dILD_error_ku100.m`
Computes the ILD error for the ku100 dataset and saves the result as a .mat file in the `mat_data` folder.

## Plot Scripts in `plots` Folder

  - `plot_HRTFs_ku100.m`
Plots the HRTFs for the ku100 dataset.

  - `plot_HRTFs_HUTUBS.m`
Plots the HRTFs for the HUTUBS-subject91 dataset.

  - `plot_mag_errors_ku100.m`
Generates a plot for the magnitude error of the ku100 dataset.

  - `plot_dILD_errors_ku100.m`
Generates a plot for the ILD error of the ku100 dataset.

## About the Original Branch
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


## About the MCA Upsampling Optimisation
This optimisation project was carried out as part of the Virtual Acoustic Reality course on the Audiocommunication and Technology curriculum at the technical university of Berlin.

Maximilian Jakob Wiedemann (485066)
Maximilian.wiedenmann@camput-tu.berlin.de

Raphaël Guillaume Gillioz (485054)
r.gillioz@campus.tu-berlin.de

Mai Pham (485086)
mai.thuc.an.pham@campus.tu-berlin.de

Erik Hansen (485059)
erik.hansen@campus.tu-berlin.de

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
[11] J. M. Arend, C. Pörschmann, S. Weinzierl, and F. Brinkmann, “Magnitude-Corrected and Time-Aligned Interpolation of Head-Related Transfer Functions,” IEEE/ACM Trans. Audio Speech Lang. Proc., vol. 31, pp. 3783–3799, 2023.   
[12] B. Bernschütz, C. Pörschmann, S. Spors, and S. Weinzierl, “SOFiA Sound Field Analysis Toolbox,” in Proceedings of the International Conference on Spatial Audio - ICSA 2011, 2011, pp. 8–16.  
[13] F. Brinkmann and S. Weinzierl, “AKtools - an open software toolbox for signal acquisition, processing, and inspection in acoustics,” in Proceedings of the 142nd AES Convention, Berlin, Germany, 2017, pp. 1–6.  
[14] P. Majdak, Y. Iwaya, T. Carpentier, R. Nicol, M. Parmentier, A. Roginska, Y. Suzuki, K. Watanabe, H. Wierstorf, H. Ziegelwanger, and M. Noisternig, “Spatially Oriented Format for Acoustics: A Data Exchange Format Representing Head-Related Transfer Functions,” in Proceedings of the 134th AES Convention, Rome, Italy, 2013, pp. 1–11.  
[15] P. Søndergaard and P. Majdak, "The Auditory Modeling Toolbox," in The Technology of Binaural Listening, edited by J. Blauert, Berlin Heidelberg: Springer-Verlag, pp. 33-56, 2013.  
[16] H. Wierstorf and S. Spors, "Sound Field Synthesis Toolbox," in Proceedings of the 132th AES Convention, Budapest, Hungary, 2012, pp. 1–4.  
[17] Jaroslaw Tuszynski (2021). Triangle/Ray Intersection (https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection), MATLAB Central File Exchange. Retrieved October 29, 2021.  
[18] T. Lübeck, J. M. Arend, and C. Pörschmann, “Spatial Upsampling of Sparse Spherical Microphone Array Signals,” IEEE/ACM Trans. Audio Speech Lang. Proc., vol. 31, pp. 1163–1174, 2023.  
[19] J. M. Arend, H. R. Liesefeld, and C. Pörschmann, “On the influence of non-individual binaural cues and the impact of level normalization on auditory distance estimation of nearby sound sources,” Acta Acust., vol. 5, no. 10, pp. 1–21, 2021.  
[20] J. M. Arend, M. Ramírez, H. R. Liesefeld, and C. Pörschmann, “Do near-field cues enhance the plausibility of non-individual binaural rendering in a dynamic multimodal virtual acoustic scene?,” Acta Acust., vol. 5, no. 55, pp. 1–14, 2021.
[21] F. Brinkmann, M. Dinakaran, R. Pelzer, P. Grosche, D. Voss, S. Weinzierl, "A Cross-Evaluated Database of Measured and Simulated HRTFs Including 3D Head Meshes, Anthropometric Features, and Headphone Impulse Responses," J. Audio Eng. Soc., vol. 67, pp. 705-718, 2019.

