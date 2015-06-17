Optical distortion calculation and correction tool (C/C++ library)
========

## Usage:

After running `cmake`, there will be two executables which user can run:  

`polyestim length_threshold sampling_factor {input1.pgm} [input2.pgm...] {polyout_filename.txt}`  

`distcorrect {poly_fname.txt} {input.pgm} {output.pgm}`  

## Main principle

Given set of calibration images (example is provided in `data` folder), obtain correction polynomial, calculate RMSE of correction and correct calibration or any other images taken under the same camera settings.

## The algorithm  

The library is based on research done by CMLA ENS-Cachan and IMAGINE LIGM ENPC, and the reference paper titled [LENS DISTORTION CORRECTION WITH A CALIBRATION HARP](http://www.researchgate.net/publication/221121089_Lens_distortion_correction_with_a_calibration_harp). The code was adapted from originally written Matlab prototype by Zhongwei Tang.

## Author information  

The software uses some `C`-based libraries for image processing etc. For more info refer to each `lib` subdirectory individually. The main [framework](https://github.com/vicrucann/DistCorr-lib/tree/master/distortion) and [polynomial estimator](https://github.com/vicrucann/DistCorr-lib/tree/master/libDistortion) were written by Victoria Rudakova. 
