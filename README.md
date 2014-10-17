DistCorr
========

Optical distortion calculation and correction; open source C/C++ library

/* Given set of calibration images (folder data), obtain correction polynomial, 
calculate RMSE of correction and correct calibration or any other images 
taken under the same camera settings. */

Usage:
// polynomial estimation:
polyestim length_threshold sampling_factor <input1.pgm> [input2.pgm...] <polyout_filename.txt>

// correction in image(-s):
distcorrect <poly_fname.txt> <input.pgm> <output.pgm>

    Copyright (C) 2014 Victoria Rudakova <vicrucann@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.