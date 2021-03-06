/* Bicubic distortion model: obtain using Levenberg-Marquardt 
minimization.
    Copyright (C) 2014 <vicrucann@gmail.com>

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
*/

#ifndef DISTORTION_H
#define DISTORTION_H

#include "LM.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace libNumerics;

template <typename T> vector<T> bicubicDistModel(const vector<T>& completeParams, const matrix<T>& coefTerm);

/// Class for one line.
template <typename T>
class LineData {
public:
	LineData() { nPoints = 0; }

private: 
	std::vector<T> _pointX, _pointY; ///< Vectors for x and y coordinates.
	int nPoints;
	matrix<T> _coefTermX, _coefTermY; ///< Matrices for keeping constant values of polynomial.

public:
	int sizeLine() { return nPoints; }
	T x(int i) { return _pointX[i]; } ///< Access x coordinate of index \a i.
	T y(int i) { return _pointY[i]; } ///< Access y coordinate of index \ai.
	void pushPoint(const T x, const T y); ///< Add point to a line.
	void theta(const vector<T>& paramsX, const vector<T>& paramsY, T& alpha, T& beta); ///< Sines and cosines of line angle for each point in a line.
	vector<T> jacobian(const vector<T>& paramsX, const vector<T>& paramsY, const vector<int>& flagX, const vector<int>& flagY) const ;
	vector<T> residuals(const vector<T>& paramsX, const vector<T>& paramsY) const;
	T RMSE(const vector<T>& paramsX, const vector<T>& paramsY) const;
	void coefTermsCalc(int degX, int degY, T xp, T yp, T scale = 0);

private:
	void coefTermFill(matrix<T>& coefTerm, int deg, T x, T y) const;
}; // LineData

template <typename T>
class LMRectifyDistortion;

/// Class for a set of lines.
template <typename T>
class DistortedLines {
public:
	std::vector< LineData<T> > _line;
	int nLines, nGroups, nImages; 
	std::vector<int> nlines4Group;
	std::vector<int> nlines4Image;
	DistortedLines() { nLines = 0; nGroups = 0; }
	int totalPointsNumber();
	void pushMemGroup(int numLines);
	void defineImages(std::vector<int> &images){ nlines4Image = images; nImages = images.size(); }
	void pullMemoryLine(void);
	void pushPoint(int idxLine, T valX, T valY) { _line[idxLine].pushPoint(valX, valY); } ///< Add a point to a line with index \a idxLine.
	T RMSE(const vector<T>& paramsX, const vector<T>& paramsY, int degX, int degY, T xp, T yp);
	T RMSE_max(const vector<T>& paramsX, const vector<T>& paramsY, int degX, int degY,  T xp, T yp);
	vector<T> correctionLMA(vector<T>& paramsX, vector<T>& paramsY, vector<int>& flagX, vector<int>& flagY, 
		int degX, int degY, T xp, T yp);
	vector<T> verification(const vector<T>& paramsX, const vector<T>& paramsY, const vector<int>& flagX, const vector<int>& flagY, 
		int b_order, int c_order, T xp, T yp);
	void write2file(char* fname);

private:
	void estimatedThetas(const vector<T>& paramsX, const vector<T>& paramsY, vector<T>& alpha, vector<T>& beta, int degX, int degY, T xp, T yp);
	void normalization(DistortedLines<T>& normDistLines, T& scale, T xp, T yp);
}; // DistortionLines

/// Class to refine the distortion polynomial parameters.
template <typename T>
class LMRectifyDistortion : public MinLM<T>
{
public:
	LMRectifyDistortion(int oX, int oY, vector<int>& flagx, vector<int>& flagy, DistortedLines<T>& normDistLines, T scale, T xp, T yp);

private:
	int orderX, orderY;
	vector<int> flagX, flagY;
	DistortedLines<T> distLines;

public:
	virtual void modelData(const vector<T>& P, vector<T>& ymodel) const;
	virtual void modelJacobian(const vector<T>& P, matrix<T>& J) const;
}; // LMRectiryDistortion

template <typename T> void denormalization(vector<T>& denormX, vector<T>& denormY, const vector<T>& normX, const vector<T>& normY, 
	T scaleX, T scaleY, int orderX, int orderY);

/// Returns a correction polynomial for the given distortion one.
template <typename T>
vector<T> getParamsInv(const vector<T>& paramsX, const vector<T>& paramsY, int degX, int degY, int w, int h, T xp, T yp);

/// Returns a correction polynomial for the given corrected and distorted coordinates.
template <typename T>
vector<T> getParamsCorrection(vector<T>& x_corr, vector<T>& y_corr, vector<T>& x_dist, vector<T>& y_dist, int degX, int degY, T xp, T yp);

// Need to see definitions for templates...
#include "distortion.cpp"

#endif
