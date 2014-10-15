#ifndef CORRECTION_H
#define CORRECTION_H

#include "matrix.h"
#include <cmath>
using namespace libNumerics;

template <typename T>
T bicubicDistModel(const vector<T>& completeParams, const vector<T>& coefTerm) { return completeParams * coefTerm; }

template <typename T>
void undistortPixel(T& xr, T& yr, const vector<T>& paramsX, const vector<T>& paramsY, 
	const int x, const int y, T xp, T yp, const int degX, const int degY)
{
	vector<T> coefTermX(paramsX.size()), coefTermY(paramsY.size());
	T powX = 0;
	T ini_monoX = 0;
	if (x != xp) {
		powX = 1 / (x-xp); 
		ini_monoX = pow(x-xp, degX); }
	T const powY = y-yp;
	T monoX, monoY;
	int idx = 0;
	for (int ii = degX; ii >= 0; ii--) {
		monoY = 1;
		monoX = ini_monoX; 
		ini_monoX *= powX;
		for (int j = 0; j <= ii; j++){
			coefTermX[idx] = monoX * monoY; //coefTermX[idx] = pow(x-xp, ii-j) * pow(y-yp, j);
			monoY *= powY; // pow( y-yp, j)
			monoX *= powX; // pow( x-xp, ii-j)
			idx++;
		}
	}
	if (x != xp) ini_monoX = pow(x-xp, degY);
	idx = 0;
	for (int ii = degY; ii >= 0; ii--) {
		monoY = 1;
		monoX = ini_monoX; 
		ini_monoX *= powX;
		for (int j = 0; j <= ii; j++){
			coefTermY[idx] = monoX * monoY; //coefTermY[idx] = pow(x-xp, ii-j) * pow(y-yp, j);
			monoY *= powY; // pow( y-yp, j)
			monoX *= powX; // pow( x-xp, ii-j)
			idx++;
		}
	}
	xr = bicubicDistModel(paramsX, coefTermX) + xp;
	yr = bicubicDistModel(paramsY, coefTermY) + yp;
}

#endif
