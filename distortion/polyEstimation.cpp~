/* Distortion correction main functions: given set of calibration images, 
obtain correction polynomial, calculate RMSE of correction and correct
calibration or any other images taken under the same camera settings
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
*/
#include "distortion.h"
#include "gaussian_convol_on_curve.h"
#include "straight_edge_points.h"
#include "ntuple_ll.h"
#include "pgm_io.h"
#include "misc.h"

#include <algorithm>
#include <cstring>
#include <cstdio>
#include <cstdlib>

 /*----------------------------------------------------------------------------*/
/* Static parameters.
 */
static const double sigma = 1.0;
static const double th_low = 1.0;
static const double th_hi = 10.0;
static const double min_length = 100.0;

/* Static parameters for Gaussian convolution along the lines */
static const double unit_sigma = 0.8;
static const double Nsigma = 0.8;
static const bool resampling = true;
static const bool eliminate_border = false;
static const double up_factor = 1.0;

/*----------------------------------------------------------------------------*/
/*                                    Funcs                                    */
/*----------------------------------------------------------------------------*/

template <typename T>
void read_images(DistortedLines<T>& distLines, int& w, int& h, int argc, char ** argv, int beginIdx, int endIdx) {
	int w_tmp = 0, h_tmp = 0;
	ntuple_ll point_set,p;
	image_double image;
	int length_thresh;
	double down_factor, x, y;
	int total_nb_lines, total_threshed_nb_lines, threshed_nb_lines;
	ntuple_list convolved_pts;
	/* check parameters */
	if( argc < 4 ) error("use: %s length_threshold sampling_factor <input1.pgm> [input2.pgm...] <polyout_filename>", argv[0]);
	assert(beginIdx >= 3 && beginIdx < argc-1 && endIdx >= 3 && endIdx < argc-1);
	length_thresh = (int)(atoi(argv[1])); 
	down_factor = (double)(atoi(argv[2])); 
	printf("There are %d input images. The minimal length of lines is set to %d\n", endIdx+1-beginIdx, length_thresh); 
	fflush(stdout);
	/* initialize memory */
	point_set = new_ntuple_ll(1);
	total_nb_lines = 0;
	total_threshed_nb_lines = 0; 
	std::vector<int> lines4img(endIdx+1-beginIdx);
	int idx = 0;
	/* process each of the input images */  
	for(int i = beginIdx; i <= endIdx; i++) {
		threshed_nb_lines = 0;
		/* open image, compute edge points, close it */
		image = read_pgm_image_double(argv[i]);
		p = straight_edge_points(image,sigma,th_low,th_hi,min_length);
		w = image->xsize; h = image->ysize;
		if (i == beginIdx) {w_tmp = w; h_tmp = h;}
		else { assert(w == w_tmp && h == h_tmp); }
		free_image_double(image);
		/* copy the points set in the new image to the global set of points */
		for(int j=0; j<(int)p->size; j++) {
			if((int)p->list[j]->size > length_thresh) {
				add_ntuple_list(point_set, p->list[j]);
				p->list[j] = 0; }
			else
				threshed_nb_lines += 1; }
		lines4img[idx] = p->size-threshed_nb_lines;
		idx++;
		printf("Image %d: totally %d lines detected, %d eliminated, %d left.\n", i-2, p->size, threshed_nb_lines, p->size-threshed_nb_lines); 
		fflush(stdout);
		total_nb_lines += p->size;
		total_threshed_nb_lines += threshed_nb_lines;
		free_ntuple_ll(p);
	}
	distLines.pushMemGroup((int)point_set->size);
	distLines.defineImages(lines4img);
	int countL = 0;
	for(unsigned int i=0; i<point_set->size; i++) {
		/* Gaussian convolution and sub-sampling */
		convolved_pts = gaussian_convol_on_curve(unit_sigma, Nsigma, resampling, eliminate_border, up_factor, down_factor, point_set->list[i]);
		countL++;
		/* Save points to DistortionLines structure */
		for(unsigned int j = 0; j < convolved_pts->size; j++) {
			x = convolved_pts->values[j*convolved_pts->dim];
			y = convolved_pts->values[j*convolved_pts->dim+1];
			distLines.pushPoint(countL-1,x,y); }
	}
	printf("Totally there are %d lines detected and %d of them are eliminated.\n", total_nb_lines, total_threshed_nb_lines);
	fflush(stdout);
	/* free memory */
	free_ntuple_ll(point_set);
	free_ntuple_list(convolved_pts);
}

template <typename T>
vector<T> incLMA(DistortedLines<T>& distLines, const int order, const int inc_order, T xp, T yp) {
	int sizexy = (order + 1) * (order + 2) / 2;  
	vector<T> paramsX = vector<T>::zeros(sizexy); paramsX[sizexy-3] = 1;
	vector<T> paramsY = vector<T>::zeros(sizexy); paramsY[sizexy-2] = 1; 
	T rmse = distLines.RMSE(paramsX, paramsY, order, order, xp, xp);
	T rmse_max = distLines.RMSE_max(paramsX, paramsY, order, order, xp, yp);
	printf("initial RMSE / maximum RMSE: %12.6g / %12.6g \n", rmse, rmse_max);
	fflush(stdout);
	const int beginOrder = 3;
	vector<T> midParams(1);
	for (int i = beginOrder; i <= order; i=i+inc_order) {
		int sizebc = (i+1) * (i+2) / 2;
		vector<T> b = vector<T>::zeros(sizebc);
	    vector<T> c = vector<T>::zeros(sizebc);
		vector<int> flagX = vector<int>::ones(sizebc);
		vector<int> flagY = vector<int>::ones(sizebc);
		flagX[sizebc-3] = flagX[sizebc-2] = flagX[sizebc-1] = 0;
		flagY[sizebc-2] = flagY[sizebc-3] = flagY[sizebc-1] = 0;
		if (i == beginOrder) { 
			b[sizebc-3] = 1; 
			c[sizebc-2] = 1;  }
		else {
			int sizebcOld = (i-inc_order+1) * (i-inc_order+2) / 2;
			for (int k = 0; k < sizebcOld; k++) {
				b[sizebc-sizebcOld+k] = midParams[k];
				c[sizebc-sizebcOld+k] = midParams[sizebcOld+k]; }
			 }
		midParams = distLines.correctionLMA(b, c, flagX, flagY, i, i, xp, yp);
		rmse = distLines.RMSE(midParams.copy(0, sizebc-1), midParams.copyRef(sizebc, sizebc+sizebc-1), i, i, xp, yp);
		rmse_max = distLines.RMSE_max(midParams.copyRef(0, sizebc-1), midParams.copyRef(sizebc, sizebc+sizebc-1), i, i, xp, yp);
		printf("degree %i current RMSE / maximum RMSE: %12.6g / %12.6g \n", i, rmse,  rmse_max);
		fflush(stdout);
	}
	printf("\n Iterative linear minimization step: \n");
	fflush(stdout);
	T diff = 100;
	T prevRmse = 1e+16; /* infinity */
	vector<int> flagX = vector<int>::ones(sizexy); 
	vector<int> flagY = vector<int>::ones(sizexy); 
	flagX[sizexy-3] = flagX[sizexy-1] = flagX[sizexy-2] = 0;
	flagY[sizexy-2] = flagY[sizexy-1] = flagY[sizexy-3] = 0;
	vector<T> estParams(sizexy+sizexy);
	int iter = 0;
	while (diff > 0.001) {
		estParams = distLines.verification(midParams.copy(0, sizexy-1), midParams.copy(sizexy, sizexy+sizexy-1), flagX, flagY, order, order, xp, yp);
		midParams = estParams;
		rmse = distLines.RMSE(estParams.copy(0, sizexy-1), estParams.copyRef(sizexy, sizexy+sizexy-1), order, order, xp, yp);
		rmse_max = distLines.RMSE_max(estParams.copyRef(0, sizexy-1), estParams.copyRef(sizexy, sizexy+sizexy-1), order, order, xp, yp);
		printf("current RMSE / maximum RMSE: %12.6g / %12.6g \n", rmse, rmse_max );
		fflush(stdout);
		
		diff = fabs(prevRmse - rmse);
		prevRmse = rmse;
		iter++;
	}
	return estParams;
}

template <typename T>
vector<T> polyInv(const vector<T>& poly_params, const int degX, const int degY, int wi, int he, T xp, T yp) {	
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	return getParamsInv(poly_params.copyRef(0, sizex-1), poly_params.copyRef(sizex, sizex+sizey-1), degX, degY, wi, he, xp, yp);;
}

template <typename T>
void printMono(FILE* pfile, T mono, int degX, int degY) {
	if (degX != 0 && degY != 0) {
		if (mono > 0)	fprintf(pfile, "+ %.16g * x^%i * y^%i\n", mono, degX, degY);
		else			fprintf(pfile, "- %.16g * x^%i * y^%i\n", -1*mono, degX, degY); }
	else if (degX == 0 && degY == 0) {
		if (mono > 0)	fprintf(pfile, "+ %.16g\n", mono);
		else			fprintf(pfile, "- %.16g\n", -1*mono); }
	else if (degX == 0) {
		if (mono > 0)	fprintf(pfile, "+ %.16g * y^%i\n", mono, degY);
		else			fprintf(pfile, "- %.16g * y^%i\n", -1*mono, degY); }
	else {
		if (mono > 0)	fprintf(pfile, "+ %.16g * x^%i\n", mono, degX);
		else			fprintf(pfile, "- %.16g * x^%i\n", -1*mono, degX); }
}

template <typename T>
void save_poly(char* fname, vector<T>& poly_params, const int degX, const int degY) {
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	vector<T> paramsX = poly_params.copyRef(0, sizex-1);
	vector<T> paramsY = poly_params.copyRef(sizex, sizex+sizey-1);
	FILE *pfile;
	pfile = fopen(fname, "wt");
	if(pfile == NULL)
		error("cannot open file %s.\n", fname);
	int idx = 0;
	fprintf(pfile, "# polyX(x,y): \n");
	for (int i = degX; i >= 0; i--) {
		for (int j = 0; j <= i; j++) {
			printMono(pfile, paramsX[idx], i-j, j);
			idx++; 	} }
	idx = 0;
	fprintf(pfile, "# polyY(x,y): \n");
	for (int i = degY; i >= 0; i--) {
		for (int j = 0; j <= i; j++) {
			printMono(pfile, paramsY[idx], i-j, j);
			idx++; } }
	fclose(pfile);
}

/*----------------------------------------------------------------------------*/
/*                                    Main                                    */
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv)
{
	printf("INFO parameters use: \npolyestim length_threshold sampling_factor <input1.pgm> [input2.pgm...] <polyout_filename.txt> \n");
	printf("=======================\n \n");
	fflush(stdout);
	/* Detect the lines from the input images. Push the data into DistortedLines */ 
	int idx_begin = 3, idx_end = argc-2; 
	int w, h;	
	DistortedLines<double> distLines;
	read_images(distLines, w, h, argc, argv, idx_begin, idx_end);
	//distLines.write2file("in/points.txt");
		
	/* Distortion correction: calcualte the polynomial parameters */
	printf("\n Incremental LMA for distortion polynomial... \n");
	fflush(stdout);
	double xp = (double)w/2+0.2, yp = (double)h/2+0.2; /* +0.2 - to avoid integers */
	const int order = 11;
	const int inc = 2; /* increment; only odd orders will be taken */
	vector<double> poly_params = incLMA <double> (distLines, order, inc, xp, yp);
	printf("\n Incremental LMA done. The polynomial was obtained. \n");
	fflush(stdout);

	/* Get an inverse polynomial */ 
	printf("\n Distortion correction polynomial... "); fflush(stdout);
	vector<double> poly_params_inv = polyInv<double>(poly_params, order, order, w, h, xp, yp);
	printf("done. \n");
	fflush(stdout);

	/* Save the correction polynomial to output file */
	printf("\n Saving polynomial to file... "); fflush(stdout);
	save_poly(argv[argc-1], poly_params_inv, order, order);
	printf("done. \n");

	return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
