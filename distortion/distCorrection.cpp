/* Distortion correction main functions: given set of calibration images, 
obtain correction polynomial, calculate RMSE of correction and correct
calibration or any other images taken under the same camera settings
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
#include "matrix.h"
#include "correction.h"
#include "pgm_io.h"
#include "misc.h"
#include "spline.h"

#include "distortion.h"
#include "gaussian_convol_on_curve.h"
#include "straight_edge_points.h"
#include "ntuple_ll.h"
#include "pgm_io.h"

#include <algorithm>
#include <ctime>
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

using namespace libNumerics;

/*----------------------------------------------------------------------------*/
/*                                    Funcs                                    */
/*----------------------------------------------------------------------------*/
/* Measure RMSE of lines on the given image */
template <typename T>
T image_RMSE(int length_thresh, double down_factor, image_double image)
{
	printf("\n Measuring RMSE of lines on the image... \n");
	DistortedLines<T> distLines;
	int w, h;
	ntuple_ll point_set,p;
	double x, y;
	int total_nb_lines, total_threshed_nb_lines, threshed_nb_lines;
	ntuple_list convolved_pts;
	point_set = new_ntuple_ll(1);
	total_nb_lines = 0;
	total_threshed_nb_lines = 0; 
	threshed_nb_lines = 0;

	p = straight_edge_points(image,sigma,th_low,th_hi,min_length);
	w = image->xsize; h = image->ysize;
	/* copy the points set in the new image to the global set of points */
		for(int j=0; j<(int)p->size; j++) {
			if((int)p->list[j]->size > length_thresh) {
				add_ntuple_list(point_set, p->list[j]);
				p->list[j] = 0; }
			else
				threshed_nb_lines += 1; }
	printf("For image there are totally %d lines detected and %d of them are eliminated.\n", p->size, threshed_nb_lines);
	free_ntuple_ll(p);

	distLines.pushMemGroup((int)point_set->size);
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
	free_ntuple_ll(point_set);
	free_ntuple_list(convolved_pts);

	int order = 3;
	T xp = (T)w/2+0.2, yp = (T)h/2+0.2;
	int sizexy = (order + 1) * (order + 2) / 2;  
	vector<T> paramsX = vector<T>::zeros(sizexy); paramsX[sizexy-3] = 1;
	vector<T> paramsY = vector<T>::zeros(sizexy); paramsY[sizexy-2] = 1; 
	T rmse = distLines.RMSE(paramsX, paramsY, order, order, xp, yp);
	T rmse_max = distLines.RMSE_max(paramsX, paramsY, order, order, xp, yp);
	printf("RMSE / maximum RMSE: %12.6g / %12.6g \n", rmse, rmse_max);
	return rmse;
}

/* Given an image and a correction polynomial. Apply it to every pixel and save result to output folder */
template <typename T>
image_double correct_image(char* imname, char* out_name, int spline_order, const vector<T>& poly_params_inv, 
	const int degX, const int degY) 
{	
	/* read the data */
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	vector<T> paramsX = poly_params_inv.copyRef(0, sizex-1);
	vector<T> paramsY = poly_params_inv.copyRef(sizex, sizex+sizey-1);
	printf("\n Undistorted image is being calculated... \n");
	image_double image_in = read_pgm_image_double(imname);
	int wi = image_in->xsize, he = image_in->ysize;
	T xp = (T)wi/2+0.2, yp = (T)he/2+0.2;
	image_double image_out = new_image_double(wi, he);
	std::cout << "[" << std::flush;
    prepare_spline(image_in, spline_order);
	for (int y = 0; y < he; y++) { 
		for (int x = 0; x < wi; x++) {
			/* do the correction for every pixel */
			T p1=0, p2=0;
			undistortPixel(p1, p2, paramsX, paramsY, x, y, xp, yp, degX, degY);
			double clr = interpolate_image_double(image_in, spline_order, p1, p2);
			if (clr < 0) clr = 0; 
			if (clr > 255) clr = 255;
			image_out->data[x+y*wi] = clr;
		}
		/* output progress */
		T percent = ((T)y / (T)he)*100;
		if (!(y % (int)(0.2*he))) std::cout << (int)percent+1 << "%" << std::flush;
		else if (!(y % (int)(0.04*he))) std::cout << "." << std::flush;
	}
	std::cout << "] " << std::flush;
	/* Save the result image to file */
	write_pgm_image_double(image_out, out_name);
	//free_image_double(image_out);
	free_image_double(image_in);
	printf(" done \n");
	return image_out;
}

/* Pop out one character from char array */
char popchar( char *c, int idx, int size) {
	char res = c[idx];
	for (int i = idx; i < size; i++)
		if (i > idx) c[i-1] = c[i];
	return res;
}

/* Difines the degee values for X and Y polinomials from file */
void read_degree(FILE *pfile, int& degX, int& degY) {
	degX = 0; degY = 0;
	char buffer[250], sign[2];
	double coef = 0;
	const char* grid = "#", *star = "*" , *plus = "+", *minu = "-";
	bool flagX = false, flagY = false;
	while (!feof(pfile)) {
		/* expects each monimial to have a form of: "+/- coef * x^deg1 * y^deg2 " */
		/* deg1 and/or deg2 can be zeros, leaving just a coefficent value */
		fscanf(pfile, "%s", sign);
		/* if sign "#" is met - it's a comment, we may skip it. */
		if (strcmp(sign, grid) != 0) {
			char mono1[5], mono2[5];
			fscanf(pfile, "%lf", &coef);
			long currPos = ftell (pfile);
			fscanf(pfile, "%s", sign);
			if (strcmp(sign, star) != 0) {
				if (!feof(pfile)) {
					fseek(pfile, currPos, SEEK_SET);
					assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	} }
			else {
				fscanf(pfile, "%s", mono1);
				currPos = ftell (pfile);
				fscanf(pfile, "%s", sign);
				if (strcmp(sign, star) != 0) {
					if (!feof(pfile)) {
						fseek(pfile, currPos, SEEK_SET);
						assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	}
				}
				else fscanf(pfile, "%s", mono2); }
			/* pop-out "x^" and "y^" so that to calculate degrees. */	
			popchar(mono1, 0, 5); popchar(mono1, 0, 5);
			popchar(mono2, 0, 5); popchar(mono2, 0, 5);
			int tmpdeg = std::max(atoi(mono1), atoi(mono2));
			/* there are two polynomials in file: degree of X is defined by flagX, second - flagY. */
			/* assumed first poly is for X direction, second is for Y. */
			if (flagX) { if (tmpdeg > degX) degX = tmpdeg; }
			else { if (tmpdeg > degY) degY = tmpdeg; }		
		}
		else {
			fscanf(pfile, "%s", buffer);
			if (degX == 0 && degY == 0) flagX = true;
			else {flagX = false; flagY = true;}
		}
	}
	fseek(pfile, 0, SEEK_SET);
}

/* Finds an index of coefTerm[] vector based on given x and y degrees. */
int coefIdx(int degree, int x, int y) {
	int a1 = degree + 1;
	int n = std::abs(x+y-degree)+1;
	int an = x+y+1;
	int Sn = n * (a1 + an);
	Sn /= 2;
	return Sn-an+y; }

/* Reads the poly coefficients and insert them into vector of coefficients - coefTerm[]. */
/* Also returns the degrees for each polynomial. */
template <typename T>
vector<T> read_poly(char* fname, int& degX, int& degY) {
	FILE *pfile;
	pfile = fopen(fname, "r");
	if(pfile == NULL) error("unable to open file %s.\n", fname);
	/* get the degrees for each polynomial - need it to know the size of vectors paramsX and paramsY. */
	read_degree(pfile, degX, degY);
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	vector<T> paramsX(sizex), paramsY(sizey);
	char buffer[500], sign[2];
	const char* grid = "#", *star = "*", *plus = "+", *minu = "-";
	bool flagX = false, flagY = false;
	while (!feof(pfile)) {
		/* reading is done by the same manner as in "read_degree(pfile, degX, degY);" */
		fscanf(pfile, "%s", sign);
		if (strcmp(sign, grid) != 0 && (strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 ) ) {
			char mono1[5], mono2[5];
			double coef = 0;
			fscanf(pfile, "%lf", &coef);
			long currPos = ftell (pfile);
			/* save the coef sign. */
			if (strcmp(sign, minu) == 0) coef *= -1;
			fscanf(pfile, "%s", sign);
			if (strcmp(sign, star) != 0) {
				if (!feof(pfile)) {
					fseek(pfile, currPos, SEEK_SET);
					assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	} 	}
			else {
				fscanf(pfile, "%s", mono1);
				currPos = ftell (pfile);
				fscanf(pfile, "%s", sign);
				if (strcmp(sign, star) != 0) {
					if (!feof(pfile)) {
						fseek(pfile, currPos, SEEK_SET);
						assert(strcmp(sign, plus) == 0 || strcmp(sign, minu) == 0 || strcmp(sign, grid) == 0); 	} 	}
				else fscanf(pfile, "%s", mono2); }
			char x = popchar(mono1, 0, 5); popchar(mono1, 0, 5);
			char y = popchar(mono2, 0, 5); popchar(mono2, 0, 5);
			/* see which degree belongs to which variable. */
			int tmpdegX = 0, tmpdegY = 0;
			if ( x == 'x')		tmpdegX = atoi(mono1);
			else if ( x == 'y')	tmpdegY = atoi(mono1);
			if ( y == 'y') 	tmpdegY = atoi(mono2);
			/* save the coef to accoring paramsX/Y; fist poly in file belongs to paramsX, second - paramsY. */
			if (flagX && !flagY) { 
				int idx = coefIdx(degX, tmpdegX, tmpdegY);
				paramsX[idx] = coef; }
			else {
				int idx = coefIdx(degY, tmpdegX, tmpdegY);
				paramsY[idx] = coef; }		
		}
		else {
			fscanf(pfile, "%s", buffer);
			if (!flagX) flagX = true;
			else flagY = true; }
	}
	fclose(pfile);
	vector<T> poly_params(sizex+sizey);
	/* copy the paramsX and paramsY to one vector. */
	/* later we will be able to separate them since the degrees for each poly are known. */
	for (int k = 0; k < sizex; k++) poly_params[k] = paramsX[k];
	for (int k = 0; k < sizey; k++) poly_params[k+sizex] = paramsY[k];
	return poly_params;
}

/*----------------------------------------------------------------------------*/
/*                                    Main                                    */
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv)
{ 
	printf("INFO parameters use: \ndistcorrect length_threshold sampling_factor <poly_fname.txt> <input.pgm> <output.pgm> \n");
	printf("=======================\n \n");
	
	int degX, degY;
	printf(" Reading the correction polynomial from file... \n");
	vector<double> poly_params_inv = read_poly<double>(argv[3], degX, degY);
	printf("done. \n");
	int spline_order = 5;
	image_double img = correct_image<double>(argv[4], argv[5], spline_order, poly_params_inv, degX, degY);
	//image_RMSE<double>(int(atoi(argv[1])), double(atoi(argv[2])), img);
	
	free_image_double(img);
	return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
