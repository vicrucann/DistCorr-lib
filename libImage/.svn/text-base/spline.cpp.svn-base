/* Spline interpolation.
    Copyright (C) 2007 Lionel Moisan <Lionel.Moisan@parisdescartes.fr>
    Copyright (C) 2010 Pascal Monasse <monasse@imagine.enpc.fr>

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

#include "spline.h"
//#include "libIO/nan.h"
#include <cmath>
#include <cfloat>
#include <cstring>

static double initcausal(double* c, int step, int n, double z)
{
    double zk,z2k,iz,sum;

    zk = z; iz = 1./z;
    z2k = pow(z,(double)n-1.);
    sum = c[0] + z2k * c[step*(n-1)];
    z2k = z2k*z2k*iz;
    for(int k = 1; k <= n-2; k++) {
        sum += (zk+z2k)*c[step*k];
        zk *= z;
        z2k *= iz;
    }
    return (sum/(1.-zk*zk));
}

static double initanticausal(double* c, int step, int n, double z)
{
    return (z/(z*z-1.)) * (z*c[step*(n-2)]+c[step*(n-1)]);
}

static void invspline1D(double* c, int step, int size, double* z, int npoles)
{
    /* normalization */
    double lambda=1.0;
    for(int k = npoles-1; k >= 0; k--)
        lambda *= (1.-z[k])*(1.-1./z[k]);
    for(int n = size-1; n >= 0; n--)
        c[step*n] *= static_cast<double>(lambda);

    for(int k=0 ; k < npoles; k++) { // Loop on poles
        /* forward recursion */
        c[0] = static_cast<double>(initcausal(c, step, size, z[k]));
        for(int n=1; n < size; n++)
            c[step*n] += static_cast<double>(z[k]*c[step*(n-1)]);
        /* backward recursion */
        c[step*(size-1)] = static_cast<double>(initanticausal(c, step, size, z[k]));
        for(int n=size-2; n >= 0; n--)
            c[step*n] = static_cast<double>(z[k]*(c[step*(n+1)]-c[step*n]));
    }
}

/// Put in array \a z the poles of the spline of given \a order.
static bool fill_poles(double* z, int order)
{
    switch(order) {
    case 1:
        break;
        // case 2: z[0]=-0.17157288;  /* sqrt(8)-3 */
        // break;
    case 3: z[0]=-0.26794919;  /* sqrt(3)-2 */ 
        break;
        // case 4: z[0]=-0.361341; z[1]=-0.0137254;
        // break;
    case 5: z[0]=-0.430575; z[1]=-0.0430963;
        break;
        // case 6: z[0]=-0.488295; z[1]=-0.0816793; z[2]=-0.00141415;
        // break;
    case 7: z[0]=-0.53528; z[1]=-0.122555; z[2]=-0.00914869;
        break;
        // case 8: z[0]=-0.574687; z[1]=-0.163035; z[2]=-0.0236323;
        // z[3]=-0.000153821;
        // break;
    case 9: z[0]=-0.607997; z[1]=-0.201751; z[2]=-0.0432226; z[3]=-0.00212131;
        break;
        // case 10: z[0]=-0.636551; z[1]=-0.238183; z[2]=-0.065727;
        // z[3]=-0.00752819; z[4]=-0.0000169828;
        // break;
    case 11: z[0]=-0.661266; z[1]=-0.27218; z[2]=-0.0897596; z[3]=-0.0166696; 
        z[4]=-0.000510558;
        break;
    default:
        return false;
    }
    return true;
}

/// Prepare image for cardinal spline interpolation.
bool prepare_spline(image_double& im, int order)
{
    if(order < 3)
        return true;

    // Init poles of associated z-filter
    double z[5];
    if(! fill_poles(z, order))
        return false;
    int npoles = order/2;

	for(int y = 0; y < im->ysize; y++) // Filter on lines
		invspline1D(im->data+y*im->xsize, 1, im->xsize, z, npoles);
	for(int x = 0; x < im->xsize; x++) // Filter on columns
		invspline1D(im->data+x, 1*im->xsize, im->ysize, z, npoles);
    return true;
}

/* c[] = values of interpolation function at ...,t-2,t-1,t,t+1,... */

/* coefficients for cubic interpolant (Keys' function) */
static void keys(double* c, double t, double a)
{
    double t2 = t*t;
    double at = a*t;
    c[0] = a*t2*(1.0f-t);
    c[1] = (2.0f*a+3.0f - (a+2.0f)*t)*t2 - at;
    c[2] = ((a+2.0f)*t - a-3.0f)*t2 + 1.0f;
    c[3] = a*(t-2.0f)*t2 + at;
}

/* coefficients for cubic spline */
static void spline3(double* c, double t)
{
    double tmp = 1.0f-t;
    c[0] = 0.1666666666f *t*t*t;
    c[1] = 0.6666666666f -0.5f*tmp*tmp*(1.0f+t);
    c[2] = 0.6666666666f -0.5f*t*t*(2.0f-t);
    c[3] = 0.1666666666f *tmp*tmp*tmp;
}

/* pre-computation for spline of order >3 */
static void init_splinen(double* a, int n)
{
    a[0] = 1.;
    for(int k=2; k <= n; k++)
        a[0]/=(double)k;
    for(int k=1; k <= n+1; k++)
        a[k] = - a[k-1] *(n+2-k)/k;
}

/* fast integral power function */
static float ipow(float x, int n)
{
    float res;
    for(res = 1.; n; n>>=1) {
        if(n&1)
            res *= x;
        x *= x;
    }
    return res;
}

/* coefficients for spline of order >3 */
static void splinen(double* c, double t, double* a, int n)
{
    memset((void*)c, 0, (n+1)*sizeof(double));
    for(int k=0; k <= n+1; k++) { 
        double xn = ipow(t+(double)k, n);
        for(int i=k; i <= n; i++) 
            c[i] += a[i-k]*xn;
    }
}

/// Spline interpolation of given \a order of image \a im at point (x,y).
/// \a out must be an array of size the number of components.
/// Supported orders: 0(replication), 1(bilinear), -3(Keys's bicubic), 3, 5, 7,
/// 9, 11.
/// \a paramKeys is Keys's parameter, only used for order -3.
/// Success means a valid order and pixel in image.
bool interpolate_spline( image_double& im, int order, double x, double y, double& out, double paramKeys)
{
	double  cx[12],cy[12];

	/* CHECK ORDER */
	if(order != 0 && order != 1 && order != -3 &&
		order != 3 && order != 5 && order != 7 && order != 9 && order != 11)
		return false;

	double ak[13];
	if(order > 3)
		init_splinen(ak, order);

	bool bInside = false;
	/* INTERPOLATION */
	if(order == 0) { /* zero order interpolation (pixel replication) */
		int xi = (int)floor((double)x);
		int yi = (int)floor((double)y);
		bInside = valid_image_double(im, xi, yi);//im.valid(xi, yi);
		if(bInside) {
			double p = im->data[xi+yi*im->xsize]; //im.pixel(xi, yi);
			out = p;
		}
	} else { /* higher order interpolations */
		bInside = (x>=0.0f && x<=(double)im->xsize && y>=0.0f && y<=(double)im->ysize);
		if(bInside) {
			x -= 0.5f; y -= 0.5f;
			int xi = (x<0)? -1: (int)x;
			int yi = (y<0)? -1: (int)y;
			double ux = x - (double)xi;
			double uy = y - (double)yi;
			switch(order)  {
			case 1: /* first order interpolation (bilinear) */
				cx[0] = ux; cx[1] = 1.0f-ux;
				cy[0] = uy; cy[1] = 1.0f-uy;
				break;
			case -3: /* third order interpolation (bicubic Keys' function) */
				keys(cx, ux, paramKeys);
				keys(cy, uy, paramKeys);
				break;
			case 3: /* spline of order 3 */
				spline3(cx, ux);
				spline3(cy, uy);
				break;
			default: /* spline of order >3 */
				splinen(cx, ux, ak, order);
				splinen(cy, uy, ak, order);
				break;
			}
			int n2 = (order==-3)? 2: (order+1)/2;
			int n1 = 1-n2;
			/* this test saves computation time */
			if(valid_image_double(im, xi+n1, yi+n1) && valid_image_double(im, xi+n2,yi+n2)) {
				out = 0.0f;
				for(int dy = n1; dy <= n2; dy++) {
					int adrs = (xi+n1) + (yi+dy) * im->xsize;
					for(int dx = n1; dx <= n2; dx++) {
						double f = im->data[adrs];
						out += cy[n2-dy]*cx[n2-dx] * f;
						adrs++;
					}
				}

			} else
				out = 0.0f;
			for(int dy = n1; dy <= n2; dy++)
				for(int dx = n1; dx <= n2; dx++) {
					double v = 0.0f; // the image is not infinite, there is no data outside
					out += cy[n2-dy]*cx[n2-dx]*v;
				}

		}
	}
	return bInside;
}
