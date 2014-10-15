#ifdef DISTORTION_H

using namespace libNumerics;

template <typename T>
vector<T> bicubicDistModel(const vector<T>& completeParams, const matrix<T>& coefTerm) 
{ 
	return coefTerm.t() * completeParams; 
}

/// One point with coordinates \a x and \a y is added to the line.
template <typename T>
void LineData<T>::pushPoint(const T x, const T y) 
{ 
	_pointX.push_back(x); 
	_pointY.push_back(y); 
	nPoints++; 
}

/// Outputs sine and cosine for the line.
template <typename T>
void LineData<T>::theta(const vector<T>& paramsX, const vector<T>& paramsY, T& alpha, T& beta) 
{ 
	vector<T> ux = bicubicDistModel(paramsX, _coefTermX); 
	vector<T> uy = bicubicDistModel(paramsY, _coefTermY);
	T Ax = mean(ux); 
	T Ay = mean(uy);
	T Axy = ux * uy / nPoints;
	T Vxx = (ux-Ax).qnorm() / nPoints;
	T Vyy = (uy-Ay).qnorm() / nPoints;
	T Vxy = Axy - Ax*Ay;
	T theta = 0.5 * atan2(-2*Vxy, Vxx-Vyy);
	alpha = std::sin(theta);
	beta = std::cos(theta);
}

/// Jacobian of the line.
template <typename T>
vector<T> LineData<T>::jacobian(const vector<T>& paramsX, const vector<T>& paramsY, const vector<int>& flagX, const vector<int>& flagY) const 
{
	T nPoints_inv = 1/ ((T)nPoints);
	vector<T> ux = bicubicDistModel(paramsX, _coefTermX);
	vector<T> uy = bicubicDistModel(paramsY, _coefTermY);

	T Ax = mean(ux); 
	T Ay = mean(uy);
	T Axy = ux * uy / nPoints;
	vector<T> x(ux-Ax);
	vector<T> y(uy-Ay);
	T Vxx = x.qnorm() / nPoints; 
	T Vyy = y.qnorm() / nPoints;
	T Vxy = Axy - Ax*Ay; //T Vxy = x*y/nPoints;
	T Vxx_yy = Vxx-Vyy;
	T koef = std::sqrt(Vxx_yy*Vxx_yy + 4*Vxy*Vxy);
	T rmse = std::sqrt(0.5*(Vxx + Vyy - koef)); 

	int sizeX = paramsX.size(); 
	int sizeY = paramsY.size();
	const matrix<T>& DxDb = _coefTermX;
	matrix<T> DxDc = matrix<T>::zeros(sizeX, nPoints);
	const matrix<T>& DyDc = _coefTermY;
	matrix<T> DyDb = matrix<T>::zeros(sizeY, nPoints); 

	vector<T> DaxDb(sizeX); 
	vector<T> DaxDc = vector<T>::zeros(sizeY); // dXdC = 0
	for (int k = 0; k < sizeX; k++) 
		DaxDb[k] = sum(DxDb.rowRef(k)) * nPoints_inv;
	vector<T> DayDb = vector<T>::zeros(sizeX); 
	vector<T> DayDc(sizeY);
	for (int k = 0; k < sizeY; k++) 
		DayDc[k] = sum(DyDc.rowRef(k)) * nPoints_inv;
	vector<T> DaxyDb = (DxDb * uy) * nPoints_inv;
	vector<T> DaxyDc = (DyDc * ux) * nPoints_inv;

	vector<T> DvxxDb(sizeX), DvyyDc(sizeY); 
	vector<T> DvxxDc = vector<T>::zeros(sizeY); 
	vector<T> DvyyDb = vector<T>::zeros(sizeX); 	
	vector<T> xn2 = 2*nPoints_inv * x;
	vector<T> yn2 = 2*nPoints_inv * y;
	for (int i = 0; i < sizeX; i++) 
		DvxxDb[i] = xn2 * (DxDb.rowRef(i) - DaxDb[i]);
	for (int i = 0; i < sizeY; i++) 
		DvyyDc[i] = yn2 * (DyDc.rowRef(i) - DayDc[i]);

	vector<T> DvxyDb = DaxyDb - DaxDb*Ay;
	vector<T> DvxyDc = DaxyDc - DayDc*Ax;	
	int jacidx = 0;
	vector<T> jac = vector<T>::zeros(sizeX+sizeY);
	for (int i = 0; i < sizeX; i++) {
		if (flagX[i] == 1) jac[jacidx] = 0.25/rmse * (DvxxDb[i] - 1/koef * (4*Vxy*DvxyDb[i] + Vxx_yy*DvxxDb[i]) ); 
		else jac[jacidx] = 0;
		jacidx++; }
	for (int i = 0; i < sizeY; i++) {
		if (flagY[i] == 1)  jac[jacidx] = 0.25/rmse * (DvyyDc[i] - 1/koef * (4*Vxy*DvxyDc[i] - Vxx_yy*DvyyDc[i]) );
		else jac[jacidx] = 0;
		jacidx++; }
	return jac;
}

/// Residual for each point in a line.
template <typename T>
vector<T> LineData<T>::residuals(const vector<T>& paramsX, const vector<T>& paramsY) const 
{
	vector<T> ux = bicubicDistModel(paramsX, _coefTermX); 
	vector<T> uy = bicubicDistModel(paramsY, _coefTermY);
	T Ax = mean(ux); 
	T Ay = mean(uy);
	T Axy = ux * uy / nPoints;
	vector<T> x(ux-Ax);
	vector<T> y(uy-Ay);
	T Vxx = x.qnorm() / nPoints;
	T Vyy = y.qnorm() / nPoints;
	T Vxy = Axy - Ax*Ay;
	T theta = 0.5 * atan2(-2*Vxy, Vxx-Vyy);
	return sin(theta) * x + cos(theta) * y; 
} 

template <typename T>
T LineData<T>::RMSE(const vector<T>& paramsX, const vector<T>& paramsY) const 
{
	vector<T> ux = bicubicDistModel(paramsX, _coefTermX);
	vector<T> uy = bicubicDistModel(paramsY, _coefTermY);
	T Ax = mean(ux); 
	T Ay = mean(uy);
	T Axy = ux * uy / nPoints;
	T Vxx = (ux-Ax).qnorm() / nPoints; 
	T Vyy = (uy-Ay).qnorm() / nPoints;
	T Vxy = Axy - Ax*Ay;
	T Vxx_yy = Vxx-Vyy;
	return std::sqrt(0.5*(Vxx + Vyy - std::sqrt(Vxx_yy*Vxx_yy + 4*Vxy*Vxy) )); 
}

/// Calcualtes the constant coefficients of polynomials. 
/// Must be done in advance for some functions.
template <typename T>
void LineData<T>::coefTermsCalc(int degX, int degY, T xp, T yp, T scale) 
{
	int sizex = (degX + 1) * (degX + 2) / 2;
	int sizey = (degY + 1) * (degY + 2) / 2;
	_coefTermX = matrix<T>(sizex, nPoints);
	_coefTermY = matrix<T>(sizey, nPoints);
	if (scale == 0) {
		coefTermFill(_coefTermX, degX, xp, yp);
		coefTermFill(_coefTermY, degY, xp, yp); }
	else {
		coefTermFill(_coefTermX, degX, 0, 0);
		coefTermFill(_coefTermY, degY, 0, 0); }
}

template <typename T>
void LineData<T>::coefTermFill(matrix<T>& coefTerm, int deg, T x, T y) const 
{
	int idx = 0;
	for (int ii = deg; ii >= 0; ii--) {
		for (int j = 0; j <= ii; j++){
			for (int k = 0; k < nPoints; k++) 
				coefTerm(idx, k) = pow(_pointX[k]-x, ii-j) * pow(_pointY[k]-y, j);
			idx++;
		}
	}
}

/// Number of all the points of all the lines.
template <typename T>
int DistortedLines<T>::totalPointsNumber() 
{
	int res = 0;
	for (int i = 0; i < nLines; i++)
		res += _line[i].sizeLine();
	return res;  
}

/// Allocate memory for \a numLines lines.
template <typename T>
void DistortedLines<T>::pushMemGroup(int numLines) 
{
	nGroups++;
	nlines4Group.push_back(numLines);
	nLines += numLines;
	_line.resize(nLines); 
}

/// Delete a line from a group.
template <typename T>
void DistortedLines<T>::pullMemoryLine(void) 
{
	nlines4Group[nlines4Group.size()-1]--;
	nLines--;
	_line.resize(nLines); 
}

template <typename T>
T DistortedLines<T>::RMSE(const vector<T>& paramsX, const vector<T>& paramsY, int degX, int degY, T xp, T yp) 
{
	T rmse = 0;
	for (int i = 0; i < nLines; i++) {
		_line[i].coefTermsCalc(degX, degY, xp, yp);
		vector<T> residuals = _line[i].residuals(paramsX, paramsY);
		rmse += residuals.qnorm(); }
	return std::sqrt(rmse / totalPointsNumber()); 
}

template <typename T>
T DistortedLines<T>::RMSE_max(const vector<T>& paramsX, const vector<T>& paramsY, int degX, int degY,  T xp, T yp) 
{
	T rmseMax = 0;
	for (int i = 0; i < nLines; i++) {
		_line[i].coefTermsCalc(degX, degY, xp, yp);
		vector<T> residuals = _line[i].residuals(paramsX, paramsY);
		for (int j = 0; j < _line[i].sizeLine(); j++)
			rmseMax = std::max(fabs(residuals(j)), rmseMax);
	}
	return rmseMax; 
}

/// LMA for lines correction.
/// Returns the distortion polynomial coefficients.
template <typename T>
vector<T> DistortedLines<T>::correctionLMA(vector<T>& paramsX, vector<T>& paramsY, vector<int>& flagX, vector<int>& flagY, int degX, int degY, T xp, T yp) 
{
	int sizex = paramsX.size();
	int sizey = paramsY.size();
	T scale = 0;
	DistortedLines<T> normDistLines; 
	normalization(normDistLines, scale, xp, yp);
	vector<T> denormX(sizex), denormY(sizey);
	denormalization(denormX, denormY, paramsX, paramsY, 1/scale, 1/scale, degX, degY);
	int sizeP = sizex + sizey;
	vector<T> P(sizeP);
	for (int i = 0; i < sizex; i++) P[i] = denormX[i]; // initial distribution
	for (int i = 0; i < sizey; i++) P[i+sizex] = denormY[i];
	vector<T> ydata = vector<T>::zeros(nLines);
	int maxIters = 1500; // thresholds
	T tolFun = 0.01;
	LMRectifyDistortion<T> lm(degX, degY, flagX, flagY, normDistLines, scale, xp, yp);  // LMA
	lm.minimize(P, ydata, tolFun, maxIters);
	vector<T> estDenormX(sizex), estDenormY(sizey);
	denormalization(estDenormX, estDenormY, P.copyRef(0, sizex-1), P.copyRef(sizex, sizex+sizey-1), scale, scale, degX, degY);
	vector<T> newP(sizeP);
	for (int i = 0; i < sizex; i++) newP[i] = estDenormX[i];
	for (int i = 0; i < sizey; i++) newP[i+sizex] = estDenormY[i];
	return newP;
}

/// Linear minimization.
template <typename T>
vector<T> DistortedLines<T>::verification(const vector<T>& paramsX, const vector<T>& paramsY, const vector<int>& flagX, const vector<int>& flagY, 
	int b_order, int c_order, T xp, T yp) 
{	
	int sizex = flagX.size();
	int sizey = flagY.size();
	vector<T> paramsB = vector<T>::zeros(sizex); 
	vector<T> paramsC = vector<T>::zeros(sizey);
	paramsB[sizex-3] = paramsX[sizex-3]; paramsB[sizex-1] = paramsX[sizex-1];
	paramsC[sizey-2] = paramsY[sizey-2]; paramsC[sizey-1] = paramsY[sizey-1];

	vector<T> alpha(nLines), beta(nLines); 
	estimatedThetas(paramsX, paramsY, alpha, beta, b_order, c_order, xp, yp); // calculates sin(theta) and cos(theta) -> alpha and beta
	T scale = 1;
	DistortedLines<T> normDistLines; 
	normalization(normDistLines, scale, xp, yp);
	int idx = 0;
	for (int ii = b_order; ii >= 0; ii--) {
		for (int jj = 1; jj <= ii+1; jj++) {
			paramsB[idx] = paramsB[idx] * pow(scale, ii);
			idx++; } }
	idx = 0;
	for (int ii = c_order; ii >= 0; ii--) {
		for (int jj = 1; jj <= ii+1; jj++) {
			paramsC[idx] = paramsC[idx] * pow(scale, ii);
			idx++; } } 
	int nSamples = 0;
	int coefMatSize = nLines + sum(flagX) + sum(flagY);
	matrix<T> coefMat = matrix<T>::zeros( coefMatSize, coefMatSize); 
	vector<T> bb = vector<T>::zeros(coefMat.nrow());

	typedef std::vector <matrix<T> > Cell;
	Cell coefTermB(nLines), coefTermC(nLines); 
	int coefMatIdx = 0;
	int lineIdx = 0;
	for (int i = 0; i < nLines; i++) {
		int prev_nb_lines = 0;
		if (i != 0) prev_nb_lines = i;

		LineData<T> one_line = normDistLines._line[prev_nb_lines];
		int nb_samples = one_line.sizeLine();
		nSamples += nb_samples;

		coefTermB[i] = matrix<T>::zeros(sizex, nb_samples); 
		coefTermC[i] = matrix<T>::zeros(sizey, nb_samples);
		idx = 0;
		for (int ii = b_order; ii >= 0; ii--) {
			for (int jj = 1; jj <= ii+1; jj++) {
				for (int k = 0; k < one_line.sizeLine(); k++)
					coefTermB[i](idx, k) = pow(one_line.x(k), ii-(jj-1)) * pow(one_line.y(k), jj-1);
				idx++; } }
		idx = 0;
		for (int ii = c_order; ii >= 0; ii--) {
			for (int jj = 1; jj <= ii+1; jj++) {
				for (int k = 0; k < one_line.sizeLine(); k++)
					coefTermC[i](idx, k) = pow(one_line.x(k), ii-(jj-1)) * pow(one_line.y(k), jj-1); 
				idx++; } }

		vector<T> tmp1 = vector<T>::zeros(nLines);
		tmp1[coefMatIdx] = -nb_samples;

		bb[coefMatIdx] = -nb_samples*alpha[i]*xp - nb_samples*beta[i]*yp;

		vector<T> tmp2 = vector<T>::zeros(sum(flagX));
		idx = 0;
		for (int jj = 0; jj < sizex; jj++) {
			if (flagX[jj] == 1) {
				tmp2[idx] = alpha[i] * sum(coefTermB[i].rowRef(jj));
				idx++; }
			else bb[coefMatIdx] -= paramsB[jj] * alpha[i] * sum(coefTermB[i].rowRef(jj)); }

		vector<T> tmp3 = vector<T>::zeros(sum(flagY));
		idx = 0;
		for (int jj = 0; jj < sizey; jj++) {
			if (flagY[jj] == 1) {
				tmp3[idx] = beta[i] * sum(coefTermC[i].rowRef(jj));
				idx++; }
			else bb[coefMatIdx] -= paramsC[jj] * beta[i] * sum(coefTermC[i].rowRef(jj)); }

		for (int k = 0; k < tmp1.size(); k++) coefMat(coefMatIdx, k) = tmp1[k];
		for (int k = 0; k < tmp2.size(); k++) coefMat(coefMatIdx, tmp1.size() + k) = tmp2[k];
		for (int k = 0; k < tmp3.size(); k++) coefMat(coefMatIdx, tmp1.size()+tmp2.size()+k) = tmp3[k];
		coefMatIdx++;
		lineIdx++;
	}

	matrix<T> tcoefTermB = matrix<T>::zeros(sizex, nSamples);
	matrix<T> tcoefTermB_bis = matrix<T>::zeros(sizex, nSamples);
	matrix<T> tcoefTermB_ter = matrix<T>::zeros(sizex, nSamples);

	matrix<T> tcoefTermC = matrix<T>::zeros(sizey, nSamples);
	matrix<T> tcoefTermC_bis = matrix<T>::zeros(sizey, nSamples);
	matrix<T> tcoefTermC_ter = matrix<T>::zeros(sizey, nSamples);

	matrix<T> bb_b = matrix<T>::zeros(sizex, nSamples); 
	matrix<T> bb_b_bis = matrix<T>::zeros(sizex, nSamples);
	matrix<T> bb_c = matrix<T>::zeros(sizey, nSamples);
	matrix<T> bb_c_bis = matrix<T>::zeros(sizey, nSamples);

	int idxb = 0, idxc = 0;
	for (int i = 0; i < nLines; i++) {
		const matrix<T>& tmpb = coefTermB[i];
		for (int k = 0; k < tmpb.ncol(); k++) {
			tcoefTermB.paste(0, idxb+k, tmpb.col(k));
			tcoefTermB_bis.paste(0, idxb+k, (-alpha[i]*alpha[i])*tmpb.col(k) );
			tcoefTermB_ter.paste(0, idxb+k, (-alpha[i]*beta[i])*tmpb.col(k) );
			bb_b.paste(0, idxb+k, (-xp*alpha[i]*alpha[i])*tmpb.col(k) );
			bb_b_bis.paste(0, idxb+k, (-yp*alpha[i]*beta[i])*tmpb.col(k) );  }
		idxb += tmpb.ncol();

		const matrix<T>& tmpc = coefTermC[i];
		for (int k = 0; k < tmpc.ncol(); k++) {
			tcoefTermC.paste(0, idxc+k, tmpc.col(k));
			tcoefTermC_bis.paste(0, idxc+k, (-alpha[i]*beta[i])*tmpc.col(k) );
			tcoefTermC_ter.paste(0, idxc+k, (-beta[i]*beta[i])*tmpc.col(k) );
			bb_c.paste(0, idxc+k, (-xp*alpha[i]*beta[i])*tmpc.col(k) );
			bb_c_bis.paste(0, idxc+k, (-yp*beta[i]*beta[i])*tmpc.col(k) ); }
		idxc += tmpc.ncol();
	}

	matrix<T> bb1 = -bb_b - bb_b_bis;
	int num_flag_gamma = nLines, num_flag_b = sum(flagX);
	for (int i = 0; i < sizex; i++) {
		if (flagX[i] == 1) {
			bb[coefMatIdx] = sum( bb1.rowRef(i) );
			const vectorRef<T> v_total_coef_term_b_bis = tcoefTermB_bis.rowRef(i);
			idx = 0;
			for (int jj = 0; jj < sizex; jj++) {
				if (flagX[jj] == 1) {
					coefMat(coefMatIdx, num_flag_gamma+idx) = tcoefTermB.rowRef(jj) * v_total_coef_term_b_bis;
					idx++; }
				else 
					bb[coefMatIdx] -= paramsB[jj] * tcoefTermB.rowRef(jj) * v_total_coef_term_b_bis; }

			const vectorRef<T> v_total_coef_term_b_ter = tcoefTermB_ter.rowRef(i);
			idx = 0;
			for (int jj = 0; jj < sizey; jj++) {
				if (flagY[jj] == 1) {
					coefMat(coefMatIdx, num_flag_gamma+num_flag_b+idx) = tcoefTermC.rowRef(jj) * v_total_coef_term_b_ter;
					idx++; }
				else 
					bb[coefMatIdx] -= paramsC[jj] * tcoefTermC.rowRef(jj) * v_total_coef_term_b_ter; }
			coefMatIdx++;
		} 	
	}

	matrix<T> bb2 = -bb_c - bb_c_bis;
	int num_flag_c = sum(flagY);
	for (int i = 0; i < sizey; i++) {
		if (flagY[i] == 1) {
			bb[coefMatIdx] = sum( bb2.rowRef(i) );
			const vectorRef<T> v_total_coef_term_c_bis = tcoefTermC_bis.rowRef(i);
			idx = 0;
			for (int jj = 0; jj < sizex; jj++) {
				if (flagX[jj] == 1) {
					coefMat(coefMatIdx, num_flag_gamma+idx) = tcoefTermB.rowRef(jj) * v_total_coef_term_c_bis;
					idx++;	}
				else 
					bb[coefMatIdx] -= paramsB[jj] * tcoefTermB.rowRef(jj) * v_total_coef_term_c_bis; }

			const vectorRef<T> v_total_coef_term_c_ter = tcoefTermC_ter.rowRef(i);
			idx = 0;
			for (int jj = 0; jj < sizey; jj++) {
				if (flagY[jj] == 1) {
					coefMat(coefMatIdx, num_flag_gamma+num_flag_c+idx) = tcoefTermC.rowRef(jj) * v_total_coef_term_c_ter;
					idx++; }
				else 
					bb[coefMatIdx] -= paramsC[jj] * tcoefTermC.rowRef(jj) * v_total_coef_term_c_ter; }
			coefMatIdx++;
		} 	
	}

	for (int i = nLines; i < coefMat.nrow(); i++) {
		for (int j = 0; j < nLines; j++)
			coefMat(i,j) = coefMat(j, i); }

	vector<T> sol(bb.size());
	solveLU(coefMat, bb, sol);

	vector<T> correction_params = sol.copyRef(nLines, sol.size()-1);
	idx = 0;
	int idx1 = 0;
	for (int ii = b_order; ii >= 0; ii--) {
		for (int jj = 1; jj <= ii+1; jj++) {
			if (flagX[idx] == 1) {
				paramsB[idx] = correction_params[idx1];
				idx1++; }
			idx++;
		} 	}

	idx = 0;
	for (int ii = c_order; ii >= 0; ii--) {
		for (int jj = 1; jj <= ii+1; jj++) {
			if (flagY[idx] == 1) {
				paramsC[idx] = correction_params[idx1];
				idx1++; }
			idx++;
		} 	}

	vector<T> estim_params_b = vector<T>::zeros(sizex);
	vector<T> estim_params_c = vector<T>::zeros(sizey);
	idx = 0;
	for (int ii = b_order; ii >= 0; ii--) {
		for (int jj = 1; jj <= ii+1; jj++) {
			estim_params_b[idx] = paramsB[idx] / pow(scale, ii);
			idx++; } }

	idx = 0;
	for (int ii = c_order; ii >= 0; ii--) {
		for (int jj = 1; jj <= ii+1; jj++) {
			estim_params_c[idx] = paramsC[idx] / pow(scale, ii);
			idx++; } }

	vector<T> estim_params(sizex + sizey);
	for (int k = 0; k < sizex; k++) estim_params[k] = estim_params_b[k];
	for (int k = 0; k < sizey; k++) estim_params[k+sizex] = estim_params_c[k];
	return estim_params;
}

/// Returns vectors of sines and cosines of angle for each line.
template <typename T>
void DistortedLines<T>::estimatedThetas(const vector<T>& paramsX, const vector<T>& paramsY, vector<T>& alpha, vector<T>& beta, int degX, int degY, T xp, T yp) 
{
	for (int i = 0; i < nLines; i++) {
		_line[i].coefTermsCalc(degX, degY, xp, yp);
		_line[i].theta(paramsX, paramsY, alpha[i], beta[i]); }
}

template <typename T>
void DistortedLines<T>::normalization(DistortedLines<T>& normDistLines, T& scale, T xp, T yp) 
{
	T dist = 0;
	for (int i = 0; i < nLines; i++) {
		for (int j = 0; j < _line[i].sizeLine(); j++)
			dist += std::sqrt( pow(_line[i].x(j) - xp, 2) + pow(_line[i].y(j) - yp, 2) ); }
	scale = dist / totalPointsNumber();
	normDistLines.pushMemGroup(nLines);
	for (int i = 0; i < nLines; i++) {
		for (int j = 0; j < _line[i].sizeLine(); j++)
			normDistLines.pushPoint(i, (_line[i].x(j) - xp) / scale , (_line[i].y(j) - yp) / scale ); }
}

template <typename T>
LMRectifyDistortion<T>::LMRectifyDistortion(int oX, int oY, vector<int>& flagx, vector<int>& flagy, DistortedLines<T>& normDistLines, T scale, T xp, T yp) 
{
	orderX = oX; orderY = oY;
	flagX = flagx; flagY = flagy;
	distLines = normDistLines;
	for (int i = 0; i < distLines.nLines; i++) distLines._line[i].coefTermsCalc(orderX, orderY, xp, yp, scale);
}

template <typename T>
void LMRectifyDistortion<T>::modelData(const vector<T>& P, vector<T>& ymodel) const {
	int sizeP = P.size();
	int sizeX = flagX.size(); int sizeY = flagY.size();
	assert( sizeX + sizeY == sizeP );
	ymodel = 0;
	for (int i = 0; i < distLines.nLines; i++)
		ymodel[i] = distLines._line[i].RMSE(P.copyRef(0, sizeX-1), P.copyRef(sizeX, sizeP-1));
}

template <typename T>
void LMRectifyDistortion<T>::modelJacobian(const vector<T>& P, matrix<T>& J) const  {
	int sizeP = P.size();
	int sizeX = flagX.size(); int sizeY = flagY.size();
	assert( sizeX + sizeY == sizeP );
	J = 0;
	for (int i = 0; i < distLines.nLines; i++) {
		vector<T> jacv = distLines._line[i].jacobian(P.copyRef(0, sizeX-1), P.copyRef(sizeX, sizeP-1), flagX, flagY);
		J.paste(i, 0, jacv.t() );
	}
}

template <typename T>
void denormalization(vector<T>& denormX, vector<T>& denormY, const vector<T>& normX, const vector<T>& normY, 
	T scaleX, T scaleY, int orderX, int orderY) 
{
	int idx = 0;
	for (int ii = orderX; ii >= 0; ii--) {
		for (int j = 1; j <= ii+1; j++) {
			denormX[idx] = normX[idx] / (pow(scaleX, ii-(j-1)) * pow(scaleY, j-1));
			idx++; }
	}
	idx = 0;
	for (int ii = orderY; ii >= 0; ii--) {
		for (int j = 1; j <= ii+1; j++) {
			denormY[idx] = normY[idx] / (pow(scaleX, ii-(j-1)) * pow(scaleY, j-1));
			idx++; }
	}
}

/// Returns a correction polynomial for the given corrected and distorted coordinates.
template <typename T>
vector<T> getParamsCorrection(vector<T>& x_corr, vector<T>& y_corr, vector<T>& x_dist, vector<T>& y_dist, int degX, int degY, T xp, T yp)
{
	int sizex = (degX+1)*(degX+2)/2;
	int sizey = (degY+1)*(degY+2)/2;
	int lenxy = x_corr.size();
	vector<T> x_dist_rad(x_dist-xp), y_dist_rad(y_dist-yp);
	vector<T> x_corr_rad(x_corr-xp), y_corr_rad(y_corr-yp);
	T norm_xy = 0;
	for (int i = 0; i < lenxy; i++) norm_xy += x_dist_rad[i]*x_dist_rad[i] + y_dist_rad[i]*y_dist_rad[i];
	norm_xy = std::sqrt(norm_xy);
	x_dist_rad  /= norm_xy; y_dist_rad  /= norm_xy;

	matrix<T> coefTermX(sizex, lenxy), coefTermY(sizey, lenxy);
	int idx = 0;
	for (int ii = degX; ii >= 0; ii--) {
		for (int j = 0; j <= ii; j++)  {
			for (int k = 0; k < lenxy; k++) 
				coefTermX(idx, k) = pow(x_dist_rad[k], ii-j) * pow(y_dist_rad[k], j);
			idx++; 
		} 	}
	idx = 0;
	for (int ii = degY; ii >= 0; ii--) {
		for (int j = 0; j <= ii; j++)  {
			for (int k = 0; k < lenxy; k++) 
				coefTermY(idx, k) = pow(x_dist_rad[k], ii-j) * pow(y_dist_rad[k], j);
			idx++; 
		} 	}

	int sizexy = sizex+sizey;
	matrix<T> coef_mat = matrix<T>::zeros(sizexy, sizexy);
	vector<T> m = vector<T>::zeros(sizexy);
	int count_i = 0;
	for (int i = 0; i < sizex; i++) {
		int count_j = 0;
		for (int j = 0; j < sizex; j++) {
			coef_mat(count_i, count_j) = coefTermX.rowRef(i) * coefTermX.rowRef(j);
			count_j++; }
		m[count_i] = x_corr_rad * coefTermX.rowRef(i);
		count_i++; }

	count_i = 0;
	for (int i = 0; i < sizey; i++) {
		int count_j = 0;
		for (int j = 0; j < sizey; j++) {
			coef_mat(sizex+count_i, sizex+count_j) = coefTermY.rowRef(i) * coefTermY.rowRef(j);
			count_j++; }
		m[sizex+count_i] = y_corr_rad * coefTermY.rowRef(i);
		count_i++; 	}

	matrix<T> normalization_mat1 = matrix<T>::zeros(sizexy, sizexy);
	matrix<T> inv_normalization_mat1 = matrix<T>::zeros(sizexy, sizexy);
	for (int i = 0; i < sizexy; i++) {
		if (i < sizex)
			normalization_mat1(i,i) = coef_mat(0, i);
		else
			normalization_mat1(i,i) = coef_mat(sizex, i); 
		inv_normalization_mat1(i,i) = 1 / normalization_mat1(i,i);
	}

	matrix<T> normalized_coef_mat = coef_mat * inv_normalization_mat1;
	matrix<T> normalization_mat2 = matrix<T>::zeros(sizexy, sizexy);
	for (int i = 0; i < sizexy; i++) normalization_mat2(i, i) = normalized_coef_mat(0, 0) / normalized_coef_mat(i, i);
	normalized_coef_mat =  normalization_mat2*normalized_coef_mat;

	vector<T> paramsInv(sizex+sizey);
	solveLU(normalized_coef_mat, normalization_mat2 * m, paramsInv);
	paramsInv = inv_normalization_mat1 * paramsInv;
	vector<T> denorm_paramsInv(sizex+sizey);
	vectorRef<T> denormX = denorm_paramsInv.copyRef(0, sizex-1);
	vectorRef<T> denormY = denorm_paramsInv.copyRef(sizex, sizex+sizey-1);
	denormalization(denormX, denormY,  // copyRef is not const?!
		paramsInv.copyRef(0, sizex),
		paramsInv.copyRef(sizex, sizex+sizey-1),
		norm_xy, norm_xy, degX, degY);

	return denorm_paramsInv;
}

/// Returns a correction polynomial for the given distortion one.
template <typename T>
vector<T> getParamsInv(const vector<T>& paramsX, const vector<T>& paramsY, int degX, int degY, int w, int h, T xp, T yp)
{
	int sizex = paramsX.size();
	int sizey = paramsY.size();
	// grid of coordinates with step 
	const int stepX = 0.007*w, stepY = 0.007*h; 
	int divresultx = w/stepX;
	int divresulty = h/stepY;
	int lenxy = divresultx * divresulty;
	vector<T> x_corr(lenxy), y_corr(lenxy);
	int idx = 0;
	for (int i = 0; i < divresulty; i++) {
		for (int j = 0; j < divresultx; j++) {
			x_corr[idx] = j*stepX;
			y_corr[idx] = i*stepY;
			idx++; 	} }
	matrix<T> coefTermX(sizex, lenxy), coefTermY(sizey, lenxy);
	idx = 0;
	for (int ii = degX; ii >= 0; ii--) {
		for (int j = 0; j <= ii; j++)  {
			for (int k = 0; k < lenxy; k++) 
				coefTermX(idx, k) = pow( x_corr[k]-xp, ii-j) * pow(y_corr[k]-yp, j);
			idx++; 
		} 	}
	idx = 0;
	for (int ii = degY; ii >= 0; ii--) {
		for (int j = 0; j <= ii; j++)  {
			for (int k = 0; k < lenxy; k++) 
				coefTermY(idx, k) = pow(x_corr[k]-xp, ii-j) * pow(y_corr[k]-yp, j);
			idx++; 	
		} 	}
	vector<T> x_dist = bicubicDistModel(paramsX, coefTermX) + xp; 
	vector<T> y_dist = bicubicDistModel(paramsY, coefTermY) + yp;

	return getParamsCorrection(x_corr, y_corr, x_dist, y_dist, degX, degY, xp, yp);

	//vector<T> x_dist_rad(x_dist-xp), y_dist_rad(y_dist-yp);
	//vector<T> x_corr_rad(x_corr-xp), y_corr_rad(y_corr-yp);
	//T norm_xy = 0;
	//for (int i = 0; i < lenxy; i++) norm_xy += x_dist_rad[i]*x_dist_rad[i] + y_dist_rad[i]*y_dist_rad[i];
	//norm_xy = std::sqrt(norm_xy);
	//x_dist_rad  /= norm_xy; y_dist_rad  /= norm_xy;

	//idx = 0;
	//for (int ii = degX; ii >= 0; ii--) {
	//	for (int j = 0; j <= ii; j++)  {
	//		for (int k = 0; k < lenxy; k++) 
	//			coefTermX(idx, k) = pow(x_dist_rad[k], ii-j) * pow(y_dist_rad[k], j);
	//		idx++; 
	//	} 	}
	//idx = 0;
	//for (int ii = degY; ii >= 0; ii--) {
	//	for (int j = 0; j <= ii; j++)  {
	//		for (int k = 0; k < lenxy; k++) 
	//			coefTermY(idx, k) = pow(x_dist_rad[k], ii-j) * pow(y_dist_rad[k], j);
	//		idx++; 
	//	} 	}

	//int sizexy = sizex+sizey;
	//matrix<T> coef_mat = matrix<T>::zeros(sizexy, sizexy);
	//vector<T> m = vector<T>::zeros(sizexy);
	//int count_i = 0;
	//for (int i = 0; i < sizex; i++) {
	//	int count_j = 0;
	//	for (int j = 0; j < sizex; j++) {
	//		coef_mat(count_i, count_j) = coefTermX.rowRef(i) * coefTermX.rowRef(j);
	//		count_j++; }
	//	m[count_i] = x_corr_rad * coefTermX.rowRef(i);
	//	count_i++; }

	//count_i = 0;
	//for (int i = 0; i < sizey; i++) {
	//	int count_j = 0;
	//	for (int j = 0; j < sizey; j++) {
	//		coef_mat(sizex+count_i, sizex+count_j) = coefTermY.rowRef(i) * coefTermY.rowRef(j);
	//		count_j++; }
	//	m[sizex+count_i] = y_corr_rad * coefTermY.rowRef(i);
	//	count_i++; 	}

	//matrix<T> normalization_mat1 = matrix<T>::zeros(sizexy, sizexy);
	//matrix<T> inv_normalization_mat1 = matrix<T>::zeros(sizexy, sizexy);
	//for (int i = 0; i < sizexy; i++) {
	//	if (i < sizex)
	//		normalization_mat1(i,i) = coef_mat(0, i);
	//	else
	//		normalization_mat1(i,i) = coef_mat(sizex, i); 
	//	inv_normalization_mat1(i,i) = 1 / normalization_mat1(i,i);
	//}

	//matrix<T> normalized_coef_mat = coef_mat * inv_normalization_mat1;
	//matrix<T> normalization_mat2 = matrix<T>::zeros(sizexy, sizexy);
	//for (int i = 0; i < sizexy; i++) normalization_mat2(i, i) = normalized_coef_mat(0, 0) / normalized_coef_mat(i, i);
	//normalized_coef_mat =  normalization_mat2*normalized_coef_mat;

	//vector<T> paramsInv(sizex+sizey);
	//solveLU(normalized_coef_mat, normalization_mat2 * m, paramsInv);
	//paramsInv = inv_normalization_mat1 * paramsInv;
	//vector<T> denorm_paramsInv(sizex+sizey);
	//vectorRef<T> denormX = denorm_paramsInv.copyRef(0, sizex-1);
	//vectorRef<T> denormY = denorm_paramsInv.copyRef(sizex, sizex+sizey-1);
	//denormalization(denormX, denormY,  // copyRef is not const?!
	//	paramsInv.copyRef(0, sizex),
	//	paramsInv.copyRef(sizex, sizex+sizey-1),
	//	norm_xy, norm_xy, degX, degY);
	//return denorm_paramsInv;
}

#endif
