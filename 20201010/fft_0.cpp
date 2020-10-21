#include <cmath>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <utility>

using namespace std;

// Private function prototypes
static size_t reverseBits(size_t val, int width);
void transform(vector<double> &real, vector<double> &imag);
void inverseTransform(vector<double> &real, vector<double> &imag);
void transformRadix2(vector<double> &real, vector<double> &imag);
void transformBluestein(vector<double> &real, vector<double> &imag);
static size_t reverseBits(size_t val, int width);
void convolve(const vector<double> &xvec, const vector<double> &yvec, vector<double> &outvec);
void convolve(
		const vector<double> &xreal, const vector<double> &ximag,
		const vector<double> &yreal, const vector<double> &yimag,
		vector<double> &outreal, vector<double> &outimag);

int main() {
    int N = 52;
    vector<double> real(N);
    vector<double> imag(N);
    for (int i = 0; i < N; i++) 
    {
        real[i] = sin(6*M_PI*i/N)+sin(9*M_PI*i/N)+sin(12*M_PI*i/N);
        imag[i] = 0;
    }
    transform(real, imag);
    cout << "\nFFT DONE" << endl;
    ofstream wFile("fft.csv");
    for (int i = 0; i < N; i++) {
        string csv = to_string(sqrt(real[i]*real[i]+imag[i]*imag[i])) + "\n";
        wFile << csv.c_str();    
    }
    wFile.close();
    return 0;
}

void transform(vector<double> &real, vector<double> &imag) {
	size_t n = real.size();
	if (n != imag.size())
		throw invalid_argument("Mismatched lengths");
	if (n == 0)
		return;
	else if ((n & (n - 1)) == 0)  // Is power of 2
		transformRadix2(real, imag);
	else  // More complicated algorithm for arbitrary sizes
		transformBluestein(real, imag);
}


void inverseTransform(vector<double> &real, vector<double> &imag) {
	transform(imag, real);
}


void transformRadix2(vector<double> &real, vector<double> &imag) {
	// Length variables
	size_t n = real.size();
	if (n != imag.size())
		throw invalid_argument("Mismatched lengths");
	int levels = 0;  // Compute levels = floor(log2(n))
	for (size_t temp = n; temp > 1U; temp >>= 1)
		levels++;
	if (static_cast<size_t>(1U) << levels != n)
		throw domain_error("Length is not a power of 2");
	
	// Trigonometric tables
	vector<double> cosTable(n / 2);
	vector<double> sinTable(n / 2);
	for (size_t i = 0; i < n / 2; i++) {
		cosTable[i] = cos(2 * M_PI * i / n);
		sinTable[i] = sin(2 * M_PI * i / n);
	}
	
	// Bit-reversed addressing permutation
	for (size_t i = 0; i < n; i++) {
		size_t j = reverseBits(i, levels);
		if (j > i) {
			swap(real[i], real[j]);
			swap(imag[i], imag[j]);
		}
	}
	
	// Cooley-Tukey decimation-in-time radix-2 FFT
	for (size_t size = 2; size <= n; size *= 2) {
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (size_t i = 0; i < n; i += size) {
			for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				size_t l = j + halfsize;
				double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
				double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
				real[l] = real[j] - tpre;
				imag[l] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
}


void transformBluestein(vector<double> &real, vector<double> &imag) {
	// Find a power-of-2 convolution length m such that m >= n * 2 + 1
	size_t n = real.size();
	if (n != imag.size())
		throw invalid_argument("Mismatched lengths");
	size_t m = 1;
	while (m / 2 <= n) {
		if (m > SIZE_MAX / 2)
			throw length_error("Vector too large");
		m *= 2;
	}
	
	// Trigonometric tables
	vector<double> cosTable(n), sinTable(n);
	for (size_t i = 0; i < n; i++) {
		uintmax_t temp = static_cast<uintmax_t>(i) * i;
		temp %= static_cast<uintmax_t>(n) * 2;
		double angle = M_PI * temp / n;
		cosTable[i] = cos(angle);
		sinTable[i] = sin(angle);
	}
	
	// Temporary vectors and preprocessing
	vector<double> areal(m), aimag(m);
	for (size_t i = 0; i < n; i++) {
		areal[i] =  real[i] * cosTable[i] + imag[i] * sinTable[i];
		aimag[i] = -real[i] * sinTable[i] + imag[i] * cosTable[i];
	}
	vector<double> breal(m), bimag(m);
	breal[0] = cosTable[0];
	bimag[0] = sinTable[0];
	for (size_t i = 1; i < n; i++) {
		breal[i] = breal[m - i] = cosTable[i];
		bimag[i] = bimag[m - i] = sinTable[i];
	}
	
	// Convolution
	vector<double> creal(m), cimag(m);
	convolve(areal, aimag, breal, bimag, creal, cimag);
	
	// Postprocessing
	for (size_t i = 0; i < n; i++) {
		real[i] =  creal[i] * cosTable[i] + cimag[i] * sinTable[i];
		imag[i] = -creal[i] * sinTable[i] + cimag[i] * cosTable[i];
	}
}


void convolve(const vector<double> &xvec, const vector<double> &yvec, vector<double> &outvec) {
	size_t n = xvec.size();
	if (n != yvec.size() || n != outvec.size())
		throw invalid_argument("Mismatched lengths");
	vector<double> outimag(n);
	convolve(xvec, vector<double>(n), yvec, vector<double>(n), outvec, outimag);
}


void convolve(
		const vector<double> &xreal, const vector<double> &ximag,
		const vector<double> &yreal, const vector<double> &yimag,
		vector<double> &outreal, vector<double> &outimag) {
	
	size_t n = xreal.size();
	if (n != ximag.size() || n != yreal.size() || n != yimag.size()
			|| n != outreal.size() || n != outimag.size())
		throw invalid_argument("Mismatched lengths");
	
	vector<double> xr = xreal;
	vector<double> xi = ximag;
	vector<double> yr = yreal;
	vector<double> yi = yimag;
	transform(xr, xi);
	transform(yr, yi);
	
	for (size_t i = 0; i < n; i++) {
		double temp = xr[i] * yr[i] - xi[i] * yi[i];
		xi[i] = xi[i] * yr[i] + xr[i] * yi[i];
		xr[i] = temp;
	}
	inverseTransform(xr, xi);
	
	for (size_t i = 0; i < n; i++) {  // Scaling (because this FFT implementation omits it)
		outreal[i] = xr[i] / n;
		outimag[i] = xi[i] / n;
	}
}


static size_t reverseBits(size_t val, int width) {
	size_t result = 0;
	for (int i = 0; i < width; i++, val >>= 1)
		result = (result << 1) | (val & 1U);
	return result;
}

