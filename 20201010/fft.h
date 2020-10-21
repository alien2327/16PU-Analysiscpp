#ifndef FFT_H
#define FFT_H
#include <iostream>
#include <complex>
#include <vector>
#include <list>
#include <fft_hardcord.h>
#define FFT_PRIME_THRESHOLD 256
using Complex = std::complex<double>;

template<typename T> void fft(T &array, int size);
template<typename T> void ifft(T &array, int size);
template<typename T> void fft2d(T &array, int row, int column);
template<typename T> void ifft2d(T &array, int row, int column);
template<typename T> void dft(T &array, int size);
template<typename T> void idft(T &array, int size);

void fft(std::vector<Complex> &array);
void ifft(std::vector<Complex> &array);
void dft(std::vector<Complex> &array);
void idft(std::vector<Complex> &array);

/*internal functions*/
template<typename T>
void dft(
    T &array,
    std::vector<Complex> &temp,
    std::vector<Complex> &omega,
    int size,
    int n
);

template<typename T> void fft2(
    T &array,
    std::vector<Complex> &temp,
    std::vector<Complex> &omega,
    int size
);

template<typename T> void fft_prime(T &array,int size,int n,int sign);
template<typename T> void fft(T &array, int size, int sign);
template<typename T> void fft2d(T &array, int row, int column, int sign);

double fft_cos(int n, int m);
double fft_sin(int n, int m);
int primePrimitiveRoot(int p);
int powerMod(int a, int b, int m);
int eulerPhi(int n);
int primitiveRoot(int n);
int primitiveRoot(int n, int phi);

std::vector<std::vector<int>> factorInteger(int size);
std::vector<std::vector<int>> tally(std::vector<int> list);
void print(std::vector<std::vector<int>> &list);

void fft_pf(std::vector<Complex> &array);
template<typename T> void fft_pf(T &array, int size);
template<typename T> void fft_pf(T &array, int size, int sign);
template<typename T> void fft_pf(
    T &array,
    std::vector<Complex> temp,
    std::vector<Complex> omega,
    int size,
    int sign,
    std::vector<std::vector<int>> list,
    int index
);

void fft_rec(std::vector<Complex> &array);
template<typename T> void fft_rec(T &array, int size);
template<typename T> void fft_rec(T &array, int size, int sign);
template<typename T> void fft_rec(
    T &array,
    std::vector<Complex> &temp,
    std::vector<Complex> &omega,
    int size,
    int n, 
    int sign
);
#endif