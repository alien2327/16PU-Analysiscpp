#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include "FftRealPair.hpp"

using namespace std;
using cdbl = complex<double>;
 
const double PI = acos(-1);
 
inline unsigned bitreverse(const unsigned n, const unsigned k) {
    unsigned r, i;
    for (r = 0, i = 0; i < k; ++i)
        r |= ((n >> i) & 1) << (k - i - 1);
    return r;
}

void fft(vector<cdbl> &a, int size=0, bool is_reverse=false) {
    unsigned n, k;
    if (size == 0) {
        n = a.size();
        k = __builtin_ctz(n);
    } else {
        n = size;
        k = __builtin_ctz(n);
    }
    unsigned s, i, j;
    for (i = 0; i < n; i++) {
        j = bitreverse(i, k);
        if (i < j)
            swap(a[i], a[j]);
    }
    for (s = 2; s <= n; s *= 2) {
        double t = 2*PI/s * (is_reverse? -1 : 1);
        cdbl ws(cos(t), sin(t));
        for (i = 0; i < n; i += s) {
            cdbl w(1);
            for (j = 0; j < s/2; j++) {
                cdbl tmp = a[i + j + s/2] * w;
                a[i + j + s/2] = a[i + j] - tmp;
                a[i + j] += tmp;
                w *= ws;
            }
        }
    }
    if (is_reverse)
        for (i = 0; i < n; i++)
            a[i] /= n;
}

void do_fft(vector<cdbl> &a) {
    int n = a.size();
    fft(a, n-1);
    fft(a);
}

int main() {
    int N = 4096;
    vector<cdbl> v(N);
    for (int i = 0; i < N; i++) v[i] = cdbl(10*cos(2*M_PI*i/3*N)+5*sin(2*M_PI*i/N), 0);
    fft(v);
    cout << "\nFFT DONE" << endl;
    ofstream wFile("fft.csv");
    for (int i = 0; i < N; i++) {
        string csv = to_string(abs(v[i])/N) + "\n";
        wFile << csv.c_str();    
    }
    wFile.close();
    return 0;
}