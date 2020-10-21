#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>

using namespace std;

double mean(double *set) {
    double res = 0.0;
    for (int i = 0; i < 10; i++) {
        res += set[i];
    }
    res /= 10;
    return res;
}

double variance(double *set) {
    double res = 0.0;
    double dev[10];
    for (int i = 0; i < 10; i++) {
        dev[i] = pow(set[i]-mean(set), 2);
    }
    res = mean(dev);
    return res;
}

double stadev(double *set) {
    double res = sqrt(variance(set));
    return res;
}

double covar(double *setA, double *setB) {
    double res = 0.0;
    double cov[10];
    for (int i = 0; i < 10; i++) {
        cov[i] = (setA[i] - mean(setA)) * (setB[i] - mean(setB));
    }
    res = mean(cov);
    return res;
}

double corel(double *setA, double *setB) {
    double res = covar(setA, setB)/(stadev(setA)*stadev(setB));
    return res;
}

int main() {
    double a[10] = {4.5,12.3,4.0,6.1,2.6,3.4,12.0,1.0,0.0,6.0};
    double b[10] = {3.7,2.0,4.0,61.0,51.0,2.3,1.2,1.0,8.8,0.0};

    double mean_a = mean(a);
    double var_a = variance(a);
    double std_a = stadev(a);

    double mean_b = mean(b);
    double var_b = variance(b);
    double std_b = stadev(b);

    cout << "Mean value of array a : " << mean_a << endl;
    cout << "Variance value of array a : " << var_a << endl;
    cout << "Std value of array a : " << std_a << endl;

    cout << "Mean value of array b : " << mean_b << endl;
    cout << "Variance value of array b : " << var_b << endl;
    cout << "Std value of array b : " << std_b << endl;

    double cov_ab = covar(a, b);
    double cor_ab = corel(a, b);

    cout << "Covariance of a and b : " << cov_ab << endl;
    cout << "Corelation of a and b : " << cor_ab << endl;
    return 0;
}