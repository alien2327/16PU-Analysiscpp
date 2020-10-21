#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>

using namespace std;

double (*mat_sum(double (*matA)[3], double (*matB)[3]))[3]{
    double (*result)[3] = new double[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i][j] = matA[i][j] + matB[i][j];
        }
    }
    return result;
}

double (*mat_sub(double (*matA)[3], double (*matB)[3]))[3]{
    double (*result)[3] = new double[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i][j] = matA[i][j] + matB[i][j];
        }
    }
    return result;
}

double (*mat_scalar(double a, double (*mat)[3]))[3]{
    double (*result)[3] = new double[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i][j] = a * mat[i][j];
        }
    }
    return result;
}

double (*mat_dot(double (*matA)[3], double (*matB)[3]))[3]{
    double (*result)[3] = new double[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                result[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }
    return result;
}

int main() {
    double (*a)[3] = new double[3][3];
    double (*b)[3] = new double[3][3];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i][j] = pow(2.1 * (i + 2), 1.4 * (j));
            b[i][j] = 1.2 * (i + j) * (j - i);
        }
    }
    cout << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << b[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    double (*res_dot)[3] = mat_dot(a, b);
    double (*res_sum)[3] = mat_sum(a, b);
    double (*res_sub)[3] = mat_sub(a, b);
    double (*res_scalar)[3] = mat_scalar(3, b);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << res_dot[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << res_sum[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << res_sub[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << res_scalar[i][j] << " ";
        }
        cout << endl;
    }
    return 0;
}