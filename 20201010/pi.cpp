#include <iostream>
#include <cmath>

using namespace std; 

void Iterate(int iter);
void DecimalApprox(int approx);

int main() {
    char input;
    bool loop = true;
    int num = 0;

    cout << "\nThis Program calculates the number Pi.\n\n";

    do {
        cout << "Would you like to calculate it through I>teration or ";
        cout << "A>pproximation: ";
        cin >> input;

        switch(input) {
            case 'a':
            case 'A':
                loop = false;
                input = 'a';
                break;

            case 'i':
            case 'I':
                loop = false;
                input = 'i';
                break;

            default:
                cout << "Incorrect Selection. Try again... \n";
                break;
        }
    }while(loop == true);

    if(input == 'i') {
        cout << "How many iterations?: ";
        cin >> num;

        Iterate(num);
    }

    if(input == 'a') {
        cout << "How many decimals of approximations?: ";
        cin >> num;

        DecimalApprox(num);
    }

    return 0;
}

void Iterate(int num) {
    double piVal = 0.0;
    int sign = 1;
    for (int i = 0; i < num; i++) {
        piVal += sign / ((2.0 * i) + 1.0);
        sign *= -1;
    }
    piVal = 4*piVal;
    cout << piVal << endl;
}

void DecimalApprox(int num) {
    bool loop = true;
    double piVal = 0.0;
    double prePi = piVal;
    int acc = 1;

    do {
        double piVal = 0.0;
        int sign = 1;

        for (int i = 0; i < acc; i++) {
            piVal += sign / ((2.0 * i) + 1.0);
            sign *= -1;
        }
        piVal = 4*piVal;
        if (abs(piVal - prePi) <= pow(10, -1.0 * num)) {
            cout << piVal << endl;
            loop = false;
        }
        prePi = piVal;
        acc += 1;

    } while (loop == true);
}