#ifndef FFT_HARDCORD_H
#define FFT_HARDCORD_H
#include <vector>
using Complex = std::complex<double>;
#define FFT_HARDCORD_SIZE 4
inline void fft_hardcord(std::vector<Complex> &array,std::vector<Complex> &temp,int n,int sign){

int m = array.size()/n;
switch(m){
      case 1:
         break;
      case 2:
        for(int i=0;i<n;i++){
            double r0 = array[i*2+0].real() + array[i*2+1].real()
            ;
            double r1 = array[i*2+0].real() + (-array[i*2+1].real())
            ;
            double i0 = 
             + (array[i*2+0].imag()) + (array[i*2+1].imag());
            double i1 = 
             + (array[i*2+0].imag()) + (-array[i*2+1].imag());
            array[i*2+0].real(r0);
            array[i*2+0].imag(i0);
            array[i*2+1].real(r1);
            array[i*2+1].imag(i1);
        }
        break;
      case 3:
        for(int i=0;i<n;i++){
            double r0 = array[i*3+0].real() + array[i*3+1].real() + array[i*3+2].real()
            ;
            double r1 = array[i*3+0].real() + (-5.00000000000000e-01)*(array[i*3+1].real()) + (-5.00000000000000e-01)*(array[i*3+2].real())
             + (-sign)*(8.66025403784439e-01)*(array[i*3+1].imag()) + (-sign)*(-8.66025403784438e-01)*(array[i*3+2].imag());
            double r2 = array[i*3+0].real() + (-5.00000000000000e-01)*(array[i*3+1].real()) + (-5.00000000000000e-01)*(array[i*3+2].real())
             + (-sign)*(-8.66025403784438e-01)*(array[i*3+1].imag()) + (-sign)*(8.66025403784439e-01)*(array[i*3+2].imag());
            double i0 = 
             + (array[i*3+0].imag()) + (array[i*3+1].imag()) + (array[i*3+2].imag());
            double i1 = sign*(8.66025403784439e-01)*(array[i*3+1].real()) + sign*(-8.66025403784438e-01)*(array[i*3+2].real())
             + (array[i*3+0].imag()) + (-5.00000000000000e-01)*(array[i*3+1].imag()) + (-5.00000000000000e-01)*(array[i*3+2].imag());
            double i2 = sign*(-8.66025403784438e-01)*(array[i*3+1].real()) + sign*(8.66025403784439e-01)*(array[i*3+2].real())
             + (array[i*3+0].imag()) + (-5.00000000000000e-01)*(array[i*3+1].imag()) + (-5.00000000000000e-01)*(array[i*3+2].imag());
            array[i*3+0].real(r0);
            array[i*3+0].imag(i0);
            array[i*3+1].real(r1);
            array[i*3+1].imag(i1);
            array[i*3+2].real(r2);
            array[i*3+2].imag(i2);
        }
        break;
      case 4:
        for(int i=0;i<n;i++){
            double r0 = array[i*4+0].real() + array[i*4+1].real() + array[i*4+2].real() + array[i*4+3].real()
            ;
            double r1 = array[i*4+0].real() + (-array[i*4+2].real())
             + (-sign)*(array[i*4+1].imag()) + sign*(array[i*4+3].imag());
            double r2 = array[i*4+0].real() + (-array[i*4+1].real()) + array[i*4+2].real() + (-array[i*4+3].real())
            ;
            double r3 = array[i*4+0].real() + (-array[i*4+2].real())
             + sign*(array[i*4+1].imag()) + (-sign)*(array[i*4+3].imag());
            double i0 = 
             + (array[i*4+0].imag()) + (array[i*4+1].imag()) + (array[i*4+2].imag()) + (array[i*4+3].imag());
            double i1 = sign*(array[i*4+1].real()) + (-sign)*(array[i*4+3].real())
             + (array[i*4+0].imag()) + (-array[i*4+2].imag());
            double i2 = 
             + (array[i*4+0].imag()) + (-array[i*4+1].imag()) + (array[i*4+2].imag()) + (-array[i*4+3].imag());
            double i3 = (-sign)*(array[i*4+1].real()) + sign*(array[i*4+3].real())
             + (array[i*4+0].imag()) + (-array[i*4+2].imag());
            array[i*4+0].real(r0);
            array[i*4+0].imag(i0);
            array[i*4+1].real(r1);
            array[i*4+1].imag(i1);
            array[i*4+2].real(r2);
            array[i*4+2].imag(i2);
            array[i*4+3].real(r3);
            array[i*4+3].imag(i3);
        }
        break;
}

};
#endif