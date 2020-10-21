#include <iostream>
#include <time.h>

using namespace std;

int main() {
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [20];

    time (&rawtime);
    timeinfo = localtime (&rawtime);

    strftime (buffer,20,"%G_%m_%e_%H_%M_%S",timeinfo);
    puts (buffer);
    return 0;
}