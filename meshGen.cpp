#include "meshGen.h"


using namespace std;

void createLinspace(double* Vect, double mini, double maxi, int n){
    double delta = (maxi - mini)/double(n-1);
    for(int i = 0; i < n; i++){
        Vect[i] = mini + double(i)*delta;
        }
}


void heightFunction(int n, double* x, double a, double b, double h1, double* h){
    for (int i =0; i < n; i++){
        h[i] = (a * x[i] * x[i] + b * x[i] + h1);
    }
}

