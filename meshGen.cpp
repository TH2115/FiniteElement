
/*
* This source generates the mesh and create the height function
*
* Written by Trieu Ho
* 18/03/18
*/

#include "meshGen.h"
using namespace std;

/**
* @brief create a linearly spaces array
* @param [out] Vect     linearly spaced vector
* @param [in] mini      lower bound of vector
* @param [in] maxi      upper bound of vector
* @param [in] n         number of elements
*/
void createLinspace(double *Vect, double mini, double maxi, int n) {
  double delta = (maxi - mini) / double(n - 1);
  for (int i = 0; i < n; i++) {
    Vect[i] = mini + double(i) * delta;
  }
}

/**
* @brief height function which describes the plate geometry h =(a*x*x + b*x +
* h1)/2
* @param [in] n         number of elements
* @param [in] x         linearly spaced array of x positions
* @param [in] a         geometric constant
* @param [in] b         geometric constant
* @param [in] h1        plate side length 1
* @param [out] h     plate height
*/
void heightFunction(int n, double *x, double a, double b, double h1,
                    double *h) {
  for (int i = 0; i < n; i++) {
    h[i] = (a * x[i] * x[i] + b * x[i] + h1);
  }
}
