/*
* This source checks if a given square matrix is positive definite using
* cholesky factorisation
*
* Written by Trieu Ho
* 18/03/18
*/

#include "checkPositiveDef.h"
#include "cblas.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

#define F77NAME(x) x##_
extern "C" {
// lapack routine to do cholesky factorisation
void F77NAME(dpotrf)(const char &uplo, const int &n, double *a, const int &lda,
                     int &info);
}

/**
* @brief Checks is a matrix is symmetric positive definite
* @param [in] n         length of square matrix
* @param [in] D         square matrix of interest
* @param [out]          returns true if D is symmetric positive definite
*/
bool checkPositiveDefinite(int n, double *D) {
  const int lda = n;
  int info = 0;
  // check if matrix is symmetric and positive definite
  // if cholesky factorisation sucessful --> symmetric and postive definite

  F77NAME(dpotrf)('U', n, D, lda, info);
  // if info = 0 --> sucessful factorisation

  if (info == 0) {
    return true;
  } else {
    cout << "D is not symmetric and positive definite." << endl;
    return false;
  }
}
