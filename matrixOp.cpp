/*
* This source contains functions that allow matrix operations
*
* Written by Trieu Ho
* 18/03/18
*/

#include "MatrixOp.h"
#include <iostream>

using namespace std;

/**
* @brief scales a matrix by a constant
* @param [out] A       matrix of interest
* @param [in] n        number of rows
* @param [in] k        number of columns
* @param [in] c        constant scalar
*/
void MatrixScale(double *A, int n, int k, double c) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      A[i * k + j] = A[i * k + j] * c;
    }
  }
}

/**
* @brief adds two matrices of same dimensions together
* @param [out] A       matrix of interest
* @param [in] n        number of rows
* @param [in] k        number of columns
* @param [in] B        secondary matrix
*/
void MatrixAdd(double *A, int n, int k, double *B) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      A[i * k + j] = A[i * k + j] + B[i * k + j];
    }
  }
}

/**
* @brief Inverts a 2x2 matrix and finds the determinant
* @param [in] A         matrix of interest
* @param [out] Ain      inverse of matrix
* @param [out] detJ     determinant of matrix
*/
void MatrixInv2by2(double *A, double *Ainv, double &detJ) {

  detJ = A[0] * A[3] - A[1] * A[2];
  if (detJ != 0) {
    Ainv[0] = (1 / detJ) * A[3];
    Ainv[1] = -(1 / detJ) * A[1];
    Ainv[2] = -(1 / detJ) * A[2];
    Ainv[3] = (1 / detJ) * A[0];
  } else {
    cout << "Singular Jacobian matrix" << endl;
  }
}
