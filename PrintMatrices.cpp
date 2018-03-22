/*
* This source prints out matrices and vectors
*
* Written by Trieu Ho
* 18/03/18
*/

#include "PrintMatrices.h"
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace std;

/**
* @brief Performs matrix print outs of a double type matrix.
* @param [in] n     number of matrix rows
* @param [in] k     number of matrix columns
* @param [in] H     printed matrix
*/
void PrintMatrix(int n, int k, double *H) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      cout << setw(1) << setprecision(6) << H[i * k + j] << " ";
    }
    cout << endl;
  }
}

/**
* @brief Performs matrix print outs of a integer type matrix.
* @param [in] n     number of matrix rows
* @param [in] k     number of matrix columns
* @param [in] H     printed matrix
*/
void PrintMatrixInt(int n, int k, int *H) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      cout << setw(1) << H[i * k + j] << " ";
    }
    cout << endl;
  }
}

/**
* @brief Performs matrix print outs of a boolean type matrix.
* @param [in] n     number of matrix rows
* @param [in] k     number of matrix columns
* @param [in] H     printed matrix
*/
void PrintMatrixBool(int n, int k, bool *H) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      cout << setw(1) << H[i * k + j] << " ";
    }
    cout << endl;
  }
}

/**
* @brief Performs vector print outs of a double type vector.
* @param [in] n     number of elements in vector
* @param [in] u     printed vector
*/
void PrintVector(int n, double *u) {
  for (int i = 0; i < n; ++i) {
    cout << u[i] << endl;
  }
}
