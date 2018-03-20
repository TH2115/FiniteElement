#ifndef MATRIXOP_H_INCLUDED
#define MATRIXOP_H_INCLUDED

void MatrixScale(double* A, int n, int k , double c);

void MatrixAdd(double* A, int n, int k , double* B);

void MatrixInv2by2(double* A, double* Ainv, double& detJ);


#endif // MATRIXOP_H_INCLUDED
