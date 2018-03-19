#ifndef MATRIXOP_H_INCLUDED
#define MATRIXOP_H_INCLUDED

MatrixScale(double* A, int n, int k , double c);

MatrixAdd(double* A, int n, int k , double* B);

MatrixInv2by2(double* A, double* Ainv, double& detJ);


#endif // MATRIXOP_H_INCLUDED
