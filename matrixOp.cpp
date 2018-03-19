#include <MatrixOp.h>


void MatrixScale(double* A, int n, int k , double c){
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            A[i*k+j] =  A[i*k+j] * c;
        }
    }
}


void MatrixAdd(double* A, int n, int k , double* B){
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            A[i*k+j] =  A[i*k+j] + B[i*k+j];
        }
    }
}


void MatrixInv2by2(double* A, double* Ainv, double& detJ){

    detJ = A[0] * A[3] - A[1] * A[2];
    if (detJ != 0){
    Ainv[0] = (1/detJ) * A[3];
    Ainv[1] = -(1/detJ) * A[1];
    Ainv[2] = -(1/detJ) * A[2];
    Ainv[3] = (1/detJ) * A[0];
    }
    else {
    cout << "Singular Jacobian matrix" << endl;
    }
}

