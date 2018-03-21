#include "checkPositiveDef.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "cblas.h"


using namespace std;

#define F77NAME(x) x##_
extern "C" {
    // lapack routine to do cholesky factorisation
    void F77NAME(dpotrf)(const char& uplo, const int& n, double* a,
                        const int& lda, int& info);
}


bool checkPositiveDefinite(int n, double* D){
    const int lda = n;
    int info = 0;
    // check if matrix is symmetric and positive definite
    // if cholesky factorisation sucessful --> symmetric and postive definite

    F77NAME(dpotrf)('U', n, D, lda, info);
    // if info = 0 --> sucessful factorisation

    if (info == 0){
        return true;
    }
    else{
        cout << "D is not symmetric and positive definite." << endl;
        return false;
    }
}

