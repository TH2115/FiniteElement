#ifndef LAPACKROUTINES_H_INCLUDED
#define LAPACKROUTINES_H_INCLUDED


#define F77NAME(x) x##_
extern "C" {
    // lapack routine to do cholesky factorisation
    void F77NAME(dpotrf)(const char& uplo, const int& n, double* a,
                        const int& lda, int& info);
    //lapack for LU facorisation
    void F77NAME(dgetrf)(const int& m, const int& n, double* a, const int& lda, int& info);

    //lapack to solve Ax = B
    void F77NAME(dgesv)(const int& n, const int& nrhs, double* a, const int& lda, int* ipiv, double* b, const int& ldb, int& info);
}




#endif // LAPACKROUTINES_H_INCLUDED
