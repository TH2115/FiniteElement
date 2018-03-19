#include "PrintMatrices.h"

// Print a matrix
void PrintMatrix(int n, int k, double* H) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            cout << setw(5) << H[i*k+j] << " ";
        }
        cout << endl;
    }
}

// Print a matrix ith integer
void PrintMatrixInt(int n, int k, int* H) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            cout << setw(5) << H[i*k+j] << " ";
        }
        cout << endl;
    }
}

// Print a matrix
void PrintMatrixBool(int n, int k, bool* H) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < k; ++j) {
            cout << setw(5) << H[i*k+j] << " ";
        }
        cout << endl;
    }
}

// Print a vector
void PrintVector(int n, double* u) {
    for (int i = 0; i < n; ++i) {
        cout << u[i] << endl;
    }
}
