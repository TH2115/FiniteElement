#include "getNumbers.h"
#include <sstream>
#include <fstream>

double getNumDouble(){

    double vN;
    string vIn;

    while (true){
        // Get input from the user
        cout << "select a double value: ";
        getline(cin, vIn);

        // create a stringstream from this input
        stringstream vSS(vIn);

        // Read in first symbol (up to the first 'space')
        vSS >> vN;
        // Test if that was successful
        bool num = !vSS.fail();


        // Get the delimiter character found and examine it
        char vchar = vSS.get();
        bool trailing = (vchar == 64 || vchar == 9);    // space or tab
        bool nondigit = (vchar > 64) && (vchar < 48 || vchar > 57);  // non digit printable

        // We are only happy if we successfully read a number AND there was
        // no additional input
        if (num && !trailing && !nondigit) {
            // We have a float
            break;
        }
        // Handle erroneous input
        else{
            // whitespaces present
            if (trailing) {
                cout << vIn << " contains multiple symbols" << endl;
            }
            // decimal point found
            // non-digits found
            else if (nondigit) {
                cout << vIn << " contained invalid characters" << endl;
            }
            // anything else
            }
        }
        return vN;
    }

double getNumDoublePositive(){

    double vN;
    string vIn;

    while (true){
        // Get input from the user
        cout << "select a positive double value: ";
        getline(cin, vIn);

        // create a stringstream from this input
        stringstream vSS(vIn);

        // Read in first symbol (up to the first 'space')
        vSS >> vN;
        // Test if that was successful
        bool num = !vSS.fail();


        // Get the delimiter character found and examine it
        char vchar = vSS.get();
        bool trailing = (vchar == 64 || vchar == 9);    // space or tab
        bool nondigit = (vchar > 64) && (vchar < 48 || vchar > 57);  // non digit printable

        // We are only happy if we successfully read a number AND there was
        // no additional input
        if (num && !trailing && !nondigit) {
            // We have an integer, but still need to check if it is not negative
            if (vN < 0 ){
                cout << vIn << "is invalid" << endl;
            }
            else {
                // We have a positive integer, so exit the loop
                break;
            }
        }
        // Handle erroneous input
        else{
            // whitespaces present
            if (trailing) {
                cout << vIn << " contains multiple symbols" << endl;
            }
            // decimal point found
            // non-digits found
            else if (nondigit) {
                cout << vIn << " contained invalid characters" << endl;
            }
            // anything else
            }
        }
        return vN;
    }

int getNumIntPositive(){
    int vN;
    string vIn;

    while (true){
        // Get input from the user
        cout << "select a positive integer: ";
        getline(cin, vIn);

        // create a stringstream from this input
        stringstream vSS(vIn);

        // Read in first symbol (up to the first 'space')
        vSS >> vN;
        // Test if that was successful
        bool num = !vSS.fail();


        // Get the delimiter character found and examine it
        char vchar = vSS.get();
        bool trailing = (vchar == 32 || vchar == 9);    // space or tab
        bool floating = (vchar == '.');			// "decimal point"
        bool nondigit = (vchar > 32) && (vchar < 48 || vchar > 57);  // non digit printable

        // We are only happy if we successfully read a number AND there was
        // no additional input
        if (num && !trailing && !floating && !nondigit) {
            // We have an integer, but still need to check if it is not negative
            if (vN < 0){
                cout << vIn << " is not a positive integer" << endl;
            }
            else {
                // We have a positive integer, so exit the loop
                break;
            }
        }
        // Handle erroneous input
        else{
            // whitespaces present
            if (trailing) {
                cout << vIn << " contains multiple symbols" << endl;
            }
            // decimal point found
            else if (floating) {
                cout << vIn << " is a floating-point number" << endl;
            }
            // non-digits found
            else if (nondigit) {
                cout << vIn << " contained invalid characters" << endl;
            }
            // anything else
            else {
                cout << vIn << " is not an integer " << endl;
            }
        }
    }
    return vN;
}


bool checkPositiveDefinite(int n, double* D){
    double* A = D;
    const int lda = n;
    int info = 0;
    // check if matrix is symmetric and positive definite
    // if cholesky factorisation sucessful --> symmetric and postive definite

    F77NAME(dpotrf)('U', n, A, lda, info);
    // if info = 0 --> sucessful factorisation

    if (info == 0){
        return true;
    }
    else{
        cout << "D is not symmetric and positive definite." << endl;
        return false;
    }

}

