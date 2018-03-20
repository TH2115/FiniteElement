#ifndef MESHGEN_H_INCLUDED
#define MESHGEN_H_INCLUDED

void createLinspace(double* Vect, double mini, double maxi, int n);

void heightFunction(int n, double* x, double a, double b, double h1, double* h);


#endif // MESHGEN_H_INCLUDED
