#ifndef VTKO_H_INCLUDED
#define VTKO_H_INCLUDED
#include <string>

void WriteVTK(double *Coord, int *ElemNodes, int n_node, int n_elem, double *T,
              std::string filename);

#endif // VTKO_H_INCLUDED
