#ifndef ELEMENTCALC_H_INCLUDED
#define ELEMENTCALC_H_INCLUDED

void FillMatrixY(int n_elx, int n_ely, double *Y, double *h);

void FillMatrixCoord(int n_elx, int n_ely, double *x, double *Y, double *coord);

void FillMatrixTopo(int n_elx, int n_ely, int *NodeTopo);

void FillMatrixElemNode(int n_elx, int n_ely, int *ElemNode, int *NodeTopo);

void FillMatrixElemXY(int n_elem, int n_node_elem, int *ElemNode, double *coord,
                      double *ElemX, double *ElemY);

#endif // ELEMENTCALC_H_INCLUDED
