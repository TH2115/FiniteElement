#ifndef GLOBALCALC_H_INCLUDED
#define GLOBALCALC_H_INCLUDED

void FillMatrixTop(int n_elem, int *ElemNode, int *Top);

void FillMatrixglobDof(int n_elem, int *ElemNode, int *n_NodeDof,
                       int n_node_elem, int n_node, int *globDof, int &nDof);

void FillMatrixKp(int n_elem_upper, int n_elem_lower, int gaussorder,
                  int *ElemNode, double *Coord, int n_node_elem, int n_node,
                  int *globDof, int *n_NodeDof, double *K, double *D,
                  double t_p, int n_Dof, double *GP, double W, int *Top);

void FillMatrixK(int n_elem, int gaussorder, int *ElemNode, double *Coord,
                 int n_node_elem, int n_node, int *globDof, int *n_NodeDof,
                 double *K, double *D, double t_p, int n_Dof, double *GP,
                 double W);

void FillMatrixFluxNodes(int *fluxNodes, int *NodeTopo, int n_elx, int n_ely,
                         char side);

void FillMatrixNBC(double *n_bc, int nFluxNodes, int *fluxNodes, double q);

void FillMatrixf(double *f, int nbe, double *n_bc, double *coord, double *GP,
                 int gauss, double t_p, double W);

void FillMatrixBC(double *BC, int *NodeTopo, double T0, int n_elx, int n_ely,
                  int *TempNode, char side);

void MatrixPartition(double *K_EE, double *K_FF, double *K_EF, double *K,
                     int nDof, int sideDimT, int *TempNode, double *T,
                     double *f, double *T_E, double *f_F, bool *mask_E);

void GlobalReconstruct(double *T_E, double *T_F, bool *mask_E, int nDof,
                       double *T);

#endif // GLOBALCALC_H_INCLUDED
