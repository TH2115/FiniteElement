/*
* This is the main source file for solving a finite element heat problem
* This runs in parallel with two processes using MPI
*
* Written by Trieu Ho
* 18/03/18
*/

#include "ElementCalc.h"
#include "MatrixOp.h"
#include "PrintMatrices.h"
#include "VTKO.h"
#include "cblas.h"
#include "checkPositiveDef.h"
#include "globalCalcp.h"
#include "lapackRoutines.h"
#include "meshGen.h"
#include "triple.h"
#include <boost/program_options.hpp>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/time.h>

// An alias to reduce typing
using namespace std;
namespace po = boost::program_options;

#define F77NAME(x) x##_
extern "C" {
// lapack routine to do cholesky factorisation
void F77NAME(dpotrf)(const char &uplo, const int &n, double *a, const int &lda,
                     int &info);
// lapack for LU facorisation
void F77NAME(dgetrf)(const int &m, const int &n, double *a, const int &lda,
                     int &info);

// lapack to solve Ax = B
void F77NAME(dgesv)(const int &n, const int &nrhs, double *a, const int &lda,
                    int *ipiv, double *b, const int &ldb, int &info);
}

int main(int argc, char *argv[]) {
  struct timeval start, endt;
  gettimeofday(&start, NULL);
  //////// using Boost library to get command line arguments//////
  //////// Default values are those for caes 1///////////////
  po::options_description opts("Test case value definitions.");
  opts.add_options()("case", po::value<int>()->default_value(1), "Case number")(
      "a", po::value<double>()->default_value(0.0),
      "Geometric parameter describing the height of the plate, a*x*x + b*x + "
      "h1.")("h1", po::value<double>()->default_value(1.0),
             "Height of the left side of the plate [m].")(
      "h2", po::value<double>()->default_value(1.0),
      "Height of the right side of the plate [m].")(
      "L", po::value<double>()->default_value(2.0), "Length of the plate [m].")(
      "t_p", po::value<double>()->default_value(0.2),
      "Thickness of plate [m].")("Kxx",
                                 po::value<double>()->default_value(250.0),
                                 "Thermal conductivity [W/mk].")(
      "Kyy", po::value<double>()->default_value(250.0),
      "Thermal conductivity [W/mk].")("Kxy",
                                      po::value<double>()->default_value(0.0),
                                      "Thermal conductivity [W/mk].")(
      "Nelx", po::value<int>()->default_value(10),
      "Number of elements in the x direction.")(
      "Nely", po::value<int>()->default_value(5),
      "Number of elements in the y direction.")(
      "T_BC_side", po::value<char>()->default_value('L'),
      "Side in which the temperature is applied.")(
      "T_BC", po::value<double>()->default_value(10.0),
      "Boundary temperature.")("q_BC_side",
                               po::value<char>()->default_value('R'),
                               "Side in which the heat flux is applied.")(
      "q_BC", po::value<double>()->default_value(2500.0),
      "Boundary heat flux.")("help", "Print help message.");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, opts), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << "Test case value definition." << endl;
    cout << opts << endl;
    return 0;
  }
  int CaseNum = vm["case"].as<int>();
  double a = vm["a"].as<double>();
  double h1 = vm["h1"].as<double>();
  double h2 = vm["h2"].as<double>();
  double L = vm["L"].as<double>();
  double t_p = vm["t_p"].as<double>();
  double Kxx = vm["Kxx"].as<double>();
  double Kyy = vm["Kyy"].as<double>();
  double Kxy = vm["Kxy"].as<double>();
  int N_elx = vm["Nelx"].as<int>();
  int N_ely = vm["Nely"].as<int>();
  char T_BC_side = vm["T_BC_side"].as<char>();
  double T_BC = vm["T_BC"].as<double>();
  char q_BC_side = vm["q_BC_side"].as<char>();
  double q_BC = vm["q_BC"].as<double>();

  /////// creating matrix D /////////////
  double *D = new double[2 * 2]();

  D[0] = Kxx;
  D[1] = Kxy;
  D[2] = Kxy;
  D[3] = Kyy;
  //// creating D_test matrix to test if positive definite
  double *D_test = new double[2 * 2]();
  bool pass;
  D_test[0] = Kxx;
  D_test[1] = Kxy;
  D_test[2] = Kxy;
  D_test[3] = Kyy;
  pass = checkPositiveDefinite(2, D_test);
  delete[] D_test;

  if (pass == false) {
    cout << "D matrix is not positive definite!" << endl;
  }

  double b; // constant in polynomial describing beam height
  b = -a * L + (h2 - h1) / L;

  ////////////////FE VALUES///////////////
  int gaussorder = 2;
  int N_node_elem = 4; // number of nodes in each element
  int N_NodeDof[4] = {1, 1, 1,
                      1}; // number of DoF per node (1 = temperature only)
  int N_eDof = 0;
  for (int i = 0; i < 4; i++) {
    N_eDof = N_eDof + N_NodeDof[i]; // total number of DoF per element
  }

  int N_elem = N_elx * N_ely;             // total number of elements
  int N_node = (N_elx + 1) * (N_ely + 1); // total number of nodes

  ///////////////////calculation of nodel coordinate matrix /////
  double *x = new double[N_elx + 1]();
  createLinspace(x, 0, L,
                 N_elx +
                     1); // fill x with a linearly space numbers between 0 and L
  double *h = new double[N_elx + 1]();
  heightFunction(N_elx + 1, x, a, b, h1, h);

  double *Y = new double[(N_ely + 1) * (N_elx + 1)]();

  FillMatrixY(N_elx, N_ely, Y, h);

  // cordinates of nodes //
  double *coord = new double[N_node * 2]();
  FillMatrixCoord(N_elx, N_ely, x, Y, coord);
  delete[] Y;

  /////////////Calculation of topology matrix NodeTopo//////////

  int *NodeTopo = new int[(N_elx + 1) * (N_ely + 1)]();
  FillMatrixTopo(N_elx, N_ely, NodeTopo);

  ////////// calculation of topology matrix ElemNode //////
  int *ElemNode = new int[N_elem * 5]();
  FillMatrixElemNode(N_elx, N_ely, ElemNode, NodeTopo);

  double *ElemX = new double[N_elem * N_node_elem]();
  double *ElemY = new double[N_elem * N_node_elem]();
  FillMatrixElemXY(N_elem, N_node_elem, ElemNode, coord, ElemX, ElemY);

  int *globDof = new int[N_node * 2]();
  int *Top = new int[N_elem * N_node_elem]();
  int nDof = 0;
  FillMatrixTop(N_elem, ElemNode, Top);

  FillMatrixglobDof(N_elem, ElemNode, N_NodeDof, N_node_elem, N_node, globDof,
                    nDof);

  /////////////////////// Assembly of global stiffness matrix K ////////

  double *K = new double[nDof * nDof]();
  ///////// gauss points and weights /////////////
  // weights
  double W = 1.0;
  // gauss points
  double GP[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};

  FillMatrixK(N_elem, gaussorder, ElemNode, coord, N_node_elem, N_node, globDof,
              N_NodeDof, K, D, t_p, nDof, GP, W);
  // PrintMatrix(nDof,nDof,K);

  /////////////////////// Begin solving///////////////

  // compute nodal boundary flux vector

  ///////Define edges --- RIGHT EDGE//////////////
  ////////////////////// USE N_ely for right/left, N_elx for top/bottom
  int sideDimFlux;

  if (q_BC_side == 'R' || q_BC_side == 'L') {
    sideDimFlux = N_ely + 1;
  } else if (q_BC_side == 'T' || q_BC_side == 'B') {
    sideDimFlux = N_elx + 1;
  } else {
    cout << "input error" << endl;
  }

  int *fluxNodes = new int[sideDimFlux]();
  int nFluxNodes = sideDimFlux;
  FillMatrixFluxNodes(fluxNodes, NodeTopo, N_elx, N_ely, q_BC_side);
  ////// Defining load /////////////

  double *n_bc = new double[4 * (nFluxNodes - 1)]();

  FillMatrixNBC(n_bc, nFluxNodes, fluxNodes, q_BC);

  int nbe = nFluxNodes - 1; // number of elements with flux load

  double f[nDof] = {0};

  FillMatrixf(f, nbe, n_bc, coord, GP, gaussorder, t_p, W);

  delete[] x;
  delete[] h;
  delete[] Top;

  ///////////// Appply boundary condiitons ////////
  int sideDimT;

  if (T_BC_side == 'R' || T_BC_side == 'L') {
    sideDimT = N_ely + 1;
  } else if (T_BC_side == 'T' || T_BC_side == 'B') {
    sideDimT = N_elx + 1;
  }

  double *BC = new double[2 * (sideDimT)]();
  int *TempNode = new int[sideDimT];
  FillMatrixBC(BC, NodeTopo, T_BC, N_elx, N_ely, TempNode, T_BC_side);

  //////////////// Assembling global force vector /////////////

  double OrgDof[nDof] = {0}; // original DoF number
  double T[nDof] = {0};      // initialis nodal temperature vector
  int rDof = nDof;           // reduced DoF
  int ind[sideDimT];

  for (int i = 0; i < sideDimT; i++) {
    ind[i] = BC[2 * i];
    OrgDof[ind[i]] = -1;
    T[ind[i]] = BC[2 * i + 1];
  }

  rDof = rDof - sideDimT;
  int RedDof[rDof] = {0};
  int counter1 = 0;

  for (int i = 0; i < nDof; i++) {
    if (OrgDof[i] == 0) {
      OrgDof[i] = counter1;
      RedDof[counter1] = i;
      counter1 = counter1 + 1;
    }
  }

  delete[] BC;
  delete[] ElemX;
  delete[] ElemY;

  ////////////////Partition matrix ////////////
  double *K_EE = new double[sideDimT * sideDimT]();
  double *K_FF = new double[(nDof - sideDimT) * (nDof - sideDimT)]();
  double *K_EF = new double[(nDof - sideDimT) * sideDimT]();
  double *T_E = new double[sideDimT]();
  double *f_F = new double[(nDof - sideDimT)]();
  bool mask_E[nDof] = {false};
  double *T_F = new double[nDof - sideDimT]();
  double *f_E = new double[sideDimT]();

  MatrixPartition(K_EE, K_FF, K_EF, K, nDof, sideDimT, TempNode, T, f, T_E, f_F,
                  mask_E);

  delete[] TempNode;

  ////////////////////// Solve for D_F ///////////////

  cblas_dgemv(CblasRowMajor, CblasTrans, sideDimT, nDof - sideDimT, 1.0, K_EF,
              nDof - sideDimT, T_E, 1, 0.0, T_F, 1);

  MatrixScale(T_F, 1, nDof - sideDimT, -1);

  MatrixAdd(T_F, 1, nDof - sideDimT, f_F);

  int *ipiv1 = new int[nDof - sideDimT];
  int info1;
  F77NAME(dgesv)
  (nDof - sideDimT, 1, K_FF, nDof - sideDimT, ipiv1, T_F, nDof - sideDimT,
   info1);

  delete[] ipiv1;

  /////////// reconstruction of global tempterature ////////////////////

  GlobalReconstruct(T_E, T_F, mask_E, nDof, T);

  ////////////// compute reaction f_E /////////////

  double *K_EEdotT_E = new double[sideDimT]();
  cblas_dgemv(CblasRowMajor, CblasNoTrans, sideDimT, sideDimT, 1.0, K_EE,
              sideDimT, T_E, 1, 0.0, K_EEdotT_E, 1);

  double *K_EFdotT_F = new double[sideDimT]();
  cblas_dgemv(CblasRowMajor, CblasNoTrans, sideDimT, nDof - sideDimT, 1.0, K_EF,
              nDof - sideDimT, T_F, 1, 0.0, K_EFdotT_F, 1);

  MatrixAdd(f_E, 1, sideDimT, K_EEdotT_E);

  MatrixAdd(f_E, 1, sideDimT, K_EFdotT_F);

  GlobalReconstruct(f_E, f_F, mask_E, nDof, f);

  string filename = "case";
  string CaseN = to_string(CaseNum);

  filename.append(CaseN);

  WriteVTK(coord, ElemNode, N_node, N_elem, T, filename);

  delete[] K_EFdotT_F;
  delete[] K_EEdotT_E;
  delete[] f_E;
  delete[] f_F;
  delete[] T_E;
  delete[] T_F;
  delete[] K_EE;
  delete[] K_EF;
  delete[] K_FF;
  delete[] ElemNode;
  delete[] NodeTopo;

  gettimeofday(&endt, NULL);

  cout << "Time: "
       << ((endt.tv_sec - start.tv_sec) * 1000000u + endt.tv_usec -
           start.tv_usec) /
              1.e6
       << endl;
}
