/*
* This source fills in all of the neccesary matrices to find the local stiffness
* matrices of each elemnent
*
* Written by Trieu Ho
* 18/03/18
*/

#include "ElementCalc.h"
#include "meshGen.h"

using namespace std;

/**
* @brief Fills in the height distribution matrix
* @param [in] n_elx     number of elements in the x domain
* @param [in] n_ely     number of elements in the y domain
* @param [in] h         height of elements matrix
* @param [out] Y        height distribution of elements
*/
void FillMatrixY(int n_elx, int n_ely, double *Y, double *h) {
  double *tempVect = new double[n_ely + 1];
  for (int i = 0; i < n_elx + 1; i++) {
    createLinspace(tempVect, -h[i] / double(2), h[i] / double(2), n_ely + 1);
    for (int j = 0; j < n_ely + 1; j++) {
      Y[i * (n_ely + 1) + j] = tempVect[j];
    }
  }
  delete[] tempVect;
}

/**
* @brief Fills in the coordinate matrix
* @param [in] n_elx     number of elements in the x domain
* @param [in] n_ely     number of elements in the y domain
* @param [in] x         position of elements in the x domain
* @param [in] Y         height distribution of elements
* @param [out] coord    coordinate matrix
*/
void FillMatrixCoord(int n_elx, int n_ely, double *x, double *Y,
                     double *coord) {
  for (int i = 0; i < n_elx + 1; i++) {
    for (int j = 0; j < n_ely + 1; j++) {
      coord[2 * i * (n_ely + 1) + (2 * j)] = x[i];
      coord[2 * i * (n_ely + 1) + (2 * j + 1)] = Y[i * (n_ely + 1) + j];
    }
  }
}

/**
* @brief Fills in the topology matrix
* @param [in] n_elx         number of elements in the x domain
* @param [in] n_ely         number of elements in the y domain
* @param [out] NodeTopo     topology matrix
*/
void FillMatrixTopo(int n_elx, int n_ely, int *NodeTopo) {
  for (int j = 0; j < n_ely + 1; j++) {
    for (int i = 0; i < n_elx + 1; i++) {
      NodeTopo[i + j * (n_elx + 1)] = j + i * (n_ely + 1);
    }
  }
}

/**
* @brief Fills in the element nodes matrix
* @param [in] n_elx         number of elements in the x domain
* @param [in] n_ely         number of elements in the y domain
* @param [in] NodeTopo      topology matrix
* @param [out] ElemNode     element node matrix
*/
void FillMatrixElemNode(int n_elx, int n_ely, int *ElemNode, int *NodeTopo) {
  int elemnr = 0;
  for (int i = 0; i < n_elx; i++) {
    for (int j = 0; j < n_ely; j++) {
      ElemNode[elemnr * 5 + 0] = elemnr;
      ElemNode[elemnr * 5 + 4] = NodeTopo[(1 + j) * (n_elx + 1) + i];
      ElemNode[elemnr * 5 + 3] = NodeTopo[(1 + j) * (n_elx + 1) + (i + 1)];
      ElemNode[elemnr * 5 + 2] = NodeTopo[(j) * (n_elx + 1) + (i + 1)];
      ElemNode[elemnr * 5 + 1] = NodeTopo[(j) * (n_elx + 1) + i];
      elemnr = elemnr + 1;
    }
  }
}

/**
* @brief Fills in the element x and y coorindates
* @param [in] n_elem        number of elements in domain
* @param [in] n_node_elem   number of nodes in each element
* @param [in] ElemNode      element node matrix
* @param [in] coord         coordinate matrix
* @param [out] ElemX        element x nodal coordinates
* @param [out] ElemY        element y nodal coordinates
*/
void FillMatrixElemXY(int n_elem, int n_node_elem, int *ElemNode, double *coord,
                      double *ElemX, double *ElemY) {
  int eNodes;
  double eCoordx;
  double eCoordy;
  for (int i = 0; i < n_elem; i++) {
    for (int j = 1; j < 5; j++) {

      eNodes = ElemNode[i * 5 + j]; // element nodes
      eCoordx = coord[2 * (eNodes)];
      eCoordy = coord[2 * eNodes + 1];

      ElemX[i * 4 + (j - 1)] = eCoordx;
      ElemY[i * 4 + (j - 1)] = eCoordy;
    }
  }
}
