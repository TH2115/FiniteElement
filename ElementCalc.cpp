#include <ElementCalc.h>




// filling  matrix Y
void FillMatrixY(int n_elx, int n_ely, double* Y, double* h) {
    double* tempVect = new double[n_ely+1];
    for(int i = 0; i < n_elx+1 ; i++){
        createLinspace(tempVect ,-h[i]/double(2), h[i]/double(2), n_ely + 1);
        for(int j = 0; j < n_ely+1; j++){
            Y[i*(n_ely+1) + j] = tempVect[j];
        }
    }
    delete[] tempVect;
}

// filling  coordinate matrix
void FillMatrixCoord(int n_elx, int n_ely, double* x, double* Y, double* coord) {

    for(int i = 0; i < n_elx+1 ; i++){
        for(int j = 0; j < n_ely+1; j++){
            coord[2*i*(n_ely+1)+(2*j)] = x[i];
            coord[2*i*(n_ely+1)+(2*j+1)] = Y[i*(n_ely+1) + j];
        }
    }
}

// filling topology matrix
void FillMatrixTopo(int n_elx, int n_ely, int* NodeTopo) {

    for(int j = 0; j < n_ely+1; j++){
        for(int i = 0; i < n_elx+1 ; i++){
            NodeTopo[i + j*(n_elx+1)] = j + i*(n_ely+1);
        }
    }
}

// filling element node matrix
void FillMatrixElemNode(int n_elx, int n_ely, int* ElemNode, int* NodeTopo) {

    int elemnr = 0;
    for(int i = 0; i < n_elx ; i++){
        for(int j = 0; j < n_ely; j++){
            ElemNode[elemnr*5 + 0] = elemnr;
            ElemNode[elemnr*5 + 4] = NodeTopo[ (1+j)*(n_elx+1) + i];
            ElemNode[elemnr*5 + 3] = NodeTopo[(1+j)*(n_elx+1) + (i + 1) ];
            ElemNode[elemnr*5 + 2] = NodeTopo[(j)*(n_elx+1) + (i + 1) ];
            ElemNode[elemnr*5 + 1] = NodeTopo[(j)*(n_elx+1) + i];
            elemnr = elemnr + 1;

        }
    }
}


// filling element node matrix
void FillMatrixElemXY(int n_elem, int n_node_elem, int* ElemNode, double* coord, double* ElemX, double* ElemY ) {

//    int* eNodes = new int[n_node_elem];
//    double* eCoord = new double[n_node_elem * 2];
    int eNodes;
    double eCoordx;
    double eCoordy;
    for (int i = 0; i < n_elem; i++){
        for (int j = 1; j < 5; j++){

            eNodes = ElemNode[i*5 + j]; //element nodes
            eCoordx = coord[2*(eNodes)];
            eCoordy = coord[2*eNodes + 1];

            ElemX[i*4 + (j - 1)] = eCoordx;
            ElemY[i*4 + (j - 1)] = eCoordy;

        }
    }
}


