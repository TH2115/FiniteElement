#include "globalCalcp.h"
#include "MatrixOp.h"
#include "triple.h"
#include <map>
#include "cblas.h"
#include <iostream>
#include <cmath>
#include "PrintMatrices.h"


using namespace std;




// filling Top matrix
void FillMatrixTop(int n_elem, int* ElemNode, int* Top) {

    for(int i = 0; i < n_elem ; i++){
        for(int j = 1; j < 5; j++){
            Top[4*i + j-1] = ElemNode[5*i + j];
        }
    }
}


// filling in the globDof matrix
void FillMatrixglobDof(int n_elem, int* ElemNode, int* n_NodeDof, int n_node_elem, int n_node, int* globDof, int &nDof){

    int* globNodes = new int[n_elem * 4]; //global node numbers of element nodes
    int nNode;
    for (int i = 0; i < n_elem; i++){
        globNodes[4*i + 0] = ElemNode[5*i + 1];
        globNodes[4*i + 1] = ElemNode[5*i + 2];
        globNodes[4*i + 2] = ElemNode[5*i + 3];
        globNodes[4*i + 3] = ElemNode[5*i + 4];
        //loop over element nodes and create first column
        for (int k = 0; k < n_node_elem; k++){
            nNode = ElemNode[(n_node_elem + 1)*i + k + 1];
            // if the already existing ndof of the present node is less
            // than the present elements ndof then replace the ndof for that node
            if (globDof[2*nNode + 0] < n_NodeDof[k]){
                globDof[2*nNode + 0] = n_NodeDof[k];
            }
        }
    }
    int eDof;
    for (int j = 0; j < n_node; j++){
        eDof = globDof[2*j];
        globDof[2*j + 1] = nDof;
        nDof = nDof +1;
    }
}



void FillMatrixK(int n_elem_upper,int n_elem_lower, int gaussorder, int* ElemNode, double* Coord, int n_node_elem,int n_node, int* globDof, int* n_NodeDof, double* K, double* D, double t_p, int n_Dof, double* GP, double W, int* Top){
    map<triple,double>data;
    int gauss = gaussorder; // gauss order
    int eNodes[1]; // node elements
    double eCoordx;
    double eCoordy;
    double eCoord[8]; //node coordinates
    int* gDofNode = new int();
    double eta;
    double xi;



    double* J = new double[2 * 2];
    double DetJ;//for determinant of jacobian
    double* Jinv = new double[2 * 2];
    double* C = new double[2 * 4];
    double* B = new double[2 * 4];
    double* E = new double[4 * 4];
    double scale;
   // double* transB = new double[2 * 4];
    triple point;


    //// assemble K /////

    for(int i = n_elem_lower ; i < n_elem_upper; i++){
    //data for element i
        double Ke[n_node_elem*n_node_elem] = {0};
        int gDof[n_node_elem*(n_node_elem-1)] = {0};// used to contruct scatter matrix

        for (int j = 1; j < 5; j++){
            eNodes[j-1] = ElemNode[i*5 + j]; //element nodes
            eCoordx = Coord[2*(eNodes[j-1])];
            eCoordy = Coord[2*eNodes[j-1] + 1];
            eCoord[2*(j-1)] = eCoordx;
            eCoord[2*(j-1) + 1] = eCoordy;
        }


        for (int j = 0; j < n_node_elem; j++){
            //global dof for node j
            for(int k = 1; k < n_NodeDof[j] + 1  ; k++){
                gDofNode[k - 1] = globDof[2*eNodes[j] + k];
                gDof[j*(n_NodeDof[j])+ k - 1] = gDofNode[k-1];
            }
        }


        // local stiffness matrix K_e found by gauss integration //

        for (int n = 0; n < gauss; n++){
            for (int m = 0; m < gauss; m++){
                eta = GP[n];
                xi = GP[m];

                // gradient of shape funcions
                double GN[8] = {-0.25 *(1 - eta), 0.25*(1 - eta), 0.25*(1 + eta), -0.25*(1+eta), -0.25*(1 - xi), -0.25*(1 + xi), 0.25*(1+ xi), 0.25*(1 - xi)} ;

                //jacobian
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 4, 1.0, GN, 4, eCoord, 2, 0.0, J, 2);

                //determinant of jacobian and inverse
                MatrixInv2by2(J, Jinv, DetJ);

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 4, 2, 1.0, Jinv , 2, GN, 4, 0.0, B, 4);

                // Ke = Ke + B'*D*B*th*DetJ*Wi*Wj
                cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 4, 2, 2, 1.0, B , 4, D, 2, 0.0, C, 2);

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, 4, 2, 1.0, C , 2, B, 4, 0.0, E, 4);

                scale = (t_p)*(DetJ)*(W)*(W);

                MatrixScale( E, 4, 4 , scale);

                MatrixAdd(Ke, 4 , 4 , E);

            }
        }
        //inserting the global stiffness matix into the global system
        //stiffness matrix

        for (int j = 0; j < n_node_elem; j++){ // column
            int TopX = Top[i*n_node_elem+j]; // find topology location
            for (int k = 0; k < n_node_elem; k++){ // row
                int TopY = Top[i*n_node_elem + k];// find topology location
                K[TopY*n_node+TopX] += Ke[j*n_node_elem+k]; // replaced nDof with nnode
            }
        }

    }
    ////clearing memory
    delete gDofNode;
    delete[] J;
    delete[] Jinv;
    delete[] C;
    delete[] B;
    delete[] E;

}



void FillMatrixFluxNodes(int* fluxNodes, int* NodeTopo,int n_elx, int n_ely, char side ){

    if (side == 'R'){
        for(int i = 0; i < n_ely + 1; i++){
            fluxNodes[i] = NodeTopo[i*(n_elx + 1) + (n_elx)];
        }
    }

    else if (side == 'L'){
        for(int i = 0; i < n_ely + 1; i++){
            fluxNodes[i] = NodeTopo[i*(n_elx + 1)];
        }
    }

    else if (side == 'B'){
        for(int j = 0; j < n_elx + 1; j++){
            fluxNodes[j] = NodeTopo[j];
        }
    }

    else if (side =='T'){
        for(int j = 0; j < n_elx + 1; j++){
            fluxNodes[j] = NodeTopo[j + (n_elx + 1) * (n_ely)];
        }
    }

    else{
        cout << "Invalid side Input." << endl;
    }

}


void FillMatrixNBC(double* n_bc, int nFluxNodes, int* fluxNodes, double q){

    for (int i = 0; i < nFluxNodes - 1; i++ ){
        n_bc[i] = double(fluxNodes[i]);
        n_bc[nFluxNodes- 1 + i] = (fluxNodes[i + 1]);
        n_bc[2*(nFluxNodes-1) + i] = q;
        n_bc[3*(nFluxNodes-1) + i] = q;

    }

}


void FillMatrixf(double* f, int nbe, double* n_bc, double* coord, double* GP, int gauss, double t_p, double W){

    int node1;
    int node2;
    double n_bce[gauss];
    double x1;
    double y1;
    double x2;
    double y2;
    double N[gauss];
    double flux;



    for (int i = 0; i < nbe; i++){
        double fq[2] = {0,0}; // inital the nodal source vector
        double leng = 0;
        double detJ =0;
        double xi =0;
        node1 = n_bc[i]; //node 1
        node2 = n_bc[nbe + i]; // node 2
        n_bce[0] = n_bc[nbe*2 + i];
        n_bce[1] = n_bc[nbe*3 + i]; //flux value at an edge


        x1 = coord[2*node1]; //x coord of the first node
        y1 = coord[2*node1 + 1]; //y coord of the first node
        x2 = coord[2*node2]; // x coord of the 2nd node
        y2 = coord[2*node2 + 1]; // y coord of the 2nd node


        leng = leng +  sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
        detJ = leng /2.0 ; //1D jacobian


        for (int j = 0; j < gauss; j++){

            xi = GP[j];
            //1D shape functions in parent domain
            N[0] = 0.5 * (1.0 - xi);
            N[1] = 0.5 * (1.0 + xi);

            flux = cblas_ddot(2, N,1,n_bce,1); //calculating flux



            //calculating the nodal flux
            // fq = fq + w[i]* N * flux * detJ * t_p
            MatrixScale(N, 1, 2 , W);


            MatrixScale(N, 2, 1 , flux);

            MatrixScale(N, 2, 1 , detJ);


            MatrixScale(N, 2, 1 , t_p);

            MatrixAdd(fq, 2, 1 , N);

        }

        MatrixScale(fq, 2, 1 , -1.0); // fq = -fq

        f[node1] = f[node1] + fq[0];
        f[node2] = f[node2] + fq[1];
    }
}


void FillMatrixBC(double* BC, int* NodeTopo,  double T0, int n_elx, int n_ely,int* TempNode, char side ){

    if (side == 'R'){
        for(int i = 0; i < n_ely + 1; i++){
            TempNode[i] = NodeTopo[i*(n_elx + 1) + (n_elx)];
            BC[2*i] = TempNode[i];
            BC[2*i + 1] = T0;
        }
    }

    else if (side == 'L'){
        for(int i = 0; i < n_ely + 1; i++){
            TempNode[i]= NodeTopo[i*(n_elx + 1)];
            BC[2*i] = TempNode[i];
            BC[2*i + 1] = T0;
        }
    }

    else if (side == 'B'){
        for(int j = 0; j < n_elx + 1; j++){
            TempNode[j]= NodeTopo[j];
            BC[2*j] = TempNode[j];
            BC[2*j + 1] = T0;
        }
    }

    else if (side =='T'){
        for(int j = 0; j < n_elx + 1; j++){
            TempNode[j] = NodeTopo[j + (n_elx+1) * (n_ely)];
            BC[2*j] = TempNode[j];
            BC[2*j + 1] = T0;
        }
    }

    else{
        cout << "Invalid side Input." << endl;
    }


}


void MatrixPartition(double* K_EE, double* K_FF, double* K_EF, double* K, int nDof, int sideDimT, int* TempNode, double* T, double* f, double* T_E, double* f_F, bool* mask_E){




    for(int i = 0; i < sideDimT ; i++){
        mask_E[TempNode[i]] = true;
    }

    int counter1 = 0;
    int counter2 = 0;

    for (int i = 0; i < nDof; i++){
        if (mask_E[i] == true){
            T_E[counter1] = T[i];
            counter1 = counter1 + 1;
        }
        else {
            f_F[counter2] = f[i];
            counter2 = counter2 + 1;
        }
    }

    //PrintMatrix(1,6,T_E);
    //PrintMatrix(1,nDof - 6, f_F);
    bool mask_E1;
    bool mask_E2;
    int counter3 = 0;
    int counter4 = 0;
    int counter5 = 0;

    for (int i = 0; i < nDof; i++){
        for(int j = 0; j < nDof; j++){

            mask_E1 = mask_E[i];
            mask_E2 = mask_E[j];

            if ((mask_E1 == true) & (mask_E2 == true)){
                K_EE[counter3] = K[i*nDof + j];
                counter3 = counter3+1;

            }

            else if ((mask_E1 == false) & (mask_E2 == false)){
                K_FF[counter4] = K[i*nDof + j];
                counter4 = counter4 + 1;

            }

            else if ((mask_E1 == true) & (mask_E2 == false)){
                K_EF[counter5] = K[i*nDof + j];
                counter5 = counter5 + 1;
            }
        }
    }

}

void GlobalReconstruct(double* T_E, double* T_F, bool* mask_E, int nDof, double* T){

    bool mask_E1;
    int counter1 = 0;
    int counter2 = 0;

    for (int i = 0; i < nDof; i++){
        mask_E1 = mask_E[i];
        if (mask_E1 == true){
            T[i] = T_E[counter1];
            counter1 = counter1 + 1;
        }
        else {
            T[i] = T_F[counter2];
            counter2 = counter2 + 1;
        }
    }
}




