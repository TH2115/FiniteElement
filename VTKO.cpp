#include "VTKO.h"
#include<iostream>
#include<iomanip>
#include<fstream>
#include <string>




void WriteVTK(double* Coord, int* ElemNodes, int n_node,int n_elem, double* T, std::string filename) {
	filename.append(".vtk");

	std::ofstream vOut(filename, std::fstream::out | std::fstream::trunc);

	vOut << "# vtk DataFile Version 4.0" << std::endl;
	vOut << "vtk output" << std::endl;
	vOut << "ASCII" << std::endl;
	vOut << "DATASET UNSTRUCTURED_GRID" << std::endl;



	vOut << "POINTS " << n_node << " double" << std::endl;

	double x,y;
	for(int i = 0; i < n_node; i++) {
		x = Coord[2*i];
		y = -Coord[2*i + 1];
		vOut << x << " " << y << " " << "0.0 ";
	}

	vOut << std::endl;

	vOut << "CELLS " << (int)n_elem << " " << int(n_elem*5) << std::endl;
	int n;

	for(int i = 0; i < n_elem; i++) {
		vOut << (int)4 << " ";

		for (int j = 0; j < 4;j++) {
			n = ElemNodes[5*i + j +1];
			vOut << n << " ";
		}
		vOut << std::endl;
	}

	vOut << "CELL_TYPES " << n_elem << std::endl;
	for(int i = 0; i < n_elem; i++) {
		vOut << (int)9 << std::endl;
	}
	vOut << "POINT_DATA " << n_node << std::endl;
	vOut << "FIELD FieldData 1" <<  std::endl;
	vOut << "disp 1 " << n_node << " double" << std::endl;
	for(int i = 0; i < n_node; i++) {
		vOut << T[i] << " ";
	}
	vOut << std::endl;
	vOut.close();
}
