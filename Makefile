default: myprog

all: myprog

main.o: main.cpp 
	g++ -Wall -fexceptions -g  -c main.cpp -o main.o

ElementCalc.o: ElementCalc.cpp ElementCalc.h 
	g++ -Wall -fexceptions -g  -c ElementCalc.cpp -o ElementCalc.o

getNumbers.o: getNumbers.cpp getNumbers.h
	g++ -Wall -fexceptions -g  -c getNumbers.cpp -o getNumbers.o

globalCalc.o: globalCalc.cpp globalCalc.h
	g++ -Wall -fexceptions -g  -c globalCalc.cpp -o globalCalc.o

matrixOp.o: matrixOp.cpp MatrixOp.h
	g++ -Wall -fexceptions -g  -c matrixOp.cpp -o matrixOp.o

meshGen.o: meshGen.cpp meshGen.h
	g++ -Wall -fexceptions -g  -c meshGen.cpp -o meshGen.o
 
PrintMatrices.o: PrintMatrices.cpp PrintMatrices.h
	g++ -Wall -fexceptions -g  -c PrintMatrices.cpp -o PrintMatrices.o

VTKO.o: VTKO.cpp VTKO.h
	g++ -Wall -fexceptions -g  -c VTKO.cpp -o VTKO.o

myprog: ElementCalc.o getNumbers.o globalCalc.o main.o matrixOp.o meshGen.o PrintMatrices.o VTKO.o
	g++  -o myprog ElementCalc.o getNumbers.o globalCalc.o main.o matrixOp.o meshGen.o PrintMatrices.o VTKO.o  -lboost_program_options -lblas -llapack 




globalCalcp.o: globalCalcp.cpp globalCalcp.h
	g++ -Wall -fexceptions -g  -c globalCalcp.cpp -o globalCalcp.o

mainp.o: mainp.cpp 
	mpicxx -Wall -fexceptions -g  -c mainp.cpp -o mainp.o

myprogp: ElementCalc.o getNumbers.o globalCalcp.o mainp.o matrixOp.o meshGen.o PrintMatrices.o VTKO.o
	mpicxx  -o myprogp ElementCalc.o getNumbers.o globalCalcp.o mainp.o matrixOp.o meshGen.o PrintMatrices.o VTKO.o  -lboost_program_options -lblas -llapack 



C1:
	 ./myprog --case=1 --a=0 --h1=1.0 --h2=1.0 --L=2.0 --t_p=0.2 --Kxx=250.0 --Kyy=250.0 --Kxy=0.0 --Nelx=10 --Nely=5 --T_BC_side='L' --T_BC=10.0 --q_BC_side='R' --q_BC=2500.0

C2:
	./myprog --case=2 --a=0 --h1=1.0 --h2=1.0 --L=2.0 --t_p=0.2 --Kxx=250.0 --Kyy=250.0 --Kxy=0.0 --Nelx=10 --Nely=5 --T_BC_side='B' --T_BC=10.0 --q_BC_side='T' --q_BC=2500.0

C3:	
	./myprog --case=3 --a=0.25 --h1=1.0 --h2=1.3 --L=3.0 --t_p=0.2 --Kxx=250.0 --Kyy=250.0 --Kxy=0.0 --Nelx=15 --Nely=8 --T_BC_side='L' --T_BC=-20.0 --q_BC_side='B' --q_BC=-5000.0

C1p:
	mpiexec -np 2 ./myprogp --case=4 --a=0 --h1=1.0 --h2=1.0 --L=2.0 --t_p=0.2 --Kxx=250.0 --Kyy=250.0 --Kxy=0.0 --Nelx=10 --Nely=5 --T_BC_side='L' --T_BC=10.0 --q_BC_side='R' --q_BC=2500.0

C2p:
	mpiexec -np 2 ./myprogp --case=5 --a=0 --h1=1.0 --h2=1.0 --L=2.0 --t_p=0.2 --Kxx=250.0 --Kyy=250.0 --Kxy=0.0 --Nelx=10 --Nely=5 --T_BC_side='B' --T_BC=10.0 --q_BC_side='T' --q_BC=2500.0

C3p:	
	mpiexec -np 2 ./myprogp --case=6 --a=0.25 --h1=1.0 --h2=1.3 --L=3.0 --t_p=0.2 --Kxx=250.0 --Kyy=250.0 --Kxy=0.0 --Nelx=15 --Nely=8 --T_BC_side='L' --T_BC=-20.0 --q_BC_side='B' --q_BC=-5000.0



.PHONY: clean run
clean:
	rm -f *.o -f *.vtk
