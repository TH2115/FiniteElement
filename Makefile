default: myprog

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

myprog: ElementCalc.o getNumbers.o globalCalc.o main.o matrixOp.o meshGen.o PrintMatrices.o 
	g++  -o myprog ElementCalc.o getNumbers.o globalCalc.o main.o matrixOp.o meshGen.o PrintMatrices.o   -lboost_program_options -lblas -llapack 



.PHONY: clean run
run:
	./myprog --a=0 --h1=1.0 --h2=1.0 --t_p=0.2 --Kxx=250.0 --Kyy=250.0 --Kxy=0.0 --Nelx=10 --Nely=5 --T_BC_side='L' --T_BC=10.0 --q_BC_side='R' --q_BC=2500.0

clean:
	rm -f *.o
