#Instructions on how to run the finite element solver

1. Go to the directory ~/FINAL in the terminal.
2. To run the serial code, first type 'make myprog' followed by 'make Cx' where x = 1,2 or 3.
3. To run the parallel code, first type 'make myprogp' followed by 'make Cxp' where x = 1,2 or 3.

following the above steps will result in casex.vtk ouput files, where x = 1,2,3 is for the
the repsective case number in the serial code and x = 4,5,6 is for parallel.

These .vtk files can be viewed using paraview.