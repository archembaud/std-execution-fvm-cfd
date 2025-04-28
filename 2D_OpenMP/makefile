COMPILER=g++
OPTIONS=-O3 -ltbb -fopenmp
2domp:
	${COMPILER} FVM.cpp ${OPTIONS} -c
	${COMPILER} 2D_OpenMP.cpp ${OPTIONS} -c
	${COMPILER} FVM.o 2D_OpenMP.o ${OPTIONS} -o 2D_OpenMP.exe

