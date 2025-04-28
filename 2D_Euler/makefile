COMPILER=clang++
OPTIONS=-O3 -ltbb 
2deuler:
	${COMPILER} FVM.cpp ${OPTIONS} -c
	${COMPILER} 2D_Euler.cpp ${OPTIONS} -c
	${COMPILER} FVM.o 2D_Euler.o ${OPTIONS} -o 2D_Euler.exe

