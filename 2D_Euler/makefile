COMPILER=g++
# Override this using:
# make COMPILER=icx (to use the Intel Compiler)
# make COMPILER=clang++  (clang)
# make COMPILER=nvc++ (NVIDIA compiler)
OPTIONS=-O3 -ltbb -lstdc++ 
2deuler:
	${COMPILER} FVM.cpp ${OPTIONS} -c
	${COMPILER} 2D_Euler.cpp ${OPTIONS} -c
	${COMPILER} FVM.o 2D_Euler.o ${OPTIONS} -o 2D_Euler.exe

