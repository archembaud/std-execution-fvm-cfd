COMPILER=g++
# Override this using:
# make COMPILER=icx (to use the Intel Compiler)
# make COMPILER=clang++  (clang)
# make COMPILER=nvc++ (NVIDIA compiler)
OPTIONS=-O3 -ltbb -lstdc++ 
1deuler:
	${COMPILER} FVM.cpp ${OPTIONS} -c
	${COMPILER} 1D_Euler.cpp ${OPTIONS} -c
	${COMPILER} FVM.o 1D_Euler.o ${OPTIONS} -o 1D_Euler.exe

clean:
	rm -f *.o *.exe

