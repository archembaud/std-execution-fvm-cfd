COMPILER=g++
# Override this using:
# make COMPILER=icx (to use the Intel Compiler)
# make COMPILER=clang++  (clang)
# make COMPILER=nvc++ (NVIDIA compiler)
OPTIONS=-O3 -ltbb -lstdc++ -fopenmp 
2domp:
	${COMPILER} FVM.cpp ${OPTIONS} -c
	${COMPILER} 2D_OpenMP.cpp ${OPTIONS} -c
	${COMPILER} FVM.o 2D_OpenMP.o ${OPTIONS} -o 2D_OpenMP.exe
clean:
	rm -f *.o *.exe

2dclang:
	clang++ -O3 -ltbb -lstdc++ -fopenmp=libomp -I/usr/lib/gcc/x86_64-linux-gnu/13/include/ FVM.cpp -c
	clang++ -O3 -ltbb -lstdc++ -fopenmp=libomp -I/usr/lib/gcc/x86_64-linux-gnu/13/include/ 2D_OpenMP.cpp -c
	clang++ FVM.o 2D_OpenMP.o -O3 -ltbb -lstdc++ -fopenmp=libomp -I/usr/lib/gcc/x86_64-linux-gnu/13/include/ -o 2D_OpenMP.exe

