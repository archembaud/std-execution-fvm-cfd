COMPILER=g++
OPTIONS=-O3 -ltbb

all: simple 1deuler 2deuler omp
	${info Compile Complete}

simple:
	${COMPILER} ./Simple_Vector_Demo/demo.cpp ${OPTIONS} -o ./Simple_Vector_Demo/demo.exe

1deuler:
	${info Compiling 1D Euler Solver}
	$(MAKE) -C 1D_Euler $@
	${info Done}

2deuler:
	${info Compiling 2D Euler Solver}
	$(MAKE) -C 2D_Euler $@
	${info Done}

omp:
	${info Compiling 1D Euler Solver (OpenMP)}
	$(MAKE) -C 1D_OpenMP $@
	${info Done}
