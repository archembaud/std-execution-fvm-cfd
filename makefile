COMPILER=g++
OPTIONS=-O3 -ltbb

all: simple 1deuler
	${info Compile Complete}

simple:
	${COMPILER} ./Simple_Vector_Demo/demo.cpp ${OPTIONS} -o ./Simple_Vector_Demo/demo.exe

1deuler:
	${info Compiling 1D Euler Solver}
	$(MAKE) -C 1D_Euler $@
	${info Done}
