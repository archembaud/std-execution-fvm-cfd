COMPILER=g++
OPTIONS=-O3

all: simple
	${info Compile Complete}

simple:
	${COMPILER} ./Simple_Vector_Demo/demo.cpp -ltbb ${OPTIONS} -o ./Simple_Vector_Demo/demo.exe
