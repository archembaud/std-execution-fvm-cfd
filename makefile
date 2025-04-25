COMPILER=g++

all: simple
	${info Compile Complete}

simple:
	${COMPILER} ./Simple_Vector_Demo/demo.cpp -ltbb -o ./Simple_Vector_Demo/test.exe
