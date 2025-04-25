# std-execution-fvm-cfd

A Finite Volume Method (FVM) solver for the compressible Euler Equations, parallelized using std::execution policies.

## Background

Conventional High Performance Computing (HPC) implementation of shared memory parallel computing makes use of the OpenMP libraries for sharing compute load across multiple physical cores in a system. This has allowed for very fine grained control of how variables in memory are shared, how work is shared, how synchronization occurs - with succesful outcomes and very performant code.

However, the introduction of the std:execution in C++17 makes is possible for us to parallelize these numerical algorithms without the use of OpenMP. This repository demonstrates this in practice for a vector-split Finite Volume Method - specifically, the SHLL method.

## Building and Running

All codes are built using:

```bash
make
```

You can check the code is operating correctly on your system using the simple demonstration provided:

```bash
./Simple_Vector_Demo/test.exe
```

which should return:

```bash
Number of available threads: 28
Value of result[0] is 2
Value of result[1] is 4
Value of result[2] is 6
Value of result[3] is 8
Value of result[4] is 10
Value of result[0] is 2
Value of result[1] is 8
Value of result[2] is 18
Value of result[3] is 32
Value of result[4] is 50
```

**NOTE** The number of available threads will depend on your own system.

