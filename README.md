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

## Performance

### First order, maximum optimization, g++

| Number of Cells | Number of Time Steps | Timing (Run A), s | Timing (Run B), s |
|----------------| ---------------| ----------------| ---------------| 
| 2048          | 3277           | 0.114             | 0.119        |
| 4096          | 6554           | 0.269             | 0.242        |
| 8192          | 13108           | 0.597            | 0.561        |
| 16384         | 26215         |  1.619            |  1.523        |
| 32768         | 52429         |  5.377           |   5.357      |
| 65536         | 104858        |  17.131           |  17.044      |

Table 1: Time required using 1st order with maximum optimization and g++ (13.3.0) with std::execution code with the maximum number of threads (28 on i7-14700HX)

### First order, maximum optimization, clang++

| Number of Cells | Number of Time Steps | Timing (Run A), s | Timing (Run B), s |
|----------------| ---------------| ----------------| ---------------| 
| 2048          | 3277           |  0.112          |   0.119     |
| 4096          | 6554           |  0.248          |   0.246     |
| 8192          | 13108           |  0.593          |  0.606       |
| 16384         | 26215         |   1.694          |   1.538s     |
| 32768         | 52429         |   5.847        |   5.923      |
| 65536         | 104858        |   17.483          |     16.918      |

Table 2: Time required using 1st order with maximum optimization and clang++ (18.1.3) with std::execution code with the maximum number of threads (28 on i7-14700HX)

### First order, maximum optimization, g++, OpenMP Code (Short-lived Threads)

| Number of Cells | Number of Time Steps | Timing (Run A), s | Timing (Run B), s |
|----------------| ---------------| ----------------| ---------------| 
| 2048          | 3277           |   0.073         |  0.086      |
| 4096          | 6554           |   0.169         |   0.169     |
| 8192          | 13108           |  0.369           |   0.372      |
| 16384         | 26215         |  0.955           |   1.020       |
| 32768         | 52429         |  6.937           |   6.879     |
| 65536         | 104858        |  28.161           |  27.609     |

Table 3: Time required using 1st order with maximum optimization and g++ (13.3.0) with Openmp code with the maximum number of threads (28 on i7-14700HX). The OpenMP implementation uses short-lived threads (as opposed to persistent threads). Performance seems comparible to that of persistent OpenMP threads.
