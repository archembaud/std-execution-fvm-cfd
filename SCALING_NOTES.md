# Weak Scaling Tests

Here is the test schedule:

	1 thread: 512 x 512 cells
	2 threads: 724 x 724 cells
	4 threads: 1024 x 1024 cells
	8 threads: 1448 x 1448 cells

Metrics I used to use to evaluate performance:

Time / cell / timestep (seconds / cell * timestep)

Recommended metrics


## About the CPUs being testing

### Intel Core i7-14700HX

Total Cores: 20 (8 P-cores, or Performance-cores, 12 E-cores, or Efficient Cores)






## Run on Intel(R) Core(TM) i7-14700HX

### std::execution::par_unseq

Computing the number of threads here is non-trivial as the number of threads, scheduling strategy, and SIMD width are not exposed via any standard API.
These can be estimated, however, using the stategy demonstrated in Thread_Counter, and implemented in the First Run computation of conserved quantities. 
Seeing as the number of threads is supposed to be fixed for weak scaling, this is not going to be easy to bring meaning to.
This test uses a constant CFL - meaning the number of timesteps increases with the problem size.
This slightly ruins the weak scaling test, as we see increased overheads due to time stepping.

| Number of Cells | Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-----------------|-------------------|------------------|-----------------------------| -------------------------|
| 512 x 512       |   [3,4,5,5]       |     1.12         |          2.3437e+05         |       5.2033e-09         |
| 724 x 724       |     [3,4,5,5]     |     3.265        |                             |                          |
| 1024 x 1024     |    [5,6,7,5]      | 9.830, 9.812, 12.829, 10.070 | 
| 1448 x 1448     | [16, 16, 24, 18]  | 38.293, 40.294, 40.947, 38.563 | 

#### Using taskset

Attempt to use taskset to limit the number of threads used by std::execution:

    taskset -c 0-0 ./count

Which is confirmed as working using the thread counter demonstration.

| Number of Cells | Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-----------------|-------------------|------------------|-----------------------------| -------------------------|
| 512 x 512       |      1       | 5.714, 5.766, 8.864, 5.796        |                  |             |
| 724 x 724       |      2       | 14.838, 12.501, 12.360, 15.486           |                             |                          |
| 1024 x 1024     |     4        | 19.991, 16.524, 20.096, 16.087 | 
| 1448 x 1448     |     8        | 51.639, 57.229, 57.459, 56.095 |
| 1448 x 1448*     |     16        | 41.805, 39.387, 39.475, 38.443 |

#### Using taskset with fixed number of time steps

Using 1000 time steps for all resolutions.

| Number of Cells | Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-----------------|-------------------|------------------|-----------------------------| -------------------------|
| 512 x 512       |      1       |   7.013, 7.300, 6.967, 7.307      |                  |             |
| 724 x 724       |      2       |   10.346, 10.974, 10.841, 10.390         |                             |                          |
| 1024 x 1024     |     4        |  9.853, 12.426, 9.408, 10.366| 
| 1448 x 1448     |     8        |  22.966, 20.706, 21.093, 24.421|
| 1448 x 1448*     |     16        | 15.480, 18.701, 18.347, 17.638 |


### OpenMP

This test uses a constant CFL - meaning the number of timesteps increases with the problem size.
This slightly ruins the weak scaling test, as we see increased overheads due to time stepping.

| Number of Cells | Number of Threads* | Number of time steps | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-----------------|-------------------|------------------| ------------------|-----------------------------| -------------------------|
| 512 x 512       |     [1]           |     820          | 5.196, 5.261, 5.130, 7.792  |         |               |
| 724 x 724       |     [2]           |     1159         |  13.973, 16.893, 13.346, 17.436     |                             |                          |
| 1024 x 1024     |     [4]           |     1639         |  18.64, 17.453, 20.996, 16.512
| 1448 x 1448     |     [8]           |     2317         |  60+33.458,  60+25.034, 60+32.965, 60+35.139 |
| 1448 x 1448*    |     [16]           |     2317         |  60+4.451,  60+3.069, 60+3.210, 60+3.542 |

#### Using taskset with OpenMP and fixed number of time steps

Using taskset together with manually setting the number of threads to level the playing field.

| Number of Cells | Number of Threads* | Number of time steps | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-----------------|-------------------|------------------| ------------------|-----------------------------| -------------------------|
| 512 x 512       |     [1]           |     1000          | 6.265, 6.443, 6.585, 6.582  |         |               |
| 724 x 724       |     [2]           |     1000         | 12.932, 10.947, 9.747, 13.152     |                             |                          |
| 1024 x 1024     |     [4]           |     1000         | 11.272, 8.678, 8.740, 12.580 | 
| 1448 x 1448     |     [8]           |     1000         | 31.070, 31.925, 33.183, 32.148  |
| 1448 x 1448*    |     [16]          |     1000         | 26.554, 28.713, 25.102 |
