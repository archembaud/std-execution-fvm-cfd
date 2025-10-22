# Weak Scaling using Intel (ICX) compiler

## Run on the AMD Ryzen 9 9950X

### std::execution::par_unseq using taskset and fixed number of timesteps (weak scaling)

Using 1000 time steps for all resolutions.

| Number of Cells | Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-----------------|-------------------|------------------|-----------------------------| -------------------------|
| 512 x 512       |      1       |   9.555, 9.557, 9.581, 9.629    |                  |             |
| 724 x 724       |      2       |   10.032, 9.751, 9.785, 9.762     |         |                  |
| 1024 x 1024     |     4        |   12.241, 12.286, 12.348, 12.259| 
| 1448 x 1448     |     8        |   19.749, 19.751, 19.812, 19.899|
| 1448 x 1448*     |     16      |   18.037, 18.048, 18.088, 18.036|

### OpenMP using taskset and fixed number of timesteps (weak scaling)

Using 1000 time steps for all resolutions.

| Number of Cells | Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-----------------|-------------------|------------------|-----------------------------| -------------------------|
| 512 x 512       |      1       |  10.386, 10.402, 10.401, 10.394  |                  |             |
| 724 x 724       |      2       |  10.458, 10.521, 10.476, 10.480  |                             |                          |
| 1024 x 1024     |     4        |  12.297, 12.329, 12.330, 12.396| 
| 1448 x 1448     |     8        |  20.069, 20.275, 20.107, 20.103|
| 1448 x 1448*     |     16        | 15.904, 15.980, 15.932, 15.928|