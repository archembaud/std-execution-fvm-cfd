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

Table A

### OpenMP using taskset and fixed number of timesteps (weak scaling)

Using 1000 time steps for all resolutions.

| Number of Cells | Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-----------------|-------------------|------------------|-----------------------------| -------------------------|
| 512 x 512       |      1       |  10.386, 10.402, 10.401, 10.394  |                  |             |
| 724 x 724       |      2       |  10.458, 10.521, 10.476, 10.480  |                             |                          |
| 1024 x 1024     |     4        |  12.297, 12.329, 12.330, 12.396| 
| 1448 x 1448     |     8        |  20.069, 20.275, 20.107, 20.103|
| 1448 x 1448*     |     16        | 15.904, 15.980, 15.932, 15.928|

Table B

# Strong Scaling using Intel (ICX) compiler

## Run on the AMD Ryzen 9 9950X

### STD::EXECUTION (ICX)

#### 512 x 512 

These codes have also been adjusted to run 1000 time steps.

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 9.548, 9.611, 9.602, 9.596   |                  |             |
|      2       | 4.824, 4.853, 4.816, 4.835    |         |                  |
|     4        | 2.513, 2.466, 2.458, 2.477   | 
|     8        | 1.317, 1.305, 1.294, 1.295|
|     16        | 1.178, 1.183, 1.171, 1.163  |

Table C

#### 1024 x 1024

These codes have also been adjusted to run 1000 time steps.

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       |  43.243, 43.450, 43.346, 43.631  |                  |             |
|      2       |  22.288, 22.323, 22.422, 22.363  |         |                  |
|     4        |  12.172, 12.172, 12.129, 12.141  | 
|     8        |  8.342, 8.348, 8.365, 8.381|
|     16        | 5.881, 5.844, 5.854, 5.873  |

Table D

#### 2048 x 2048

These codes have also been adjusted to run 1000 time steps. This also includes time required for memory management, problem initialization, result saving (not parallel).

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 180+1.003, 180+0.927, 180+2.022, 180+1.963   |                  |             |
|      2       | 60+33.448, 60+34.448, 60+34.683,    |         |                  |
|     4        | 52.554, 53.172, 52.871, 52.569     | 
|     8        | 41.712, 41.575, 41.835, 41.614 |
|     16       | 40.856, 40.820, 40.818, 40.839 |

Table E

### OpenMP (ICX)

#### 512 x 512 

These codes have also been adjusted to run 1000 time steps.

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 10.353, 10.354, 10.327, 10.329  |                  |             |
|      2       | 5.205,  5.184, 5.210, 5.223   |         |                  |
|     4        | 2.664, 2.649, 2.653, 2.738  | 
|     8        | 1.395, 1.391, 1.365, 1.369       |
|     16       | 0.788, 0.776, 0.776, 0.766   |

Table F

#### 1024 x 1024

These codes have also been adjusted to run 1000 time steps.

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 46.183, 46.408, 46.400, 46.701 |                  |             |
|      2       | 23.481, 23.663, 23.511, 23.372   |         |                  |
|     4        | 12.317, 12.157, 12.425, 12.317   | 
|     8        | 8.439, 8.404, 8.341, 8.376
|     16        | 3.638, 3.710, 3.682, 3.687

Table G

#### 2048 x 2048

These codes have also been adjusted to run 1000 time steps. This also includes time required for memory management, problem initialization, result saving (not parallel).

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 180+11.713, 180+12.121, 180+11.992, 180+12.224
|      2       | 60+39.315, 60+39.547, 60+38.240, 60+38.097
|     4        | 54.203, 53.624, 54.091, 53.368
|     8        | 41.975, 41.882, 41.968, 41.816
|     16       | 39.969, 40.047, 39.962, 40.086

Table H

## Weak scaling on Intel Core Ultra 9 285K using ICX

### std::execution::par_unseq using taskset and fixed number of timesteps (weak scaling)

Using 1000 time steps for all resolutions.

| Number of Cells | Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-----------------|-------------------|------------------|-----------------------------| -------------------------|
| 512 x 512       |      1       |  9.733, 9.617, 9.791, 9.672    |                  |             |
| 724 x 724       |      2       |  10.271, 10.193, 10.255, 10.223  |         |                  |
| 1024 x 1024     |     4        |  11.031, 11.046, 11.105, 11.137 | 
| 1448 x 1448     |     8        |  13.268, 13.064, 13.128, 13.124 |
| 1448 x 1448*     |     16      |  11.616, 11.570, 11.585, 11.608 |

Table I

### OpenMP using taskset and fixed number of timesteps (weak scaling)

Using 1000 time steps for all resolutions.

| Number of Cells | Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-----------------|-------------------|------------------|-----------------------------| -------------------------|
| 512 x 512       |      1       |  10.542, 10.266, 10.368, 10.384  |                  |             |
| 724 x 724       |      2       |  11.099, 11.061, 10.951, 11.049 |                             |                          |
| 1024 x 1024     |     4        |  11.815, 11.818, 11.856, 11.819 | 
| 1448 x 1448     |     8        |  13.198, 13.133, 13.158, 13.110 |
| 1448 x 1448*     |     16      |  21.292, 21.218, 21.207, 21.265 |

Table J

## Run on the Intel Core Ultra 9 285K

### STD::EXECUTION (ICX)

#### 512 x 512 

These codes have also been adjusted to run 1000 time steps.

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 9.684, 9.682, 9.676, 9.663
|      2       | 4.957, 4.981, 5.002, 5.002
|     4        | 2.703, 2.709, 2.711, 2.712
|     8        | 1.446, 1.445, 1.443, 1.444
|     16       | 1.137, 1.142, 1.145, 1.131

Table K

#### 1024 x 1024

These codes have also been adjusted to run 1000 time steps.

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       |  41.214, 41.030, 41.087, 41.144
|      2       | 
|     4        | 
|     8        | 
|     16        |

Table L

#### 2048 x 2048

These codes have also been adjusted to run 1000 time steps. This also includes time required for memory management, problem initialization, result saving (not parallel).

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 
|      2       | 
|     4        | 
|     8        | 
|     16       | 

Table M

### OpenMP (ICX)

#### 512 x 512 

These codes have also been adjusted to run 1000 time steps.

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 
|      2       | 
|     4        |
|     8        |
|     16       | 

Table N

#### 1024 x 1024

These codes have also been adjusted to run 1000 time steps.

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 
|      2       | 
|     4        | 
|     8        | 
|     16        |

Table O

#### 2048 x 2048

These codes have also been adjusted to run 1000 time steps. This also includes time required for memory management, problem initialization, result saving (not parallel).

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 
|      2       | 
|     4        | 
|     8        | 
|     16       | 

Table P


