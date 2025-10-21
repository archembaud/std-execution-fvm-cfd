# Strong Scaling Notes

Strong scaling tests also require understanding of how to set the number of threads used by std::execution parallel execution models.

Hence, we need to run these again.

## Run on the AMD Ryzen 9 9950X

### STD::EXECUTION (G++)

#### 512 x 512 

These codes have also been adjusted to run 1000 time steps.

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       |  6.076, 6.045, 6.048, 6.032    |                  |             |
|      2       |  3.058, 3.055, 3.068, 3.054      |         |                  |
|     4        |  1.585, 1.563, 1.577, 1.572| 
|     8        |  0.855, 0.843, 0.842, 0.869|
|     16        | 1.051, 1.048, 1.044, 1.043 |

#### 1024 x 1024

These codes have also been adjusted to run 1000 time steps. This also includes time required for memory management, problem initialization, result saving (not parallel).

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 29.072, 29.147, 29.157, 29.330     |                  |             |
|      2       | 15.063, 15.300, 15.366, 15.329     |         |                  |
|     4        | 9.121, 9.060, 9.257, 9.175     | 
|     8        | 7.952, 7.921, 7.978, 7.911  |
|     16       | 5.758, 5.781, 5.767, 5.761 |

Time required for initialization, saving, freeing etc: 1 thread: 0.091, 2 threads: 0.095, 16 threads: 0.095.

#### 2048 x 2048

These codes have also been adjusted to run 1000 time steps. This also includes time required for memory management, problem initialization, result saving (not parallel).

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 120+5.284, 120+4.937, 120+5.658,120+5.190   |                  |             |
|      2       | 60+8.181, 60+8.494, 60+8.571,60+7.752     |         |                  |
|     4        | 43.241, 43.408, 44.120, 43.690     | 
|     8        | 41.233, 41.117, 41.146, 41.283  |
|     16       | 40.659, 40.534, 40.583, 40.679 |

### OpenMP (G++)

#### 512 x 512 

These codes have also been adjusted to run 1000 time steps.

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 5.422, 5.400, 5.382, 5.382   |                  |             |
|      2       | 2.734, 2.738, 2.729, 2.716      |         |                  |
|     4        | 1.402, 1.387, 1.390, 1.393 | 
|     8        | 0.747, 0.728, 0.713, 0.728 |
|     16       | 0.432, 0.428, 0.426, 0.426 |

#### 1024 x 1024

These codes have also been adjusted to run 1000 time steps. This also includes time required for memory management, problem initialization, result saving (not parallel).

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       |  26.414, 26.802, 26.724, 26.767   |                  |             |
|      2       |  14.130, 14.104, 14.089, 14.135   |         |                  |
|     4        |  8.739, 8.936, 8.884, 8.814    | 
|     8        |  7.787, 7.873, 7.828, 7.891 |
|     16       |  2.813, 2.845, 2.787, 2.838|

Initialisation times and saving time: 1 core: 0.024, 2 cores: 0.024 ... 16 cores: 0.023
This is hardly worth mentioning.

#### 2048 x 2048

These codes have also been adjusted to run 1000 time steps. This also includes time required for memory management, problem initialization, result saving (not parallel).

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 60+56.529, 60+56.217, 60+55.571, 60+56.766  |                  |             |
|      2       | 60+2.939, 60+3.380, 60+4.469, 60+2.424    |         |                  |
|     4        | 42.890, 43.452, 42.603, 43.702 | 
|     8        | 41.297, 41.038, 41.339, 41.098 |
|     16       | 39.487, 39.455, 39.539, 39.647 |