# Strong Scaling Notes

Strong scaling tests also require understanding of how to set the number of threads used by std::execution parallel execution models.

Hence, we need to run these again.

## Run on the AMD Ryzen 9 9950X

### STD::EXECUTION

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

#### 2048 x 2048

These codes have also been adjusted to run 1000 time steps. This also includes time required for memory management, problem initialization, result saving (not parallel).

| Number of Threads* | Average Time (s) | Throughput (cells / second) |  Time / cell / timestep  |
|-------------------|------------------|-----------------------------| -------------------------|
|      1       | 120+5.284, 120+4.937, 120+5.658,120+5.190   |                  |             |
|      2       | 60+8.181, 60+8.494, 60+8.571,60+7.752     |         |                  |
|     4        | 43.241, 43.408, 44.120, 43.690     | 
|     8        | 41.233, 41.117, 41.146, 41.283  |
|     16       | 40.659, 40.534, 40.583, 40.679 |
