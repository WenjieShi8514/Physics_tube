# Contents
Mathematica scripts coded by Charlie Duclut at Institut Curie France (2024). MATLAB scripts coded by Wenjie Shi at Ocean University of China (2024).

## Mathematica scripts


## MATLAB scripts
### Header 1:
    NotochordCellCortexThicknessFitting.m
Realize a Tri-Gaussian distribution fitting and extract characteristic parameters from cortex fluorescence intensity distribution in basal-lateral domain, based on the Z-stack confocal images:

    Fit the characteristic parameters of cells with no adjacent cells on either side
    calculate the cell A-P axis planar cell polarity (PCP) intensity ratio
    Fit the characteristic parameters of cells with adjacent cells on both sides
    Includes three sub-functions: CellOverlapDataSplitting.m, MultiGaussian_fit.m, MultiGaussianTransfertoEquivalentTriGaussian.m

### Header 2:
    CellOverlapDataSplitting.m
It is a sub-function of NotochordCellCortexThicknessFitting.m, aiming at split the overlap fluorescence signal from two neighbour cells:

    Determine the adjacency of the current cell with the previous and next cells
    Calculate the fluorescence intensity ratio of the previous and next cells based on the average fluorescence intensity of the basal domain
    Use different strategies to split the original data based on various conditions
    
### Header 3:
    MultiGaussian_fit.m
It is a sub-function of NotochordCellCortexThicknessFitting.m, aiming at cell grayscale data fitting with Multi-Gaussian functions:
    Identify the boundary point between the basal domain and lateral domain as the segmentation basis, and extract the middle section
    Fit the lateral domain section with a mono-Gaussian function, and use spline regression for the basal domain, determine the number of Gaussian functions based on the extrema
    Perform multi-Gaussian fittings for the basal domain using the number of Gaussians determined by the spline regression.  
    Fix the mean values of each Gaussian in the basal domain's multi-Gaussian fit results and perform a global multi-Gaussian fit for the basal-lateral domain

### Header 3:
    MultiGaussianTransfertoEquivalentTriGaussian.m, aiming at Convert multi-Gaussian fitting results to equivalent Tri-Gaussian distribution
It is a sub-function of NotochordCellCortexThicknessFitting.m, aiming at cell grayscale data fitting with Multi-Gaussian functions:








