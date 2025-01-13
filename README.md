# Contents
Mathematica scripts coded by Charlie Duclut at Institut Curie France (2024). MATLAB scripts coded by Wenjie Shi at Ocean University of China (2024).

## Mathematica scripts
### Header 1:
    ciona_asymmetric_lumen_dynamics.nb
In this notebook, the dynamical equations describing the growth of a lumen made of two (symmetric or asymmetric) spherical caps are derived:

    Equations and dimensionless form of symmetric spherical cap shape
    Full equations and dimensionless form of non-symmetric spherical cap shape (no symmetry breaking, asymmetric surface tension, asymmetric actin turnover)

## MATLAB scripts
### Header 2:
    NotochordCellCortexThicknessFitting.m
Realize a Tri-Gaussian distribution fitting and extract characteristic parameters from cortex fluorescence intensity distribution in basal-lateral domain, based on the Z-stack confocal images:

    Fit the characteristic parameters of cells with no adjacent cells on either side
    calculate the cell A-P axis planar cell polarity (PCP) intensity ratio
    Fit the characteristic parameters of cells with adjacent cells on both sides
    Includes three sub-functions: CellOverlapDataSplitting.m, MultiGaussian_fit.m, MultiGaussianTransfertoEquivalentTriGaussian.m

### Header 3:
    CellOverlapDataSplitting.m
It is a sub-function of NotochordCellCortexThicknessFitting.m, aiming at split the overlap fluorescence signal from two neighbour cells:

    Determine the adjacency of the current cell with the previous and next cells
    Calculate the fluorescence intensity ratio of the previous and next cells based on the average fluorescence intensity of the basal domain
    Use different strategies to split the original data based on various conditions
    
### Header 4:
    MultiGaussian_fit.m
It is a sub-function of NotochordCellCortexThicknessFitting.m, aiming at cell grayscale data fitting with Multi-Gaussian functions:

    Identify the boundary point between the basal domain and lateral domain as the segmentation basis, and extract the middle section
    Fit the lateral domain section with a mono-Gaussian function, and use spline regression for the basal domain, determine the number of Gaussian functions based on the extrema
    Perform multi-Gaussian fittings for the basal domain using the number of Gaussians determined by the spline regression.  
    Fix the mean values of each Gaussian in the basal domain's multi-Gaussian fit results and perform a global multi-Gaussian fit for the basal-lateral domain

### Header 5:
    MultiGaussianTransfertoEquivalentTriGaussian.m
It is a sub-function of NotochordCellCortexThicknessFitting.m, aiming at converting multi-Gaussian fitting results output by MultiGaussian_fit.m to equivalent Tri-Gaussian distribution:

    Calculate the areas of the fitted multi-Gaussian distribution in the basal domain anterior/posterior-lateral domain regions
    Fix the area of the basal domain region and fit the equivalent mono-Gaussian distribution for the basal domain region to obtain its characteristic parameters
    Recalculate the areas of the anterior/posterior-lateral domain regions and fine-tune the parameters of the three-Gaussian distribution to ensure area conservation across the basal domain and the anterior/posterior-lateral domain regions

### Header 6:
    RelativeLateralEdgeSignalCalculation.m
Calculate the acto-myosin fluorescence signal at later-lumen boundary from Tri-Gaussian distribution using NotochordCellCortexThicknessFitting.m:

    Read the mix_all data file (xlsx) that stores the fitted Tri-Gaussian distribution characteristic parameters and the protein file (e.g., lifeact) (xlsx) that stores file names corresponding to the characteristic parameters
    Exclude cells with excessive fitting errors or cells where characteristic parameters significantly deviate from the average
    Calculate the average thickness and average length of the lateral domain cortex before lumen formation
    Calculate the cortex thickness and length of the anterior- and posterior-lateral domains of each cell at different stages after lumen formation
    Calculate the contractile force exerted on each lumen by the contractile rings of notochord cells on both sides

### Header 7:
    RelativeLumenGeometryunderZO1Overexpression.m
Calculate the relative lumen TD, LD and volume after ZO1 dominant negative form overexpression

    Assign the values of Transverse Diameter (radial) and Longitudinal Radius (axial) from experiments
    Calculate the mean and standard deviation of the lumen volume for each individual in the ZO1 WT overexpression group and the control group
    Calculate the mean and standard deviation of volume and length for the control group and the experimental group
    Perform statistical tests

### Header 8:
    LinearSstabilityAnalysisofSphericalCapCase.m
Numerically calculate of the steady-state lumen and its stability as a function of model parameters

    Estimate an initial value of the model parameters
    Calculate the steady-state solution
    Compute the eigenvalues of the linear stability matrix to determine the stability of each given parameter value on an one-dimensional parameter plane

### Header 9:
    PhaseDiagramofSphericalCapCase.m
Numerically calculate of the state diagram of the steady-state lumen as a function of the given two-dimensional model parameters plane

    Estimate an initial value of the model parameters
    Give a precise range of two chosen parameters
    Calculate the steady-state solution
    Compute the eigenvalues of the linear stability matrix to determine the stability of each given parameter value on a two-dimensional parameter plane

### Header 10:
    FullSolutionofSphericalCapLumenDynamics.m
Numerically solve the differential-algebraic equations and predict lumen dynamics

    Give a series of values for the dimensionless parameters. Those parameters may be fixed in time or may be time-dependent
    Give an initial condition for the lumen geometry and the ion concentration difference
    Solve the equations using the custom-made first order DAEs solver





