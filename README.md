# Semi-parametric estimation of an isotropic Spatio-Temporal Hawkes process for car accidents on a road network
This repository includes data and codes to efficiently estimate a semi-parametric spatio-temporal Hawkes process on a road network. The model is applied to a dataset of road accidents that occurred within the Great Ring Road (GRA) surrounding the urban area of Rome. In particular, the model's specification can account for both spatial and temporal patterns characterizing such phenomena, and the road network topology is needed to make valid inference on the data generative process.

The repository contains:
- a folder, named ExampleRome_M1M2_42019, containing the codes to estimate the model on an example dataset concerning the road accidents that occurred in Municipio I and Municipio II of the City of Rome in April 2019.
- a subfolder ExampleRome_M1M2_42019/Data/ including the shapefiles and the raw data in .csv format
- a C++ file, named _auxiliary.cpp_, including auxiliary functions that are needed to efficiently compute the model ingredients
- Three additional R scripts that should be run sequentially: 
    - 0.data_preparation.R
    - 1.kernel_computation.R
    - 2.model_estimation.R

In what follows, we provide a general overview of the considered semi-parametric spatio-temporal Hawkes Process.

## The periodic spatio-temporal Hawkes process

We seek to model the occurrence of traffic collisions over a spatio-temporal domain $\mathcal{Q}=\mathcal{D}\times\mathcal{T}$, where $\mathcal{D}\subseteq\mathbb{R}^2$ denotes the spatial dimension and $\mathcal{T}=[0, T]$ the temporal dimension.
We assume that the number of car-crashes $N(B\times[t_1,t_2])$, where $B\subset\mathcal{D}$ and $[t_1,t_2]\subset\mathcal{T}$, is the result of a simple and (locally) finite spatio-temporal point process. It can be defined through the suitable specification of the conditional intensity function $\lambda_c(\boldsymbol{s}, t)$.

In particular, we express these two components as in the semi-parametric spatio-temporal periodic Hawkes Process:
$$\lambda_c(\boldsymbol{s}, t)=\mu_0\cdot\mu_{s}(\boldsymbol{s})\cdot\mu_t\left( t\right) + A\cdot\int_0^t\int_\mathcal{D} g_s(\boldsymbol{s}-\boldsymbol{u})\cdot g_t(t-\tau) N\left( d \boldsymbol{u}\times d \tau\right),$$
where $\mu_s(\cdot), \mu_t(\cdot)$ are the spatial and temporal background intensities such that their average value over $\mathcal{D}$ and $\mathcal{T}$ is $1$, $g_s(\cdot), g_t(\cdot)$ are the spatial and temporal excitation functions such that their integral over $\mathcal{D}$ and $\mathcal{T}$ is $1$, and $\mu_0,\, A>0$ are two real-valued parameters that regulate the overall level of the background and the excitation.
The spatial excitation function depends on the Euclidean distance between the primary event and nearby locations:
$$g_s\left(\boldsymbol{s}'-\boldsymbol{s}\right)=g_s\left( ||\boldsymbol{s}'-\boldsymbol{s}||\right)=g_s\left(\sqrt{(x'-x)^2+(y'-y)^2}\right),$$
so that $g_s(\cdot):\mathbb{R}^+\rightarrow [0,+\infty)$.

Let $\boldsymbol{x}_{i}$ be a $((k+1)\times 1)$ vector of covariates available on each event, we can express:
$$\log\left( A_i\right) = \boldsymbol{x}_i^\top\cdot\boldsymbol{\beta},\quad  i=1,\dots, n,$$
as in a Generalized Linear Model, with $\boldsymbol{\beta}$ as a vector of intercept and coefficients.

## Estimation
The perform estimation we need to evaluate the relative impact of background and excitation in each event. We here propose a semi-parametric estimation procedure where the various functions' shapes are estimated non-parametrically through weighted kernel smoothing and the coefficients $\mu_0$ and $A$ are estimated through maximum likelihood. The two estimations can be unified in an alternate procedure inside an EM-type estimation algorithm.
