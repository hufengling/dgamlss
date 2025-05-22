# dGAMLSS
--------
**Maintainer**: Fengling Hu, fengling.hu@pennmedicine.upenn.edu

## Table of content
- [1. Installation](#id-section1)
- [2. Background](#id-section2)
- [3. Software](#id-section3)
- [4. Citation](#id-section4)

<div id='id-section1'/>

## 1. Installation
The R package can be installed via devtools by running the following code

```
# install.packages("devtools")
devtools::install_github("hufengling/dgamlss")
```

Then, you can load this package via

```
library(dgamlss)
```

## 2. Background
There is growing interest in estimating population reference ranges across age and sex to better identify atypical clinically-relevant measurements throughout the lifespan. For this task, the World Health Organization recommends using Generalized Additive Models for Location, Scale, and Shape (GAMLSS) which can model non-linear growth trajectories under complex distributions that address the heterogeneity in human populations.
Fitting GAMLSS models requires large, generalizable sample sizes, especially for accurate estimation of extreme quantiles, but obtaining such multi-site data can be challenging due to privacy concerns and practical considerations. In settings where patient data cannot be shared, privacy-preserving distributed algorithms for federated learning can be used, but no such algorithm exists for GAMLSS.
We propose distributed GAMLSS (dGAMLSS), a distributed algorithm which can fit GAMLSS models across multiple sites without sharing patient-level data. We demonstrate the effectiveness of dGAMLSS in constructing population reference charts across clinical, genomics, and neuroimaging settings.

<img width="468" alt="image" src="https://github.com/user-attachments/assets/67798915-b7b3-4ecd-a2cf-675c98b853bb" />

This Figure provides an overview of the dGAMLSS algorithm which enables fitting of GAMLSS models via privacy-preserving federated learning. Fitting occurs via an outer cycle over the parameters and a nested inner cycle of Newton-Raphson updates for parameter-specific coefficients. Broadly, the algorithm is as follows. Step 1: Initialize all coefficients for all parameters and set the outer loop to the μ parameter. Step 2: For the given parameter, send current coefficients and receive site-level matrices. Step 3: Calculate coefficient updates for the given parameter. If the updated coefficients are identical to the previous coefficients, advance to the next parameter. Step 4: Repeat Steps 2 and 3 until all coefficients converge for all parameters.

<div id='id-section3'/>

## 3. Software

The `dgamlss` package provides a number of functions for end-users, based on the `gamlss` R package introduced by Rigby and Stasinopoulos (2005), <doi:10.1111/j.1467-9876.2005.00510.x>

`dgamlss_coordinating()` and `dgamlss_coordinating_penalized()` provide examples for how to use other functions in the package when fixed basis splines or penalized basis splines are desired, respectively. These functions represent how a central, or "coordinating" site, would aggregate information across a number of remote, or "local" sites.

Other functions in the package are used to 1) define model formulations, including spline basis and penalty matrix, 2) summarize local site data such that privacy is preserved but information can be communicated, 3) aggregate local site data to fit the dGAMLSS model, and 4) conduct inference and plot dGAMLSS reference charts. A number of behind-the-scenes helper functions are also included.

## 4. Citation
If you are using DeepComBat for harmonization, please cite the following article:

dGAMLSS: An exact, distributed algorithm to fit Generalized Additive Models for Location, Scale, and Shape for privacy-preserving population reference charts
Fengling Hu, Jiayi Tong, Margaret Gardner, Lifespan Brain Chart Consortium, Andrew A. Chen, Richard A.I. Bethlehem, Jakob Seidlitz, Hongzhe Li, Aaron Alexander-Bloch, Yong Chen, Russell T. Shinohara
bioRxiv 2024.12.20.629834; doi: https://doi.org/10.1101/2024.12.20.629834
