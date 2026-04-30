# soilmapping-in-python
A simple Python repo with some tools useful in the world of digital soil mapping, taken from other work done in the R language.

## Overview
This repo exists because most soil sampling scientists use R, not Python. I use Python, so I wanted to make scientific work available to a wider audience.
All of the work in this repo is either based off existing works. Several common data science libraries are used as well.
I like to use the jupyter lab environment, but the code can be adapted to other Python runtimes as well.
The versions listed below are what was used, and I am in the process of updating where newer versions are available.


## PyPI Packages Required
- numpy==1.26.4
- matplotlib==3.9.0
- tqdm==4.66.5
- joblib==1.4.2
- scipy==1.13.1
- scikit-learn==1.6.1
- clhs==1.0.2
- pyreadr==0.5.4
- jupyterlab==4.2.1

## Works Included
- https://github.com/newdale/opendsm by Daniel Saurette: this is the core inspiration for this repo and opt_samp.py.
- https://doi.org/10.1016/j.geoderma.2023.116553 by Saurette et al., 2023: this is the main paper which provides the test case for this code, and is also the scientific description for the above repo.
- https://doi.org/10.7717/peerj.6451 by Malone et al., 2019: this paper provides algorithmic details used in the quantile, cov_matrix, and h_matrix code chunks.
- The packages from PyPI are also works in their own right. Check out their PyPI pages!
