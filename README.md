
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Code for *Look-Ahead Screening Rules for the Lasso*

<!-- badges: start -->

[![R-CMD-check](https://github.com/jolars/LookAheadScreening/workflows/R-CMD-check/badge.svg)](https://github.com/jolars/HessianScreening/actions)
<!-- badges: end -->

## Results

The results from the simulations, which were run on a dedicated HPC
cluster, are stored in the [results folder](results/). The figures in
the paper, generated from these results, are stored in
[`figures/`](figures/).

## Reproducing the Results

The results from our paper were run through a singularity container.
Check the releases for pre-built singularity containers that you can
download and use.

To reproduce the results, **always** use the singularity container. To
run an experiment from the singularity container, call

``` shell
singularity run --bind results:/Project/results container.sif <script>
```

where `<script>` should be a name of a script in the [experiments
folder](experiments/), such as `experiments/simulateddata.R`.

### Re-building the Singularity Container

If you want to re-build the singularity container from scratch (or
simply want to clone the repo to your local drive), you can do so via
the following steps.

1.  Clone the repository to your local hard drive. On linux, using SSH
    authentication, run
    
    ``` shell
    git clone git@github.com:jolars/LookAheadScreening.git
    ```

2.  Navigate to the root of the repo and build the singularity container
    by calling
    
    ``` shell
    cd LookAheadScreening
    sudo singularity build container.sif Singularity
    ```

Then proceed as in [Reproducing the Results](#reproducing-the-results)
to run the experiments.

### Running Experiments without Singularity (Not Recommended\!)

Alternatively, you may also reproduce the results by cloning this
repository and starting R in the root directory of this folder (which
will activate the renv repository) and then run

``` r
renv::restore()
```

to restore the project library. Then build the R package (see below) and
run the simulations directly by running the scripts in the experiments
folder. This is **not recommended**, however, since it, unlike the
Singularity container approach, does not exactly reproduce the software
environment used when these simulations where originally run and may
result in discrepancies due to differences in for instance operating
systems, compilers, and BLAS/LAPACK implementations.

## R Package

If you want to build and experiment with the package, you can do so by
calling

``` shell
 R CMD INSTALL  .
```

provided you have `cd`ed to the root folder of this repository. First
ensure, however, that you have enabled the renv project library by
calling `renv::restore()` (see the section above).

## Data

The datasets used in these simulations are stored in the [data
folder](data/). Scripts to retrieve these datasets from their original
sources can be found in [`data-raw/`](data-raw/).
