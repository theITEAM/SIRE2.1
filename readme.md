# Susceptibility Infectivity and Recoverability Estimation (SIRE2.1)

## Introduction

In the era of rapid expansion of the human population with increasing demands on food security, effective solutions that reduce the incidence and impact of infectious diseases in plants and livestock are urgently needed. Even within a species hosts differ widely in their response to infection and therefore also in their relative contribution to the spread of infection within and across populations. Three key epidemiological host traits affect infectious disease spread: susceptibility (propensity to acquire infection), infectivity (propensity to pass on infection to others) and recoverability (propensity to recover quickly). Disease control strategies aimed at reducing disease spread may, in principle, target improvement in any one of these three traits.

SIRE allows for estimation of different effects that can potentially change the suceptibility, infectivity or transition time (between compartmental states) of individuals. These effects can either be from: 1) additive genetic effects, 2) environmental effects, 3) fixed effects, and/or 4) single nucleotide polymorphisms (SNPs). 

SIRE implements a Bayesian algorithm which makes use of temporal data (consisting of any combination of recorded infection times, recovery times or disease status measurements) from multiple epidemics.

Version 2.1 allow for the specification of different types of sequential compartmental model (e.g. SI, SIR, SEIR etc...) and provides flexibility to apply different effects to different aspects of the model.

## Getting started

The software is compiled using "make".

MPI is loaded using "module load mpi/openmpi-x86_64"

The code is run using: "mpirun -n 20 ./sire examples/example1.xml 0"

Here -n 20 sets the number of cores (and MCMC chains) which are run on and the number 0 sets the random seed (which can be changed).
