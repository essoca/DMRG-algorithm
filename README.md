# DMRG-algorithm
Principal component analysis for 1D quantum many-body systems

This repo contains sample code in MATLAB of the Density Matrix Renormalization Group or [DMRG algorithm](https://en.wikipedia.org/wiki/Density_matrix_renormalization_group), which uses an idea similar to the statistical [principal component analysis](https://en.wikipedia.org/wiki/Principal_component_analysis) to investigate the low-energy physics of one-dimensional quantum many-body systems. The code is organized as follows:

- OBCdmrg: Implements the ground-state DMRG (at zero temperature) with open boundary conditions.
- t-dmrg: Implements time-dependent DMRG at zero temperature.
- LowTdmrg: Extends t-dmrg to the evolution in imaginary time to investigate finite-temperature physics.
