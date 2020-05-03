# Phase-stability-analyzer

Set of routines to perform electrochemical stability analysis of solid-state electrolytes for Lithium-ion batteries. A complete enumeration of possible decomposition reactions is found with their associated potentials, as described in https://arxiv.org/pdf/1901.02251.pdf .

In this implementation, materials and their total energies are queried from the Materials Project ( https://materialsproject.org/ ), the data post-processed and error estimates for the potential of the instability reactions are provided. An example notebook "Example.ipynb" shows example of a possible usage.
