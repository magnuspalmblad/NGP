# NGP
This repository is a collection of simulation and analysis algorithms for single-molecule, next-generation proteomics (NGP). Specifically, the scripts can be used to study the influence of choice of labeled amino acids on the number of unique and the number of proteotypic reads for a given proteome.

The "simulations.R" script is has no inputs, but these parameters are defined in the beginning of the script:

| Parameter        | Type    | Meaning                                                            |
| ---------------- |---------| -------------------------------------------------------------------|
| N_SIMULATIONS    | integer | The number of simulations of 1-20 randomly chosen amino acids      |
| MAX_LENGTH       | integer | The maximum peptide read length considered                         |
| MISSED_CLEAVAGES | integer | The maximum number of missed cleavage sites                        |
| ENZYME           | string  | The proteolytic enzyme simulated (e.g. "trypsin")                  |
| N_PROTEINS       | integer | The number of proteins from the database considered (-1 for all)   |
| IDEAL            | boolean | whether the sequencing is assumed to be ideal or containing errors |

The outputs are matrices containing the selection of amino acids, the number of proteins lacking a proteotypic partial read and the number of unique possible reads. 

The seed for the random number generator is set by set.seed(). Running simulations with the same seed but different parameters will repeat the simulations for the exact same selection of amino acids.

The "plots.R" makes pretty visualizations from the simulations using ggplot2. It also uses multidimensional linear regression to determine the relative discrimination power of the different amino acids. This R script was used to produce the figures in the JPR technical note describing the generic simulations of single-molecule/next-generation proteomics.

The scripts are written in R for maximum readability. For deeper simulations, refactoring them in C/C++ will likely be beneficial.
