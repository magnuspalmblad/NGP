# NGP
This repository is a collection of next-generation proteomics (NGP) simulation and analysis algorithms.

The "simulations.R" script is has no inputs, but these parameters are defined in the beginning of the script:

| Parameter        | Type    | Meaning                                                            |
| ---------------- |---------| -------------------------------------------------------------------|
| MAX_LENGTH       | integer | The maximum peptide read length considered                         |
| MISSED_CLEAVAGES | integer | The maximum number of missed cleavage sites                        |
| ENZYME           | string  | The proteolytic enzyme simulated (e.g. "trypsin")                  |
| N_PROTEINS       | integer | The number of proteins from the database considered (-1 for all)   |
| IDEAL            | boolean | whether the sequencing is assumed to be ideal or containing errors |

The seed for the random number generator is set by set.seed(). Running simulations with the same seed but different parameters will repeat the simulations for the exact same selection of amino acids.

The "plots.R" makes pretty visualizations from the simulations using ggplot2. It also uses multidimensional linear regression to determine the relative discrimination power of the different amino acids. This R script was used to produce the figures in the technical note describing the generic simulations of single-molecule/next-gen proteomics.
