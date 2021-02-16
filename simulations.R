library(plyr)
library(UniProt.ws)
library(cleaver)

## simulation parameters
N_SIMULATIONS <- 20
MAX_LENGTH <- 50
MISSED_CLEAVAGES <- 1
ENZYME <- "trypsin"
N_PROTEINS <- 20322 # use -1 for all proteins in database
IDEAL <- TRUE
set.seed(1)

## select species (e.g. 9606 = Homo sapiens)
up <- UniProt.ws(taxId = 9606)

## download all sequences (e.g. NeXtProt)
nextprot <- keys(up, keytype = "NEXTPROT")
if (N_PROTEINS == -1)
  N_PROTEINS <- length(nextprot)
proteins <-
  select(up,
         nextprot[1:N_PROTEINS],
         keytype = "NEXTPROT",
         columns = c("SEQUENCE"))
names(proteins) <- nextprot[1:N_PROTEINS]

## in-silco digestion of the protein sequences
peptides <-
  cleave(proteins$SEQUENCE,
         enzym = ENZYME,
         missedCleavages = 0:MISSED_CLEAVAGES)

## calculate average % of proteins with at least one unique partial read
results <- matrix(0, nrow = N_SIMULATIONS, ncol = 20)
choice <- matrix("", nrow = N_SIMULATIONS, ncol = 20)
unique_reads <- matrix(0, nrow = N_SIMULATIONS, ncol = 20)
for (run in 1:N_SIMULATIONS) {
  # run N_SIMULATIONS simulations
  for (aa in 1:20) {
    # number of readable amino acids
    
    ## generate partial reads
    temp <- peptides
    digest <- list()
    AMINO_ACIDS <- "ARNDCQEGHILKMFPSTWYVUOBJZX"
    READ_AA <- "ARNDCQEGHILKMFPSTWYV"
    rank <- sort(runif(20), index.return = TRUE)$ix[1:aa] # select amino acids
    for (i in 1:20)
      if (!is.element(i, rank)) substr(READ_AA, i, i) <- "i"
    visible <- gsub("i", "", READ_AA) # "visible", i.e. read aa's
    for (i in 1:20)
      if (is.element(i, rank)) substr(AMINO_ACIDS, i, i) <- "r"
    invisible <- gsub("r", "", AMINO_ACIDS) # "invisible", i.e. not read aa's
    invisible <- paste("[", invisible, "]", sep = "")
    i <- 1 # protein index
    
    for (x in temp) {
      x <- gsub("[[:space:]]", "", x) # remove whitespace (including \n)
      x <- gsub(invisible, "-", x) # replace invisible aa's with "-"
      digest[[i]] <- x # add vector of possible peptide reads for protein i
      if (IDEAL == FALSE) {
        for (x1 in x) {
          for (j in 1:nchar(x1)) # simulate non-ideal situation, allowing - to be --
          {
            if (substr(x1, j, j) == "-") {
              x_alt <-
                paste(substr(x1, 1, j - 1),
                      "--",
                      substr(x1, j + 1, nchar(x1)),
                      sep = "")
              if (!(is.element(x_alt, digest[[i]])))
                digest[[i]] <- append(digest[[i]], x_alt)
            }
          }
        }
      }
      i <- i + 1
    }
    
    ## find unique reads (reads only appearing once)
    counts <- count(unlist(digest)) # get read frequencies
    singles <- counts[counts$freq == 1,][1] # get reads appearing only once
    singles <- singles$x[nchar(singles[,]) <= MAX_LENGTH] # remove long peptides
    
    ## count number of globally unique (proteotypic) reads per protein
    N <- c()
    for (i in 1:N_PROTEINS)
      N[i] <- length(intersect(digest[[i]], singles))
    
    ## count number of proteins with zero unique reads and store in results
    results[run, aa] <- results[run, aa] + sum(N == 0) # proteins with no proteotypic read
    choice[run, aa] <- visible # which amino acids were read
    unique_reads[run, aa] <-  dim(counts)[1] # how many unique peptide reads
    print(results)
    print(choice)
    print(unique_reads)
    print(run)
  }
}
