set -e

############################################################################
## Two datasets, fit all batches in one go, regress out Pantry phenotypes ##
############################################################################

## Process annotations and determine number of gene batches
latent-rna setup
## Bin the gene regions using coverage data from all datasets
latent-rna binning
## Prepare the coverage data for each dataset
latent-rna coverage -d dset1
latent-rna coverage -d dset2
## Fit gene models using all datasets
latent-rna fit
## Apply the models to the coverage data to get latent phenotypes
latent-rna transform -d dset1
latent-rna transform -d dset2

#################################################
## Process batches and/or datasets in parallel ##
#################################################

# latent-rna setup
# N_BATCHES=$(cat info/n_batches.txt)
# batches=$(seq 0 $((N_BATCHES-1)))
# parallel --jobs 3 latent-rna binning -b {} ::: $batches
# parallel --jobs 3 latent-rna coverage -d dset1 -b {} ::: $batches
# parallel --jobs 3 latent-rna coverage -d dset2 -b {} ::: $batches
# parallel --jobs 3 latent-rna fit -b {} ::: $batches
# parallel --jobs 3 latent-rna transform -d {} ::: dset1 dset2
