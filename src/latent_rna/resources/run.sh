set -e

n_batches=3  # This will depend on batch size (in config) and number of protein-coding genes
batches=$(seq 0 $((n_batches-1)))

############################################################################
## Two datasets, fit all batches in one go, regress out Pantry phenotypes ##
############################################################################

## Bin the gene regions using coverage data from all datasets
latent-rna binning
## Prepare the coverage data for each dataset
latent-rna prepare -d dset1
latent-rna prepare -d dset2
## Fit gene models using all datasets
latent-rna fit
## Apply the models to the coverage data to get latent phenotypes
latent-rna transform -d dset1
latent-rna transform -d dset2

#################################################
## Process batches and/or datasets in parallel ##
#################################################

# parallel --jobs 3 latent-rna binning -b {} ::: $batches
# parallel --jobs 3 latent-rna prepare -d dset1 -b {} ::: $batches
# parallel --jobs 3 latent-rna prepare -d dset2 -b {} ::: $batches
# parallel --jobs 3 latent-rna fit -b {} ::: $batches
# parallel --jobs 3 latent-rna transform -d {} ::: dset1 dset2
