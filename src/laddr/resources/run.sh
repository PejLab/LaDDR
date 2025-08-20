set -e

############################################################################
## Two datasets, fit all batches in one go, regress out Pantry phenotypes ##
############################################################################

## Process annotations and determine number of gene batches (do not run if using existing info/ directory)
laddr setup
## Bin the gene regions using coverage data from all datasets (do not run if using existing gene_bins/ directory)
laddr binning
## Prepare the coverage data for each dataset
laddr coverage -d dset1
laddr coverage -d dset2
## Fit gene models using all datasets (do not run if using existing models/ directory)
laddr fit
## Apply the models to the coverage data to get latent phenotypes
laddr transform -d dset1
laddr transform -d dset2

#################################################
## Process batches and/or datasets in parallel ##
#################################################

# laddr setup
# N_BATCHES=$(cat info/n_batches.txt)
# batches=$(seq 0 $((N_BATCHES-1)))
# parallel --jobs 3 laddr binning -b {} ::: $batches
# parallel --jobs 3 laddr coverage -d dset1 -b {} ::: $batches
# parallel --jobs 3 laddr coverage -d dset2 -b {} ::: $batches
# parallel --jobs 3 laddr fit -b {} ::: $batches
# parallel --jobs 3 laddr transform -d {} ::: dset1 dset2

########################
## Extract model data ##
########################

# laddr inspect -g ENSG00000008128 -d dset1
# laddr inspect -g ENSG00000008128 -d dset2
