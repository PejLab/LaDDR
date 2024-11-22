set -e

batch_size=20
n_batches=3

# seq 0 $((n_batches-1)) > data_bash_script/batches.txt

## Single dataset, fit all batches in one go, no regressing out Pantry phenotypes
latent-rna binning
latent-rna prepare
latent-rna fit
latent-rna transform
