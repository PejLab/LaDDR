set -e

script=../latent_RNA.py
gtf=data_input/Homo_sapiens.GRCh38.106.chr1_0-2Mb.gtf
bamdir=data_input/bam

mkdir -p data_bash_script/ex1
mkdir -p data_bash_script/ex2
mkdir -p data_bash_script/ex3
mkdir -p data_bash_script/ex4
mkdir -p data_bash_script/ex5

## head -n 8 data_input/samples.txt > data_input/samples.dset1.txt
## tail -n +9 data_input/samples.txt > data_input/samples.dset2.txt
## echo -e 'dset1\ndset2' > data_input/datasets.txt

###################
## Get gene bins ##
###################

python ../scripts/collapse_annotation.py $gtf data_bash_script/collapsed.gtf --collapse_only
gzip data_bash_script/collapsed.gtf

python ../scripts/get_gene_bins.py \
    --gtf data_bash_script/collapsed.gtf.gz \
    --chromosomes data_input/chr_lengths.genome \
    --output data_bash_script/gene_bins.bed.gz

######################
## Get bin coverage ##
######################

mkdir -p data_bash_script/covg_sample
cat data_input/samples.txt | while read sample; do
    echo $sample
    bedtools coverage -split -sorted -counts \
        -a data_bash_script/gene_bins.bed.gz \
        -b $bamdir/$sample.bam \
        -g data_input/chr_lengths.genome \
        | cut -f7 \
        > data_bash_script/covg_sample/$sample.txt
done

awk '{print "data_bash_script/covg_sample/"$1".txt"}' data_input/samples.txt > data_bash_script/covg_sample_files.txt
awk '{print "data_bash_script/covg_sample/"$1".txt"}' data_input/samples.dset1.txt > data_bash_script/covg_sample_files.dset1.txt
awk '{print "data_bash_script/covg_sample/"$1".txt"}' data_input/samples.dset2.txt > data_bash_script/covg_sample_files.dset2.txt

#########################
# Fit and apply models ##
#########################

## Single dataset, fit all batches in one go, no regressing out Pantry phenotypes
python $script prepare -i data_bash_script/covg_sample_files.txt -r data_bash_script/gene_bins.bed.gz -d data_bash_script/ex1/covg_batch/ --batch-size=20
python $script fit -d data_bash_script/ex1/covg_batch/ -m data_bash_script/ex1/models/
python $script transform -d data_bash_script/ex1/covg_batch/ -m data_bash_script/ex1/models/ -o data_bash_script/ex1/latent_phenos.tsv.gz

## Single dataset, fit all batches in one go, regress out Pantry phenotypes
python $script prepare -i data_bash_script/covg_sample_files.txt -r data_bash_script/gene_bins.bed.gz --pheno-paths-file data_input/pheno_files.txt -d data_bash_script/ex2/covg_batch/ --batch-size=20
python $script fit -d data_bash_script/ex2/covg_batch/ -m data_bash_script/ex2/models/
python $script transform -d data_bash_script/ex2/covg_batch/ -m data_bash_script/ex2/models/ -o data_bash_script/ex2/latent_phenos.tsv.gz

## Single dataset, fit batches separately
python $script prepare -i data_bash_script/covg_sample_files.txt -r data_bash_script/gene_bins.bed.gz --pheno-paths-file data_input/pheno_files.txt -d data_bash_script/ex3/covg_batch/ --batch-size=20
parallel -j3 python $script fit -d data_bash_script/ex3/covg_batch/ -b {} -m data_bash_script/ex3/models/ :::: data_bash_script/ex3/covg_batch/batches.txt
python $script transform -d data_bash_script/ex3/covg_batch/ -m data_bash_script/ex3/models/ -o data_bash_script/ex3/latent_phenos.tsv.gz

## Multiple datasets (prepare separately, train together, transform separately)
python $script prepare -i data_bash_script/covg_sample_files.dset1.txt -r data_bash_script/gene_bins.bed.gz --pheno-paths-file data_input/pheno_files.txt -d data_bash_script/ex4/covg_batch_dset1/ --batch-size=20
python $script prepare -i data_bash_script/covg_sample_files.dset2.txt -r data_bash_script/gene_bins.bed.gz --pheno-paths-file data_input/pheno_files.txt -d data_bash_script/ex4/covg_batch_dset2/ --batch-size=20
python $script fit -d data_bash_script/ex4/covg_batch_dset1/ data_bash_script/ex4/covg_batch_dset2/ -b 0 -m data_bash_script/ex4/models/
python $script fit -d data_bash_script/ex4/covg_batch_dset1/ data_bash_script/ex4/covg_batch_dset2/ -b 1 -m data_bash_script/ex4/models/
python $script fit -d data_bash_script/ex4/covg_batch_dset1/ data_bash_script/ex4/covg_batch_dset2/ -b 2 -m data_bash_script/ex4/models/
python $script transform -d data_bash_script/ex4/covg_batch_dset1/ -m data_bash_script/ex4/models/ -o data_bash_script/ex4/latent_phenos.dset1.tsv.gz
python $script transform -d data_bash_script/ex4/covg_batch_dset2/ -m data_bash_script/ex4/models/ -o data_bash_script/ex4/latent_phenos.dset2.tsv.gz

## Many datasets (prepare separately, train together, transform separately)
# Same as above but looping through datasets and passing a file listing the directories:
cat data_input/datasets.txt | while read dset; do
    python $script prepare -i data_bash_script/covg_sample_files.${dset}.txt -r data_bash_script/gene_bins.bed.gz --pheno-paths-file data_input/pheno_files.txt -d data_bash_script/ex5/covg_batch_${dset}/ --batch-size=20
done
awk '{print "data_bash_script/ex5/covg_batch_"$1}' data_input/datasets.txt > data_bash_script/ex5/covg_batch_dirs.txt
cat data_bash_script/ex5/covg_batch_dset1/batches.txt | while read batch; do
    python $script fit --dir-file data_bash_script/ex5/covg_batch_dirs.txt -b $batch -m data_bash_script/ex5/models/
done
cat data_input/datasets.txt | while read dset; do
    python $script transform -d data_bash_script/ex5/covg_batch_${dset}/ -m data_bash_script/ex5/models/ -o data_bash_script/ex5/latent_phenos.${dset}.tsv.gz
done
