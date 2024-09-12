set -e

script=../latent_RNA.py
gtf=data_input/Homo_sapiens.GRCh38.106.chr1_0-2Mb.gtf
bamdir=data_input/bam

batch_size=20
n_batches=3

mkdir -p data_bash_script/ex1
mkdir -p data_bash_script/ex2
mkdir -p data_bash_script/ex3
mkdir -p data_bash_script/ex4
mkdir -p data_bash_script/ex5

## head -n 8 data_input/samples.txt > data_input/samples.dset1.txt
## tail -n +9 data_input/samples.txt > data_input/samples.dset2.txt
## echo -e 'dset1\ndset2' > data_input/datasets.txt

seq 0 $((n_batches-1)) > data_bash_script/batches.txt

###################
## Get gene bins ##
###################

python ../scripts/collapse_annotation.py $gtf data_bash_script/collapsed.gtf --collapse_only
gzip data_bash_script/collapsed.gtf

python ../scripts/get_gene_bins.py \
    --gtf data_bash_script/collapsed.gtf.gz \
    --chromosomes data_input/chr_lengths.genome \
    --binning-method bin-width \
    --bin-width-coding 16 \
    --bin-width-noncoding 128 \
    --batch-size $batch_size \
    --outdir data_bash_script/gene_bins

######################
## Get bin coverage ##
######################

mkdir -p data_bash_script/covg_bigwig

## Generate bigWig files using deeptools
mkdir -p data_bash_script/covg_bigwig
cat data_input/samples.txt | while read sample; do
    echo $sample
    # samtools index $bamdir/$sample.bam
    bamCoverage -b $bamdir/$sample.bam -o data_bash_script/covg_bigwig/$sample.bw -of bigwig --binSize 1 -p 4
done

## bigWig to mean covg per bin using deeptools
mkdir -p data_bash_script/covg_binned
bw_files=$(cat data_input/samples.txt | sed 's/^/data_bash_script\/covg_bigwig\//' | sed 's/$/.bw/' | tr '\n' ' ')
labels=$(cat data_input/samples.txt | tr '\n' ' ')
cat data_bash_script/batches.txt | while read batch; do
    multiBigwigSummary BED-file \
        --bwfiles $bw_files \
        --labels $labels \
        --BED data_bash_script/gene_bins/batch_$batch.bed.gz \
        -o data_bash_script/covg_binned/batch_$batch.npz
done

## Same as above but for two separate datasets
cat data_input/datasets.txt | while read dset; do
    mkdir -p data_bash_script/covg_binned_${dset}
    bw_files=$(cat data_input/samples.${dset}.txt | sed 's/^/data_bash_script\/covg_bigwig\//' | sed 's/$/.bw/' | tr '\n' ' ')
    labels=$(cat data_input/samples.${dset}.txt | tr '\n' ' ')
    cat data_bash_script/batches.txt | while read batch; do
        multiBigwigSummary BED-file \
            --bwfiles $bw_files \
            --labels $labels \
            --BED data_bash_script/gene_bins/batch_${batch}.bed.gz \
            -o data_bash_script/covg_binned_${dset}/batch_${batch}.npz
    done
done

##########################
## Fit and apply models ##
##########################

## Single dataset, fit all batches in one go, no regressing out Pantry phenotypes
python $script prepare -i data_bash_script/covg_binned/ --bins-dir data_bash_script/gene_bins/ --n-batches 3 -o data_bash_script/ex1/covg_norm/
python $script fit -d data_bash_script/ex1/covg_norm/ -m data_bash_script/ex1/models/ --n-batches 3
python $script transform -d data_bash_script/ex1/covg_norm/ -m data_bash_script/ex1/models/ --n-batches 3 -o data_bash_script/ex1/latent_phenos.tsv.gz

## Single dataset, fit all batches in one go, regress out Pantry phenotypes
python $script prepare -i data_bash_script/covg_binned/ --bins-dir data_bash_script/gene_bins/ --n-batches 3 --pheno-paths-file data_input/pheno_files.txt -o data_bash_script/ex2/covg_norm/
python $script fit -d data_bash_script/ex2/covg_norm/ -m data_bash_script/ex2/models/ --n-batches 3
python $script transform -d data_bash_script/ex2/covg_norm/ -m data_bash_script/ex2/models/ --n-batches 3 -o data_bash_script/ex2/latent_phenos.tsv.gz

## Single dataset, fit batches separately
parallel -j$n_batches python $script prepare -i data_bash_script/covg_binned/ --bins-dir data_bash_script/gene_bins/ -b {} --pheno-paths-file data_input/pheno_files.txt -o data_bash_script/ex3/covg_norm/ :::: data_bash_script/batches.txt
parallel -j$n_batches python $script fit -d data_bash_script/ex3/covg_norm/ -b {} -m data_bash_script/ex3/models/ :::: data_bash_script/batches.txt
python $script transform -d data_bash_script/ex3/covg_norm/ -m data_bash_script/ex3/models/ --n-batches 3 -o data_bash_script/ex3/latent_phenos.tsv.gz

## Multiple datasets (prepare separately, train together, transform separately)
parallel -j$n_batches python $script prepare -i data_bash_script/covg_binned_dset1/ --bins-dir data_bash_script/gene_bins/ -b {} --pheno-paths-file data_input/pheno_files.txt -o data_bash_script/ex4/covg_norm_dset1/ :::: data_bash_script/batches.txt
parallel -j$n_batches python $script prepare -i data_bash_script/covg_binned_dset2/ --bins-dir data_bash_script/gene_bins/ -b {} --pheno-paths-file data_input/pheno_files.txt -o data_bash_script/ex4/covg_norm_dset2/ :::: data_bash_script/batches.txt
python $script fit -d data_bash_script/ex4/covg_norm_dset1/ data_bash_script/ex4/covg_norm_dset2/ --n-batches 3 -m data_bash_script/ex4/models/
python $script transform -d data_bash_script/ex4/covg_norm_dset1/ -m data_bash_script/ex4/models/ --n-batches 3 -o data_bash_script/ex4/latent_phenos.dset1.tsv.gz
python $script transform -d data_bash_script/ex4/covg_norm_dset2/ -m data_bash_script/ex4/models/ --n-batches 3 -o data_bash_script/ex4/latent_phenos.dset2.tsv.gz

## Many datasets (prepare separately, train together, transform separately)
Same as above but looping through datasets and passing a file listing the directories:
cat data_input/datasets.txt | while read dset; do
    python $script prepare -i data_bash_script/covg_binned_${dset}/ --bins-dir data_bash_script/gene_bins/ --n-batches 3 --pheno-paths-file data_input/pheno_files.txt -o data_bash_script/ex5/covg_norm_${dset}/
done
awk '{print "data_bash_script/ex5/covg_norm_"$1}' data_input/datasets.txt > data_bash_script/ex5/covg_norm_dirs.txt
cat data_bash_script/batches.txt | while read batch; do
    python $script fit --norm-covg-dir-file data_bash_script/ex5/covg_norm_dirs.txt -b $batch -m data_bash_script/ex5/models/
done
cat data_input/datasets.txt | while read dset; do
    python $script transform -d data_bash_script/ex5/covg_norm_${dset}/ -m data_bash_script/ex5/models/ --n-batches 3 -o data_bash_script/ex5/latent_phenos.${dset}.tsv.gz
done
