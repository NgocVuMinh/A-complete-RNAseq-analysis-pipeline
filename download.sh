#!bin/bash
# This script is used to retrive indices, references and tools
# Usage: ./download.sh


# Download reference genome:
mkdir -p ref
cd ref
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
gzip -d grch38_genome.tar.gz
mv grch38_genome gh38_index

# Download annotation of the reference genome
wget ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz  
gzip -d Homo_sapiens.GRCh38.84.gtf.gz
mv Homo_sapiens.GRCh38.84.gtf gh38.gtf

# Or you could use the reference transcriptome instead:
# wget https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz


### Converting GTF to BED12
## The following part of the script is adapted from https://gist.github.com/gireeshkbogu/f478ad8495dca56545746cd391615b93
## The credit goes to author: https://gist.github.com/gireeshkbogu 

# First, download UCSC scripts
wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.arm64/gtfToGenePred
wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.arm64/genePredToBed
# wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.arm64/bedToBigBed

# Second, download chromosome sizes and filter out unnecessary chromosomes
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
grep -v chrM hg38.chrom.sizes| grep -v _hap | grep -v Un_gl |grep -v random > hg38.chrom.filtered.sizes
rm hg38.chrom.sizes

# Third, make them executable
chmod +x gtfToGenePred genePredToBed # bedToBigBed

# Convert Gtf to genePred
./gtfToGenePred gh38.gtf gh38.genePred

# Convert genPred to bed12
./genePredToBed gh38.genePred gh38.bed12

# Convert sorted bed12 to bigBed (useful for trackhubs)
# ./bedToBigBed gh38.sorted.bed hg19.chrom.filtered.sizes gh38.bb