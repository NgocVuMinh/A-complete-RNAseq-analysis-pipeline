# A complete RNAseq analysis pipeline

### [Ngoc M. Vu](https://github.com/NgocVuMinh)
This repository contains command line scripts and R scripts used to process RNA sequencing data. This also serves as the supplementary materials in support of my undergraduate thesis titled "RNAseq analysis of breast cancer cells co-cultured with adipose-derived stem cells" and an in-submission research paper (DOI will be updated soon). 

Two sample types were included in this analysis:
* CT - control samples: cancer cells cultured in standard condition
* CC - co-cultured samples: cancer cells co-cultured with adipose-derived stem cells

Details of the experimental design will be updated.
The number of sample types does not affect the command line scripts, but adjustments to the R scripts must be if more than 2 sample types are applied.

## Pipeline overview
<div style="text-align:center;">
  <img src="https://github.com/NgocVuMinh/RNAseq-analysis-pipeline-from-scratch/blob/main/overview.png" width="450" />
</div>

## Bulk RNA sequencing data

The sample FASTQ files will be provided once our paper is published. The original data was provided by the Department of Bioinformatics, Institute of Biotechnology, Vietnam Academy of Science and Technology (VAST).


## Sequence processing pipeline

All software was executed on a local MacOS machine, with Conda 23.7.3 and Python 3.11.

The first pipeline processes the raw RNAseq data (FASTQ files). Three MAIN steps along with the required packages and softwares were as follows:
1. **Quality control - Assessing the quality of reads in each FASTQ file**: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used. [MultiQC](https://multiqc.info/docs/getting_started/installation/) was then used to merge the quality control output from all samples. After quality assessment, trimming low-quality bases at the 5'- and 3'- ends of the sequenced reads might be needed. In that case, [Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html) is the tool to go with.

2. **Read alignment - Mapping the sequenced reads to the reference genome**: [HISAT2](https://daehwankimlab.github.io/hisat2/manual/)
This step requires building the genomic or transcriptomic index from scratch, or alternatively, downloading the prebuilt index from the HISAT2 manual:  [*GRCh38* Homo sapiens - genome](https://daehwankimlab.github.io/hisat2/download/#h-sapiens) (see more at [readme.md](https://github.com/NgocVuMinh/RNAseq-analysis-pipeline-from-scratch/tree/main/ref/readme.md)). HISAT2 produces SAM files as the output, which should then be converted to BAM format using [Samtools](https://www.htslib.org/doc/samtools.html). The mapping results (including gene coverage and read distribution) could be assessed using [Qualimap](http://qualimap.conesalab.org/) and [RSeQC](https://rseqc.sourceforge.net/#read-distribution-py).

3. **Gene quantification**: using BAM files, [featureCounts](https://subread.sourceforge.net/featureCounts.html) identifies the gene annotation and quantify the number of reads for each gene. The annotation file for the genome *GRCh38.p14* was used. This step produces individual output file for each sample. 

### Installation and set-up

Firstly, install all softwares and tools. The list of softwares used in this pipeline is provided in the table below. You may click on the attached URLs for installation instructions:

| **Softwares** | **Version** | **Purpose** |
|--------------|--------------|--------------|
| **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** | 0.12.1 | Assessing the quality of sequenced reads in each sample |
| **[Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)** | 4.8 | Trimming and filtering low-quality reads |
| **[HISAT2](https://daehwankimlab.github.io/hisat2/manual/)** | 2.2.1 | Mapping reads to the reference genome or transcriptome |
| **[Samtools](https://www.htslib.org/doc/samtools.html)** | 1.19.2 | Converting SAM files to BAM files, index and sort BAM files |
| **[Qualimap (RNAseq)](http://qualimap.conesalab.org/)** | 2.3 | Assessing quality of reads and mapping results |
| **[RSeQC](https://rseqc.sourceforge.net/#read-distribution-py)** | 5.0.1 | Assessing mapping results |
| **[MultiQC](https://multiqc.info/docs/getting_started/installation/)** | 1.21 | Merging outputs from FastQC, Qualimap and RSeQC |
| **[featureCounts](https://subread.sourceforge.net/featureCounts.html)** | 2.0.6 | Gene annotation and quantification |

Then download and process all indices, references and other tools:

```
./download.sh
```
### Run the entire pipeline
```
./run_processing.sh
```

Finally, merge all output files from featureCounts. The merged count file (featureCounts_merged.txt) was then used for gene expression analysis in R.

```
cd featureCounts/

for file in *clean.txt; do tail -n +2 "$file" > "${file%.txt}2.txt"; done # remove the first line in the output file

Rscript ./Rscripts/merge_counts.R
```


## Gene expression and functional analysis in R

All scripts were run in R-Studio, R version 4.3.1.

The Bioconductor’s *DESeq2* package was used to normalize count data, estimate dispersions, and identify differentially expressed genes. Genes with count lower than 10 were filtered out. A gene is considered differentially expressed in the experimental group compared to the control group if its *adjusted p-value* ≤ 0.05 and | *log2foldchange* | ≥ 1.

Gene Ontology (GO) terms and Kyoto Encyclopedia of Genes and Genomes (KEGG) enrichment analysis were implemented to assess biological processes and pathways. The R Bioconductor’s *clusterProfiler* package was used to retrieve the enriched GO categories and KEGG pathways associated with the differentially expressed geneset. Functions from the *Pathview* package and *enrichplot* package were used to visualize the enriched pathways.

Visualization methods:
* Volcano plot
* Heatmap
* Category netplot (cnetplot)
* GSEA plot
* etc.

R packages and libraries used:

```
library(DESeq2)
library(tidyverse)
library(airway)
library(apeglm)
library(pheatmap)
library(RColorBrewer)
library(stringr) 
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(gage)
library(pathview)
library(gageData)
library(enrichplot)
```

R scripts:

```
count_processing.R      # running DESeq2
visualization.R         # visualization
heatmap.R               # heatmap
gene_ontology.R         # GO analysis
kegg.R                  # KEGG analysis
```

## Tutorials and other references

The ComplexHeatmap package documentation:

* https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html

Differential gene expression workshop (Archived):

* https://hbctraining.github.io/DGE_workshop/lessons/09_functional_analysis.html

Part II: Enrichment analysis, The Biological Knowledge Mining Book - Bioinformatics Group @ SMU:

* https://yulab-smu.top/biomedical-knowledge-mining-book/clusterProfiler-dplyr.html

* https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html 

Youtube Tutorials:

* https://youtube.com/playlist?list=PLi1VnGoeDGjvHvl83QySD2oAQYFHPRYso&si=39DEmrLNNskfZqzQ 

* https://youtube.com/playlist?list=PLJefJsd1yfhYa97VUXkQ4T90YH6voYjCn&si=neiWOYiqvyvelmcr 