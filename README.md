# geneIS-nf

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

## Introduction

**geneIS-nf** is a bioinformatics analysis pipeline for assigning initiation sites, mapped from NS-seq data (see [`iniseq-nf`](https://github.com/pavrilab/inisite-nf) and [`classifyIS-nf`](https://github.com/pavrilab/classifyIS-nf)) their nearest gene in a reproducible way.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

## Pipeline summary

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install the [`pandas`](https://pandas.pydata.org/docs/index.html), [`numpy`](https://numpy.org/), [`scipy`](https://www.scipy.org/) and [`matplotlib`](https://matplotlib.org/) and [`BioPython`](https://github.com/biopython/biopython) Python packages and the [`argparser`](https://cran.r-project.org/web/packages/argparser/index.html), [`BiocManager`](https://github.com/Bioconductor/BiocManager), [`GenomicFeatures`](http://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) and [`ChIPseeker`](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html) R packages

iii. Clone repository with 
```bash
nextflow pull pavrilab/geneIS-nf
```

iv. Start running your own analysis!
```bash
nextflow run pavrilab/geneIS-nf --masterTable IS.master.tsv --txDb annotation.sql --email e@mail.com --xCol WT_col --yCol KD_col
```

## Main arguments
#### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. For example `-profile cbe` invokes the execution of processes using the [`slurm`](https://slurm.schedmd.com/documentation.html) workload manager. If no profile is given the pipeline will be executed locally.

#### `--masterTable`
Mastertable containing at least quantification results (see [`classifyIS-nf`](https://github.com/pavrilab/classifyIS-nf))

#### `--txDb`
SQL file containing a gene annoation generated with the GenomicFeatures R-package

#### `--email`
Email address to fetch EntrezIDs with BioPython

#### `--xCol`
column in the mastertable holding SNS-seq read quantification results for WT

#### `--yCol`
column in the mastertable holding SNS-seq read quantification results for KD

## Generic arguments
#### `--foldChange`
adds diagonal lines in distance of foldChange to plot

#### `--axMin`
Lower bound for axis values of x- and y-axis (Default: 0)

#### `--axMax`
Upper bound for axis values of x- and y-axis (Default: 8)

#### `--filePrefix`
Prefix for the result files name

#### `--outputDir`
Folder to which results will be written (is created if not existing)

## Results
The pipeline generates five result files:

1.  The `*.master.tsv` file is a tab-separated file with the basic layout of a BEDfile containing the genomic coordinates of the peaks their names and quantification results for either condition and the assigned class
2.  The `*.chipseeker.mapped.tsv` file is a tab-separated file containing the annotation results
3.  The two PDF files `*.density.pdf` and `*.class.pdf` are the visualization of the quantification results. The density plot show the distribution of peaks in the quantification space and the class plot visualizes the class assignment results

## Credits
The pipeline was developed by [Daniel Malzl](mailto:daniel.malzl@gmx.at) for use at the [IMP](https://www.imp.ac.at/), Vienna.

Many thanks to others who have helped out along the way too, including (but not limited to): [@t-neumann](https://github.com/t-neumann), [@pditommaso](https://github.com/pditommaso).

## Citations
### Pipeline tools
* [Nextflow](https://www.ncbi.nlm.nih.gov/pubmed/28398311/)
  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.
  
### Python Packages
* [BiocManager](https://github.com/Bioconductor/BiocManager)

* [GenomicFeatures](http://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
  > Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M, Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” PLoS Computational Biology, 9. doi: 10.1371/journal.pcbi.1003118
  
* [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html)
  > Yu G, Wang L, He Q (2015). “ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization.” Bioinformatics, 31(14), 2382-2383. doi: 10.1093/bioinformatics/btv145

* [BioPython](https://github.com/biopython/biopython)

* [pandas](https://pandas.pydata.org/docs/index.html)
  > Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)
  
* [numpy](https://numpy.org/)
  > Stéfan van der Walt, S. Chris Colbert and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011). doi: 10.1109/MCSE.2011.37

* [scipy](https://www.scipy.org/)
  > Virtanen, P., Gommers, R., Oliphant, T.E. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods 17, 261–272 (2020). doi: 10.1038/s41592-019-0686-2
  
* [matplotlib](https://matplotlib.org/)
  > John D. Hunter. Matplotlib: A 2D Graphics Environment, Computing in Science & Engineering, 9, 90-95 (2007). doi: 10.1109/MCSE.2007.55
