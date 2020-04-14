# CovProfile

CovProfile: profiling the viral genome and gene expressions of
SARS-COV-2

## Getting Started
To get CovProfile, run the following command:
```
    git clone https://github.com/yonghanyu/CovProfile.git
```

To see the options and usage for CovProfile:
```
    cd CovProfile
    python3 CovProfile.py --help
```

We recommend to use CovProfile with SARS-COV-2 reference file and gene annotation file from NCBI database

CovProfile currently support PCR captured Nanopore sequenced reads. The program require PCR primer information to run. Prepare your primer_information file with tsv format wtih the following header
```
name    pool    seq    length
``` 

### Prerequisites
python version 3.5+ is required with package pysam, numpy, scipy and BCBio. To install them with pip, run the following command:
```
pip3 install pysam numpy scipy BCBio 
```
Several software is required for CovProfile to perform alignment, processing bam file and call SNP:
* [Minimap2](https://github.com/lh3/minimap2) 
* [Samtools](http://samtools.sourceforge.net/)
* [Bcftools](https://samtools.github.io/bcftools/)
* [Blat](https://genome.ucsc.edu/FAQ/FAQblat.html)

##Cite CovProfile
If you use CovProfile in your work, please cite with doi: doi.org/10.1101/2020.04.05.026146


## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details