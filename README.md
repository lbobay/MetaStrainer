# MetaStrainer

MetaStrainer is a python tool used for recovering strain-level genotypes of bacterial genomes from shotgun metagenomic datasets.

## Requirements

### Programs
- samtools
- bowtie2

### Python libraries
- emcee
- numpy
- biopython

## Installation


## Running

Basic usage:
```bash
python MetaStrainer_master.py -1 Example_trimmed_R1.fastq.gz -2 Example_trimmed_R2.fastq.gz -r Genome_fullgenome.gbff -f region -o outputfolder 
```

Strain threshold is at 99.5% whole-genome average nucleotide identity (ANI) between the generated strains.
Use a lower or higher threshold with the ``` -s/--StrainThreshold``` option.

Advanced usage:

## Citation
