# MetaStrainer

MetaStrainer is a python tool used for recovering strain-level genotypes of bacterial genomes from shotgun metagenomic datasets.

## Dependencies

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

MetaStrainer expects a quality controlled Fastq read pair files.

The reference should be a full Genbank annotation file (.gbff) having DNA sequences included. 

Advanced usage:

To improve mapping, additional flanking regions are extended to upstream and downstream of the genes. The recommended setting is equal to the sequenced length of a read (e.g. 150bp). This could be set using the ```-f/--READLENGTH``` option based on the sequenced library.


Strain threshold is at 99.5% whole-genome average nucleotide identity (ANI) between the generated strains.
Use a lower or higher threshold with the ``` -s/--StrainThreshold``` option.



## Citation
