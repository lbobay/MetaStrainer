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


## Important Notes

The community composition of the sample should be characterized before running MetaStrainer using appropriate taxonomic classification methods. If more than one closely related species, under the same genus, are present it strongly advised to use the right species genome that is closest to the ones close to the samples. While MetaStrainer can stop processing if a distant species is used, a closely related but wrong species would still map and cause problems with strain inferences due to issues such as cross-mapping. 

Some unique environments are not suffeciently represented in general databases. It may be useful to prepare a custom reference database for proper taxonomic characterization of non-model organisms or unique niches in host-specific environments.


## Citation

Sharaf H, Bobay L-M. MetaStrainer: Accurate reconstruction of bacterial strain genotypes from short-read metagenomic samples. bioRxiv. 2026. doi:10.64898/2026.03.02.709061. https://www.biorxiv.org/content/10.64898/2026.03.02.709061v1
