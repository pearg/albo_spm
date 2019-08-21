# Albopictus sex prediction model

## Overview 

This program reads a set of FASTQ files from a paired-end ddRAD sequencing 
Aedes albopictus sample (nlaIII and mluCI restriction enzymes) and predicts 
the sex of the sample.

## Installation

Make sure that the following are accessable in your $PATH:

* samtools
* bedtools
* bowtie2

The following R libraries are also required:

* optparse
* randomForest
* ps

Clone the github repository:

```
git clone https://github.com/pearg/albo_spm.git
```

## Usage

### Help message

Usage information can be displayed using the `-h` or `--help` argument:

```
$ ./albo_spm.R -h
Usage: ./albo_spm.R [OPTIONS] --r1 [R1_FASTQ] --r2 [R2_FASTQ] --output [OUTPUT_FILENAME]
Predict sex of Aedes albopictus from ddRAD sequencing (nlaIII and mluCI restriction enzymes).


Options:
	--r1=R1_FASTQ
		R1 FASTQ file [REQUIRED]

	--r2=R2_FASTQ
		R2 FASTQ file [REQUIRED]

	--output=OUTPUT_FILENAME
		Output filename [REQUIRED]

	--sample_name=SAMPLE_NAME
		Sample name

	--ctrl_depth
		Output control depth of coverage to output

	--save_rdata_file=RDATA_FILENAME
		RData filename to save objects for debugging

	--threads=N
		Number of threads to use [default: 2]

	--version
		Print version and exit

	-h, --help
		Show this help message and exit

```

### Example usage

Example:

```
albo_spm.R \
    --r1 sample_01.R1.fq \
    --r2 sample_01.R2.fq \
    --output sample_01.sex.tsv \
    --sample_name sample_01 \
    --threads 5
``` 

Example output:
```
$ cat sample_01.sex.tsv

sample	class	probability
sample_01	F	0.955
```
