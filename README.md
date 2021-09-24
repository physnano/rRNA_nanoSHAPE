# rRNA_nanoSHAPE
Code accompanying "Direct detection of RNA modifications and structure using single molecule nanopore sequencing"

# nanoSHAPE analysis:

In vitro nanoSHAPE analysis can be performed on a single target using the commands below. First the multi-fast5 files must be converted to single read fast5 files prior to processing with tombo. Tombo is then used to resquiggle the reads, detect modifications, and generate per-read stats file and optional plot. Comprehensive documentation regarding tombo can be found [here](https://nanoporetech.github.io/tombo/). Finally, the <em>extract_array.py</em> script is used to generate an output numpy array which can be easily processed using python for per-read modification visualization/plotting. 

## Step 1: convert the multi fast5 files to single fast5 files

	multi_to_single_fast5 -i <directory of multi fast5s> --recursive --save_path <output path> --flowcell <flowcell spec> --kit <kit spec> --fast5_out --cpu_threads_per_Caller <num threads>

## Step 2: Resquiggle reads using tombo

	tombo resquiggle <path to the directory containing single fast5 files>  <path to fasta reference> --include-event-stdev --num-most-common-errors 5 --ignore-read-locks --overwrite --processes <num processes>

## Step 3: Detect modifications using tombo sample compare

	tombo detect_modifications model_sample_compare --fast5-basedirs <path to single fast5 files MODIFIED> --control-fast5-basedirs <path to single fast5 files CONTROL/UNMODIFIED> --statistics-file-basename <basename for stats output> --minimum-test-reads <num min test reads, typically 100 for in vitro> --multiprocess-region-size 250

## Step 4: Generate per-read stats file

	tombo detect_modifications model_sample_compare --fast5-basedirs <path to single fast5 files MODIFIED> --control-fast5-basedirs <path to single fast5 files CONTROL/UNMODIFIED> --statistics-file-basename <basename for stats output> --per-read-statistics-basename <basename for per-read stats output> --minimum-test-reads <num min test reads, typically 100 for in vitro> --multiprocess-region-size 250

## Step 5: Tombo plot of per-read stats (modifications)

	tombo plot per_read --genome-locations <genomic location of target> --per-read-statistics-filename <per-read statistics filename generated in step4> --genome-fasta <FASTA reference file> --fast5-basedirs <path to single fast5 files MODIFIED> --num-reads <number of reads> --num-bases <num bases> --pdf-filename <output pdf>

## Step 6: Custom python script for extracting modifications sites as array

	python extract_array.py --fa <FASTA reference file> --f5 <path to single fast5 files to process, either MODIFIED or UNMODIFIED> --out <numpy array out path with extension ".npy"> --rna-name <RNA name from FASTA reference> --sort-reads --current (or --dwell, --stdev)




