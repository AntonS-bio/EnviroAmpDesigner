# EnviroAmpDesigner

This tool designs primers for multiplex PCR of environmental or similarly contaminated samples to detect and genotype target organism (e.g. *Salmonella* Typhi or *Salmonella* Typhimurium). While tool can design primer of any length, they are intended for Oxford Nanopore devices. The tool has been tested for bacteria, and may struggle with large genomes such as fungal or protists.

<br/>
EnviroAmpDesigner requires three core inputs:

A) User generated VCF files where each file is a sample of target organism mapped to user chosen reference sequence

B) User defined genotypes for each VCF file (ex. GenoTyphi for *Salmonella* Typhi or MLST for invasive non-typhoidal *Salmonella*). 

C) User chosen list of assemblies of organism that primers **should not target** - this can be easily obtained via "datasets" tool from NCBI and for *Salmonella* Typhi can be set of non-Typhi *Salmonella* and other Enterobacterales.

<br/>
EnviroAmpDesigner works in four stages:

First, it uses VCFs from (A) and genotypes data from (B) to identify which SNPs are uniquely associated with genotypes

Second, it compares target organism reference sequence around SNPs identified in previous step to the assemblies in (C) to identify homologues in non-target organism.

Third, EnviroAmpDesigner looks for sites in homologues identified in previous step where target and non-target organisms differ by one or more nucleotides. 

Fourth, it generates primers pairs where at least one 3' primer end is placed at sites identified at previous step. Mismatches at primer 3' end have high impact on primer efficiency and such placement improves detection of target organism in contaminated samples.

## Intended Workflow

Design of multiplex panel is an iterative process that consist of design of primers *in silico*, testing of primers in the lab, evaluation of result, design of new primers *in silico* to cover the targets of failed primers, testing new primers in the lab, etc. To illustrate this process in detail, this is how we use it:
1. Identify genotypes of interest in S. Typhi (21 target genotypes)
2. Run EnviroAmpDesigner specifying all 21 target genotypes in "hierarchy_file" as well as all of their descendant genotypes
3. Take produced list of primer and select from them those of roughly equal length and melting temperatures.
4. Test the primers in the lab using both pure culture, environmental and pure culture spiked with off-target organisms [Add J's paper reference later]
5. Analyse amplification results
6. Keep primers that worked
7. Select replacement primer pairs
  - for genotype target that have alternative primers in (3), use those alternatives or
  - for targets without alternative in (3) rerun EnviroAmpDesigner specifying only the these targets in "hierarchy_file" and relaxing search parameters using combination of "snp_specificity", "snp_sensitivity", "max_matching_negative_genomes" and "flank_len_to_check"
8. repeat from (4). However, to minimise number of iterations, for each genotype target in (7) test 3-5 primer pairs in each iteration instead of one pair per target per iteration.

## Installation

At present, the tool can be downloaded via GitHub, but it will soon be available as a bioconda package. 
Dependencies:
--mafft
--blast
--minimap2
--primer3-py (python wrapper for Primer3)

## Running
The tool requires two input arguments and can be run as follows:
```
python run.py -m Amplicon -c config.json
or
python run.py -m SNP -c config.json
```
"SNP" option only identified and reports the number of SNPs that differentiate genotypes, whereas Amplicon also attempts to design amplicons for these SNPs.

Due to large number of options and to improve reproducibility most inputs are specified via a JSON file (config.json above) an example of which is in this repository "sample_files" directory.

### JSON input file
The fields in JSON config are:

  
  "name_stubs": a optional list of suffixes by which VCF files differ from samples names in file with genotypes. For example, if VCF names are Sample123.sorted.vcf, Sample456.sorted.vcf, etc. and file with genotypes has "Sample123" and "Sample456", you should add [".sorted"] to this field.
  
  "max_cpus": optional field - number of CPUs to use
  
  "reference_fasta": FASTA sequence to which the VCFs where mapped<br/>
  
  "repeats_bed_file": optional list of regions (as .BED file) on FASTA which should be excluded from analysis. Mainly used to exclude repeat regions.
  
  "hierarchy_file": list of genotypes in JSON format showing the how the genotypes are related. See example in "sample_files" directory.
  
  "meta_data_file": The file where first column contains the samples names and some column (specified in "genotypes_column" field below) contains genotype names.
  
  "existing_primers": Optional list of existing primers, these will be used to check that new primers don't interfere with existing ones.
  
  
  "vcf_dir": directory with VCF files for target organism.
  
  "negative_genomes": directory with assemblies of organism that primers **should not target**. See (C) at the top of this README.
  
  "use_negative_genomes_subdir": True/False - instructs EnviroAmpDesigner to check subdirectories of "negative_genoes". Useful if genomes were downloaded using NCBI datasets.
  
  "temp_blast_db": directory for temporary files
  
  
  "delimiter": separator (usually "," or "\t") for columns in "meta_data_file"
  
  "genotype_column": name of column in "meta_data_file" which contains the genotype values. Eg. if column in file is "Final_Genotypes", specify "Final_Genotypes" in this field"
  
  
  "snp_specificity": Number between 0 and 100. Sometimes, there are no SNPs that perfectly separate genotypes or perhaps very few such SNPs, but without viable primers. This allows to relax specificity of SNPs that define genotypes. 
  
  "snp_sensitivity": Number between 0 and 100, Same as "snp_specificity", only sensitivity.
  
  "gts_with_few_snps": by default, EnviroAmpDesigner attempts to capture multiple genotypes with a single primer by targeting closely located genotype defining SNPs, but for some genotypes this may produce suboptimal results. This options forces the tool to use all SNPs 
  for genotypes specified here. Eg. if the tool should use all SNPs for genotypes 3.2.1, 2.4.9, and 4.3.1 the field should be  ["3.2.1", "2.4.9", "4.3.1"],
  
  "flank_len_to_check": Integer >0, but ideally >500. This is the length of region upstream and downstream of genotype defining SNPs that will be checked for homologues among off-target organisms. 
  
  "min_amplicon_length": Integer >0, minimum length of amplicon
  
  "blast_e_value": value betweem 0 and 1, e-value cut-off to use when looking for homologues among off-target organisms
  
  "blast_word_size": Integer >11, BLASTn minimum word size when looking for homologues among off-target organisms
  
  "max_matching_negative_genomes": Number between >=0. When EnviroAmpDesigner is looking for nucleotides that distinguish target and off-target organisms, sometimes there isn't nucleotide that perfectly separates them perfectly. This specifies how many off-target organisms can have the same nucleotide as target organisms at position X for position X to still be valid site for 3' end of primers. Relaxing this potentially make primers less discriminating, but increases number of possible primers due to higher number of place the 3' end can be position.
  
  
  "output_dir": Directory for outputs.
  
  "genotype_snps": List of SNPs that were identified as unique to some genotypes.
  
  "snps_bed": BED files with identified SNPs.
  
  "genotypes_data": For testing only.
  
  "species_data": For testing only.
  
  "multi_gt_intervals": BED file specifying intervals which contain multiple genotype defining SNPs or SNPs defining genotypes listed in "gts_with_few_snps". To minimise the number of primers required, EnviroAmpDesigner will design primers only for these regions.
  
  "msa_dir": Directory for Multiple Sequence Alignment files output
  
  "genoptype_snps_vcf": VCF file with genotype defining SNPs
  
  "gt_and_species_snps_vcf": VCF file with genotype and target organism defining SNPs
  
  "PRIMER_OPT_SIZE": Integer >0, but ideally >19, optimal size of primer, these parameters are for Primer3, please check "https://primer3.org/manual.html"
  
  "PRIMER_OPT_TM": Decimal number between 0 and 100, but ideally between 50 and 65, target primer melting temperature
  
  "PRIMER_MIN_TM": Decimal number between 0 and 100, but ideally few degrees below "PRIMER_OPT_TM", min primer melting temperature
  
  "PRIMER_MAX_TM": Decimal number between 0 and 100, but ideally few degrees above "PRIMER_OPT_TM", max primer melting temperature
