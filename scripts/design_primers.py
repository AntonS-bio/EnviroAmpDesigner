#from typing import Tuple, List, Dict
import pickle
from os.path import exists, expanduser, join
from os import makedirs, remove, listdir
from shutil import which
from sys import exit
from data_classes import SNP, Genotype, Genotypes, InputConfiguration
from inputs_validation import ValidateFiles
from load_vcfs import VCFutilities
from identify_genotype_snps import GenotypeSnpIdentifier
import name_converters
from snp_optimiser import SnpOptimiser
from identify_species_snps import IdentifySpeciesSnps
from primers_generator import PrimersGenerator
import metadata_utils
import argparse
import warnings

def _check_inputs(config_data: InputConfiguration):

    ### START Input validation tests ####
    print("Validating input files")
    file_validator=ValidateFiles()
    if not file_validator.validate_fasta(config_data.reference_fasta):
        exit(1)

    if not exists(config_data.existing_primers):
        warnings.warn(f'Existing primers file {config_data.existing_primers} does not exist.')
        exit(1)

    if config_data.specific_target_snps!="":
        if not file_validator.validate_specific_target_snps(config_data.specific_target_snps):
            exit(1)
        
    if config_data.repeats_bed_file!="":
        if not file_validator.validate_bed(config_data.repeats_bed_file):
            exit(1)
        if not file_validator.contigs_in_fasta(config_data.repeats_bed_file, config_data.reference_fasta):
            exit(1)
    
    if not file_validator.validate_vcf(config_data.vcf_dir):
        exit(1)
    
    if not file_validator.validate_hierarchy(config_data.hierarchy_file, config_data):
        exit(1)
    
    if not file_validator.validate_negative_genomes(config_data.negative_genomes):
        exit(1)
        
    if not metadata_utils.load_metadata(config_data):
        exit(1)
    print("Input validation complete. ")
    #### END Input validation tests ####

def _parse_arguments():
    parser = argparse.ArgumentParser(description='Generate list of SNPs that uniquely identify one or more genotypes')
    parser.add_argument('-c','--config_file', metavar='', type=str,
                        help='Config file, see sample.json for example.', required=True)
    parser.add_argument('-m','--mode', metavar='', type=str, choices=['SNP', 'Amplicon'],
                        help='"SNP" to only get genotype defining SNPs, "Amplicon" to also generate amplicons', required=True)

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        exit(0)
    return args

def _check_tools() -> bool:
    for command, name in zip(["mafft", "makeblastdb", "blastn", "minimap2"], ["MAFFT", "NCBI BLAST", "NCBI BLAST", "minimap2"]):
        if which(command) is None:
            print(f'Missing {name} program. It is available via Conda.')
            exit(0)
    return True

def _setup_analysis(config_data: InputConfiguration):
    #create directory for output is does not exit
    if exists(config_data.msa_dir):
        warnings.warn("MSA directory exists, removing all .fasta and .fna to save memory")
        [remove( join(config_data.msa_dir,f) ) for f in listdir(config_data.msa_dir) if f.split(".")[-1]=="fna" or f.split(".")[-1]=="fasta"]
    for value in [config_data.output_dir, config_data.msa_dir, config_data.temp_blast_db]:
        if not exists(value):
            makedirs(value)
    
    #### START Identify genotype defining SNPs ####
    for stub in config_data.name_stubs:
        name_converters.name_stubs.add(stub)

def _load_specific_target_snps(config_data: InputConfiguration, genotypes: Genotypes):
    """
    Use "specific_target_snps" field from config to load the reference sequence positions 
    at which a specific SNP needs to be targeted. As example are genes with AMR mutations which are 
    not genotype specific.

    :param config_data: the configuration data for the run
    :type config_data: InputConfiguration

    :param genotypes: List of genotypes that are targeted by the run
    :type config_data: Genotypes
    """
    with open(config_data.specific_target_snps) as snps_file:
        for line in snps_file:
            contig, position, name=line.strip().split()
            extra_genotype=Genotype(name)
            #The ref and alt base are irrelevant as specific position is targeted for amplicon capture
            target_snp=SNP(ref_contig_id=contig, ref_base="T", alt_base="C", position=int(position),passes_filters=True)
            target_snp.is_genotype_snp=True
            target_snp.sensitivity=1
            target_snp.specificity=1
            #the actual depth is irrelevant
            extra_genotype.add_genotype_allele(target_snp, target_snp.alt_base, depth=100)
            genotypes.genotypes.append(extra_genotype)
            config_data.gts_with_few_snps.append(name)

def _identify_genotype_SNPs(config_data: InputConfiguration):

    snp_identifier=GenotypeSnpIdentifier(config_data)

    genotypes: Genotypes = snp_identifier.identify_snps()

    #genotypes.genotypes_to_snp_matrix().to_csv(config_data.genotype_snps, sep="\t", index=False)

    with open(config_data.genotype_snps, "w") as snps_file:
        for genotype in genotypes.genotypes:
            for snp in genotype.defining_snps:
                if snp.passes_filters:
                    fourth_col=f'SNP:{snp.ref_base}/{snp.alt_base}/GT:{genotype.name}/SP:{snp.specificity:.2f}/SE:{snp.sensitivity:.2f}'
                    snps_file.write("\t".join( [str(f) for f in [snp.ref_contig_id, snp.position, snp.position+1, fourth_col] ])+"\n")
    return genotypes

def _optimise_snps(config_data: InputConfiguration, genotypes: Genotypes):
    snp_opimiser=SnpOptimiser()
    max_iterval_len=config_data.max_amplicon_len
    gts_with_few_snps=config_data.gts_with_few_snps+ [f for f in genotypes.genotypes if len(f.defining_snps)<=10]
    amplicon_intervals=snp_opimiser.optimise(max_iterval_len,genotypes, rare_gts=gts_with_few_snps)
    if len(amplicon_intervals)==0:
        print("No intervals with multiple genotypes were identified and none of the genotypes are listed are rare. Add genotypes to 'gts_with_few_snps' in config. Exiting.")
        exit()
    
    with open(config_data.multi_gt_intervals, "w") as output_bed:
        for interval in amplicon_intervals:
            interval_start=min([f.position for f in interval["snps"]])
            interval_end=max([f.position for f in interval["snps"]])
            contig_id=interval["snps"][0].ref_contig_id
            gts="_".join(interval["genotypes"])
            output_bed.write('\t'.join([contig_id, str(interval_start),str(interval_end),gts])+"\n")


def main():

    _check_tools()
    args=_parse_arguments()
    
    run_mode=args.mode
    config_file=expanduser(args.config_file)
    if not exists(config_file):
        print(f'File {config_file} not found, please check spelling.')
        exit(0)

    try:
        config_data = InputConfiguration( config_file )
    except IOError as error:
        print(error)
        exit(1)

    _check_inputs(config_data)

    _setup_analysis(config_data)

    genotypes: Genotypes = _identify_genotype_SNPs(config_data)

    #For amplification of AMR genes, the examination of the target species is not required
    #however, we need to identify the SNP that uniquely identify the target species

    _load_specific_target_snps(config_data, genotypes)
   
    ### DEBUG command
    # with open(config_data.genotypes_data, "wb") as output:
    #     pickle.dump(genotypes, output)
    ### DEBUG command

    ### DEBUG command
    # with open(config_data.genotypes_data, "rb") as pickled_file:
    #     genotypes: Genotypes = pickle.load(pickled_file)
    ### DEBUG command

    _optimise_snps(config_data, genotypes)
    # # ### END optimise the selection of SNPs ####

    # # #### START MSA generation section ####

    # vcf_utils=VCFutilities()
    # vcf_utils.output_genotypes_vcf(genotypes, config_data.gt_snps_vcf)
    if run_mode!="Amplicon":
        exit(0)

    snp_identifier=IdentifySpeciesSnps(ref_fasta=config_data.reference_fasta,
                                    msa_dir=config_data.msa_dir,
                                    negative_genomes_dir=config_data.negative_genomes,
                                    temp_blast_db_dir=config_data.temp_blast_db,
                                    amplicons_bed=config_data.multi_gt_intervals)


    species_genotype: Genotype=snp_identifier.generate_flanking_amplicons()

    # ### DEBUG command
    # with open(config_data.output_dir+"/species_gt.pkl", "wb") as output:
    #     pickle.dump(species_genotype, output)
    # with open(config_data.output_dir+"/species_gt.pkl", "rb") as pickled_file:
    #     species_genotype: Genotype = pickle.load(pickled_file)
    # ### DEBUG command
    flanking_amplicons=snp_identifier.get_bifurcating_snps(species_genotype)

   
    # ### DEBUG command
    # with open(config_data.species_data, "wb") as output:
    #     pickle.dump(flanking_amplicons, output)
    # ### DEBUG command


    # #### END MSA generation section ####

    ### DEBUG command
    # with open(config_data.species_data, "rb") as pickled_file:
    #     species: Genotype = pickle.load(pickled_file)
    # with open(config_data.genotypes_data, "rb") as pickled_file:
    #     genotypes: Genotypes = pickle.load(pickled_file)

    ### DEBUG command

    target_gts= [genotype.name for genotype in genotypes.genotypes]
    genotypes.genotypes.append(flanking_amplicons)

    genotypes.get_duplicate_snps()

    vcf_utils=VCFutilities()
    vcf_utils.output_species_vcf(genotypes, config_data.gt_species_snps_vcf)

    #### START Generate primers ####
    generator=PrimersGenerator(config_data)
    generator.genotypes = genotypes
    #generator._for_testing_load_gts(config_data.species_data, config_data.genotypes_data)
    generator.find_candidate_primers(target_gts)

    with open(config_data.output_dir+"primers.tsv","w") as output_file:
        header="\t".join(["Name", "Forward Species SNPs", "Reverse Species SNPs",
                        "Penalty", "Contig", "Start","End","Length",
                        "Forward","Forward Tm", "Forward GC",
                        "Reverse","Reverse Tm", "Reverse GC"])+"\n"
        output_file.write(header)
        for pair in generator.new_primer_pairs:
            output_file.write(pair.to_string()+"\n")
    exit()

main()

