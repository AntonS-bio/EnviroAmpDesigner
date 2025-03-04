### Checks all the inputs are valid before executing the program
### this avoids program crashing mid-way due to bad inputs.

from os.path import exists
from os import listdir
import warnings
from tqdm import tqdm
from Bio import SeqIO
from typing import Set
from hierarchy_utils import HierarchyUtilities
from data_classes import InputConfiguration

class ValidateFiles:

    def __init__(self) -> None:
        self._validated_files: Set[str]=set()
    
    @property
    def validated_files(self) -> Set[str]:
        return self._validated_files


    def validate_many(self, file_list, files_type) -> True:
        if len(file_list)==0:
            raise ValueError(f'Cannot validate an empty list')
        with tqdm(total=len(file_list)) as progress_meter:
            for file in file_list:
                if files_type=="fasta":
                    self.validate_fasta(file)
                elif files_type=="vcf":
                    self.validate_vcf(file)
                elif files_type=="bed":
                    self.validate_bed(file)
                else:
                    raise ValueError(f'Unknown files type: {files_type}')
                progress_meter.update(1)
        return True

    def validate_bed(self, bed_file_name: str, **kwargs) -> True:
        """Class for validate input bedfile
        :key min_col_number: minimum required number of column in bed file, default: 3, int
        """        
        if not exists(bed_file_name):
            warnings.warn(f'File or directory {bed_file_name} does not exist')
            return False
        min_column_count=kwargs.get("min_col_number",3)
        with open(bed_file_name) as bed_file:
            for line_counter, line in enumerate(bed_file):
                if line.find("\t")==-1:
                    warnings.warn(f'No tab-delimited found in file {bed_file_name} on line {line_counter}:\n {line}')
                    return False
                line_values=line.strip().split("\t")
                if len(line_values)<min_column_count:
                    raise ValueError(f'Min number of tab-delimited columns is {min_column_count}, but only {len(line_values)} were found on line {line_counter}:\n {line} ')
                if not line_values[1].isdigit() or not line_values[2].isdigit():
                    raise ValueError(f'Expecting numeric values in column 1 and 2, but none numeric values (posibly decimal) values were found in {bed_file_name} on line {line_counter}:\n {line}')
                if int(line_values[2])<=int(line_values[1]):
                    raise ValueError(f'Value in column 2 must be greater than value in column 1. Check {bed_file_name} on line {line_counter}:\n {line}')
            if 'line_counter' not in vars():
                warnings.warn(f'{bed_file_name} is empty')
                return False
        self._validated_files.add(bed_file_name)
        return True

    def validate_fasta(self, fasta_file_name: str) -> True:
        if not exists(fasta_file_name):
            warnings.warn(f'File or directory {fasta_file_name} does not exist')
            return False
        with open(fasta_file_name) as fasta_file:
            first_fifty_char=fasta_file.readline()[0:50]
            if len(first_fifty_char)==0:
                warnings.warn(f'Fasta file {fasta_file_name} is empty')
                return False
            if first_fifty_char[0]!=">":
                warnings.warn(f'Fasta file must have ">" on first line in {fasta_file_name}\n {first_fifty_char}')
                return False
        self._validated_files.add(fasta_file_name)
        return True
    
    def fasta_has_dashes(self, fasta_file_name: str) -> bool:
        with open(fasta_file_name) as fasta_file:
            for line in fasta_file:
                if line[0]!=">":
                    if line.find("-")>-1:
                        warnings.warn(f'Fasta file {fasta_file_name} has "-". This is likely to cause problems')
                        return True
        return False
    
    def contigs_in_fasta(self, bed_file_name: str, fasta_file_name: str) -> bool:
        ### Assume that bed and fasta files have already been validated
        bed_contigs=set()
        with open(bed_file_name) as bed_file:
            for line in (bed_file):
                bed_contigs.add(line.split("\t")[0])

        for record in SeqIO.parse(fasta_file_name,"fasta"):
            if record.id in bed_contigs:
                bed_contigs.remove(record.id)
            if len(bed_contigs)==0:
                break
        if len(bed_contigs)!=0:
            missing_contigs="\n".join( list(bed_contigs)[0:min(len(bed_contigs),10)] )
            warnings.warn(f'Some bed file {bed_file_name} contig IDs not found in fasta file \n Ex. {fasta_file_name} (first 10):\n {missing_contigs} ')
            return False
        return True
            
    def contigs_in_vcf(self, bed_file_name: str, vcf_file_name: str) -> bool:
        '''Checks that at least one bedfile contig is present in the VCF file
        Does NOT check if ALL bed file contigs are in the VCF file because some contigs
        may not have SNPs'''
        bed_contigs=set()
        with open(bed_file_name) as bed_file:
            for line in (bed_file):
                bed_contigs.add(line.split("\t")[0])
        if len(bed_contigs)==0:
            warnings.warn(f'Bed file {bed_file_name} is empty')
            return True
        
        vcf_contigs=set()
        with open(vcf_file_name) as vcf_file:
            for line in vcf_file:
                if line[0]!="#":
                    contig_id=line.split("\t")[0]
                    vcf_contigs.add(contig_id)
                    if contig_id in bed_contigs:
                        return True
        if len(vcf_contigs)==0:
            warnings.warn(f'VCF file {vcf_file_name} had no variants (i.e. no lines that do not start with #)')
            return True
        else:
            warnings.warn("\n"+f'None of the contigs in VCF file {vcf_file_name} are present in bedfile {bed_file_name}')
            return False

    def validate_negative_genomes(self, genomes_dir: str) ->bool:
        if not exists(genomes_dir):
            warnings.warn(f'Negative genomes directory {genomes_dir} does not exist. Please check spelling.')
            return False
        negative_genomes=[file for file in listdir(genomes_dir) if file.split(".")[-1]=="fna" or file.split(".")[-1]=="fasta"]
        if len(negative_genomes)==0:
            warnings.warn(f'Negative genomes directory {genomes_dir} has no files ending in .fna or .fasta')
            return False
        return True

    def validate_vcf(self, vcfs_dir: str) -> bool:
        if not exists(vcfs_dir):
            warnings.warn(f'VCF directory {vcfs_dir} does not exist. Please check spelling.')
            return False
        vcf_files=[file for file in listdir(vcfs_dir) if file.split(".")[-1]=="vcf"]
        if len(vcf_files)==0:
            warnings.warn(f'Directory {vcfs_dir} have no VCF files')
            return False
        vcf_file_name=vcfs_dir+vcf_files[0]
        if not exists(vcf_file_name):
            warnings.warn(f'VCF file or directory {vcf_file_name} does not exist')
            return False
        vcf_file= open(vcf_file_name)
        first_two_char=vcf_file.readline()[0:2]
        if first_two_char!="##":
            vcf_file.close()
            warnings.warn(f'First line of VCF file {vcf_file_name} does not start with ##')
            return False
        for line in vcf_file:
            if len(line)>=6:
                if line[0:2]=="##":
                    pass #waiting to find line with #CHROM in it
                if line[0:6]=="#CHROM" and len(line.split("\t"))<=9:
                    vcf_file.close()
                    warnings.warn(f'The header VCF line (#CHROM...) should have 9 or more lines, but has fewer in file {vcf_file_name}')
                    return False
        self._validated_files.add(vcf_file_name)                    
        return True

    def validate_hierarchy(self, hierarchy_file_name: str, config_data: InputConfiguration) -> bool:
        if not exists(hierarchy_file_name):
            warnings.warn(f'Hierarchy file {hierarchy_file_name} does not exist')
            return False
        with open(hierarchy_file_name) as input_file:
            for line in input_file:
                if line.strip()=="":
                    warnings.warn(f'Genotype hierarchy file {hierarchy_file_name} has no lines or empty lines')
                    return False
        self._validated_files.add(hierarchy_file_name)

        hierarchy_data=HierarchyUtilities().load_hierarchy(config_data.hierarchy_file)
        self.rare_genotypes_in_hierarchy(config_data,hierarchy_data=hierarchy_data)
        
        return True

    def rare_genotypes_in_hierarchy(self, config_data: InputConfiguration, hierarchy_data: HierarchyUtilities) -> bool:
        rare_gts_no_in_hierarchy=[]
        for genotype in config_data.gts_with_few_snps:
            if genotype not in hierarchy_data:
                rare_gts_no_in_hierarchy.append(genotype)
        if len(rare_gts_no_in_hierarchy)>0:
            mismatched_gts=",".join(rare_gts_no_in_hierarchy)
            raise ValueError(f'Config file has genotype {mismatched_gts} as rare genotype, but its absent in first column of hierarchy file')

    def validate_specific_target_snps(self, specific_snps_file_name:str) -> bool:
        """
        Check that the file with list of specific SNPs has three columns per line,
        that separator is Tab, that value in last column is unique
        and that second column contains integers

        :param specific_snps_file_name: File name for file with specific SNPs to target
        :type specific_snps_file_name: str
        """
        lines_counter=1
        names=set()
        if not exists(specific_snps_file_name):
            warnings.warn(f'File with specific target SNPs is specified, but does not exist. Check the address: {specific_snps_file_name} does not exist.')
            return False
        with open(specific_snps_file_name) as snps_file:
            for line in snps_file:
                values=line.strip().split("\t")
                if len(values)!=3:
                    if lines_counter==0:
                        warnings.warn(f'SNPs file {specific_snps_file_name} is either empty or does not use Tab as delimiter')
                        warnings.warn(f'This is not a problem, unless you meant to specify something in it.')
                    if line.strip()!="":
                        warnings.warn(f'SNPs file {specific_snps_file_name} line {str(lines_counter)} has information, but contains more or fewer than 3 values.')
                        warnings.warn(f'Please check that it uses Tab as delimiter and only has three values per line.')
                        return False
                if not values[1].isdigit():
                    warnings.warn(f'Line {str(lines_counter)} in file {specific_snps_file_name} has non-integer value. Position in contig can only be integer.')
                    return False
                if values[2] in names:
                    warnings.warn(f'Some names (third column) in file {specific_snps_file_name} are not unique. Only unique values are permitted.')
                    return False
                names.add(values[2])
                lines_counter+=1
        return True