import primer3
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
import warnings
from os.path import exists
from typing import Dict, List, Tuple
root_dir="/home/lshas17/"
from itertools import product
from sys import path
path.insert(1, f'{root_dir}/HandyAmpliconTool/scripts/')
from data_classes import BlastResult, PrimerPair, Primer, InputConfiguration, Genotype, Genotypes, SNP
from run_blast import BlastRunner
from inputs_validation import ValidateFiles
from load_vcfs import VCFutilities
import pickle


class PrimersGenerator():
    def __init__(self, config: InputConfiguration) -> None:
        self.ref_seq: Dict[str, str]={}
        self.config=config
        if not exists(config.reference_fasta):
            raise IOError(f'Reference FASTA {config.reference_fasta} does not exist')
        for record in SeqIO.parse(config.reference_fasta, "fasta"):
            self.ref_seq[record.id]=str(record.seq)
        self.existing_primers: List[Primer]=self._load_existing_primers()
        self.new_primer_pairs: List[PrimerPair]=[]

    def _load_existing_primers(self) -> List[Primer]:
        result: List[Primer] = []
        if self.config.existing_primers!="":
            with open(self.config.existing_primers, "r") as file:
                for line in file:
                    values = line.strip().split("\t")
                    if len(values)==2:
                        result.append(self._sequence_to_primer(values[0], False) )
                        result.append(self._sequence_to_primer(values[1], True) )
        return result

    def _get_seq_coordinates_in_ref(self, sequence) -> int:
        seq=Seq(sequence)
        seq_rc=str(seq.reverse_complement())
        for ref_contig in self.ref_seq.values():
            direct_start=ref_contig.find(str(seq))
            if direct_start!=-1:
                return direct_start
            else:
                rc_start=ref_contig.find(seq_rc)
                if rc_start!=-1:
                    return rc_start
        return -1

    def _sequence_to_primer(self, primer_seq:str, is_reverse:bool) -> Primer:
        """Generates a Primer object from primer sequence.
        :param primer_seq: sequence of the primer
        :type primer_seq: str
        :param is_reverse: Indicates if primer sequence is for reverse primer
        :type is_reverse: bool        
        :return: primer object 
        :rtype: Primer
        """
        primer=Primer(primer_seq, gc_fraction(primer_seq), primer3.calc_tm(primer_seq), is_reverse )
        primer.ref_start=self._get_seq_coordinates_in_ref(primer_seq)
        return primer

    def _id_has_direction_info(self, id_value: str) -> Tuple[bool, str]:
        if str.lower(id_value).find("forward")>-1:
            return (True, "forward")
        elif str.lower(id_value).find("reverse")>-1:
            return (True, "reverse")
        else:
            return (False, "")
        
    # def load_existing_primer_from_fasta(self, existing_primers_fasta: str):
    #     self.existing_primers.clear()
    #     for record in SeqIO.parse(existing_primers_fasta, "fasta"):
    #         direction_info=self._id_has_direction_info(record.id)
    #         if direction_info[0]:
    #             primer=self._sequence_to_primer( str(record.seq), is_reverse=direction_info[1]=="reverse" )
    #         else: #this is worst case, add primer as both forward and reverse
    #             primer=self._sequence_to_primer( str(record.seq), is_reverse=True )
    #             self.existing_primers.append(primer)
    #             primer=self._sequence_to_primer( str(record.seq), is_reverse=False)
    #         self.existing_primers.append(primer)

    # def load_existing_primer_from_bed(self, existing_primers_bed: str):
    #     self.existing_primers.clear()        
    #     if existing_primers_bed!="":
    #         file_validator=ValidateFiles()
    #         file_validator.validate_bed(existing_primers_bed, min_col_number=6)
    #         with open(existing_primers_bed) as bed_file:
    #             for i, line in enumerate(bed_file):
    #                 values=line.strip().split("\t")
    #                 if values[0] not in self.ref_seq:
    #                     raise ValueError(f'Contig {values[0]} from existing amplicons file {existing_primers_bed} not found in reference fasta {self.config.reference_fasta}')
    #                 if int(values[1])+1>len(self.ref_seq[values[0]]):
    #                     raise ValueError(f'Amplicon {i+1} in existing amplicons file {existing_primers_bed} has coordinate {values[0]} {int(values[1])} greater than lenght of reference sequence: {self.config.reference_fasta}')
    #                 primer_seq=self.ref_seq[values[0]][int(values[1]):int(values[2])+1] # +1 because python excludes the last item of slice
    #                 is_reverse = values[5]=="-"
    #                 if is_reverse:
    #                     primer_seq=str(Seq(primer_seq).reverse_complement())
    #                 primer=self._sequence_to_primer( primer_seq, is_reverse )
    #                 self.existing_primers.append(primer)

    def add_global_primer_args(self, for_seq: str, rev_seq:str, template: str):
        global_args={
            'PRIMER_OPT_SIZE': self.config.primer_opt_size,
            'PRIMER_OPT_TM': self.config.primer_opt_tm,
            'PRIMER_MIN_TM': self.config.primer_min_tm,
            'PRIMER_MAX_TM': self.config.primer_max_tm,
            'PRIMER_PRODUCT_SIZE_RANGE': f'100-{len(template)}'}
        if for_seq!=None:
            global_args['SEQUENCE_PRIMER']=for_seq
        if rev_seq!=None:
            global_args['SEQUENCE_PRIMER_REVCOMP']=rev_seq
        return global_args
    
    def forms_homodimers(self, primer_sequence: str) -> bool:
        """Checks if sequence forms homodimer at above 10C
        :param primer_sequence: sequence of the primers
        :type primer_sequence: str
        :return: True if sequence forms homodimer
        :rtype: bool
        """
        thermoresult_forward=primer3.bindings.calc_homodimer(primer_sequence)
        return (thermoresult_forward.structure_found and thermoresult_forward.tm>10)
        
    def count_heterodimers(self, primer_sequence: str, orientation: str) -> int:
        """Checks if sequence forms heterodimers at above 10C
        with existing primers
        :param primer_sequence: sequence of the primers
        :type primer_sequence: str
        :param orientation: Forward, Reverse or Unknown
        :type orientation: str
        :return: number of heterodimers sequences forms with existing primers
        :rtype: int
        """
        total_structures=0
        for primer in self.existing_primers:
            if (orientation=="Forward" and not primer.is_reverse) or \
                (orientation=="Reverse" and primer.is_reverse) or \
                    (orientation=="Unknown"):
                temp_result=primer3.bindings.calc_heterodimer(primer_sequence, primer.seq)
                total_structures+=1 if temp_result.structure_found and temp_result.tm>10 else 0
        return total_structures
    
    def process_p3_output(self, output: Dict, target: SNP, id_suffix: str, template_start: int, template_contig:str) -> List[PrimerPair]:
        """Converts the primer3 output format into list of PrimerPair objects
        with existing primers
        :param output: dictionary object generated by primer3
        :type output: Dict
        :param id_prefix: suffix to assign to primer pair name
        :type id_prefix: str
        :return: number of heterodimers sequences forms with existing primers
        :rtype: int
        """
        primer_pairs=len(output["PRIMER_PAIR"])
        result: List[PrimerPair]=[]
        for i in range(0,primer_pairs):
            primer_name=id_suffix
            forward=Primer(seq=output[f'PRIMER_LEFT_{str(i)}_SEQUENCE'],
                        g_c=output[f'PRIMER_LEFT_{str(i)}_GC_PERCENT'],
                        t_m=output[f'PRIMER_LEFT_{str(i)}_TM'], is_reverse=False)
            forward.ref_start=template_start+int(output[f'PRIMER_LEFT_{str(i)}'][0])
            reverse=Primer(seq=output[f'PRIMER_RIGHT_{str(i)}_SEQUENCE'],
                        g_c=output[f'PRIMER_RIGHT_{str(i)}_GC_PERCENT'],
                        t_m=output[f'PRIMER_RIGHT_{str(i)}_TM'], is_reverse=True)
            reverse.ref_start=template_start+int(output[f'PRIMER_RIGHT_{str(i)}'][0])-len(reverse.seq)+1
            if forward.ref_end > target.position or reverse.ref_start < target.position:
                continue #pair doesn't include the target SNP
            new_pair=PrimerPair(name_suffix=primer_name, forward=forward, reverse=reverse)
            new_pair.ref_contig=template_contig
            new_pair.penalty=output[f'PRIMER_PAIR_{str(i)}_PENALTY']
            amplicon_length=new_pair.reverse.ref_end-new_pair.forward.ref_start
            if not self.forms_homodimers(forward.seq) and \
                    not self.forms_homodimers(reverse.seq) and \
                    self.count_heterodimers(forward.seq, reverse.seq)==0 and \
                    amplicon_length>self.config.min_amplicon_length and \
                    amplicon_length<self.config.max_amplicon_len:
                result.append(new_pair)
        return result

    # def _for_testing_load_gts(self, species_pkl, gts_pkl):
    #     with open(species_pkl, "rb") as pickled_file:
    #         species: Genotype = pickle.load(pickled_file)

    #     with open(gts_pkl, "rb") as pickled_file:
    #         self.genotypes: Genotypes = pickle.load(pickled_file)

    #     self.genotypes.genotypes.append(species)

    def _snps_within_interval(self, snps: List[SNP], ref_contig: int, interval_start:int, interval_end: int) -> List[SNP]:
        result=[snp for snp in snps if snp.ref_contig_id==ref_contig and snp.position > interval_start and snp.position < interval_end]
        return result
    
    def _get_ref_sequence(self, seq_id, start, end) -> str:
        if seq_id not in self.ref_seq:
            raise ValueError(f'Sequence {seq_id} not found in reference FASTA {self.config.reference_fasta}')
        return str(self.ref_seq[seq_id][start:end])
    
    def _both_primers_given(self, left_snps: List[SNP], right_snps: List[SNP], target_snp:SNP) -> List[PrimerPair]:
        results : List[PrimerPair] = []
        for left_snp, right_snp in product(left_snps, right_snps):
            #case 1
            distance_between_snps=right_snp.position-left_snp.position
            #just because both primers can be fixed, doesn't mean it's an optimal pairing.
            #Check left and right primers separately as well.
            results += self._left_primers_given(left_snps=[left_snp], target_snp=target_snp)
            results += self._right_primers_given(right_snps=[right_snp], target_snp=target_snp)
            #check the fixed primer pair
            if distance_between_snps > self.config.max_amplicon_len or distance_between_snps < InputConfiguration.min_amplicon_length: #distance between SNPs is too long:
                continue
            template_contig=left_snp.ref_contig_id
            forward=self._get_ref_sequence( template_contig, left_snp.position-1, left_snp.position+20  )
            reverse=self._get_ref_sequence( template_contig, right_snp.position-20,right_snp.position )
            reverse=str(Seq(reverse).reverse_complement())
            template_start=left_snp.position-1
            seq_args={"SEQUENCE_ID":f'{self.target_gt}_{str(target_snp.position)}',
                    "SEQUENCE_TEMPLATE": self._get_ref_sequence( template_contig , template_start, right_snp.position ),
                    "SEQUENCE_INCLUDED_REGION": [0, right_snp.position-left_snp.position+1 ]}
            #add "fixed" to name to show that this pair has both sides fixed
            global_args=self.add_global_primer_args(for_seq=forward,rev_seq=reverse, template=seq_args["SEQUENCE_TEMPLATE"])
            primers=primer3.bindings.design_primers(seq_args=seq_args, global_args=global_args)
            results += self.process_p3_output(primers, target_snp, '_both_fixed', template_start, template_contig)
        return results

    def _right_primers_given(self, right_snps:List[SNP], target_snp: SNP) ->  List[PrimerPair]:
        for right_snp in right_snps:
            template_contig=right_snp.ref_contig_id
            template_start=right_snp.position-self.config.max_amplicon_len
            reverse=self._get_ref_sequence( template_contig, right_snp.position-20, right_snp.position )
            reverse=str(Seq(reverse).reverse_complement())
            seq_args={"SEQUENCE_ID":f'{self.target_gt}_{str(target_snp.position)}',
                    "SEQUENCE_TEMPLATE": self._get_ref_sequence( template_contig, template_start , right_snp.position),
                    "SEQUENCE_INCLUDED_REGION": [0, self.config.max_amplicon_len]}
            global_args=self.add_global_primer_args(for_seq=None,rev_seq=reverse,  template=seq_args["SEQUENCE_TEMPLATE"])
            primers=primer3.bindings.design_primers(seq_args=seq_args, global_args=global_args)
            return self.process_p3_output(primers, target_snp, '_right_fixed', template_start, template_contig)

    def _left_primers_given(self, left_snps: List[SNP], target_snp: SNP) ->  List[PrimerPair]:
        for left_snp in left_snps:
            template_contig=left_snp.ref_contig_id
            template_start=left_snp.position-20
            forward=self._get_ref_sequence( template_contig, left_snp.position-1, left_snp.position+20  )
            seq_args={"SEQUENCE_ID":f'{self.target_gt}_{str(target_snp.position)}',
                    "SEQUENCE_TEMPLATE":    self._get_ref_sequence(template_contig, template_start, template_start+self.config.max_amplicon_len),
                    "SEQUENCE_INCLUDED_REGION": [0, self.config.max_amplicon_len ]}
            global_args=self.add_global_primer_args(for_seq=forward, rev_seq=None, template=seq_args["SEQUENCE_TEMPLATE"])
            primers=primer3.bindings.design_primers(seq_args=seq_args, global_args=global_args)
            return self.process_p3_output(primers, target_snp, '_left_fixed', template_start, template_contig)

    def _remove_duplicate_primer_pairs(self, all_pairs:List[PrimerPair]) -> None:
            """Sometimes the same pair of F/R sequences is generated multiple times 
            because there are multiple species SNPs around same genotype SNP
            :param all_pairs: list of primer pairs to de-duplicate
            :type all_pairs: List[PrimerPair]
            """
            existing_sequences=set()
            i=1
            while i<len(all_pairs):
                sequences=all_pairs[i].forward.seq+"_"+all_pairs[i].reverse.seq
                if sequences not in existing_sequences:
                    existing_sequences.add(sequences)
                    i+=1
                else:
                    del all_pairs[i]

    def _remove_primers_in_repeat_regions(self, all_pairs:List[PrimerPair]) -> None:
        """Removes primer pairs that overlap with repeat intevals from BED file in config
        :param all_pairs: list of primer pairs to check against repeat regions
        :type all_pairs: List[PrimerPair]
        """
        if self.config.repeats_bed_file=="":
            return None
        
        invalid_pairs:List[PrimerPair]=[]
        for pair in all_pairs:
            if ( pair.ref_contig, pair.primers[0].ref_start) in VCFutilities.repeat_coordinates or \
                ( pair.ref_contig, pair.primers[1].ref_end) in VCFutilities.repeat_coordinates:
                    invalid_pairs.append(pair)
        for pair in invalid_pairs:
            all_pairs.remove(pair)

    def _primers_list_to_string(self, primers: List[Primer]) -> str:
        result=""
        for i, primer in enumerate(primers):
            result=result+f'>Primer_{i+1}'+"\n"+primer.seq+"\n"
        return result
    
    def _primer_pairs_list_to_string(self, primer_pairs: List[PrimerPair]) -> str:
        result=""
        for pair in primer_pairs:
            result=result+f'>{pair.uuid}_Forward'+"\n"+pair.forward.seq+"\n"
            result=result+f'>{pair.uuid}_Reverse'+"\n"+pair.reverse.seq+"\n"
        return result

    def _new_primer_header_to_object(self, fasta_header: str) -> PrimerPair:
        if fasta_header not in [f.uuid for f in self.new_primer_pairs]:
            raise ValueError("Unable to find primer pair uuid among newly identified candidate primers.")
        return [f for f in self.new_primer_pairs if f.uuid==fasta_header][0]

    def _remove_interfering_primers(self, primer_pairs: List[PrimerPair]) -> None:
        """Uses BLAST to primers that map in between existing primers. Will not check primers against themselves.
        :param primer_pairs: list of new primer pairs to check against existing
        :type primer_pairs: List[PrimerPair]
        """
        existing_primers_string=self._primers_list_to_string(self.existing_primers)
        primers_to_check_string=self._primer_pairs_list_to_string( primer_pairs )
        blast_runner=BlastRunner(word_size=5, e_value=0.1)
        blast_runner.db_from_file(self.config.reference_fasta, self.config.temp_blast_db)
        existing_primers_hits:List[BlastResult] = blast_runner.run_from_multi_sequence_string(existing_primers_string, self.config.temp_blast_db)
        existing_primers_hits=sorted(existing_primers_hits, key = lambda x: (x.sseqid, x.sstart))
        new_primers_hits:List[BlastResult] = blast_runner.run_from_multi_sequence_string(primers_to_check_string, self.config.temp_blast_db)
        new_primers_hits=sorted(new_primers_hits, key = lambda x: (x.sseqid, x.sstart))
        
        #This is suboptimal, but the only way I see to allow standard file format as input for existing primers
        start_index=0
        interfering_new_primers_headers: List[str]=[]
        for i in range(0, len(existing_primers_hits)-1): #the last hit is a corner case for complete genoes (i.e. the other primer is at the start of sequence), but that's a lot of hassle for little gain
            if existing_primers_hits[i].sseqid!=existing_primers_hits[i+1].sseqid: #hits on different contigs
                continue
            if abs(existing_primers_hits[i].sstart-existing_primers_hits[i+1].sstart)>self.config.max_amplicon_len:#hits too far apart
                continue
            if existing_primers_hits[i].is_flipped or not existing_primers_hits[i+1].is_flipped: #Primers are not in forward/reverse orientation
                continue
            new_start_index=-1
            #Hits of existing primers are on the same contig and closer than maximum permitted interval length
            for j, new_primer_hit in enumerate(new_primers_hits[start_index:]):
                if new_primer_hit.sseqid==existing_primers_hits[i].sseqid: #hit on same contig
                    if new_primer_hit.sstart>=existing_primers_hits[i].sstart and new_primer_hit.sstart<=existing_primers_hits[i+1].send: #hits falls between to existing primer hits
                        if new_start_index==-1: #for efficiency, don't rerun the whole list, because both lists are sorted
                            new_start_index=start_index+j
                        #We don't know orientation of existing primers, so assume worst case - that primer aligns in intended diretion
                        interfering_new_primers_headers.append(new_primer_hit.qseqid.replace("_Forward","").replace("_Reverse",""))
            start_index=new_start_index
        interfering_new_primers_headers=list(set(interfering_new_primers_headers))
        for intefering_primer in interfering_new_primers_headers:
            self.new_primer_pairs.remove(self._new_primer_header_to_object(intefering_primer))

    def _add_extra_gts(self, primer_pairs: List[PrimerPair], all_gt_snps:List[SNP]) -> None:
        for pair in primer_pairs:
            additional_snps=self._snps_within_interval(all_gt_snps, pair.ref_contig, pair.forward.ref_start, pair.reverse.ref_end)
            if len(additional_snps)>0:
                for snp in additional_snps:
                    for gt in [f.name for f in self.genotypes.genotypes_with_snp(snp)]:
                        if gt!=InputConfiguration.SPECIES_NAME:
                            pair.targets.add(gt)

    def find_candidate_primers(self, target_gts: List[str]) -> List[PrimerPair]:
        """
        Identifies a set of optimal primers using Primer3

        :param target_gts: List of genotypes for which to design primers
        :type is_reverse: List[str]

        :return: list of primers
        :rtype: List[PrimerPair]
        """

        self.new_primer_pairs.clear()
        # for every SNP in target_lineage, identify the nearby SNPs 
        all_species_snps=[snp for genotype in self.genotypes.genotypes for snp in genotype.defining_snps if genotype.name==InputConfiguration.SPECIES_NAME]
        all_species_snps=sorted(all_species_snps, key=lambda x: (x.ref_contig_id, x.position) )
        interval_len=self.config.flank_len_to_check

        with open(self.config.output_dir+"snps.tsv","w") as output_file:
            for genotype in target_gts:
                print(f'checking SNPs for genotype {genotype} ')
                self.target_gt=genotype
                target_genotype=self.genotypes.get_genotype(self.target_gt)
                for snp in target_genotype.defining_snps:
                    output_file.write(snp.ref_contig_id+"\t"+str(snp.position)+
                                    "\t"+str(snp.position+1)+
                                    "\t"+target_genotype.name+"_"+
                                    '{0:.2f}'.format(snp.specificity)+
                                    "-"+str(snp.sensitivity)+"\n")

                for i, snp in enumerate(target_genotype.defining_snps):
                    species_gt_snps=self._snps_within_interval(all_species_snps, snp.ref_contig_id, snp.position-interval_len, snp.position+interval_len)
                    if len(species_gt_snps)==0:
                        continue
                    #if len(species_gt_snps)>0 and len(species_gt_snps)<30: # the target SNP has at least one flanking species SNPs, but too many is indicative of problematic region
                    left_species_snps=[species_snp for species_snp in species_gt_snps if species_snp.position<snp.position]
                    right_species_snps=[species_snp for species_snp in species_gt_snps if species_snp.position>snp.position]
                    #four cases:
                    # 1 - left and right serovar SNPs anchored primers
                    # 2 - left serovar anchored primer
                    # 3 - right serovar anchored primer
                    # 4 - neither side has serovar primer, ignore this kind of SNP
                    if len(left_species_snps)>0 and len(right_species_snps)>0:
                        snp_primer_pairs = self._both_primers_given(left_snps=left_species_snps, right_snps=right_species_snps, target_snp=snp)
                    elif len(left_species_snps)>0:
                        snp_primer_pairs = self._left_primers_given(left_snps=left_species_snps, target_snp=snp)
                    elif len(right_species_snps)>0:
                        snp_primer_pairs = self._right_primers_given(right_snps=right_species_snps, target_snp=snp)
                    else:
                        continue
                    print(f'Found {len(snp_primer_pairs)} primer pairs for SNP {snp.ref_contig_id} {snp.position}')
                    for pair in snp_primer_pairs:
                        pair.targets.add(genotype)
                        #Check how many SNP conincide with the generated primers
                        pair.forward.species_snps=len(self._snps_within_interval(all_species_snps, snp.ref_contig_id, pair.forward.ref_start, pair.forward.ref_end))
                        pair.reverse.species_snps=len(self._snps_within_interval(all_species_snps, snp.ref_contig_id, pair.reverse.ref_start, pair.reverse.ref_end))
                    print(f'Adding {len(snp_primer_pairs)} primer pairs for SNP {snp.ref_contig_id} {snp.position}')
                    self.new_primer_pairs+=snp_primer_pairs
                    print(len(self.new_primer_pairs))
            print(f'{len(self.new_primer_pairs)} prior to removing duplicates')
            self._remove_duplicate_primer_pairs(self.new_primer_pairs)
            print(f'{len(self.new_primer_pairs)} left after removing duplicates')
            self._remove_interfering_primers(self.new_primer_pairs)
            print(f'{len(self.new_primer_pairs)} left after removing interfering primers')
            self._remove_primers_in_repeat_regions(self.new_primer_pairs)
            print(f'{len(self.new_primer_pairs)} left after removing primers in repeat regions')

            #Check if other genotypes are captured by selected primers
            all_genotype_snps=[snp for genotype in self.genotypes.genotypes for snp in genotype.defining_snps if genotype.name!=InputConfiguration.SPECIES_NAME and genotype.name!=self.target_gt]
            all_genotype_snps=sorted(all_genotype_snps, key=lambda x: (x.ref_contig_id, x.position) )
            self._add_extra_gts(self.new_primer_pairs,all_genotype_snps)
        return self.new_primer_pairs
