#!/usr/bin/env nextflow

/*
Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology

This file is part of GINGER.

GINGER is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GINGER is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with GINGER; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

if (file(params.PDIR_PREP_MAPPING).isDirectory()) {
    error "The publish directory already exists: \"${params.PDIR_PREP_MAPPING}\""
} else {
    print "No publish directory : \"${params.PDIR_PREP_MAPPING}\""
    if (file(params.PDIR_PREP).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP}\""
        file(params.PDIR_PREP).mkdir()
    }
    if (file(params.PDIR_PREP_MAPPING).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP_MAPPING}\""
        file(params.PDIR_PREP_MAPPING).mkdir()
    }
}

if (file(params.INPUT_GENOME).isFile())    {} else { error "No path : \"${params.INPUT_GENOME}\"" }
if (file(params.INPUT_REPOUT).isFile())    {} else { error "No path : \"${params.INPUT_REPOUT}\"" }
if (file(params.INPUT_RNASEQR1).isFile())  {} else { error "No path : \"${params.INPUT_RNASEQR1}\"" }
if (file(params.INPUT_RNASEQR2).isFile())  {} else { error "No path : \"${params.INPUT_RNASEQR2}\"" }
if (file(params.INPUT_MAPPED_READ_SAM).isFile())  {} else { error "No path : \"${params.INPUT_MAPPED_READ_SAM}\"" }

SCRATCH = params.SCRATCH

// --- Mapping of RNA-seq reads ---
process mapping {
    
    cpus params.N_THREAD
    scratch "${SCRATCH}"

    publishDir "${params.PDIR_PREP_MAPPING_MAPPING}", mode: 'copy', overwrite: false
    
    input:
    path genomeFasta from params.INPUT_GENOME
    path rnaSeqRead1 from params.INPUT_RNASEQR1
    path rnaSeqRead2 from params.INPUT_RNASEQR2
    path mappedReadInSam from params.INPUT_MAPPED_READ_SAM
    
    output:
    file "${params.OPREFIX}.gtf" into gtf
    
    shell:
    '''
    mkdir !{params.OPREFIX}_tmp_index
    mkdir tmp_mapping_!{params.OPREFIX}
    
    # Mapping
    !{params.HISAT2BUILD} -p !{params.N_THREAD} !{genomeFasta} !{params.OPREFIX}_tmp_index/!{params.OPREFIX}_index
    !{params.HISAT2} -p !{params.N_THREAD} -x !{params.OPREFIX}_tmp_index/!{params.OPREFIX}_index -1 !{rnaSeqRead1} -2 !{rnaSeqRead2} -S tmp_mapping_!{params.OPREFIX}/!{params.OPREFIX}.sam --dta --no-discordant --no-mixed
    !{params.SAMTOOLS} view -bS -@ !{params.N_THREAD} mappedReadInSam | !{params.SAMTOOLS} sort -@ !{params.N_THREAD} -o tmp_mapping_!{params.OPREFIX}/!{params.OPREFIX}.sorted.bam
    
    # Gene prediction
    !{params.STRINGTIE} -p !{params.N_THREAD} -o !{params.OPREFIX}.gtf tmp_mapping_!{params.OPREFIX}/!{params.OPREFIX}.sorted.bam
    '''
}

// --- Merge ---
process to_merge {

    publishDir "${params.PDIR_PREP_MAPPING_MERGE}", mode: 'copy', overwrite: false

    input:
    path genomeFasta from params.INPUT_GENOME
    path gtf
        
    output:
    file "${params.OPREFIX}_merge.gff3" into gff3toMerge
    
    shell:
    '''
    !{params.SEQKIT} seq -u !{genomeFasta} > !{params.OPREFIX}_tmpgenome.fa
    
    awk '$7 == "."' !{gtf} > !{params.OPREFIX}_single.gtf
    
    !{params.UTILPATH_MAPPING}/gtf_genome_to_cdna_fasta.pl !{params.OPREFIX}_single.gtf !{params.OPREFIX}_tmpgenome.fa > !{params.OPREFIX}_single.fa
    
    !{params.UTILPATH_MAPPING}/ORF_finder !{params.OPREFIX}_single.fa !{params.OPREFIX}_single_orf !{params.MIN0} false
    
    !{params.SEQKIT} fx2tab -nl !{params.OPREFIX}_single_orf.cds | awk '{print $1"\t"$2"\t"$5}' > !{params.OPREFIX}_single_orf.list
    
    !{params.UTILPATH_MAPPING}/longest_transcript !{params.OPREFIX}_single_orf.list > !{params.OPREFIX}_single_orf_longest.list
    
    !{params.UTILPATH_MAPPING}/gff_trimmer !{params.OPREFIX}_single_orf.gff3 !{params.OPREFIX}_single_orf_longest.list !{params.OPREFIX}_single_orf_longest.gff3
    
    !{params.UTILPATH_MAPPING}/strand_replace !{params.OPREFIX}_single_orf_longest.gff3 !{gtf} !{params.OPREFIX}_replace.gtf
    
    awk '$7 != "."' !{params.OPREFIX}_replace.gtf > !{params.OPREFIX}_new.gtf
    
    !{params.UTILPATH_MAPPING}/gtf_genome_to_cdna_fasta.pl !{params.OPREFIX}_new.gtf !{params.OPREFIX}_tmpgenome.fa > !{params.OPREFIX}.fa
    
    !{params.UTILPATH_MAPPING}/gtf_to_alignment_gff3.pl !{params.OPREFIX}_new.gtf > !{params.OPREFIX}.gff3
    
    !{params.TD_LONGORFS} -t !{params.OPREFIX}.fa -m !{params.MIN0} -S
    !{params.TD_PREDICT} -t !{params.OPREFIX}.fa 
    
    !{params.UTILPATH_MAPPING}/cdna_alignment_orf_to_genome_orf.pl !{params.OPREFIX}.fa.transdecoder.gff3 !{params.OPREFIX}.gff3 !{params.OPREFIX}.fa > !{params.OPREFIX}_merge.gff3
    '''
}

// --- Learn(1st) ---
process to_learn_1st {

    publishDir "${params.PDIR_PREP_MAPPING_TOLEARN1ST}", mode: 'copy', overwrite: false
    
    input:
    path genomeFasta from params.INPUT_GENOME
    path gtf


    output:
    file "${params.OPREFIX}_learn_1st.gff3" into gff3_1st

    shell:
    '''
    !{params.SEQKIT} seq -u !{genomeFasta} > !{params.OPREFIX}_tmpgenome.fa
    
    !{params.UTILPATH_MAPPING}/exon_num_filter !{gtf} !{params.OPREFIX}_multi.gtf 2 10000
    
    !{params.UTILPATH_MAPPING}/gtf_genome_to_cdna_fasta.pl !{params.OPREFIX}_multi.gtf !{params.OPREFIX}_tmpgenome.fa > !{params.OPREFIX}_multi.fa
    
    !{params.UTILPATH_MAPPING}/gtf_to_alignment_gff3.pl !{params.OPREFIX}_multi.gtf > !{params.OPREFIX}_multi.aln.gff3
    
    !{params.UTILPATH_MAPPING}/ORF_finder !{params.OPREFIX}_multi.fa !{params.OPREFIX}_multi_orf !{params.MIN1} true
    
    !{params.SEQKIT} fx2tab -nl !{params.OPREFIX}_multi_orf.cds | awk '{print $1"\t"$2"\t"$5}' > !{params.OPREFIX}_multi_orf.list
    
    !{params.UTILPATH_MAPPING}/longest_transcript !{params.OPREFIX}_multi_orf.list > !{params.OPREFIX}_multi_orf_longest.list
    
    !{params.UTILPATH_MAPPING}/gff_trimmer !{params.OPREFIX}_multi_orf.gff3 !{params.OPREFIX}_multi_orf_longest.list !{params.OPREFIX}_multi_orf_longest.gff3
    
    !{params.GFFREAD} -x !{params.OPREFIX}_multi_orf_longest.cds -g !{params.OPREFIX}_multi.fa !{params.OPREFIX}_multi_orf_longest.gff3
    
    !{params.CD_HIT} -c 1.0 -i !{params.OPREFIX}_multi_orf_longest.cds -o !{params.OPREFIX}_multi_orf_longest_nr.cds
    
    !{params.SEQKIT} seq -n !{params.OPREFIX}_multi_orf_longest_nr.cds | awk '{print $1}' > !{params.OPREFIX}_multi_orf_longest_nr.list
    
    !{params.UTILPATH_MAPPING}/gff_trimmer !{params.OPREFIX}_multi_orf_longest.gff3 !{params.OPREFIX}_multi_orf_longest_nr.list !{params.OPREFIX}_multi_orf_longest_nr.gff3
    
    !{params.UTILPATH_MAPPING}/cdna_alignment_orf_to_genome_orf.pl !{params.OPREFIX}_multi_orf_longest_nr.gff3 !{params.OPREFIX}_multi.aln.gff3 !{params.OPREFIX}_multi.fa > !{params.OPREFIX}_learn_1st.gff3
    '''
}

// --- Learn(2nd) ---
process to_learn_2nd {

    publishDir "${params.PDIR_PREP_MAPPING_TOLEARN2ND}", mode: 'copy', overwrite: false

    input:
    path genomeFasta  from params.INPUT_GENOME
    path genomeRepOut from params.INPUT_REPOUT
    path gff3_1st

    output:
    file "${params.OPREFIX}_learn_2nd.gff3" into gff3_2nd

    shell:
    '''
    awk '$11 != "Simple_repeat" && $11 != "Low_complexity" {print $5"\t"$6"\t"$7"\t"$11}' !{genomeRepOut} | sed -e '1,3d' > repeat_DB.tsv
    
    !{params.SEQKIT} fx2tab -inl !{genomeFasta} | awk '{print $1"\t"$2}' | sed -e '/^$/d' > genome_DB.tsv
    
    awk '$3 == "exon" {print $1"\t"$4"\t"$5"\t"$9}' !{gff3_1st} | sed -e '/^$/d' > annotation_DB.tsv
    
    !{params.UTILPATH_MAPPING}/repeat_checker genome_DB.tsv repeat_DB.tsv annotation_DB.tsv repeat_output.tsv
    
    awk '$2 != 0 {print $1}' repeat_output.tsv > repeat_output.list
    
    !{params.UTILPATH_MAPPING}/tag_trimmer repeat_output.list Parent | sort | uniq > repeat_output.id
    
    awk '$3 == "mRNA" {print $9}' !{gff3_1st} > !{params.OPREFIX}_all.tag
    !{params.UTILPATH_MAPPING}/tag_trimmer !{params.OPREFIX}_all.tag ID > !{params.OPREFIX}_all.id
    
    !{params.UTILPATH_MAPPING}/set_difference !{params.OPREFIX}_all.id repeat_output.id norep_output.id
    
    !{params.UTILPATH_MAPPING}/gff_trimmer !{gff3_1st} norep_output.id !{params.OPREFIX}_learn_2nd.gff3
    '''
}
