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

if (file(params.PDIR_PREP_DENOVO).isDirectory()) {
    error "The publish directory already exists: \"${params.PDIR_PREP_DENOVO}\""
} else {
    print "No publish directory : \"${params.PDIR_PREP_DENOVO}\""
    if (file(params.PDIR_PREP).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP}\""
        file(params.PDIR_PREP).mkdir()
    }
    if (file(params.PDIR_PREP_DENOVO).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP_DENOVO}\""
        file(params.PDIR_PREP_DENOVO).mkdir()
    }
}
if (file(params.INPUT_GENOME).isFile())    {} else { error "No path : \"${params.INPUT_GENOME}\"" }
if (file(params.INPUT_RNASEQR1).isFile())  {} else { error "No path : \"${params.INPUT_RNASEQR1}\"" }
if (file(params.INPUT_RNASEQR2).isFile())  {} else { error "No path : \"${params.INPUT_RNASEQR2}\"" }

SCRATCH = params.SCRATCH

// --- Trinity ---
process trinity {

    cpus params.N_THREAD
    scratch "${SCRATCH}"

    publishDir "${params.PDIR_PREP_DENOVO_TRINITY}", mode: 'copy', overwrite: false
    
    input:
    path rnaSeqRead1 from params.INPUT_RNASEQR1
    path rnaSeqRead2 from params.INPUT_RNASEQR2
    
    output:
    file "${params.OPREFIX_TRINITY}.Trinity.fasta" into trinityFasta
    file "${params.OPREFIX_TRINITY}.Trinity.fasta.gene_trans_map"

    script:
    preprocessing = ""
    if (params.SRA_FLAG == 1) {
        preprocessing = """
        awk '{{print (NR%4 == 1) ? "@1_" ++i "/1": \$0}}' ${rnaSeqRead1} > rnaseq_1_renamed.fastq
        awk '{{print (NR%4 == 1) ? "@1_" ++i "/2": \$0}}' ${rnaSeqRead2} > rnaseq_2_renamed.fastq
        """
    } else {
        preprocessing = """
        ln -s ${rnaSeqRead1} rnaseq_1_renamed.fastq
        ln -s ${rnaSeqRead2} rnaseq_2_renamed.fastq
        """
    }

    shell:
    '''
    !{preprocessing}

    export PATH=`!{params.UTILPATH_DENOVO}/getPath.pl !{params.SAMTOOLS}`:\$PATH

    # Assembling with Trinity
    mkdir !{params.OPREFIX_TRINITY}
    time !{params.TRINITY} --seqType fq --left rnaseq_1_renamed.fastq --right rnaseq_2_renamed.fastq --output !{params.OPREFIX_TRINITY} --CPU !{params.N_THREAD} --max_memory !{params.MAX_MEMORY} --full_cleanup
    '''
}

// --- Oases ---
process oases {

    publishDir "${params.PDIR_PREP_DENOVO_OASES}", mode: 'copy', overwrite: false

    input:
    path rnaSeqRead1 from params.INPUT_RNASEQR1
    path rnaSeqRead2 from params.INPUT_RNASEQR2
        
    output:
    file "${params.OPREFIX}/transcripts.fa" into oasesFasta
    
    shell:
    '''
    !{params.VELVETH} !{params.OPREFIX} 31 -fastq -short -separate !{rnaSeqRead1} !{rnaSeqRead2}
    !{params.VELVETG} !{params.OPREFIX} -read_trkg yes 
    !{params.OASES} !{params.OPREFIX}
    '''
}

// --- De novo ---
process denovo {

    cpus params.N_THREAD

    publishDir "${params.PDIR_PREP_DENOVO}", mode: 'copy', overwrite: false
    
    input:
    path genomeFasta from params.INPUT_GENOME
    path trinityFasta
    path oasesFasta

    output:
    file "${params.OPREFIX}_denovo.cd.fasta.gmap.filtered"
    file "${params.OPREFIX}_final.gff3"

    shell:
    '''
    ### sequences ID simplification ###
    touch !{params.OPREFIX}_denovo.cds
    !{params.SEQKIT} seq -i -m 300 !{trinityFasta} >> !{params.OPREFIX}_denovo.cds
    !{params.SEQKIT} seq -i -m 300 !{oasesFasta} >> !{params.OPREFIX}_denovo.cds

    ### de novo sequences compression ###
    !{params.CD_HIT_EST} -i !{params.OPREFIX}_denovo.cds -c 1.0 -o !{params.OPREFIX}_denovo.cd.fasta -T !{params.N_THREAD} -M 100000

    ### GMAP alignment ###
    !{params.GMAP_BUILD} -D . -d !{params.OPREFIX}_genome_gmap_DB !{genomeFasta} > !{params.OPREFIX}_gmap_build.log 2>&1
    !{params.GMAP} -S -t !{params.N_THREAD} -n 1 -D . -d !{params.OPREFIX}_genome_gmap_DB !{params.OPREFIX}_denovo.cd.fasta > !{params.OPREFIX}_denovo.cd.fasta.gmap 2> !{params.OPREFIX}_gmap.stderr
    
    ### filtering alignment result ###
    !{params.DENOVO_PYTHON} !{params.UTILPATH_DENOVO}/filtering.py !{params.OPREFIX}_denovo.cd.fasta.gmap 95 > !{params.OPREFIX}_denovo.cd.fasta.gmap.filtered

    ### format conversion for ORF prediction
    !{params.UTILPATH_DENOVO}/gmap_native_to_format_converter.pl !{params.OPREFIX}_denovo.cd.fasta.gmap.filtered GTF > !{params.OPREFIX}_denovo.cd.fasta.gtf
    !{params.UTILPATH_DENOVO}/cufflinks_gtf_to_alignment_gff3.pl !{params.OPREFIX}_denovo.cd.fasta.gtf > !{params.OPREFIX}_denovo.cd.fasta.gff3
    !{params.UTILPATH_DENOVO}/cufflinks_gtf_genome_to_cdna_fasta.pl !{params.OPREFIX}_denovo.cd.fasta.gtf !{genomeFasta} > !{params.OPREFIX}_denovo_genome.fa

    ### ORF prediction ###
    !{params.TD_LONGORFS} -t !{params.OPREFIX}_denovo_genome.fa -m 90
    !{params.TD_PREDICT} --no_refine_starts --single_best_only -t !{params.OPREFIX}_denovo_genome.fa
    
    ### format conversion ###
    !{params.UTILPATH_DENOVO}/cdna_alignment_orf_to_genome_orf_nolimit.pl !{params.OPREFIX}_denovo_genome.fa.transdecoder.gff3 !{params.OPREFIX}_denovo.cd.fasta.gff3 !{params.OPREFIX}_denovo_genome.fa > !{params.OPREFIX}_final.gff3
    '''
}
