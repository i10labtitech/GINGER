SRC=src
DEST=util
GXX=g++

all: mapping denovo homology abinitio merge_phase0 merge_phase1 merge_phase2 summary evaluation

clean: clean_mapping clean_denovo clean_homology clean_abinitio \
clean_merge_phase0 clean_merge_phase1 clean_merge_phase2 clean_summary clean_evaluation


mapping: gff_trimmer exon_num_filter longest_transcript \
repeat_checker strand_replace set_difference \
tag_trimmer ORF_finder

gff_trimmer: ${SRC}/gff_trimmer.cpp
	${GXX} ${SRC}/gff_trimmer.cpp -o ${DEST}/gff_trimmer -std=c++0x -O3

exon_num_filter: ${SRC}/exon_num_filter.cpp
	${GXX} ${SRC}/exon_num_filter.cpp -o ${DEST}/exon_num_filter -std=c++0x -O3

longest_transcript: ${SRC}/longest_transcript.cpp
	${GXX} ${SRC}/longest_transcript.cpp -o ${DEST}/longest_transcript -std=c++0x -O3

repeat_checker: ${SRC}/repeat_checker.cpp
	${GXX} ${SRC}/repeat_checker.cpp -o ${DEST}/repeat_checker -std=c++0x -O3

strand_replace: ${SRC}/strand_replace.cpp
	${GXX} ${SRC}/strand_replace.cpp -o ${DEST}/strand_replace -std=c++0x -O3

set_difference: ${SRC}/set_difference.cpp
	${GXX} ${SRC}/set_difference.cpp -o ${DEST}/set_difference -std=c++0x -O3

tag_trimmer: ${SRC}/tag_trimmer.cpp
	${GXX} ${SRC}/tag_trimmer.cpp -o ${DEST}/tag_trimmer -std=c++0x -O3

ORF_finder: ${SRC}/ORF_finder.cpp
	${GXX} ${SRC}/ORF_finder.cpp -o ${DEST}/ORF_finder -std=c++0x -O3

clean_mapping:
	rm -f \
$(DEST)/gff_trimmer $(DEST)/exon_num_filter \
$(DEST)/longest_transcript $(DEST)/repeat_checker \
$(DEST)/strand_replace $(DEST)/set_difference $(DEST)/tag_trimmer \
$(DEST)/ORF_finder $(DEST)/ORF_finder 


homology: fastarepair fastarepair2 gff_2_proteinfasta flameshiftfilter

fastarepair: ${SRC}/fastarepair.cpp ${SRC}/function.hpp
	${GXX} ${SRC}/fastarepair.cpp -o ${DEST}/fastarepair -std=c++0x -O3

fastarepair2: ${SRC}/fastarepair2.cpp ${SRC}/function.hpp
	${GXX} ${SRC}/fastarepair2.cpp -o ${DEST}/fastarepair2 -std=c++0x -O3

gff_2_proteinfasta: ${SRC}/gff_2_proteinfasta.cpp
	${GXX} ${SRC}/gff_2_proteinfasta.cpp -o ${DEST}/gff_2_proteinfasta -std=c++0x -O3

flameshiftfilter: ${SRC}/flameshiftgrep.cpp
	${GXX} ${SRC}/flameshiftgrep.cpp -o ${DEST}/flameshiftfilter -std=c++0x -O3

clean_homology:
	rm -f $(DEST)/fastarepair $(DEST)/fastarepair2 $(DEST)/gff_2_proteinfasta $(DEST)/flameshiftfilter 


denovo:

clean_denovo:


abinitio: simple_low_norepeatmask inframe_stopcodon_exclude makefasta

simple_low_norepeatmask: ${SRC}/simple_low_norepeatmask.cpp
	${GXX} ${SRC}/simple_low_norepeatmask.cpp -o ${DEST}/simple_low_norepeatmask -std=c++0x -O3

inframe_stopcodon_exclude: ${SRC}/inframe_stopcodon_exclude.cpp
	${GXX} ${SRC}/inframe_stopcodon_exclude.cpp -o ${DEST}/inframe_stopcodon_exclude -std=c++0x -O3

makefasta: ${SRC}/makefasta.cpp
	${GXX} ${SRC}/makefasta.cpp -o ${DEST}/makefasta -std=c++0x -O3

clean_abinitio:
	rm -f $(DEST)/simple_low_norepeatmask $(DEST)/inframe_stopcodon_exclude $(DEST)/makefasta


merge_phase0: phase0_gff_editor Row2_rename RNA-seq_reform Spaln_reform \
Augustus_reform

phase0_gff_editor: ${SRC}/190521_gff_editor.cpp ${SRC}/phase0_function.hpp
	${GXX} ${SRC}/190521_gff_editor.cpp -o ${DEST}/phase0_gff_editor -std=c++0x -O3

Row2_rename: ${SRC}/row2_rename.cpp ${SRC}/phase0_function.hpp
	${GXX} ${SRC}/row2_rename.cpp -o ${DEST}/Row2_rename -std=c++0x -O3

RNA-seq_reform: ${SRC}/rnaseq_reform.cpp ${SRC}/phase0_function.hpp
	${GXX} ${SRC}/rnaseq_reform.cpp -o ${DEST}/RNA-seq_reform -std=c++0x -O3

Spaln_reform: ${SRC}/spaln_reform.cpp ${SRC}/phase0_function.hpp
	${GXX} ${SRC}/spaln_reform.cpp -o ${DEAT}/Spaln_reform -std=c++0x -O3

Augustus_reform: ${SRC}/augustus_reform.cpp ${SRC}/phase0_function.hpp
	${GXX} ${SRC}/augustus_reform.cpp -o ${DEST}/Augustus_reform -std=c++0x -O3

clean_merge_phase0:
	rm -f \
$(DEST)/phase0_gff_editor $(DEST)/Row2_rename $(DEST)/RNA-seq_reform $(DEST)/Spaln_reform \
$(DEST)/Augustus_reform


merge_phase1: Grouping subgroup new_subgroup Searchalgo \
phase1_gff_editor initial_exon_polish

Grouping: ${SRC}/grouping.cpp
	${GXX} ${SRC}/grouping.cpp -o ${DEST}/Grouping -std=c++0x -O3

subgroup: ${SRC}/subgroup_v2.2.cpp
	${GXX} ${SRC}/subgroup_v2.2.cpp -o ${DEST}/subgroup -std=c++0x -O3

new_subgroup: ${SRC}/new_subgroup.cpp
	${GXX} ${SRC}/new_subgroup.cpp -o ${DEST}/new_subgroup -std=c++0x -O3

Searchalgo: ${SRC}/searchalgo.cpp
	${GXX} ${SRC}/searchalgo.cpp -o ${DEST}/Searchalgo -std=c++0x -O3

phase1_gff_editor: ${SRC}/phase1_gff_editor.cpp ${SRC}/phase1_function.hpp
	${GXX} ${SRC}/phase1_gff_editor.cpp -o ${DEST}/phase1_gff_editor -std=c++0x -O3

initial_exon_polish: ${SRC}/initial_exon_polish.cpp
	${GXX} ${SRC}/initial_exon_polish.cpp -o ${SRC}/initial_exon_polish -std=c++0x -O3

clean_merge_phase1:
	rm -f \
$(DEST)/Grouping $(DEST)/subgroup $(DEST)/new_subgroup $(DEST)/Searchalgo \
$(DEST)/phase1_gff_editor $(DEST)/initial_exon_polish


merge_phase2: geneadd_v191115 geneadd_v191119 grouping_v1

geneadd_v191115: ${SRC}/geneadd_v191115.cpp
	${GXX} ${SRC}/geneadd_v191115.cpp -o ${DEST}/geneadd_v191115 -std=c++0x -O3

geneadd_v191119: ${SRC}/geneadd_v191119.cpp
	${GXX} ${SRC}/geneadd_v191119.cpp -o ${DEST}/geneadd_v191119 -std=c++0x -O3

grouping_v1: ${SRC}/grouping_v1.cpp
	${GXX} ${SRC}/grouping_v1.cpp -o ${DEST}/grouping_v1 -std=c++0x -O3

clean_merge_phase2:
	rm -f \
$(DEST)/geneadd_v191115 $(DEST)/geneadd_v191119 $(DEST)/grouping_v1


summary: final_reform

final_reform: ${SRC}/final_reform.cpp
	${GXX} ${SRC}/final_reform.cpp -o ${DEST}/final_reform -std=c++0x -O3

clean_summary:
	rm -f $(DEST)/final_reform


evaluation: evaluate preevaluate

evaluate: ${SRC}/evaluation4.cpp
	${GXX} ${SRC}/evaluation4.cpp -o ${DEST}/evaluate -std=c++0x -O3

preevaluate: ${SRC}/preevaluation.cpp
	${GXX} ${SRC}/preevaluation.cpp -o ${DEST}/preevaluate -std=c++0x -O3

clean_evaluation:
	rm -f $(DEST)/evaluate $(DEST)/preevaluate
