all: mapping denovo homology abinitio merge_phase0 merge_phase1 merge_phase2 summary evaluation

clean: clean_mapping clean_denovo clean_homology clean_abinitio \
clean_merge_phase0 clean_merge_phase1 clean_merge_phase2 clean_summary clean_evaluation

mapping:
	cd src/mapping && $(MAKE) all

denovo:
	cd src/denovo && $(MAKE) all

homology:
	cd src/homology && $(MAKE) all

abinitio:
	cd src/abinitio && $(MAKE) all

merge_phase0:
	cd src/merge_phase0 && $(MAKE) all

merge_phase1:
	cd src/merge_phase1 && $(MAKE) all

merge_phase2:
	cd src/merge_phase2 && $(MAKE) all

summary:
	cd src/summary && $(MAKE) all

evaluation:
	cd src/evaluation && $(MAKE) all

clean_mapping:
	cd src/mapping && $(MAKE) clean

clean_denovo:
	cd src/denovo && $(MAKE) clean

clean_homology:
	cd src/homology && $(MAKE) clean

clean_abinitio:
	cd src/abinitio && $(MAKE) clean

clean_merge_phase0:
	cd src/merge_phase0 && $(MAKE) clean

clean_merge_phase1:
	cd src/merge_phase1 && $(MAKE) clean

clean_merge_phase2:
	cd src/merge_phase2 && $(MAKE) clean

clean_summary:
	cd src/summary && $(MAKE) clean

clean_evaluation:
	cd src/evaluation && $(MAKE) clean

###
DEST=../../util/mapping
DEST2=../../util/denovo
GXX=g++

all: gff_trimmer exon_num_filter longest_transcript \
repeat_checker strand_replace set_difference \
tag_trimmer ORF_finder install

gff_trimmer: gff_trimmer.cpp
	${GXX} gff_trimmer.cpp -o gff_trimmer -std=c++0x -O3

exon_num_filter: exon_num_filter.cpp
	${GXX} exon_num_filter.cpp -o exon_num_filter -std=c++0x -O3

longest_transcript: longest_transcript.cpp
	${GXX} longest_transcript.cpp -o longest_transcript -std=c++0x -O3

repeat_checker: repeat_checker.cpp
	${GXX} repeat_checker.cpp -o repeat_checker -std=c++0x -O3

strand_replace: strand_replace.cpp
	${GXX} strand_replace.cpp -o strand_replace -std=c++0x -O3

set_difference: set_difference.cpp
	${GXX} set_difference.cpp -o set_difference -std=c++0x -O3

tag_trimmer: tag_trimmer.cpp
	${GXX} tag_trimmer.cpp -o tag_trimmer -std=c++0x -O3

ORF_finder: ORF_finder.cpp
	${GXX} ORF_finder.cpp -o ORF_finder -std=c++0x -O3

install:
	install -s \
gff_trimmer exon_num_filter longest_transcript \
repeat_checker strand_replace set_difference \
tag_trimmer ORF_finder $(DEST); \
	install -s ORF_finder $(DEST2)

clean:
	rm -f \
gff_trimmer exon_num_filter longest_transcript \
repeat_checker strand_replace set_difference \
tag_trimmer ORF_finder \
$(DEST)/gff_trimmer $(DEST)/exon_num_filter \
$(DEST)/longest_transcript $(DEST)/repeat_checker \
$(DEST)/strand_replace $(DEST)/set_difference $(DEST)/tag_trimmer \
$(DEST)/ORF_finder $(DEST2)/ORF_finder 

###
all:
clean:

###
DEST=../../util/homology
GXX=g++

all: fastarepair fastarepair2 gff_2_proteinfasta flameshiftfilter install

fastarepair: fastarepair.cpp function.hpp
	${GXX} fastarepair.cpp -o fastarepair -std=c++0x -O3

fastarepair2: fastarepair2.cpp function.hpp
	${GXX} fastarepair2.cpp -o fastarepair2 -std=c++0x -O3

gff_2_proteinfasta: gff_2_proteinfasta.cpp
	${GXX} gff_2_proteinfasta.cpp -o gff_2_proteinfasta -std=c++0x -O3

flameshiftfilter: flameshiftgrep.cpp
	${GXX} flameshiftgrep.cpp -o flameshiftfilter -std=c++0x -O3

install:
	install -s fastarepair fastarepair2 gff_2_proteinfasta flameshiftfilter $(DEST)

clean:
	rm -f fastarepair fastarepair2 gff_2_proteinfasta flameshiftfilter \
$(DEST)/fastarepair $(DEST)/fastarepair2 $(DEST)/gff_2_proteinfasta $(DEST)/flameshiftfilter 

###
DEST=../../util/abinitio
GXX=g++

all: simple_low_norepeatmask inframe_stopcodon_exclude makefasta install

simple_low_norepeatmask: simple_low_norepeatmask.cpp
	${GXX} simple_low_norepeatmask.cpp -o simple_low_norepeatmask -std=c++0x -O3

inframe_stopcodon_exclude: inframe_stopcodon_exclude.cpp
	${GXX} inframe_stopcodon_exclude.cpp -o inframe_stopcodon_exclude -std=c++0x -O3

makefasta: makefasta.cpp
	${GXX} makefasta.cpp -o makefasta -std=c++0x -O3

install:
	install -s simple_low_norepeatmask inframe_stopcodon_exclude makefasta $(DEST)

clean:
	rm -f simple_low_norepeatmask inframe_stopcodon_exclude makefasta \
$(DEST)/simple_low_norepeatmask $(DEST)/inframe_stopcodon_exclude $(DEST)/makefasta

###
DEST=../../util/merge_phase0
GXX=g++

all: gff_editor Row2_rename RNA-seq_reform Spaln_reform \
Augustus_reform install

gff_editor: 190521_gff_editor.cpp function.hpp
	${GXX} 190521_gff_editor.cpp -o gff_editor -std=c++0x -O3

Row2_rename: row2_rename.cpp function.hpp
	${GXX} row2_rename.cpp -o Row2_rename -std=c++0x -O3

RNA-seq_reform: rnaseq_reform.cpp function.hpp
	${GXX} rnaseq_reform.cpp -o RNA-seq_reform -std=c++0x -O3

Spaln_reform: spaln_reform.cpp function.hpp
	${GXX} spaln_reform.cpp -o Spaln_reform -std=c++0x -O3

Augustus_reform: augustus_reform.cpp function.hpp
	${GXX} augustus_reform.cpp -o Augustus_reform -std=c++0x -O3

install:
	install -s \
gff_editor Row2_rename RNA-seq_reform Spaln_reform \
Augustus_reform $(DEST)

clean:
	rm -f \
gff_editor Row2_rename RNA-seq_reform Spaln_reform \
Augustus_reform \
$(DEST)gff_editor $(DEST)Row2_rename $(DEST)RNA-seq_reform $(DEST)Spaln_reform \
$(DEST)Augustus_reform

###
DEST=../../util/merge_phase1
GXX=g++

all: Grouping subgroup new_subgroup Searchalgo \
gff_editor initial_exon_polish install

Grouping: grouping.cpp
	${GXX} grouping.cpp -o Grouping -std=c++0x -O3

subgroup: subgroup_v2.2.cpp
	${GXX} subgroup_v2.2.cpp -o subgroup -std=c++0x -O3

new_subgroup: new_subgroup.cpp
	${GXX} new_subgroup.cpp -o new_subgroup -std=c++0x -O3

Searchalgo: searchalgo.cpp
	${GXX} searchalgo.cpp -o Searchalgo -std=c++0x -O3

gff_editor: gff_editor.cpp function.hpp
	${GXX} gff_editor.cpp -o gff_editor -std=c++0x -O3

initial_exon_polish: initial_exon_polish.cpp
	${GXX} initial_exon_polish.cpp -o initial_exon_polish -std=c++0x -O3

install:
	install -s \
Grouping subgroup new_subgroup Searchalgo \
gff_editor initial_exon_polish $(DEST)

clean:
	rm -f \
Grouping subgroup new_subgroup Searchalgo \
gff_editor initial_exon_polish gff_editor \
$(DEST)Grouping $(DEST)subgroup $(DEST)new_subgroup $(DEST)Searchalgo \
$(DEST)gff_editor $(DEST)initial_exon_polish

###
DEST=../../util/merge_phase2
GXX=g++

all: geneadd_v191115 geneadd_v191119 grouping_v1 install

geneadd_v191115: geneadd_v191115.cpp
	${GXX} geneadd_v191115.cpp -o geneadd_v191115 -std=c++0x -O3

geneadd_v191119: geneadd_v191119.cpp
	${GXX} geneadd_v191119.cpp -o geneadd_v191119 -std=c++0x -O3

grouping_v1: grouping_v1.cpp
	${GXX} grouping_v1.cpp -o grouping_v1 -std=c++0x -O3

install:
	install -s geneadd_v191115 geneadd_v191119 grouping_v1 $(DEST)

clean:
	rm -f \
geneadd_v191115 geneadd_v191119 grouping_v1 \
$(DEST)/geneadd_v191115 $(DEST)/geneadd_v191119 $(DEST)/grouping_v1

###
DEST=../../util/summary

all: final_reform install

final_reform: final_reform.cpp
	g++ final_reform.cpp -o final_reform -std=c++0x -O3

install:
	install -s final_reform $(DEST)
#	install final_reform $(DEST)

clean:
	rm -f final_reform $(DEST)/final_reform

###
DEST=../../util/evaluation

all: evaluate preevaluate install

evaluate: evaluation4.cpp
	g++ evaluation4.cpp -o evaluate -std=c++0x -O3

preevaluate: preevaluation.cpp
	g++ preevaluation.cpp -o preevaluate -std=c++0x -O3

install:
	install -s evaluate preevaluate $(DEST)

clean:
	rm -f evaluate preevaluate $(DEST)/evaluate $(DEST)/preevaluate
