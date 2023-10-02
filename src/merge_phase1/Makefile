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
