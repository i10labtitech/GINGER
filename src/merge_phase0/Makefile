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
