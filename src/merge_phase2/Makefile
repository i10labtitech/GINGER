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
