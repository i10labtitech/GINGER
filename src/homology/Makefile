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
