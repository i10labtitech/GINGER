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
