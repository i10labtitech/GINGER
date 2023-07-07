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
