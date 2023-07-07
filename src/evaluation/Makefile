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
