DEST=../../util/denovo

all: faspro install

faspro: faspro.cpp
	g++ faspro.cpp -o faspro -std=c++0x -O3

install:
	install -s faspro $(DEST)

clean:
	rm -f faspro $(DEST)/faspro
