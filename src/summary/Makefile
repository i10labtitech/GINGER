DEST=../../util/summary

all: final_reform install

final_reform: final_reform.cpp
	g++ final_reform.cpp -o final_reform -std=c++0x -O3

install:
	install -s final_reform $(DEST)
#	install final_reform $(DEST)

clean:
	rm -f final_reform $(DEST)/final_reform
