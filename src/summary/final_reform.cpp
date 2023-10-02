/*
Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology

This file is part of GINGER.

GINGER is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GINGER is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with GINGER; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


// ./a.out reference.fasta fout.txt

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include "function.hpp"

using namespace std;

int main(int argc, char*argv[])
{

	if(argc != 3)
	{
		cout << "error: arguments may be wrong." << endl;
		cout << "./final_reform fin.txt fout.txt" << endl;
		return 0;
	}

	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:fin file not open" << endl;
		return 0;
	}

	ofstream fout(argv[2]);
	if(fout.fail())
	{
		cout << "error:fout file not open" << endl;
		return 0;
	}

	string zero;

	string lin, geneid, rnaid, exonid;
	int count = 0;
	int count2 = 0;
	vector<string> I;

	while(getline(ifs,lin) )
	{
		I = Split(lin, '\t');
		if(I.size() > 8)
		{
			if(I[2] == "mRNA")
			{
				count++;
				count2 = 0;
				geneid = "gene" + ItoS(count);
				rnaid = geneid + ".t1";
				fout << I[0] << '\t' << I[1] << '\t' << "gene" << '\t' << I[3] << '\t'
					<< I[4] << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t'
					<< "ID=" << geneid << ";Name=ORF" << endl;
				fout << I[0] << '\t' << I[1] << '\t' << "mRNA" << '\t' << I[3] << '\t'
					<< I[4] << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t'
					<< "ID=" << rnaid << ";Parent=" << geneid << ";Name=ORF" << endl;
			}
			else if(I[2] == "CDS")
			{
				count2++;
				exonid = geneid + ".exon" + ItoS(count);
				fout << I[0] << '\t' << I[1] << '\t' << "exon" << '\t' << I[3] << '\t'
					<< I[4] << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t' 
					<< "ID=" << exonid << ";Parent=" << rnaid << ";Name=ORF" << endl;
				fout << I[0] << '\t' << I[1] << '\t' << "CDS" << '\t' << I[3] << '\t'
					<< I[4] << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t' 
					<< "ID=cds." << rnaid << ";Parent=" << rnaid << ";Name=ORF" << endl;
			}
		}
	}
	return 0;
}
