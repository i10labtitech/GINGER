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

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

vector<string> Split(const string &s, char delim);
string ItoS(int number);
string Replaced_V(string str, string replacestr, int replace);

int main(int argc, char*argv[])
{

	if(argc != 3)
	{
		cout << "error" << endl;
		cout << "./a.out [spaln.gff3] [output.gff3]" << endl;
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

	string lin;

	vector<string> I, E, C;
	string id0, id1, id2, id3, tmplin, mrnalin1, mrnalin2;
	string exonlin1, exonlin2, exonlin3, cdslin1, cdslin2, info, counts;
	int count = 0;
	int checkp = 0;

	while(getline(ifs,lin) )
	{
		I = Split(lin, '\t');
		if(I[2] == "gene")
		{
			if(checkp == 1)
			{
				for(int i = 0; i < E.size(); i++)
				{
					fout << E[i] << endl;
				}
				for(int i = 0; i < C.size(); i++)
				{
					fout << C[i] << endl;
				}
			}
			string::size_type a0 = I[8].find(";");
			id0 = I[8].substr(3, a0-3);
			string::size_type a1 = I[8].find("Name=");
			id1 = I[8].substr(a1+5);
			fout << lin << endl;
			tmplin = lin;
			count = 0;
			checkp = 0;
			E.clear();
			C.clear();
		}
		else if(I[2] == "cds")
		{
			string::size_type a2 = I[8].find("Parent=");
			string::size_type a3 = I[8].find(";", a2+1);
			id2 = I[8].substr(a2+7, a3-a2-7);
			
			info = "ID=" + id2 + ";Parent=" + id0 + ";Name=" + id1;
			mrnalin1 = Replaced_V(tmplin, "mRNA", 3);
			mrnalin2 = Replaced_V(mrnalin1, info, 9);
			if(checkp == 0)
			{
				fout << mrnalin2 << endl;
				checkp = 1;
			}

			string::size_type a4 = I[8].find(";Name=");
			id3 = I[8].substr(a4);

			count++;
			counts = ItoS(count);
			info = "ID=" + id2 + ".exon" + counts + ";Parent=" + id2 + id3;
			exonlin1 = Replaced_V(lin, "exon", 3);
			exonlin2 = Replaced_V(exonlin1, ".", 8);
			exonlin3 = Replaced_V(exonlin2, info, 9);
			E.push_back(exonlin3);

			info = "ID=cds." + id2 + ";Parent=" + id2 + id3;
			cdslin1 = Replaced_V(lin, "CDS", 3);
			cdslin2 = Replaced_V(cdslin1, info, 9);
			C.push_back(cdslin2);
		}
	}
	for(int i = 0; i < E.size(); i++)
	{
		fout << E[i] << endl;
	}
	for(int i = 0; i < C.size(); i++)
	{
		fout << C[i] << endl;
	}
	return 0;
}


vector<string> Split(const string &s, char delim)
{
	vector<string> elems;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim))
	{
		if (!item.empty())
		{
			elems.push_back(item);
		}
	}
	return elems;
}

string ItoS(int number)
{
	stringstream ss;
	ss << number;
	return ss.str();
}

string Replaced_V(string str, string replacestr, int replace)
{
	string output;
	vector<string> I;
	I = Split(str, '\t');
	for(int i = 0; i < replace - 1; i++)
	{
		output += I[i];
		output += '\t';
	}
	output += replacestr;
	for(int i = replace; i < I.size(); i++)
	{
		output+= '\t';
		output+= I[i];
	}
	return(output);
}
