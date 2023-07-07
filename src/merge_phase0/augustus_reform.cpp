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
#include <algorithm>

using namespace std;


vector<string> Split(const string &s, char delim);
int StoI(string str);
string ItoS(int number);
string Replaced_V(string str, string replacestr, int replace);


int main(int argc, char*argv[])
{

	if(argc != 3)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "./a.out [augustus,gtf] [output.gff]" << endl;
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
	vector<int> P;
	string tmplin, tmplin2, tmplin3, genelin, mrnalin, counts, sposs, eposs, strand, info, info1, info2;
	int count = 0;
	int checkp = 0;
	int tmpi, tmpi2, pos, spos, epos;

	while(getline(ifs,lin) )
	{
		I = Split(lin, '\t');
		//		if(I[2] == "gene")
		if(I[2] == "transcript" || I[2] == "mRNA")
		{
			if(checkp == 1)
			{
				spos = *min_element(P.begin(), P.end());
				sposs = ItoS(spos);
				epos = *max_element(P.begin(), P.end());
				eposs = ItoS(epos);
				//tmplin = Replaced_V(genelin, sposs, 4);
				//tmplin2 = Replaced_V(tmplin, eposs, 5);
				//fout << tmplin << endl;
				tmplin = Replaced_V(mrnalin, sposs, 4);
				tmplin2 = Replaced_V(tmplin, eposs, 5);
				fout << tmplin << endl;
				if(strand == "+")
				{
					for(int i = 0; i < E.size(); i++)
					{
						tmpi2 = i + 1;
						counts = ItoS(tmpi2);
						info = "ID=" + info2 + ".exon" + counts + ";Parent=" + info2 + ";Name=Augustus%20prediction";
						tmplin = Replaced_V(E[i], info, 9);
						tmplin2 = Replaced_V(tmplin, ".", 8);
						tmplin3 = Replaced_V(tmplin2, ".", 6);
						fout << tmplin3 << endl;
					}
					for(int i = 0; i < C.size(); i++)
					{
						fout << C[i] << endl;
					}
				}
				else if(strand == "-")
				{
					for(int i = 0; i < E.size(); i++)
					{
						tmpi = E.size() - i - 1;
						tmpi2 = i + 1;
						counts = ItoS(tmpi2);
						info = "ID=" + info2 + ".exon" + counts + ";Parent=" + info2 + ";Name=Augustus%20prediction";
						tmplin = Replaced_V(E[tmpi], info, 9);
						tmplin2 = Replaced_V(tmplin, ".", 8);
						tmplin3 = Replaced_V(tmplin2, ".", 6);
						fout << tmplin3 << endl;
					}
					for(int i = 0; i < C.size(); i++)
					{
						tmpi = C.size() - i - 1;
						fout << C[tmpi] << endl;
					}
				}
			}
			checkp = 1;
			count = 0;
			E.clear();
			C.clear();
			P.clear();
			info1 = I[0] + "-" + I[8];
			info = "ID=" + info1 + ";Name=Augustus%20prediction";
			genelin = Replaced_V(lin, info, 9);
			strand = I[6];

			info2 = I[0] + "-" + I[8];// + ".t1";
			info = "ID=" + info2 + ";Parent=" + info1 + ";Name=Augustus%20prediction";
			tmplin = Replaced_V(lin, "mRNA", 3);
			mrnalin = Replaced_V(tmplin, info, 9);
		}
		else if(I[2] == "CDS")
		{
			pos = StoI(I[3]);
			P.push_back(pos);
			pos = StoI(I[4]);
			P.push_back(pos);

			tmplin = Replaced_V(lin, "exon", 3);
			E.push_back(tmplin);
			info = "ID=cds." + info2 + ";Parent=" + info2 + ";Name=Augustus%20prediction";
			tmplin = Replaced_V(lin, info, 9);
			C.push_back(tmplin);
		}
	}
	spos = *min_element(P.begin(), P.end());
	sposs = ItoS(spos);
	epos = *max_element(P.begin(), P.end());
	eposs = ItoS(epos);
	//tmplin = Replaced_V(genelin, sposs, 4);
	//tmplin2 = Replaced_V(tmplin, eposs, 5);
	//fout << tmplin << endl;
	tmplin = Replaced_V(mrnalin, sposs, 4);
	tmplin2 = Replaced_V(tmplin, eposs, 5);
	fout << tmplin << endl;
	if(strand == "+")
	{
		for(int i = 0; i < E.size(); i++)
		{
			tmpi2 = i + 1;
			counts = ItoS(tmpi2);
			info = "ID=" + info2 + ".exon" + counts + ";Parent=" + info2 + ";Name=Augustus%20prediction";
			tmplin = Replaced_V(E[i], info, 9);
			tmplin2 = Replaced_V(tmplin, ".", 8);
			tmplin3 = Replaced_V(tmplin2, ".", 6);
			fout << tmplin3 << endl;
		}
		for(int i = 0; i < C.size(); i++)
		{
			fout << C[i] << endl;
		}
	}
	else if(strand == "-")
	{
		for(int i = 0; i < E.size(); i++)
		{
			tmpi = E.size() - i - 1;
			tmpi2 = i + 1;
			counts = ItoS(tmpi2);
			info = "ID=" + info2 + ".exon" + counts + ";Parent=" + info2 + ";Name=Augustus%20prediction";
			tmplin = Replaced_V(E[tmpi], info, 9);
			tmplin2 = Replaced_V(tmplin, ".", 8);
			tmplin3 = Replaced_V(tmplin2, ".", 6);
			fout << tmplin3 << endl;
		}
		for(int i = 0; i < C.size(); i++)
		{
			tmpi = C.size() - i - 1;
			fout << C[i] << endl;
		}
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

int StoI(string str)
{
	int number;
	istringstream iss(str);
	iss >> number;
	return number;
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
