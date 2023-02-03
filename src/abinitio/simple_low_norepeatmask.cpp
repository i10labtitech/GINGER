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
// #include "/data2/yuasa/script/abinitio/script/function.hpp"

using namespace std;

struct group
{
	int gspos;
	int gepos;
};

//split関数 使い方 Split(文字列,diliminator('\t'など)）
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

//stringをintに変える関数
int StoI(string str)
{
	int number;
	istringstream iss(str);
	iss >> number;
	return number;
}

int main(int argc, char*argv[])
{
  //パラメータ数の取得
	if(argc != 4)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "./a.out [fasta] [repeatmask.out] [fout]" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:fasta file not open" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs2(argv[2]);
	if(ifs2.fail())
	{
		cout << "error:repeatmask.out file not open" << endl;
		return 0;
	}
//出力ファイル
	ofstream fout(argv[3]);
	if(fout.fail())
	{
		cout << "error:fout file not open" << endl;
		return 0;
	}

	string zero;


	vector<string> I;
	string lin, scf, stat;
	group groups1;
	int spos, epos;
	int count = 0;
	unordered_multimap<string, group> hash;
//fileを行単位で読む
	while(getline(ifs2,lin) )
	{
		count++;
		if(count > 3)
		{
			I = Split(lin, ' ');
			stat = I[10];
			if(stat != "Simple_repeat" && stat != "Low_complexity")
			{
				scf = I[4];
				spos = StoI(I[5]);
				epos = StoI(I[6]);
				groups1 = {spos, epos};
				hash.insert(pair<string, group>(scf, groups1));
			}
		}
	}

//flleを行単位で読む
	int cutfasnum, check;
	string cutfas, fasta;
	int genelength = 0;
	int maskcount = 0;

	while(getline(ifs, lin) )
	{
		string::size_type a1 = lin.find(">");
		if(a1 != string::npos)
		{
			if(fasta != zero)
			{
				auto range = hash.equal_range(scf);
				for (auto iterator = range.first; iterator != range.second; iterator++)
				{
					auto target = *iterator;
					spos = target.second.gspos;
					epos = target.second.gepos;
					for(int i = spos - 1; i < epos; i++)
					{
						if(fasta[i] != 'N')
						{
							fasta[i] = 'N';
							maskcount++;
						}
					}
				}
				check = 0;
				genelength += fasta.length();
				cutfasnum = fasta.length() / 50;
				check = fasta.length() % 50;

				for(int i = 0; i <= cutfasnum; i++)
				{
					if(i == cutfasnum && check != 0)
					{
						cutfas = fasta.substr(i*50);
						fout << cutfas << endl;
					}
					else
					{
						cutfas = fasta.substr(i*50, 50);
						fout << cutfas << endl;
					}
				}
			}


			string::size_type a2 = lin.find(" ");
			if(a2 != string::npos)
			{
				scf = lin.substr(1,a2-1);
			}
			else
			{
				scf = lin.substr(1);
			}
			fout << ">" << scf << endl;
			fasta = zero;
		}
		else
		{
			fasta += lin;
		}
	}

	auto range = hash.equal_range(scf);
	for (auto iterator = range.first; iterator != range.second; iterator++)
	{
		auto target = *iterator;
		spos = target.second.gspos;
		epos = target.second.gepos;
		for(int i = spos - 1; i < epos; i++)
		{
			if(fasta[i] != 'N')
			{
				fasta[i] = 'N';
				maskcount++;
			}
		}
	}
	check = 0;
	genelength += fasta.length();
	cutfasnum = fasta.length() / 50;
	check = fasta.length() % 50;

	for(int i = 0; i <= cutfasnum; i++)
	{
		if(i == cutfasnum && check != 0)
		{
			cutfas = fasta.substr(i*50);
			fout << cutfas << endl;
		}
		else
		{
			cutfas = fasta.substr(i*50, 50);
			fout << cutfas << endl;
		}
	}

	double rate;
	rate = 100 * (double)maskcount / (double)genelength;
	cout << "mask(bp)" << '\t' << maskcount << endl;
	cout << "genomelength(bp)" << '\t' << genelength << endl;
	cout << "maskrate(%)" << '\t' << rate << endl;
	return 0;
}
