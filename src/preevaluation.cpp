#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <algorithm>

using namespace std;

vector<string> Split(const string &s, char delim);
int StoI(string str);
string ItoS(int number);
string GFF3Change(string str, string str2);

struct group
{
	int gspos;
	int gepos;
	string glin;
};

int main(int argc, char*argv[])
{

	if(argc != 5)
	{
		cout << "Error: wrong arguments." << endl;
		cout << "./a.out [Input GFF file] [Converted file] [Delimiter (mRNA or gene)][Target (CDS)]" << endl;
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
	string str1 = argv[4];
	string str2 = argv[3];

	string zero;

	string lin, scf, stat, seplin, olin;
	int spos, epos, max;
	group groups1;
	vector<string> I;
	vector<group> G;
	vector<int> S;

	int checkp = 0;
	string str;

	while(getline(ifs,lin) )
	{
		string::size_type a1 = lin.find("#");
		if(a1 == string::npos)
		{
			I = Split(lin,'\t');
			if(I.size() > 7)
			{
				stat = I[2];
				if(stat == str2 && checkp == 0)
				{
					checkp = 1;
					seplin = lin;
				}
				else if(stat == str2 && checkp != 0)
				{
					checkp = 1;
					sort(S.begin(), S.end());
					olin = GFF3Change(seplin, "gene");
					fout << olin << endl;
					for(int i = 0; i < S.size(); i++)
					{
						for(int j = 0; j < G.size(); j++)
						{
							if(G[j].gspos == S[i])
							{
								olin = GFF3Change(G[j].glin, "exon");
								fout << olin << endl;
							}
						}
					}
					G.clear();
					S.clear();
					seplin = lin;
				}
				else if(stat == str1 && (checkp == 1 || checkp == 2) )
				{
					spos = StoI(I[3]);
					epos = StoI(I[4]);
					groups1 = {spos, epos, lin};
					G.push_back(groups1);
					S.push_back(spos);
					checkp = 2;
				}
				else if(stat != str1 && stat != str2 && checkp == 2)
				{
					checkp = 0;
					sort(S.begin(), S.end());
					olin = GFF3Change(seplin,"gene");
					fout << olin << endl;
					for(int i = 0; i < S.size(); i++)
					{
						for(int j = 0; j < G.size(); j++)
						{
							if(G[j].gspos == S[i])
							{
								olin = GFF3Change(G[j].glin,"exon");
								fout << olin << endl;
							}
						}
					}
					G.clear();
					S.clear();
				}
			}
		}
	}
	if(checkp != 0)
	{
		sort(S.begin(), S.end());
		olin = GFF3Change(seplin,"gene");
		fout << olin << endl;
		for(int i = 0; i < S.size(); i++)
		{
			for(int j = 0; j < G.size(); j++)
			{
				if(G[j].gspos == S[i])
				{
					olin = GFF3Change(G[j].glin,"exon");
					fout << olin << endl;
				}
			}
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

string GFF3Change(string str, string str2)
{
	string str3;
	string::size_type a1 = str.find('\t');
	string::size_type a2 = str.find('\t', a1+1);
	string::size_type a3 = str.find('\t', a2+1);
	str3 = str.substr(0,a2) + '\t' + str2 + '\t' + str.substr(a3+1);
	return str3;
}
