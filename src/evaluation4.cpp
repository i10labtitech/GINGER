#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>

using namespace std;

vector<string> Split(const string &s, char delim);
int StoI(string str);
string ItoS(int number);

struct group
{
	int num1;
	int num2;
};

int main(int argc, char*argv[])
{
	const string str1 = "gene";
	const string str2 = "exon";
	const string str3 = "gene";
	const string str4 = "exon";

	if(argc != 3)
	{
		cout << "Error: wrong arguments." << endl;
		cout << "./a.out [a file containing validation data] [a file containing evaluation data]" << endl;
		return 0;
	}
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:fin file not open" << endl;
		return 0;
	}
	ifstream ifs2(argv[2]);
	if(ifs2.fail())
	{
		cout << "error:fin file not open" << endl;
		return 0;
	}

	string zero;

	string lin, scf, stat, key, key2, key3, pos, sposn, eposn, strand;
	int spos, epos, hashvalue, sensitivity_bunbo_base, sensitivity_bunbo_exon, sensitivity_bunbo_gene;
	vector<string> I;
	group groups1;

//nucleotide level
	//key:scaffoldname + "pos" + pos 	value:flag
	unordered_map<string,group> hash_base;
	//出力用のvector
	vector<string> ON;

//exon level
	//key:scaffoldname + "spos" + spos + "epos" + epos + "strand" + strand
	//value:flag
	unordered_map<string,group> hash_exon;
	//出力用のvector
	vector<string> OE;

//gene level
	//key:scaffoldname + "strand" + strand + "_" + spos + "_" + epos + …
	//value:flag
	unordered_map<string,group> hash_gene;
	//出力用のvector
	vector<string> OG;

	//for gene
	int checkp = 0;
	string str;

	//gene count
	int seikaigenenum = 0;
	int predictgenenum = 0;

//loading data
	while(getline(ifs,lin) )
	{
		string::size_type a1 = lin.find("#");
		if(a1 == string::npos)
		{
			I = Split(lin,'\t');
			scf = I[0];
			stat = I[2];
			sposn = I[3];
			eposn = I[4];
			spos = StoI(sposn);
			epos = StoI(eposn);
			strand = I[6];
			if(stat == str2)
			{
				//nucleotide level
				for(int i = spos; i <= epos; i++)
				{
					pos = ItoS(i);
					key = scf + "pos" + pos;
					auto itr = hash_base.find(key);
					if( itr == hash_base.end() )
					{
						groups1 = {1,0};
						hash_base[key] = groups1;
						ON.push_back(key);
					}
					else
					{
						hashvalue = hash_base[key].num1;
						hashvalue++;
						hash_base[key].num1 = hashvalue;
					}
				}

				//exon level
				key2 = scf + "spos" + sposn + "epos" + eposn + "strand" + strand;
				auto itr = hash_exon.find(key2);
				if( itr == hash_exon.end() )
				{
					groups1 = {1,0};
					hash_exon[key2] = groups1;
					OE.push_back(key2);
				}
				else
				{
					hashvalue = hash_exon[key2].num1;
					hashvalue++;
					hash_exon[key2].num1 = hashvalue;
				}

				//gene level
				//second order
				str = "_" + sposn + "_" + eposn;
				key3 += str;
			}

		//gene level
			//first order
			else if(stat == str1 && checkp == 0)
			{
				key3 = scf + "strand" + strand;
				checkp = 1;
				seikaigenenum++;
			}
			//second order
			else if(stat == str1 && checkp == 1)
			{
				auto itr = hash_gene.find(key3);
				if( itr == hash_gene.end() )
				{
					groups1 = {1,0};
					hash_gene[key3] = groups1;
					OG.push_back(key3);
				}
				else
				{
					hashvalue = hash_gene[key3].num1;
					hashvalue++;
					hash_gene[key3].num1 = hashvalue;
				}
				key3 = scf + "strand" + strand;
				seikaigenenum++;
			}
		}
	}
	//gene level last order
	auto itr = hash_gene.find(key3);
	if( itr == hash_gene.end() )
	{
		groups1 = {1,0};
		hash_gene[key3] = groups1;
		OG.push_back(key3);
	}
	else
	{
		hashvalue = hash_gene[key3].num1;
		hashvalue++;
		hash_gene[key3].num1 = hashvalue;
	}
	sensitivity_bunbo_base = ON.size();
	sensitivity_bunbo_exon = OE.size();
	sensitivity_bunbo_gene = OG.size();

	checkp = 0;
	//key3 = zero;

//loading data
	while(getline(ifs2,lin) )
	{
		string::size_type a1 = lin.find("#");
		if(a1 == string::npos)
		{
			I = Split(lin,'\t');
			scf = I[0];
			stat = I[2];
			sposn = I[3];
			eposn = I[4];
			spos = StoI(sposn);
			epos = StoI(eposn);
			strand = I[6];
			if(stat == str4)
			{
			//nucleotide level
				for(int i = spos; i <= epos; i++)
				{
					pos = ItoS(i);
					key = scf + "pos" + pos;
					auto itr = hash_base.find(key);
					if( itr != hash_base.end() )
					{
						hashvalue = hash_base[key].num2;
						hashvalue++;
						hash_base[key].num2 = hashvalue;
					}
					else
					{
						groups1 = {0,1};
						hash_base[key] = groups1;
						ON.push_back(key);
					}
				}
			//exon level
				key2 = scf + "spos" + sposn + "epos" + eposn + "strand" + strand;
				auto itr = hash_exon.find(key2);
				if( itr != hash_exon.end() )
				{
					hashvalue = hash_exon[key2].num2;
					hashvalue++;
					hash_exon[key2].num2 = hashvalue;
				}
				else
				{
					groups1 = {0,1};
					hash_exon[key2] = groups1;
					OE.push_back(key2);
				}

				//gene level second order
				str = "_" + sposn + "_" + eposn;
				key3 += str;
			}

		//gene level
			//first order
			else if(stat == str3 && checkp == 0)
			{
				key3 = scf + "strand" + strand;
				checkp = 1;
				predictgenenum++;
			}
			//second order
			else if(stat == str3 && checkp == 1)
			{
			//	cout << key3 << endl;
				auto itr = hash_gene.find(key3);
				if( itr != hash_gene.end() )
				{
					hashvalue = hash_gene[key3].num2;
					hashvalue++;
					hash_gene[key3].num2 = hashvalue;
				}
				else
				{
					groups1 = {0,1};
					hash_gene[key3] = groups1;
					OG.push_back(key3);
				}
				key3 = scf + "strand" + strand;
				predictgenenum++;
			}
		}
	}

	//cout << key3 << endl;
//gene level last order
	auto itr2 = hash_gene.find(key3);
	if( itr2 == hash_gene.end() )
	{
		groups1 = {0,1};
		hash_gene[key3] = groups1;
		OG.push_back(key3);
	}
	else
	{
		hashvalue = hash_gene[key3].num2;
		hashvalue++;
		hash_gene[key3].num2 = hashvalue;
	}

//sensitivity of nuc level, ppv
	int sensitivity_bunsi = 0;
	int ppv_bunbo = 0;
	int ppv_bunsi = 0;
	double nsensitivity, nppv, esensitivity, eppv, gsensitivity, gppv;
	for(int i = 0; i < ON.size(); i++)
	{
		if(hash_base[ON[i]].num1 >= 1 && hash_base[ON[i]].num2 >= 1)
		{
			sensitivity_bunsi++;
		}
		ppv_bunbo += hash_base[ON[i]].num2;
		if(hash_base[ON[i]].num2 <= hash_base[ON[i]].num1)
		{
			ppv_bunsi += hash_base[ON[i]].num2;
		}
		else
		{
			ppv_bunsi += hash_base[ON[i]].num1;
		}
	}
	nsensitivity = 100 * (double)sensitivity_bunsi / (double)sensitivity_bunbo_base;
	nppv = 100 * (double)ppv_bunsi / (double)ppv_bunbo;


	sensitivity_bunsi = 0;
	ppv_bunbo = 0;
	ppv_bunsi = 0;
//sensitivity of exon level, ppv
	for(int i = 0; i < OE.size(); i++)
	{
		if(hash_exon[OE[i]].num1 >= 1 && hash_exon[OE[i]].num2 >= 1)
		{
			sensitivity_bunsi++;
		}
		ppv_bunbo += hash_exon[OE[i]].num2;
		if(hash_exon[OE[i]].num2 <= hash_exon[OE[i]].num1)
		{
			ppv_bunsi += hash_exon[OE[i]].num2;
		}
		else
		{
			ppv_bunsi += hash_exon[OE[i]].num1;
		}
	}
	esensitivity = 100 * (double)sensitivity_bunsi / (double)sensitivity_bunbo_exon;
	eppv = 100 * (double)ppv_bunsi / (double)ppv_bunbo;
	
	sensitivity_bunsi = 0;
	ppv_bunbo = 0;
	ppv_bunsi = 0;
//sensitivity of gene level, ppv
	for(int i = 0; i < OG.size(); i++)
	{
		if(hash_gene[OG[i]].num1 >= 1 && hash_gene[OG[i]].num2 >= 1)
		{
			sensitivity_bunsi++;
		}
		ppv_bunbo += hash_gene[OG[i]].num2;
		if(hash_gene[OG[i]].num2 <= hash_gene[OG[i]].num1)
		{
			ppv_bunsi += hash_gene[OG[i]].num2;
		}
		else
		{
			ppv_bunsi += hash_gene[OG[i]].num1;
		}
	}

	gsensitivity = 100 * (double)sensitivity_bunsi / (double)sensitivity_bunbo_gene;
	gppv = 100 * (double)ppv_bunsi / (double)ppv_bunbo;

	//F-score
	double nuclfscore, exonfscore, genefscore;
	nuclfscore = 2 * nsensitivity * nppv / (nsensitivity + nppv);
	exonfscore = 2 * esensitivity * eppv / (esensitivity + eppv);
	genefscore = 2 * gsensitivity * gppv / (gsensitivity + gppv);

	cout << endl;

	cout << "#correct" << '\t' << seikaigenenum << endl;
	cout << "#predict" << '\t' << predictgenenum << endl;
	
	cout << endl;
	cout << "level" << '\t' << "sens" << '\t' << "ppv" << '\t' << "fscore" << endl;
	cout << "nucl" << '\t' << nsensitivity << '\t' << nppv << '\t' << nuclfscore << endl;
	cout << "exon" << '\t' << esensitivity << '\t' << eppv << '\t' << exonfscore << endl;
	cout << "gene" << '\t' << gsensitivity << '\t' << gppv << '\t' << genefscore << endl;

	cout << endl;
	cout << "nucl sensitivity" << '\t' << "nucl ppv" << '\t' << "nucl fscore" << '\t' << "exon sensitivity" << '\t' << "exon ppv" << '\t' << "exon fscore" << '\t' << "gene sensitivity" << '\t' << "gene ppv" << '\t' << "gene fscore" << endl;
	cout << nsensitivity << '\t' << nppv << '\t' << nuclfscore << '\t' << esensitivity << '\t' << eppv << '\t' << exonfscore << '\t' << gsensitivity << '\t' << gppv << '\t' << genefscore << endl;
	
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
