// ./a.out reference.fasta fout.txt

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "function.hpp"

using namespace std;

int main(int argc, char*argv[])
{
  //パラメータ数の取得
	if(argc != 3)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "[fasta.fa] [fout]" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:fasta file not open" << endl;
		return 0;
	}
//出力ファイル
	ofstream fout(argv[2]);
	if(fout.fail())
	{
		cout << "error:fout file not open" << endl;
		return 0;
	}

	string lin, zero;
	int alignlength;
	int checkp = 0;
	string genen, fas, fas1, fas2;
	while(getline(ifs, lin))
	{
		string::size_type a1 = lin.find(">");
		if(a1 != string::npos && checkp == 0)
		{
			checkp = 1;
			genen = lin.substr(1);
		}
		else if(a1 != string::npos && checkp == 1)
		{
			fas1 = fas;
			fas = zero;
		}
		else if(a1 == string::npos && checkp == 1)
		{
			fas += lin;
		}
	}
	fas2 = fas;

	//ifout << genen << endl;
	int count = 0;
	int fascount = 0;
	int fascount2 = 0;
	int difflength;
	for(int i = 0; i < fas1.size(); i++)
	{
		if(fas1[i] == fas2[i])
		{
			count++;
		}
		if(fas1[i] != '-')
		{
			fascount++;
		}
		if(fas2[i] != '-')
		{
			fascount2++;
		}
	}

	if(fascount2 >= fascount)
	{
		difflength = fascount2 - fascount;
	}
	else
	{
		difflength = fascount - fascount2;
	}

	alignlength = fas1.length();
	double identity;
	identity = 100* ( (double)count / (double)alignlength);
	fout << genen << '\t' << identity << '\t' << difflength << endl;
	return 0;
}
