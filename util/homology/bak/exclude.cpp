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
  //パラメータ数の取得
	if(argc != 4)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "[list] [gff] [fout]" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:fin file not open" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs2(argv[2]);
	if(ifs2.fail())
	{
		cout << "error:fin file not open" << endl;
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

	string lin, genen0, genen, fas;
	unordered_map<string,int> hash;
	vector<string> V;
	int checkp = 0;
//fileを行単位で読む
	while(getline(ifs,lin) )
	{
		hash[lin] = 0;
	}
	vector<string> I;
//fileを行単位で読む
	while(getline(ifs2,lin) )
	{
		checkp = 0;
		I = Split(lin,'\t');
		if(I.size() > 8)
		{
			if(I[2] == "gene")
			{
				string::size_type a1 = I[8].find(";");
				genen = I[8].substr(3,a1-3);
				auto itr = hash.find(genen);
				if( itr == hash.end() )
				{
					fout << lin << endl;
				}
			}
			if(I[2] == "cds")
			{
				string::size_type a1 = I[8].find("Parent=");
				string::size_type a2 = I[8].find(";",a1+1);
				genen0 = I[8].substr(a1+11, a2-a1-11);
				genen = "gene" + genen0;
				auto itr = hash.find(genen);
				if( itr == hash.end() )
				{
					fout << lin << endl;
				}
			}
		}
	}
	return 0;
}
