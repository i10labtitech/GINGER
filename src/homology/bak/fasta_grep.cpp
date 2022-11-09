// ./a.out reference.fasta fout.txt
// 対応するfastaを抽出してくる

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
	if(argc != 4)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "[fasta.fa] [gff] [fout]" << endl;
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
		cout << "error:gff file not open" << endl;
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
	string lin, fas, fasname, fasname2;

	unordered_map<string, string> hash;
	while(getline(ifs,lin))
	{
		string::size_type a1 = lin.find(">");
		if(a1 != string::npos && fasname2 == zero)
		{
			fasname = lin;
			string::size_type a2 = lin.find(" ");
			fasname2 = lin.substr(1,a2-1);
		}
		else if(a1 != string::npos)
		{
			hash[fasname2] = fas;
			fas = zero;
			fasname = lin;
			string::size_type a2 = lin.find(" ");
			fasname2 = lin.substr(1,a2-1);
		}
		else
		{
			fas += lin;
		}
	}
	hash[fasname2] = fas;
	fas = zero;

	cout << "check" << endl;

	fasname = zero;
	vector<string> I;
	int checkp;
//fileを行単位で読む
	while(getline(ifs2,lin))
	{
		string::size_type a1 = lin.find("#");
		if(a1 == string::npos)
		{
			I = Split(lin,'\t');
			if(I[2] == "gene")
			{
				checkp = 0;
			}
			else if(I[2] == "cds" && checkp == 0)
			{
				string::size_type a2 = I[8].find("Target=");
				string::size_type a3 = I[8].find(" ");
				fasname = I[8].substr(a2+7, a3-a2-7);
				fout << ">" << fasname << endl;
				fout << hash[fasname] << endl;
				checkp = 1;
			}
		}
	}
	return 0;
}
