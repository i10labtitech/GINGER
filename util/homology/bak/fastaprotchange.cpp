// ./a.out reference.fasta fout.txt
// 塩基配列をタンパク質配列に変換

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
	if(argc != 3)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "./a.out [nuclfasta.fa] [fout.fa]" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:nuclfasta.fa file not open" << endl;
		return 0;
	}
//出力ファイル
	ofstream fout(argv[2]);
	if(fout.fail())
	{
		cout << "error:fout file not open" << endl;
		return 0;
	}

	string zero;

	string lin, fas, prot;
//fasta fileをhashに格納
	while(getline(ifs,lin) )
	{
		string::size_type p1 = lin.find(">");
		if(p1 != string::npos && fas == zero)
		{
			fout << lin << endl;
		}
		else if(p1 != string::npos && fas != zero)
		{
			prot = PChange(fas);
			fout << prot << endl;
			fas = zero;
			fout << lin << endl;
		}
		else
		{
			fas += lin;
		}
	}
	prot = PChange(fas);
	fout << prot << endl;
	return 0;
}
