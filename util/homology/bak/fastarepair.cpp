//reference fastaを修正する 
//塩基配列を小文字を大文字にとか、fastaの名前を簡略化したりとか

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
		cout << "./fastarepair [prerepair.fa] [output]" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:fin file not open" << endl;
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

	string lin;

	string tmpstr;
//fileを行単位で読む
	while(getline(ifs,lin) )
	{
		string::size_type p1 = lin.find(">");
		if(p1 != string::npos)
		{
			string::size_type p2 = lin.find(" ");
			if(p2 != string::npos)
			{
				tmpstr = lin.substr(0,p2);
				fout << tmpstr << endl;
			}
			else
			{
				fout << lin << endl;
			}
		}
		else
		{
			tmpstr = stoS(lin);
			fout << tmpstr << endl;
		}
	}
	return 0;
}
