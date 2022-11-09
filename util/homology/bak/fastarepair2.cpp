//protein配列の修正
//変な文字列の削除とか、geneの名前の簡略化とか

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
		cout << "./fastarepair2 [prerepair.faa] [output]" << endl;
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

	string tmpstr, fas, cutfas, cutfas2;
	int length;
//fileを行単位で読む
	while(getline(ifs,lin) )
	{
		string::size_type p1 = lin.find(">");
		if(p1 != string::npos)
		{
			if(fas != zero)
			{
		//ありえない文字列の除去
				cutfas = Erace(fas," ");
				cutfas2 = Erace(cutfas,"*");
				cutfas = Erace(cutfas2,"B");
				cutfas2 = Erace(cutfas,"J");
				cutfas = Erace(cutfas2,"O");
				cutfas2 = Erace(cutfas,"U");
				cutfas = Erace(cutfas2,"X");
				cutfas2 = Erace(cutfas,"Z");
				cutfas = Erace(cutfas2,"#");
				fout << cutfas << endl;
			}
			string::size_type p2 = lin.find(" ");
			if(p2 != string::npos)
			{
				tmpstr = lin.substr(0, p2);
				fout << tmpstr << endl;
			}
			else
			{
				fout << lin << endl;
			}
			fas = zero;
		}
		else
		{
			fas += lin;
		}
	}
	cutfas = Erace(fas," ");
	cutfas2 = Erace(cutfas,"*");
	cutfas = Erace(cutfas2,"B");
	cutfas2 = Erace(cutfas,"J");
	cutfas = Erace(cutfas2,"O");
	cutfas2 = Erace(cutfas,"U");
	cutfas = Erace(cutfas2,"X");
	cutfas2 = Erace(cutfas,"Z");
	cutfas = Erace(cutfas2,"#");
	fout << cutfas << endl;
	return 0;
}
