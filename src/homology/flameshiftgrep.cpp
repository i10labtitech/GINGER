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
	if(argc != 3)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "[protein] [fout]" << endl;
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

	string lin, genen, fas;
	vector<string> V;
	int checkp = 0;
//fileを行単位で読む
	while(getline(ifs,lin) )
	{
		string::size_type p1 = lin.find(">");
		if(p1 != string::npos)
		{
			genen = lin.substr(1);
		}
		else
		{
			fas = lin;
			checkp = 0;
			for(int i = 0; i < fas.length() - 1; i++)
			{
				if(fas[i] == '*')
				{
					fout << genen << endl;
					checkp = 1;
					break;
				}
			}
			if(checkp == 1)
			{
				continue;
			}
		}
	}
	return 0;
}
