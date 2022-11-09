// ./a.out reference.fasta fout.txt

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include "/data/tkazuki/onihitode/function.hpp"

using namespace std;

int main(int argc, char*argv[])
{
  //パラメータ数の取得
	if(argc != 3)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "fin.txt fout.txt" << endl;
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
	vector<string> I;
//fileを行単位で読む
	while(getline(ifs,lin) )
	{
		I = Split(lin,'\t');
		fout << I[0] << '\t' << I[1] << '\t' << "CDS" << '\t' << I[3] << '\t' << I[4] << '\t' << "." << '\t' << I[6] << '\t' << "Name=" << I[8] << endl;
	}
	return 0;
}
