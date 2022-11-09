#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

vector<string> Split(const string &s, char delim);
string Replaced_V(string str, string replacestr, int replace);

int main(int argc, char*argv[])
{
  //パラメータ数の取得
	if(argc != 4)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "./a.out [input.gff3] [output.gff3] [2列目に入れたい文字列]" << endl;
		return 0;
	}
//ファイルの読み込み
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

	string row2;
	row2 = argv[3];

	string zero;

	string lin, tmplin;

//fileを行単位で読む
	while(getline(ifs,lin) )
	{
		tmplin = Replaced_V(lin, row2, 2);
		fout << tmplin << endl;
	}
	return 0;
}

//split関数 使い方 Split(文字列,diliminator('\t'など)）
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

//あるtabのやつを入れ替える
string Replaced_V(string str, string replacestr, int replace)
{
	string output;
	vector<string> I;
	I = Split(str, '\t');
	for(int i = 0; i < replace - 1; i++)
	{
		output += I[i];
		output += '\t';
	}
	output += replacestr;
	for(int i = replace; i < I.size(); i++)
	{
		output+= '\t';
		output+= I[i];
	}
	return(output);
}
