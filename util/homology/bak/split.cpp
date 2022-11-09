//linuxのsplitコマンドだと-dとか-aで連番を指定しないといけないのが
//あまりに面倒なので、c++で書いた

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

//stringをintに変える関数
int StoI(string str)
{
	int number;
	istringstream iss(str);
	iss >> number;
	return number;
}

//intをstringに変える関数
string ItoS(int number)
{
	stringstream ss;
	ss << number;
	return ss.str();
}

int main(int argc, char*argv[])
{
//パラメータ数の取得
	if(argc != 3)
	{
		cout << "error:パラメーターの数が違います" << endl;
		cout << "fin.txt [行数]" << endl;
		return 0;
	}
//ファイルの読み込み
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:fin file not open" << endl;
		return 0;
	}

	string splitns;
	splitns = argv[2];
	int splitn;
	splitn = StoI(splitns);

	string lin, filename;
	ofstream fout;
	int count = 0;
	int count2 = 0;
	int checkp = 0;
//fileを行単位で読む
	while(getline(ifs,lin) )
	{
		if(count % splitn == 0)
		{
			if(checkp == 1)
			{
				fout.close();
			}
			filename = ItoS(count2);
			fout.open(filename);
			fout << lin << endl;
			checkp = 1;
			count2++;
		}
		else
		{
			fout << lin << endl;
		}
		count++;
	}
	return 0;
}
