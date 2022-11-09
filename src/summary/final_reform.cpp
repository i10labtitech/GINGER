// ./a.out reference.fasta fout.txt

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include "/data2/tkazuki/onihitode/function.hpp"

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

	string lin, geneid, rnaid, exonid;
	int count = 0;
	int count2 = 0;
	vector<string> I;
//fileを行単位で読む
	while(getline(ifs,lin) )
	{
		I = Split(lin, '\t');
		if(I.size() > 8)
		{
			if(I[2] == "mRNA")
			{
				count++;
				count2 = 0;
				geneid = "gene" + ItoS(count);
				rnaid = geneid + ".t1";
				fout << I[0] << '\t' << I[1] << '\t' << "gene" << '\t' << I[3] << '\t'
					<< I[4] << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t'
					<< "ID=" << geneid << ";Name=ORF" << endl;
				fout << I[0] << '\t' << I[1] << '\t' << "mRNA" << '\t' << I[3] << '\t'
					<< I[4] << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t'
					<< "ID=" << rnaid << ";Parent=" << geneid << ";Name=ORF" << endl;
			}
			else if(I[2] == "CDS")
			{
				count2++;
				exonid = geneid + ".exon" + ItoS(count);
				fout << I[0] << '\t' << I[1] << '\t' << "exon" << '\t' << I[3] << '\t'
					<< I[4] << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t' 
					<< "ID=" << exonid << ";Parent=" << rnaid << ";Name=ORF" << endl;
				fout << I[0] << '\t' << I[1] << '\t' << "CDS" << '\t' << I[3] << '\t'
					<< I[4] << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t' 
					<< "ID=cds." << rnaid << ";Parent=" << rnaid << ";Name=ORF" << endl;
			}
		}
	}
	return 0;
}
