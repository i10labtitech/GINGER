// ./a.out reference.fasta fout.txt
// spalnのGFFをfastaに変換

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
		cout << "./a.out [gfffile] [fastafile] [foutfile]" << endl;
		return 0;
	}
//ソiートしたいファイルの読み込み
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:fin file not open" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs2(argv[2]);
	if(ifs.fail())
	{
		cout << "error:fin2 file not open" << endl;
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

	string lin;
	vector<vector<string> > FF;
	FF = Fasta_V(ifs2);

	vector<string> V;
	int spos, epos;
	int checkp = 0;
	string status, scf, strand, info, rna, fas, cutfas, tmpfas;
//gff fileからscname = I[8].substr(a2+9,a3-a2-9);ffold情報を得る
	while(getline(ifs,lin) )
	{
		V = Split(lin,'\t');
		if(V.size() > 7)
		{
			status = V[2];
			if(status == "gene")
			{
				scf = V[0];
				strand = V[6];
				string::size_type a1 = V[8].find("ID=");
				string::size_type a2 = V[8].find(";", a1+9);
				rna = V[8].substr(a1+3,a2-a1-3);
				if(checkp == 1)
				{
					fout << endl;
				}
				checkp = 1;
				fout << ">" << rna << endl;
			}
			if(status == "cds")
			{
				spos = StoI(V[3]);
				epos = StoI(V[4]);
				if(spos <= epos)
				{
					//fout << FF.size() << endl;
					for(int i = 0; i < FF.size(); i++)
					{
						string::size_type a3 = FF[i][0].find(scf);
						if(a3 != string::npos)
						{
							fas = FF[i][1];
							break;
						}
					}
					cutfas = fas.substr(spos-1,epos-spos+1);
					if(strand == "-")
					{
						tmpfas = cutfas;
						cutfas = Greverse(tmpfas);
					}
					fout << cutfas;
				}
			}
		}
	}
	fout << endl;
	//fout << cutfas << endl;
	return 0;
}
