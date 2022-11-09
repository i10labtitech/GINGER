//spalnの変な出力を修正する
//開始が0になるとか, 終止コドンが抜けるとか

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
		cout << "./a.out [spaln.gff] [reference.fa] [fout.gff]" << endl;
		return 0;
	}
//ソiートしたいファイルの読み込み
	ifstream ifs(argv[1]);
	if(ifs.fail())
	{
		cout << "error:gff file not open" << endl;
		return 0;
	}
//ソートしたいファイルの読み込み
	ifstream ifs2(argv[2]);
	if(ifs.fail())
	{
		cout << "error:fasta file not open" << endl;
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

//fastaの読み込み
	string lin;
	vector<vector<string> > FF;
	FF = Fasta_V(ifs2);

//gffの読み込み
	vector<string> I;
	vector<string> I2;
	vector<string> V;
	int spos, epos, faslength;
	int checkp = 0;
	string status, scf, strand, end1, end2, end3, tmpend1, tmpend2, fas;
	while(getline(ifs,lin) )
	{
		I = Split(lin,'\t');
		if(I.size() > 7)
		{
			status = I[2];
			if(status == "gene")
			{
				if(V.size() >= 1)
				{
					//cds結果 出力
					for(int i = 0; i < V.size(); i++)
					{
						if(checkp == 1 && i == V.size() - 1)
						{
							I2 = Split(V[i], '\t');
							epos = StoI(I2[4]);
							fout << I2[0] << '\t' << I2[1] << '\t' << I2[2] << '\t' << I2[3] << '\t' << epos + 3 << '\t' << I2[5] << '\t' << I2[6] << '\t' << I2[7] << '\t' << I2[8] << endl;
						}
						else if(checkp == 2 && i == V.size() - 1)
						{
							I2 = Split(V[i], '\t');
							spos = StoI(I2[3]);
							fout << I2[0] << '\t' << I2[1] << '\t' << I2[2] << '\t' << spos - 3 << '\t' << I2[4] << '\t' << I2[5] << '\t' << I2[6] << '\t' << I2[7] << '\t' << I2[8] << endl; 
						}
						else if(checkp == 3)
						{
						}
						else
						{
							fout << V[i] << endl;
						}
					}
				}
				V.clear();
				checkp = 0;
				scf = I[0];
				spos = StoI(I[3]);
				epos = StoI(I[4]);
				strand = I[6];
				for(int i = 0; i < FF.size(); i++)
				{
					string::size_type a0 = FF[i][0].find(scf);
					if(a0 != string::npos)
					{
						fas = FF[i][1];
						break;
					}
				}
				faslength = fas.length();
				if(spos < 1)
				{
					checkp = 3;
				}
				else if(strand == "+" && epos + 3 <= faslength)
				{
					end1 = fas.substr(epos - 3, 3);
					end2 = fas.substr(epos, 3);
					if(end1 != "TAA" && end1 != "TAG" && end1 != "TGA" && (end2 == "TAA" || end2 == "TAG" || end2 == "TGA") )
					{
						fout << I[0] << '\t' << I[1] << '\t' << I[2] << '\t' << I[3] << '\t' << epos + 3 << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t' << I[8] << endl; 
						checkp = 1;
					}
				}
				else if(strand == "-" && spos - 3 >= 1)
				{
					tmpend1 = fas.substr(spos - 1, 3);
					tmpend2 = fas.substr(spos - 4, 3);
					end1 = Greverse(tmpend1);
					end2 = Greverse(tmpend2);
					if(end1 != "TAA" && end1 != "TAG" && end1 != "TGA" && (end2 == "TAA" || end2 == "TAG" || end2 == "TGA") )
					{
						fout << I[0] << '\t' << I[1] << '\t' << I[2] << '\t' << spos - 3 << '\t' << I[4] << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t' << I[8] << endl; 
						checkp = 2;
					}
				}
				if(checkp == 0)
				{
					fout << lin << endl;
				}
			}
			else if(status == "mRNA")
			{
				spos = StoI(I[3]);
				epos = StoI(I[4]);
				if(checkp == 0)
				{
					fout << lin << endl;
				}
				else if(checkp == 1)
				{
					fout << I[0] << '\t' << I[1] << '\t' << I[2] << '\t' << I[3] << '\t' << epos + 3 << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t' << I[8] << endl;
				}
				else if(checkp == 2)
				{
					fout << I[0] << '\t' << I[1] << '\t' << I[2] << '\t' << spos - 3 << '\t' << I[4] << '\t' << I[5] << '\t' << I[6] << '\t' << I[7] << '\t' << I[8] << endl; 
				}
			}
			else if(status == "cds")
			{
				V.push_back(lin);
			}
		}
	}
	//cds結果 出力
	for(int i = 0; i < V.size(); i++)
	{
		if(checkp == 1 && i == V.size() - 1)
		{
			I2 = Split(V[i], '\t');
			epos = StoI(I2[4]);
			fout << I2[0] << '\t' << I2[1] << '\t' << I2[2] << '\t' << I2[3] << '\t' << epos + 3 << '\t' << I2[5] << '\t' << I2[6] << '\t' << I2[7] << '\t' << I2[8] << endl;
		}
		else if(checkp == 2 && i == V.size() - 1)
		{
			I2 = Split(V[i], '\t');
			spos = StoI(I2[3]);
			fout << I2[0] << '\t' << I2[1] << '\t' << I2[2] << '\t' << spos - 3 << '\t' << I2[4] << '\t' << I2[5] << '\t' << I2[6] << '\t' << I2[7] << '\t' << I2[8] << endl; 
		}
		else if(checkp == 3)
		{
		}
		else
		{
			fout << V[i] << endl;
		}
	}
	return 0;
}
