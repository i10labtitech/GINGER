#include <iostream>//標準データ入出力
#include <fstream>//ファイルの入出力
#include <vector>//vector(動的配列クラス)
#include <string>//string(文字列クラス)
#include <sstream>//文字列の入出力
#include <unordered_map>//unordered_map関数の導入
#include <map>//map関数の導入
#include <set>//set関数の導入
#include <stdlib.h>//絶対値を用いるための関数
#include <algorithm>
#include <iomanip>//少数点以下表示
#include <time.h>
#include "function.hpp"
using namespace std;

//関数プロトタイプ宣言
void repair(ifstream &ifs);
void calculate(ifstream &fin2,int num,int flag);
void calculate_fix(ifstream &fin1);
void longestc(ifstream &fin1);

//main関数
int main(int argc,char**argv)
{
      ifstream fin; //mapping result
      int num = 1;
      int flag=0;
      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-f"){fin.open(argv[i+1]);}
            if(ss=="-num"){num = atoi(argv[i+1]);}
            if(ss=="-mlen"){num = atoi(argv[i+1]);flag=5;}
            if(ss=="-v"){flag=1;}
            if(ss=="-r"){flag=2;}
            if(ss=="-CDSfix"){flag=3;}
            if(ss=="-l"){flag=4;}
      }

//入力ファイルが存在しなかった時の出力
      if(argc<=2){            
            cout << endl;
            cout <<"\t"<<"Gff_filter tool"<<"\n\n\n";
            cout <<"version 1.1" <<"\n";
            cout <<"updated 2018/11/09"<<"\t"<<"TAKAHIRO SHINODA"<<"\n\n";
            cout <<"min lengthで切れるように変更しました"<<"\n";
            cout <<"longest modeを修正しました"<<"\n\n";

            cout <<"---How to use---"<<"\n";
            cout << "[Repair mode]"<<"\t\t\t"<< "0. gff_editor -f rna.gff -r " <<"\n";
            cout << "[Filter mode]"<<"\t\t\t"<<"1. gff_editor -f rna.gff -num 1 -v " <<"\n";
            cout << "[Filter mode]"<<"\t\t\t"<<"1. gff_editor -f rna.gff -mlen 500" <<"\n";
            cout << "[Gene region repair mode]" <<"\t"<<"2. gff_editor -f rna.gff -CDSfix" <<"\n";
            cout << "[Longest mode]" <<"\t\t\t"<<"3. gff_editor -f rna.gff -l" <<"\n";
            cout <<"------------"<<"\n\n\n";
            return 1;
      }
      
//-----------------------------------------------------------------------------------------------------
      if(flag==2){
            cout << "START Repairing." << endl;
            repair(fin);
            cout << "END Repairing." << endl;     
      }else if(flag==3){
            cout << "START Gene region repair." << endl;
            calculate_fix(fin);
            cout << "END Gene region repair." << endl;
      }else if(flag==4){
            cout << "START Longest mode." << endl;
            longestc(fin);
            cout << "END Longest mode." << endl;
      }else{
            cout << "START filtering." << endl;
            calculate(fin,num,flag);
            cout << "END filtering." << endl;
      }
      return 0;
}

void repair(ifstream &ifs)
{
      string zero;
      int checkp = 0;
      int checkp2 = 0;
      string lin, genelin;
      vector<string> I;
      ofstream fout("repaired.gff"); //出力ファイル      

//fileを行単位で読む
      while(getline(ifs,lin) )
      {
            I = split(lin,'\t');
            if(I[2] == "gene")
            {
                  genelin = lin;checkp = 1;checkp2 = 1;     
            }
            else if(I[2] == "mRNA" && checkp == 1)
            {
                  if(checkp2 == 1)
                  {
                        fout << genelin << endl;
                        checkp2 = 0;     
                  }
                  fout << lin << endl;
                  checkp = 2;     
            }
            else if(I[2] == "mRNA" && checkp == 2)
            {fout << lin << endl;}
            else if(I[2] == "exon" && checkp == 2)
            {fout << lin << endl;}
            else if(I[2] == "CDS" && checkp == 2)
            {fout << lin << endl;}
            else if(checkp == 2 && I[2] != "exon" && I[2] != "CDS" && I[2] != "mRNA")
            {checkp = 1;}     
      }

}

void calculate(ifstream &fin,int num,int flag)
{
      string lin;
      int CDSlen=0;
      vector<string> ON,CDS;
      ofstream fout1("filtered.gff"); //出力ファイル
      
//1行目処理
      getline(fin,lin);
      CDS.push_back(lin);
      
      while(getline(fin,lin)){
            if(lin != ""){
                  ON = split(lin,'\t');
                  if(ON[2] == "mRNA" ||ON[2] == "gene"){
                        
                        if(flag==0){
                              if(CDS.size() == num+1){
                                    for(int i=0;i<CDS.size();i++){
                                          fout1 << CDS[i] << "\n";
                                    }
                              }
                        }
                        else if(flag==1){
                              if(CDS.size() != num+1){
                                    for(int i=0;i<CDS.size();i++){
                                          fout1 << CDS[i] << "\n";
                                    }
                              }
                        }
                        else if(flag==5){
                              if(CDSlen >= num){
                                    for(int i=0;i<CDS.size();i++){
                                          fout1 << CDS[i] << "\n";
                                    }
                              }
                        }
                        CDSlen=0;
                        CDS.clear();
                        CDS.push_back(lin);
                        //CDSlen+= stoi(ON[4])-stoi(ON[3])+1;
                        
                  }
                  else if(ON[2] == "CDS"){
                        CDS.push_back(lin);
                        CDSlen+= stoi(ON[4])-stoi(ON[3])+1;
                  }
            }
      }
//最終行
      if(flag==0){
            if(CDS.size() == num+1){
                  for(int i=0;i<CDS.size();i++){
                        fout1 << CDS[i] << "\n";
                  }
            }
      }
      else if(flag==1){
            if(CDS.size() != num+1){
                  for(int i=0;i<CDS.size();i++){
                        fout1 << CDS[i] << "\n";
                  }
            }
      }
      else if(flag==5){
            if(CDSlen >= num){
                  for(int i=0;i<CDS.size();i++){
                        fout1 << CDS[i] << "\n";
                  }
            }
      }
}


//mapping結果で置換
void calculate_fix(ifstream &fin1)
{
      string lin,lin3;
      int st=2147483647,ed=0;
      vector<string> ON,ON2,TMP;
      ofstream fout("Gene_region_repair.gff");
      //mappingの出力結果をvectorに格納
      getline(fin1,lin);
      //最初のmRNA
      lin3=lin;
      while(getline(fin1,lin)){
            ON = split(lin,'\t');
            if(ON[2] == "mRNA" ||ON[2] == "gene"){
                  TMP =  split(lin3,'\t');            
                  fout << TMP[0]<<"\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<st<<"\t"<<ed<<"\t"<<TMP[5]<<"\t"<<TMP[6]<<"\t"<<TMP[7]<<"\t"<<TMP[8]<< "\n";
                  for(int i=0;i<ON2.size();i++){
                        fout << ON2[i] << "\n";     
                  }
                  lin3 = lin;
                  st=2147483647,ed=0;
                  ON2.clear();
            }
            else if(ON[2] == "CDS"){
                  if(stoi(ON[3]) <= st){st = stoi(ON[3]);}
                  if(stoi(ON[4]) <= st){st = stoi(ON[4]);}
                  if(ed <= stoi(ON[3])){ed = stoi(ON[3]);}
                  if(ed <= stoi(ON[4])){ed = stoi(ON[4]);}
                  //cout << st <<"\t" << ed <<"\n";
                  ON2.push_back(lin);     
            }     
      }
      //最後の1mRNA
      TMP =  split(lin3,'\t');
      fout << TMP[0]<<"\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<st<<"\t"<<ed<<"\t"<<TMP[5]<<"\t"<<TMP[6]<<"\t"<<TMP[7]<<"\t"<<TMP[8]<< "\n";
      for(int i=0;i<ON2.size();i++){
            fout << ON2[i] << "\n";
      }      
}

//被っている遺伝子Longestの結果にしぼりこむ
void longestc(ifstream &fin1)
{      
      ofstream fout("longest.gff");
      vector<string>Vec,gff,gfftmp;
      string lin;
      while(getline(fin1,lin)){
            Vec = split(lin,'\t');
            if(Vec[2]=="mRNA"||Vec[2]=="CDS"){
                  gff.push_back(lin);
            }
      }
      gfftmp=Grouping(gff);      
      longest(gfftmp,fout);
}
