#include <iostream>//標準データ入出力
#include <fstream>//ファイルの入出力
#include <vector>//vector(動的配列クラス)
#include <string>//string(文字列クラス)
#include <sstream>//文字列の入出力
#include <stdlib.h>//絶対値を用いるための関数
#include <algorithm>
#include <iomanip>//少数点以下表示
#include <time.h>
using namespace std;


//------------------------------------ declare function prototype
vector<string> split(const string &str, char sep);
string ItoS(int number);

//------------------------------------ main function
int main(int argc,char**argv)
{
      ifstream fin; //gff result
      ofstream fout; //gff result
      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-i"){fin.open(argv[i+1]);}
            if(ss=="-o"){fout.open(argv[i+1]);}
      }

//入力ファイルが存在しなかった時の出力
      if(argc<=2){            
            cout << endl;
            cout <<"\t"<<"ORFframe repair tool"<<"\n\n\n";
            cout <<"version 1.0" <<"\n";
            cout <<"updated 2018/09/26"<<"\t"<<"TAKAHIRO SHINODA"<<"\n\n\n";
            cout <<"---How to use---"<<"\n";
            cout << "./a.out -i in.gff -o out.gff" <<"\n";
            cout <<"----------------"<<"\n\n";
            return 1;
      }
      
      //declare variable
      string lin,strand;
      vector<string> tmp;
      vector< pair<int,string> >CDS;
      int flag=0;


      //1行目の処理
      getline(fin,lin);
      fout << lin <<"\n";
      tmp = split(lin,'\t');
      if(tmp[6] =="+"){flag=0;}else if(tmp[6] =="-"){flag=1;}
      else{cout <<"strandが記載されていません\n";return 1;}
      

      //全て
      while(getline(fin,lin)){
            tmp = split(lin,'\t');      
            if(tmp[2] == "mRNA"||tmp[2] == "gene"){      
                  //前のCDSsetの処理
                  sort(CDS.begin(),CDS.end());
                  if(flag==1){reverse(CDS.begin(),CDS.end());}
                  int orf=0;
                  //cout << "hoge\n";
                  for(int i=0;i<CDS.size();i++){
                        tmp = split(CDS[i].second,'\t');
                        fout << tmp[0] <<"\t"<<tmp[1] <<"\t"<<tmp[2] <<"\t"<<tmp[3] <<"\t"<<tmp[4] <<"\t"<<tmp[5] <<"\t"<<tmp[6] <<"\t"<<orf<<"\t"<<tmp[8] <<"\n";
                        orf =(stoi(tmp[4])-stoi(tmp[3])+1+orf)%3; //次のORFfragを保存
                  }
                  CDS.clear();
                  fout << lin <<"\n";
                  tmp = split(lin,'\t');
                  if(tmp[6] =="+"){flag=0;}else if(tmp[6] =="-"){flag=1;}
                  else{cout <<"strandが記載されていません\n";return 1;}                  
            }
            else if(tmp[2] == "CDS"){
                  //cout << "hoge\n";
                  CDS.push_back(make_pair(stoi(tmp[3]),lin));}
      }

      //最終行だけ別途処理
      sort(CDS.begin(),CDS.end());
      if(flag==1){reverse(CDS.begin(),CDS.end());}
      int orf=0;
      for(int i=0;i<CDS.size();i++){
            tmp = split(CDS[i].second,'\t');
            fout << tmp[0] <<"\t"<<tmp[1] <<"\t"<<tmp[2] <<"\t"<<tmp[3] <<"\t"<<tmp[4] <<"\t"<<tmp[5] <<"\t"<<tmp[6] <<"\t"<<orf<<"\t"<<tmp[8] <<"\n";
            orf =(stoi(tmp[4])-stoi(tmp[3])+1+orf)%3; //次のORFfragを保存
      }
      return 0;
}


//--------------------------------------------------- split
vector<string> split(const string &str, char sep)
{
      vector<string> v;
      stringstream ss(str);
      string buffer;
      while( getline(ss, buffer, sep) ) {
            v.push_back(buffer);     
      }
      return v;      
}

//---------------------------------------------------convert int to string 
string ItoS(int number)
{
      stringstream ss;
      ss << number;      
      return ss.str();     
}
