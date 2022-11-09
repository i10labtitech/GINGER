#include <iostream>//標準データ入出力
#include <fstream>//ファイルの入出力
#include <vector>//vector(動的配列クラス)
#include <string>//string(文字列クラス)
#include <sstream>//文字列の入出力
#include <stdlib.h>//絶対値を用いるための関数
#include <algorithm>
#include <iomanip>//少数点以下表示
#include <time.h>
#include <map>

using namespace std;


//------------------------------------ declare function prototype
vector<string> split(const string &str, char sep);

//------------------------------------ main function
int main(int argc,char**argv){

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
    cout <<"\t"<<"convert to lemon frame format tool"<<"\n\n\n";
    cout <<"version 1.0" <<"\n";
    cout <<"updated 2019/05/14"<<"\t"<<"YUTA NAKAMURA"<<"\n\n\n";
    cout <<"---How to use---"<<"\n";
    cout << "./a.out -i in.gff -o out.gff" <<"\n";
    cout <<"----------------"<<"\n\n";
    return 1;
  }
      
  //declare variable
  string lin;
  vector<string> tmp;

  while(getline(fin,lin)){
    tmp = split(lin,'\t');      
    if(tmp[7] == "1") fout << tmp[0] <<"\t"<<tmp[1] <<"\t"<<tmp[2] <<"\t"<<tmp[3] <<"\t"<<tmp[4] <<"\t"<<tmp[5] <<"\t"<<tmp[6] <<"\t"<< 2 <<"\t"<<tmp[8] <<"\n";
    else if(tmp[7] == "2") fout << tmp[0] <<"\t"<<tmp[1] <<"\t"<<tmp[2] <<"\t"<<tmp[3] <<"\t"<<tmp[4] <<"\t"<<tmp[5] <<"\t"<<tmp[6] <<"\t"<< 1 <<"\t"<<tmp[8] <<"\n";
    else fout << tmp[0] <<"\t"<<tmp[1] <<"\t"<<tmp[2] <<"\t"<<tmp[3] <<"\t"<<tmp[4] <<"\t"<<tmp[5] <<"\t"<<tmp[6] <<"\t"<<tmp[7]<<"\t"<<tmp[8] <<"\n";
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
