/*
Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology

This file is part of GINGER.

GINGER is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GINGER is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with GINGER; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
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
