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
#include <unordered_map>
#include <map>
#include <set>
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include <time.h>
using namespace std;


struct group
{      
      int st,ed;
      string chr,tool;
      vector<string> gff;
      
      bool operator<(const group &another) const
      {     
            return chr == another.chr ? st < another.st : chr < another.chr;
      };
};

vector<string> split(const string &str, char sep);
string ItoS(int number);
void Grouping(ifstream &ifs, string &f_name);
void Output(ofstream &fout,vector<group> &Box,int &gro_num);


int main(int argc,char**argv)
{
      ifstream fin; //mapping result
      string output_filename = "tmp";
      int num = 1;
      int flag=0;
      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-f"){fin.open(argv[i+1]);}
            if(ss=="-o"){output_filename=argv[i+1];}
      }

      if(argc<=4){            
            cout << endl;
            cout <<"\t"<<"annotation_grouping tool"<<"\n\n\n";
            cout <<"version 1.1" <<"\n";
            cout <<"updated 2019/10/16"<<"\t"<<"FUMIYA KOBAYASHI"<<"\n";
            cout <<"updated 2018/10/04"<<"\t"<<"TAKAHIRO SHINODA"<<"\n\n\n";

            cout <<"---How to use---"<<"\n";
            cout << "[Group mode]"<<"\t\t\t"<< "0. Group -f rna.gff -o <output>" <<"\n";
            cout <<"------------"<<"\n\n\n";
            return 1;
      }
      
//-----------------------------------------------------------------------------------------------------
      if(flag==0){
            cout << "START Grouping." << endl;
            Grouping(fin, output_filename);
            cout << "END Grouping." << endl;
      }
      return 0;
}

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

string ItoS(int number)
{
      stringstream ss;
      ss << number;      
      return ss.str();     
}

void Grouping(ifstream &ifs, string &f_name)
{
       
      ofstream fout(f_name);
      string lin;
      group gff;
      vector<string> I,gff_vec;
      vector<group> Box;
      vector<group> Box2;
      int st=0,ed=0;
      string chr,tool;
      string strand;
      int gro_num=1;

      getline(ifs,lin);
      gff_vec.push_back(lin);
      I = split(lin,'\t');
      tool = I[1];st = stoi(I[3]);ed = stoi(I[4]);
      chr = I[0];strand = I[6];

      while(getline(ifs,lin)){
            I = split(lin,'\t');
            if(I[2] == "mRNA" || I[2] == "gene")
            {
                  gff = {st,ed,chr,tool,gff_vec};
                  if(strand=="+"){Box.push_back(gff);
                  }else{Box2.push_back(gff);}                  
                  tool = I[1];
                  st = stoi(I[3]);ed = stoi(I[4]);
                  chr = I[0];strand = I[6];
                  gff_vec.clear();gff_vec.push_back(lin);
            }    
            else if(I[2] == "CDS"){gff_vec.push_back(lin);}
      }
      gff = {st,ed,chr,tool,gff_vec};
      if(strand=="+"){Box.push_back(gff);
      }else{Box2.push_back(gff);}                  
      sort(Box.begin(),Box.end());
      sort(Box2.begin(),Box2.end());
      Output(fout,Box,gro_num);
      Output(fout,Box2,gro_num);
}

void Output(ofstream &fout,vector<group> &Box,int &gro_num)
{

      int ele_num=0,map_num=0,tmpst=Box[0].st,tmped=Box[0].ed;
      string tmpchr=Box[0].chr;

      for(int i=0;i<Box.size();i++){
            if(tmpchr == Box[i].chr)
            {
                  if(Box[i].st <= tmped ){
                        if(tmped < Box[i].ed ){tmped = Box[i].ed;}           
                        ele_num ++;
                        if(Box[i].tool == "mappingbase"){map_num++;}  
                  }else{
                              cout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num<<"\n";
                              
                              for(int j=i-ele_num;j<i;j++){
                                    for(int k=0;k<Box[j].gff.size();k++){
                                          fout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num<<"\t"<<Box[j].gff[k]<<"\n";
                                    }
                              }
                        gro_num++;ele_num=1;map_num=0;
                        tmpst = Box[i].st;
                        tmped = Box[i].ed;  
                        if(Box[i].tool == "mappingbase"){map_num ++;}
                  }
            }
            else{
                        cout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num<<"\n";

                        for(int j=i-ele_num;j<i;j++){
                              for(int k=0;k<Box[j].gff.size();k++){
                                    fout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num <<"\t"<<Box[j].gff[k]<<"\n";
                              }
                        }
                  gro_num++;
                  ele_num=1;
                  map_num=0;
                  tmpchr = Box[i].chr;
                  tmpst = Box[i].st;
                  tmped = Box[i].ed;  
                  if(Box[i].tool == "mappingbase"){map_num ++;}
            }            
      }
            cout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num<<"\n";
            for(int j=Box.size()-ele_num;j<Box.size();j++){
                  for(int k=0;k<Box[j].gff.size();k++){
                        fout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num <<"\t"<<Box[j].gff[k]<<"\n";
                  }
            }
}

