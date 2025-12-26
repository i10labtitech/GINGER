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
#include <algorithm>

using namespace std;

struct cds{
  string scaffold;
  string status;
  int stpos;
  int edpos;
  string strand;
  string frame;
  string info;
  
  bool operator<(const cds &another) const
  {
    return scaffold == another.scaffold ? stpos < another.stpos : scaffold < another.scaffold;
  };

};

vector<string> split(const string &str, char sep){
  vector<string> v;
  stringstream ss(str);
  string buffer;
  while( getline(ss, buffer, sep) ) {
    v.push_back(buffer);
  }
  return v;
}

string seq_reverse(const string &frev){
  string rev;
  for (int rv=0;rv<(int)frev.size();++rv){
    if (frev[frev.length()-1 - rv] == 'A') rev += "T";
    else if (frev[frev.length()-1 - rv] == 'T') rev += "A";
    else if (frev[frev.length()-1 - rv] == 'C') rev += "G";
    else if (frev[frev.length()-1 - rv] == 'G') rev += "C";
    else rev += "N";
  }
  return(rev);
}

void to_upper(string &str){
  for(int i = 0; i < (int)str.size(); ++i){
    if(str[i] >= 0x61 && str[i] <= 0xA7){
      char up = str[i] - 0x20;
      str[i] = up;
    }
  }
}

int main(int argc, char**argv){

  ifstream fin,fin2;
  ofstream fout;
  bool frame_flag = false;
  for(int i=0;i<argc;i++){
    string ss=argv[i];
    if(ss=="-i"){fin.open(argv[i+1]);}
    if(ss=="-g"){fin2.open(argv[i+1]);}
    if(ss=="-o"){fout.open(argv[i+1]);}
    if(ss=="-f"){frame_flag = true;}
  }

  if(argc<7){
    cout << endl;
    cout <<"\t"<<"makefasta [myliftover pipeline utility]"<<"\n\n\n";
    cout <<"version 1.1" <<"\n";
    cout <<"updated 2019/11/20"<<"\t"<<"YUTA NAKAMURA"<<"\n\n\n";
    cout <<"---How to use---"<<"\n";
    cout << "./a.out -i in.gff -g genome.fa -o out.fa" << "\n";
    cout <<"----------------"<<"\n\n";
    cout <<"---Option---"<<"\n";
    cout <<"-f : consider first exon frame" << "\n";
    cout <<"----------------"<<"\n\n";
    return 1;
  }

  string str;
  unordered_map<string,string> fasta_db;
  string scaffold_id="";
  
  while(getline(fin2,str)){
    if(str.size()==0) continue;
    if(str[0]=='>'){
      scaffold_id=str.substr(1);
      fasta_db[scaffold_id] = "";
    }else{
      to_upper(str);
      fasta_db[scaffold_id] += str;
    }
  }
  
  string lin;
  vector<string> tmp, tmp2;
  vector<vector<string>> gene;
  bool flg=false;
  
  while(getline(fin,lin)){
    if(lin[0]=='#') continue; 
    if(lin.size()==0) continue;
    tmp = split(lin,'\t');
    if(tmp[0]=="UNKNOWN")continue;
    if(tmp[2] == "mRNA" && flg == true){
      gene.push_back(tmp2);
      tmp2.clear();
    }
    tmp2.push_back(lin);
    flg=true;
  }
  gene.push_back(tmp2);

  
  for(int i=0;i<(int)gene.size();i++){
    
    cds hoge_cds;
    vector<cds> CDS;
    string id = split(split(split(gene[i][0],'\t')[8],';')[0],'=')[1];
    
    for(int j=0;j<(int)gene[i].size();j++){
      tmp = split(gene[i][j],'\t');
      if(tmp[2] == "CDS"){
	hoge_cds.scaffold=tmp[0];
	hoge_cds.status=tmp[1];
	hoge_cds.stpos=stoi(tmp[3]);
	hoge_cds.edpos=stoi(tmp[4]);
	hoge_cds.strand=tmp[6];
	hoge_cds.frame=tmp[7];
	hoge_cds.info=tmp[8];
	CDS.push_back(hoge_cds);
      }
    }

    if(CDS.size() == 0) continue;
    if(fasta_db.find(CDS[0].scaffold) == fasta_db.end()) continue; 

    int f_num = 0;
    int r_num = 0;

    for(int j=0;j<(int)CDS.size();j++){
      if(CDS[j].strand == "+") f_num++;
      else if(CDS[j].strand == "-") r_num++;
    }
    string strand;
    if(f_num >= r_num) strand="+";
    else strand="-";

    sort(CDS.begin(),CDS.end());
    if(strand == "-") reverse(CDS.begin(),CDS.end());

    int frame=0;
    if(frame_flag) frame = stoi(CDS[0].frame);
    string outseq = "";
    
    fout << ">" << id << "\n";
    for(int j=0;j<(int)CDS.size();j++){
      if(CDS[j].strand == "+") outseq += fasta_db[CDS[j].scaffold].substr(CDS[j].stpos-1,CDS[j].edpos-CDS[j].stpos+1);
      else outseq += seq_reverse(fasta_db[CDS[j].scaffold].substr(CDS[j].stpos-1,CDS[j].edpos-CDS[j].stpos+1));
    }
    fout << outseq.substr(frame) << "\n";
  }
  
  return 0;
}
