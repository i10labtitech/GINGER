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
#include <map>
#include <algorithm>
#include <utility>

using namespace std;

vector<string> split(const string &str, char sep);
bool check_inframestopcodon(const string &seq);
string seq_reverse(const string &frev);
void to_upper(string &str);

struct cds{
  string scaffold;
  int stpos;
  int edpos;
  string strand;
  int frame;

  bool operator<(const cds &another) const
  {
    return scaffold == another.scaffold ? stpos < another.stpos : scaffold < another.scaffold;
  };

};

int main(int argc, char**argv){

  string del1,del2;
  ifstream fin;
  ifstream fin_2;
  ofstream fout;
  for(int i=0;i<argc;i++){
    string ss=argv[i];
    if(ss=="-i"){fin.open(argv[i+1]);}
    if(ss=="-o"){fout.open(argv[i+1]);}
    if(ss=="-g"){fin_2.open(argv[i+1]);}
    if(ss=="-d1"){del1 = argv[i+1];}
    if(ss=="-d2"){del2 = argv[i+1];}
  }

  if(argc<=10){
    cout << endl;
    cout <<"\t"<<"filter out in-frame stopcodon genes tool"<<"\n\n\n";
    cout <<"version 1.0" <<"\n";
    cout <<"updated 2019/05/14"<<"\t"<<"YUTA NAKAMURA"<<"\n\n\n";
    cout <<"---How to use---"<<"\n";
    cout << "./a.out -i in.gff -g genome.fa -d1 (mRNA|gene) -d2 (exon|CDS) -o out.gff" <<"\n";
    cout <<"----------------"<<"\n\n";
    return 1;
  }

  //making fasta_db
  string str;
  map<string,string> fasta_db;
  string scaffold_id="";

  while(getline(fin_2,str)){
    if(str[0]=='>'){
      scaffold_id=str.substr(1);
      fasta_db[scaffold_id] = "";
    }else{
      to_upper(str);
      fasta_db[scaffold_id] += str;
    }
  }

  //making gff_index
  string lin,strand;
  vector<string> tmp, tmp2;
  vector<vector<string>> gene;
  bool flg=false;

  while(getline(fin,lin)){
    if(lin[0]=='#') continue; 
    if(lin.size()==0) continue;
    tmp = split(lin,'\t');
    if(tmp[2] == del1 && flg == true){
      gene.push_back(tmp2);
      tmp2.clear();
    }
    tmp2.push_back(lin);
    flg=true;
  }
  gene.push_back(tmp2);

  //search in-frame stopcodon cds
  for(int i=0;i<(int)gene.size();i++){

    pair<int,cds> pos_cds;
    vector<pair<int,cds>> CDS;

    for(int j=0;j<(int)gene[i].size();j++){
      tmp = split(gene[i][j],'\t');
      if(tmp[2] == del2){
	pos_cds.first=stoi(tmp[3]);
	pos_cds.second.scaffold=tmp[0];
	pos_cds.second.stpos=stoi(tmp[3]);
	pos_cds.second.edpos=stoi(tmp[4]);
	pos_cds.second.strand=tmp[6];
	pos_cds.second.frame=stoi(tmp[7]);
	CDS.push_back(pos_cds);
      }
    }
    
    if(CDS.size() == 0) continue;

    sort(CDS.begin(),CDS.end());
    if(CDS[0].second.strand == "-") reverse(CDS.begin(),CDS.end());
    
    string seq;

    if(CDS[0].second.strand == "+"){
      for(int k=0;k<(int)CDS.size();k++){
	seq+=fasta_db[CDS[k].second.scaffold].substr(CDS[k].second.stpos-1,CDS[k].second.edpos-CDS[k].second.stpos+1);
      }
    }else{
      for(int k=0;k<(int)CDS.size();k++){
	seq+=seq_reverse(fasta_db[CDS[k].second.scaffold].substr(CDS[k].second.stpos-1,CDS[k].second.edpos-CDS[k].second.stpos+1));
      }
    }

    if(check_inframestopcodon(seq.substr(CDS[0].second.frame))){
      for(int l=0;l<(int)gene[i].size();l++) fout << gene[i][l] << endl;
    }
    
  }

  return 0;
}

vector<string> split(const string &str, char sep){
  vector<string> v;
  stringstream ss(str);
  string buffer;
  while( getline(ss, buffer, sep) ) {
    v.push_back(buffer);
  }
  return v;
}

bool check_inframestopcodon(const string &seq){
  int protnum = seq.size()/3;
  for(int i=0; i < protnum-1; i++){
      if(seq[3*i] == 'T' && seq[3*i+1] == 'A' && seq[3*i+2] =='A' ) return false;
      if(seq[3*i] == 'T' && seq[3*i+1] == 'A' && seq[3*i+2] =='G' ) return false;
      if(seq[3*i] == 'T' && seq[3*i+1] == 'G' && seq[3*i+2] =='A' ) return false;
  }
  return true;
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
