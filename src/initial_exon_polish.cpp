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

struct exon_profile{
  string chr,strand,tool;
  int stpos,edpos,frame;
  double score,post_score;
  
  bool operator<(const exon_profile &another) const
  {
    return chr == another.chr ? stpos < another.stpos : chr < another.chr;
  };
  
};

vector<string> split(const string &str, char sep);
void to_upper(string &str);

int main(int argc,char **argv){
  
  ifstream fin; //base_gff
  ifstream fin2; //replace_gff
  ifstream fin3; //fasta
  ofstream fout; //output_gff
  
  for(int i=0;i<argc;i++){
    string ss=argv[i];
    if(ss=="-i") fin.open(argv[i+1]);
    if(ss=="-c") fin2.open(argv[i+1]);
    if(ss=="-f") fin3.open(argv[i+1]);
    if(ss=="-o") fout.open(argv[i+1]);
  }
  
  if(argc<=2){
    cout << endl;
    cout <<"\t"<<"initial exon polish"<<"\n\n\n";
    cout <<"version 1.1" <<"\n";
    cout <<"updated 2019/11/27"<<"\t"<<"YUTA NAKAMURA"<<"\n\n\n";
    cout <<"---How to use---"<<"\n";
    cout << "[Usage]"<<"\t"<< "./a.out -i base.gff -c replace.gff -f genome.fa -o outputfile" <<"\n";
    cout <<"----------------"<<"\n\n\n";
    return 1;
  }

  vector<string> stvec,stvec2;
  bool flg = false;
  string lin;
  
  //////////////////////////////////////////////////////
  ///construct fasta db
  //////////////////////////////////////////////////////
  
  //making fasta_db
  unordered_map<string,string> hash;
  string scaffold_id="";
  
  while(getline(fin3,lin)){
    if(lin.size()==0) continue;
    if(lin[0]=='>'){
      scaffold_id=lin.substr(1);
      hash[scaffold_id] = "";
    }else{
      to_upper(lin);
      hash[scaffold_id] += lin;
    }
  }
  
  //////////////////////////////////////////////////////
  ///construct initial exon map
  //////////////////////////////////////////////////////
  
  vector<vector<string>> gff_db;
  
  while(getline(fin2,lin)){
    if(lin[0]=='#' || lin.size()==0) continue;
    stvec = split(lin,'\t');
    if(stvec[1] == "mappingbase" || stvec[1] == "denovobase") continue;
    if(stvec[2] == "mRNA" && flg == true){
      gff_db.push_back(stvec2);
      stvec2.clear();
    }
    stvec2.push_back(lin);
    flg=true;
  }
  gff_db.push_back(stvec2);
    
  unordered_multimap<string,exon_profile> initial_exon_map; //key:scaffoldID+_+strand
  vector<string> key_vec;
  exon_profile tmp_exon_profile;
  
  for(int i=0;i<(int)gff_db.size();i++){
    if(gff_db[i].size()==1) continue;
    stvec = split(gff_db[i][1],'\t'); //initial exon
    string key;
    if(stvec[6]=="+"){
      string codon = hash[stvec[0]].substr(stoi(stvec[3])-1+stoi(stvec[7]),3);
      if(codon != "ATG") continue;
    }
    else if(stvec[6]=="-"){
      string codon = hash[stvec[0]].substr(stoi(stvec[4])-3-stoi(stvec[7]),3);
      if(codon != "CAT") continue;
    }else{
      cout << "[WARNING] Something is wrong with strand.Please check your gff file." << "\n";
      continue;
    }
    tmp_exon_profile.chr = stvec[0];
    tmp_exon_profile.tool = stvec[1];
    tmp_exon_profile.stpos = stoi(stvec[3]);
    tmp_exon_profile.edpos = stoi(stvec[4]);
    tmp_exon_profile.score = stod(stvec[5]);
    tmp_exon_profile.post_score = stod(stvec[5]);
    tmp_exon_profile.strand = stvec[6];
    tmp_exon_profile.frame = stoi(stvec[7]);
    key = stvec[0]+"_"+stvec[6];
    key_vec.push_back(key);//fot detection of exact match exon
    initial_exon_map.insert(pair<string,exon_profile>(key,tmp_exon_profile));
  }

  sort(key_vec.begin(),key_vec.end());
  key_vec.erase(unique(key_vec.begin(),key_vec.end()),key_vec.end());
  
  //////////////////////////////////////////////////////
  ///detect exact match exon
  //////////////////////////////////////////////////////

  for(int i=0,lim=key_vec.size();i<lim;i++){
    auto range = initial_exon_map.equal_range(key_vec[i]);
    for(auto itr = range.first; itr != range.second; itr++){
      auto target = *itr;
      double aug = 0;
      double snap = 0;
      double hom = 0;
      for(auto itr2 = itr; itr2 != range.second; itr2++){
	if(itr==itr2) continue;
	auto target2 = *itr2;
	if(target.second.stpos == target2.second.stpos && \
	   target.second.edpos == target2.second.edpos && \
	   target.second.frame == target2.second.frame && \
	   target.second.tool != target2.second.tool){
	  if(target2.second.tool=="AUGUSTUS" && aug < itr2->second.score) aug = itr2->second.score;
	  if(target2.second.tool=="SNAP" && snap < itr2->second.score) snap = itr2->second.score;
	  if(target2.second.tool=="homology" && hom < itr2->second.score) hom = itr2->second.score;
	}
      }
      itr->second.post_score += (aug + snap + hom);
    }
  }
  
  //////////////////////////////////////////////////////
  ///replace input gff's initial exons
  //////////////////////////////////////////////////////
  
  vector<vector<exon_profile>> exon_db;
  vector<exon_profile> exvec;
  flg = false;  
  
  while(getline(fin,lin)){
    if(lin[0]=='#' || lin.size()==0) continue;
    stvec = split(lin,'\t');
    if(stvec[2] != "mRNA" && stvec[2] != "CDS") continue;
    if(stvec[2] == "mRNA" && flg == true){
      exon_db.push_back(exvec);
      exvec.clear();
    }
    if(stvec[2] == "CDS"){
      tmp_exon_profile.chr = stvec[0];
      tmp_exon_profile.tool = stvec[1];
      tmp_exon_profile.stpos = stoi(stvec[3]);
      tmp_exon_profile.edpos = stoi(stvec[4]);
      tmp_exon_profile.score = 0;
      tmp_exon_profile.strand = stvec[6];
      tmp_exon_profile.frame = stoi(stvec[7]);
      exvec.push_back(tmp_exon_profile);
    }
    flg=true;
  }
  exon_db.push_back(exvec);

  for(int i=0;i<(int)exon_db.size();i++){
    if(exon_db[i].size()==0) continue;
    if(exon_db[i][0].strand=="+"){
      if(exon_db[i][0].stpos-1+exon_db[i][0].frame+3>(int)hash[exon_db[i][0].chr].size()) continue; //for avoiding segmental fault
      string codon = hash[exon_db[i][0].chr].substr(exon_db[i][0].stpos-1+exon_db[i][0].frame,3);
      if(codon == "ATG") continue;
    }
    else if(exon_db[i][0].strand=="-"){
      if(exon_db[i][0].edpos-3-exon_db[i][0].frame<0) continue; //for avoiding segmental fault
      string codon = hash[exon_db[i][0].chr].substr(exon_db[i][0].edpos-3-exon_db[i][0].frame,3);
      if(codon == "CAT") continue;
    }else{
      cout << "[WARNING] Something is wrong with strand.Please check your gff file." << "\n";
      continue;
    }
    string key = exon_db[i][0].chr+"_"+exon_db[i][0].strand;
    auto range = initial_exon_map.equal_range(key);
    for(auto iterator = range.first; iterator != range.second; iterator++){
      auto target = *iterator;
      if(
	 ((exon_db[i][0].strand=="+" && target.second.edpos == exon_db[i][0].edpos) || (exon_db[i][0].strand=="-" && target.second.stpos == exon_db[i][0].stpos)) && \
	 (target.second.edpos-target.second.stpos+1-target.second.frame)%3 == (exon_db[i][0].edpos-exon_db[i][0].stpos+1-exon_db[i][0].frame)%3
	 ) //filtering statement
	{
	  if(target.second.post_score < exon_db[i][0].score) continue;
	  exon_db[i][0].chr = target.second.chr;
	  exon_db[i][0].tool = target.second.tool;
	  exon_db[i][0].stpos = target.second.stpos;
	  exon_db[i][0].edpos = target.second.edpos;
	  exon_db[i][0].score = target.second.post_score;
	  exon_db[i][0].strand = target.second.strand;
	  exon_db[i][0].frame = target.second.frame;
	}
    }
  }

  //////////////////////////////////////////////////////
  ///Output
  //////////////////////////////////////////////////////
  
  for(int i=0;i<(int)exon_db.size();i++){
    sort(exon_db[i].begin(),exon_db[i].end());
    if(exon_db[i].front().strand=="+"){
      int start,end;
      start = exon_db[i].front().stpos;
      end = exon_db[i].back().edpos;
      fout << exon_db[i].front().chr << "\t" << "REPLACE" << "\t" << "mRNA" << "\t" <<  start << "\t" << end << "\t" << "." << "\t" << "+" << "\t" << "." << "\t" << "ID=mRNA_" << i+1 << "\n";
      for(int j=0;j<(int)exon_db[i].size();j++) fout << exon_db[i].front().chr << "\t" << "REPLACE" << "\t" << "CDS" << "\t" << exon_db[i][j].stpos << "\t" << exon_db[i][j].edpos << "\t" << "." << "\t" << "+" << "\t" << exon_db[i][j].frame << "\t" << "ID=mRNA_" << i+1 << ".cds" << j+1 << ";Parent=mRNA_" << i+1 << "\n";
    }
    if(exon_db[i].front().strand=="-"){
      reverse(exon_db[i].begin(),exon_db[i].end());
      int start,end;
      start = exon_db[i].back().stpos;
      end = exon_db[i].front().edpos;
      fout << exon_db[i].front().chr << "\t" << "REPLACE" << "\t" << "mRNA" << "\t" <<  start << "\t" << end << "\t" << "." << "\t" << "-" << "\t" << "." << "\t" << "ID=mRNA_" << i+1 << "\n";
      for(int j=0;j<(int)exon_db[i].size();j++) fout << exon_db[i].front().chr << "\t" << "REPLACE" << "\t" << "CDS" << "\t" << exon_db[i][j].stpos << "\t" << exon_db[i][j].edpos << "\t" << exon_db[i][j].score << "\t" << "-" << "\t" << exon_db[i][j].frame << "\t" << "ID=mRNA_" << i+1 << ".cds" << j+1 << ";Parent=mRNA_" << i+1 << "\n";
    }
  }

  return 0;
}

/////////////////////////////////////////////////////////////////////////
///sub functions
/////////////////////////////////////////////////////////////////////////

vector<string> split(const string &str, char sep){
  vector<string> v;
  stringstream ss(str);
  string buffer;
  while( getline(ss, buffer, sep) ) {
    v.push_back(buffer);
  }
  return v;
}

void to_upper(string &str){
  for(int i = 0; i < (int)str.size(); ++i){
    if(str[i] >= 0x61 && str[i] <= 0xA7){
      char up = str[i] - 0x20;
      str[i] = up;
    }
  }
} 
