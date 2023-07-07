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

/////////////////////////////////////////////////////////////////////////
///set up
/////////////////////////////////////////////////////////////////////////


vector<string> split(const string &str, char sep);
string seq_reverse(const string &frev);
string make_proteinseq(const string &seq,unordered_map<string, string> &hash);
void to_upper(string &str);
void init_prot (unordered_map<string, string> &hash);

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

/////////////////////////////////////////////////////////////////////////
///main function
/////////////////////////////////////////////////////////////////////////

int main(int argc, char**argv){

  ifstream fin;
  ifstream fin_2;
  ofstream fout;
  for(int i=0;i<argc;i++){
    string ss=argv[i];
    if(ss=="-i"){fin.open(argv[i+1]);}
    if(ss=="-o"){fout.open(argv[i+1]);}
    if(ss=="-g"){fin_2.open(argv[i+1]);}
  }

  if(argc!=7){
    cout << endl;
    cout <<"\t"<<"make protein seq from gff tool"<<"\n\n\n";
    cout <<"version 1.1" <<"\n";
    cout <<"updated 2019/10/09"<<"\t"<<"YUTA NAKAMURA"<<"\n\n\n";
    cout <<"---How to use---"<<"\n";
    cout << "./a.out -i in.gff -g genome.fa -o out.faa" <<"\n";
    cout <<"----------------"<<"\n\n";
    return 1;
  }

  //making fasta_db
  string lin;
  unordered_map<string,string> fasta_db;
  string scaffold_id="";

  while(getline(fin_2,lin)){
    if(lin.size()==0) continue;
    if(lin[0]=='>'){
      scaffold_id=lin.substr(1);
      fasta_db[scaffold_id] = "";
    }else{
      to_upper(lin);
      fasta_db[scaffold_id] += lin;
    }
  }

  //making gff_index
  vector<string> tmp,tmp2;
  vector<vector<string>> gene;
  bool flg=false;

  while(getline(fin,lin)){
    if(lin[0]=='#' || lin.size()==0) continue; 
    tmp = split(lin,'\t');
    if(tmp[2] == "mRNA" && flg == true){
      gene.push_back(tmp2);
      tmp2.clear();
    }
    tmp2.push_back(lin);
    flg=true;
  }
  gene.push_back(tmp2);

  unordered_map<string, string> hash;
  init_prot(hash);
  
  //construct protein sequences.
  for(int i=0;i<(int)gene.size();i++){
    
    cds CDS;
    vector<cds> CDS_vec;
    string geneID;
    
    for(int j=0;j<(int)gene[i].size();j++){
      tmp = split(gene[i][j],'\t');
      if(tmp[2] == "mRNA"){
	if(split(split(tmp[8],';')[0],'=').size()==1) geneID = split(split(tmp[8],';')[0],'=')[0];
	else geneID = split(split(tmp[8],';')[0],'=')[1];
      }
      if(tmp[2] == "CDS"){
	CDS.scaffold=tmp[0];
	CDS.stpos=stoi(tmp[3]);
	CDS.edpos=stoi(tmp[4]);
	CDS.strand=tmp[6];
	CDS.frame=0;
	if(tmp[7]!=".") CDS.frame=stoi(tmp[7]);
	CDS_vec.push_back(CDS);
      }
    }
    
    if(CDS_vec.size() == 0) continue;
    
    string strand = CDS_vec[0].strand;
    if(strand != "+" && strand != "-") {cout << "[ERROR] Check your gff! Something is wrong with strand!" << "\n"; return 1;}

    string scaffold = CDS_vec[0].scaffold;
    if(fasta_db.find(scaffold) == fasta_db.end()) {cout << "[ERROR] Check your gff or fasta! Something is wrong with scaffold ID!" << "\n"; return 1;}
    
    sort(CDS_vec.begin(),CDS_vec.end());
    if(strand == "-") reverse(CDS_vec.begin(),CDS_vec.end());
    
    string seq;
    int frame = CDS_vec[0].frame;
    bool hoge = false;
    
    if(strand == "+"){
      for(int k=0;k<(int)CDS_vec.size();k++){
	if(CDS_vec[k].stpos==0 || CDS_vec[k].edpos > (int)fasta_db[scaffold].size()){
	  cerr << "[WARNING] " << geneID << ":There is a invalid exon. The transrating has skipped." << "\n";
	  hoge = true;
	  continue;
	}
	seq+=fasta_db[scaffold].substr(CDS_vec[k].stpos-1,CDS_vec[k].edpos-CDS_vec[k].stpos+1);
      }
    }else{
      for(int k=0;k<(int)CDS_vec.size();k++){
	if(CDS_vec[k].stpos==0 || CDS_vec[k].edpos > (int)fasta_db[scaffold].size()){
          cerr << "[WARNING] " << geneID << ":There is a invalid exon. The transrating has skipped." << "\n";
	  hoge = true;
	  continue;
        }
	seq+=seq_reverse(fasta_db[scaffold].substr(CDS_vec[k].stpos-1,CDS_vec[k].edpos-CDS_vec[k].stpos+1));
      }
    }

    if(hoge==true) continue;
    if((int)seq.size() - frame < 3){
      cerr << "[WARNING] " << geneID << ":This gene is too short!! The transrating has skipped." << "\n";
      hoge = true;
      continue;
    }
    
    fout << ">" << geneID << "\n";
    fout << make_proteinseq(seq.substr(frame),hash) << "\n";
    
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

void init_prot (unordered_map<string, string> &hash){
    hash["AAA"] = "K"; hash["AAC"] = "N"; hash["AAG"] = "K"; hash["AAT"] = "N";
    hash["ACA"] = "T"; hash["ACC"] = "T"; hash["ACG"] = "T"; hash["ACT"] = "T";
    hash["AGA"] = "R"; hash["AGC"] = "S"; hash["AGG"] = "R"; hash["AGT"] = "S";
    hash["ATA"] = "I"; hash["ATC"] = "I"; hash["ATG"] = "M"; hash["ATT"] = "I";

    hash["CAA"] = "Q"; hash["CAC"] = "H"; hash["CAG"] = "Q"; hash["CAT"] = "H";
    hash["CCA"] = "P"; hash["CCC"] = "P"; hash["CCG"] = "P"; hash["CCT"] = "P";
    hash["CGA"] = "R"; hash["CGC"] = "R"; hash["CGG"] = "R"; hash["CGT"] = "R";
    hash["CTA"] = "L"; hash["CTC"] = "L"; hash["CTG"] = "L"; hash["CTT"] = "L";

    hash["GAA"] = "E"; hash["GAC"] = "D"; hash["GAG"] = "E"; hash["GAT"] = "D";
    hash["GCA"] = "A"; hash["GCC"] = "A"; hash["GCG"] = "A"; hash["GCT"] = "A";
    hash["GGA"] = "G"; hash["GGC"] = "G"; hash["GGG"] = "G"; hash["GGT"] = "G";
    hash["GTA"] = "V"; hash["GTC"] = "V"; hash["GTG"] = "V"; hash["GTT"] = "V";

    hash["TAA"] = "*"; hash["TAC"] = "Y"; hash["TAG"] = "*"; hash["TAT"] = "Y";
    hash["TCA"] = "S"; hash["TCC"] = "S"; hash["TCG"] = "S"; hash["TCT"] = "S";
    hash["TGA"] = "*"; hash["TGC"] = "C"; hash["TGG"] = "W"; hash["TGT"] = "C";
    hash["TTA"] = "L"; hash["TTC"] = "F"; hash["TTG"] = "L"; hash["TTT"] = "F";

    hash["CTN"] = "L"; hash["GTN"] = "V"; hash["TCN"] = "S"; hash["CCN"] = "P";
    hash["ACN"] = "T"; hash["GCN"] = "A"; hash["CGN"] = "R"; hash["GGN"] = "G";

}

string make_proteinseq(const string &seq,unordered_map<string, string> &hash){
  int protnum = seq.size()/3;
  string proteinseq;
  for(int i=0; i < protnum; i++){
    string codon;
    codon += seq[3*i];
    codon += seq[3*i+1];
    codon += seq[3*i+2];
    if(hash[codon] == "") proteinseq += "X";
    else proteinseq += hash[codon];
  }
  return proteinseq;
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
