#include<iostream>
#include<string>
#include<vector>
#include<fstream>

using namespace std;

struct gene{
  string contents;
  string ID;
};

int main(const int argc,const char* argv[]){

  if(argc < 3){
    cout << "[USAGE]gbgrep original.gb geneID.list" << endl;
    return 0;
  }
  
  vector<string> genelist;
  vector<gene> genbank_index;
  string templine;

  ifstream ifs_genelist(argv[2]);
  if(!ifs_genelist){
    cout << "[USAGE]gbgrep original.gb geneID.list" << endl;
    return 1;
  }

  while(!ifs_genelist.eof() && getline(ifs_genelist, templine)){
    genelist.push_back(templine);
  }

  ifstream ifs_genbank(argv[1]);
  if(!ifs_genbank){
    cout << "[USAGE]gbgrep original.gb geneID.list" << endl;
    return 1;
  }

  string hoge;
  string hoge_ID;
  gene hogegene; 
  
  while(!ifs_genbank.eof() && getline(ifs_genbank, templine)){
    if(templine=="//"){
      hogegene.contents=hoge;
      hogegene.ID=hoge_ID;
      genbank_index.push_back(hogegene);
      hoge="";
      hoge_ID="";
      hogegene.contents="";
      hogegene.ID="";
    }else{
    hoge += templine;
    hoge += "\n";
    if(templine.find("gene=")!=string::npos)hoge_ID=templine;
    } 
  }

  for(int i=0;i<(int)genelist.size();i++){
    for(int p=0;p<(int)genbank_index.size();p++){
      if(genbank_index[p].ID.find(genelist[i])!=string::npos){
	cout << genbank_index[p].contents << "//" << endl;
	break;
      }
    } 
  }
  return 0;
}
