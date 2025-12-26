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
#include "scoring.hpp"
using namespace std;

//------------------------------------ declare function prototype
void Combination(const vector<pair<int,string>> &CDS_list,
                 const int &gro_num,				
                 stringstream &ss,				
                 unordered_map<string,string> &hash,			
                 unordered_map<string,double> &hashcds,			
                 unordered_map<string,double> &hashintron,		
                 unordered_multimap<string,string> &hashinter,  
                 unordered_multimap<string,groupinter> &hashinter2,	
		 const int &out);

void OUTPUT(const int &p,					
            const vector<long long int> &score_list,		
            const vector<vector<int>> &inter_list,		
            const vector<vector<string>> &correct_list,		
            stringstream &ss,					
            const int &gro_num,					
            int &num);

void DP(const vector<pair<int,string>> &CDS_list,			
        const int &gro_num,						
        stringstream &ss,						
        unordered_map<string,string> &hash,			
        unordered_map<string,double> &hashcds,          
        unordered_map<string,double> &hashintron,			
        unordered_multimap<string,string> &hashinter,			
        unordered_multimap<string,groupinter> &hashinter2,		
        const string &w,                            
        vector<long long int> &score_list,				
        vector< vector<string>> &correct_list,			
        vector< vector<int>> &inter_list);

void trace_back(const vector<pair<int,string>> &CDS_list,		
		unordered_map <int,pair<long long int,int>> &maxdp,	
		int t,							
		vector <string> &correct,				
		vector <int> &inter,					
		int *intdp);

//------------------------------------ main function
int main(int argc,char**argv)
{
  ifstream fin; //gff result
  ifstream fin2; //genome fasta
  ofstream fout; //output file
  int out=0;
  
  for(int i=0;i<argc;i++){
    string ss=argv[i];
    if(ss=="-fa"){fin2.open(argv[i+1]);}
    if(ss=="-gff"){fin.open(argv[i+1]);}
    if(ss=="-o"){fout.open(argv[i+1]);}
    if(ss=="-a"){out=1;}
  }
  
  if(argc<=2){            
    cout << endl;
    cout <<"\t"<<"Search algorithm"<<"\n\n\n";
    cout <<"version 2.0" <<"\n";
    cout <<"updated 2019/11/18"<<"\t"<<"YUTA NAKAMURA"<<"\n\n\n";
    cout <<"exon scoringに対応しました\n";
    cout <<"---How to use---"<<"\n";
    cout << "[Search all combination]"<<"\t"<< "./a.out -fa genome.fa -gff rna.gff -o outputfile" <<"\n\n";
    cout <<"----------------"<<"\n\n\n";
    return 1;
  }
  
  if(!fin || !fin2){
    cout << "[ERROR] Failed in opening fasta or gff file!" << "\n";
    cout << "[ERROR] Please check the files!" << "\n";
    return 1;
  }
  
  //declare variable
  vector<string> group_db;
  string lin;
  unordered_map<string,string> hash;
  unordered_map<string,double> hashcds,hashintron;
  unordered_multimap<string,string> hashinter;
  unordered_multimap<string,groupinter> hashinter2;
  stringstream ss;
  
  while(getline(fin,lin)) group_db.push_back(lin);
  fin.close();

  cout << "Start making hash table\n";
  Genome(fin2,hash);
  cout << "End making hash table\n";

  cout << "Start CDS scoring\n";
  hashcds = Score_cds_hash(group_db); 
  hashintron = Score_intron_hash(group_db); 

  make_partial_exon_flag(group_db,hash);
  hashinter = Score_intergenic_hash(group_db); 
  hashinter2 = Score_intergenic_hash2(group_db); 
  cout << "End CDS scoring\n";

  cout << "Start CDS clustering\n";
  vector<vector<string>> group = group_split(group_db);
  int size = group.size();
  int count = 1;

  cout << "detected group : " << size << "\n";
  
  for(int i=0;i<size;i++){ 
    vector<string> Vec;
    pair<int,string> tmp;
    vector<pair<int,string>> CDS_list;
    vector<pair<int,string>> CDStmp;
    int gro_num = 0;

    if(i==(int)size*count/10){
      cout <<count*10 << "% done"<< "\n";
      count++;
    }
    
    for(int j=0;j<(int)group[i].size();j++){
      lin = group[i][j];
      Vec = Split(lin,'\t');
      gro_num = StoI(Vec[0]);
      if(Vec[5] =="CDS"){
        if(Vec[9] == "+"){tmp = make_pair(stoi(Vec[6])+stoi(Vec[7])+stoi(Vec[10]),Vec[3]+"\t"+Vec[6]+"\t"+Vec[7]+"\t"+Vec[10]+"\t"+Vec[9]+"\t"+"STbase"+"\t"+"EDbase"+"\t");
        }else{tmp = make_pair(stoi(Vec[6])+stoi(Vec[7])+stoi(Vec[10])+1,Vec[3]+"\t"+Vec[6]+"\t"+Vec[7]+"\t"+Vec[10]+"\t"+Vec[9]+"\t"+"STbase"+"\t"+"EDbase"+"\t");}
        CDStmp.push_back(tmp);
      }else if(Vec[5] =="mRNA"||Vec[5] =="gene")CDS_set(CDS_list,CDStmp,hash); 
    }
    CDS_set(CDS_list,CDStmp,hash);
    CDS_filter(CDS_list); //CDS filtering
    Combination(CDS_list,gro_num,ss,hash,hashcds,hashintron,hashinter,hashinter2,out);
  }
  
  cout << "End CDS clustering\n";
  
  fout << ss.str();
  return 0;
}

//---------------------------------------------------filtering CDS which be used for search program
void Combination(const vector<pair<int,string>> &CDS_list,\
		 const int &gro_num,				\
		 stringstream &ss,			\
		 unordered_map<string,string> &hash,	\
		 unordered_map<string,double> &hashcds,		\
		 unordered_map<string,double> &hashintron,		\
		 unordered_multimap<string,string> &hashinter,		\
		 unordered_multimap<string,groupinter> &hashinter2,	\
		 const int &out)
{
  //declare variable
  vector<string> TMP;
  int num=0;
  //decalre tmporary variable
  vector<pair<int,string >> CDS_list_p,CDS_list_m;
  //
  vector<long long int> score_list;
  vector<vector<int>> inter_list;
  vector<vector<string>> correct_list;

  //
  for(int i=0;i<(int)CDS_list.size();i++){
    TMP =Split(CDS_list[i].second,'\t');
    if(TMP[4] == "+"){CDS_list_p.push_back(CDS_list[i]);
    }else if(TMP[4] == "-"){CDS_list_m.push_back(CDS_list[i]);
    }else{cout <<"strand 情報が記載されていません。\n"; return;
    }
  }
  
  // Search optimal 
  // + strand
  if(CDS_list_p.size()!=0){DP(CDS_list_p,gro_num,ss,hash,hashcds,hashintron,hashinter,hashinter2,"+", score_list,correct_list,inter_list);}
  
  // - strand
  if(CDS_list_m.size()!=0){
    reverse(CDS_list_m.begin(),CDS_list_m.end());
    DP(CDS_list_m,gro_num,ss,hash,hashcds,hashintron,hashinter,hashinter2,"-",score_list, correct_list,inter_list);
  }
  
  // output
  if(correct_list.size()!=0){
    if(out ==0){
      int q=-1;
      if(correct_list.size()==1){q=0;}
      else{if(score_list[0] > score_list[1]){q=0;}else{q=1;}}
      OUTPUT(q,score_list,inter_list,correct_list,ss,gro_num,num);
    }else if(out ==1){
      for(int p=0;p<(int)correct_list.size();p++){
	OUTPUT(p,score_list,inter_list,correct_list,ss,gro_num,num);      
      }
    }
  }

}

void OUTPUT(const int &p,				\
	    const vector<long long int> &score_list,		\
	    const vector<vector<int>> &inter_list,		\
	    const vector<vector<string>> &correct_list,	\
	    stringstream &ss,				\
	    const int &gro_num,				\
	    int &num)
{
  //
  vector <string> tmp;
  string st,ed;
  vector<int> o;
  long long int score = score_list[p];
  
  //
  o.push_back(0);
  for(int j=0;j<(int)inter_list[p].size();j++){if(inter_list[p][j] != -1){o.push_back(j);}}
  o.push_back(inter_list[p].size());

  for(int j=0;j<(int)o.size()-1;j++){
    num++;
    tmp = Split(correct_list[p][o[j+1]-1],'\t');
    if(tmp[4]=="+"){ed = tmp[2];}else{ed = tmp[1];}
    tmp = Split(correct_list[p][o[j]],'\t');
    if(tmp[4]=="+"){st = tmp[1];}else{st = tmp[2];}
    if(tmp[4]=="+"){
      ss <<tmp[0]<<"\tcandidate\tgene\t"<<st<<"\t"<<ed<<"\t"<< score <<"\t"<<tmp[4]<<"\t"<<"."<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<"\n";
      ss <<tmp[0]<<"\tcandidate\tmRNA\t"<<st<<"\t"<<ed<<"\t"<< score <<"\t"<<tmp[4]<<"\t"<<"."<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1;Parent="<<"group_num_"<<gro_num<<"_"<<num<<"\n";
    }
    else{
      ss <<tmp[0]<<"\tcandidate\tgene\t"<<ed<<"\t"<<st<<"\t"<< score <<"\t"<<tmp[4]<<"\t"<<"."<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<"\n";
      ss <<tmp[0]<<"\tcandidate\tmRNA\t"<<ed<<"\t"<<st<<"\t"<< score <<"\t"<<tmp[4]<<"\t"<<"."<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1;Parent="<<"group_num_"<<gro_num<<"_"<<num<<"\n";
    }
    int n=0;
    for(int k=o[j];k<o[j+1];k++){
      n++;
      tmp = Split(correct_list[p][k],'\t');
      ss <<tmp[0]<<"\tcandidate\texon\t"<<tmp[1]<<"\t"<<tmp[2]<<"\t.\t"<<tmp[4]<<"\t"<<tmp[3]<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1.exon"<<n<<";Parent=group_num_"<<gro_num<<"_"<<num<<".mrna1\n";
    }
    n=0;
    for(int k=o[j];k<o[j+1];k++){
      n++;
      tmp = Split(correct_list[p][k],'\t');
      ss <<tmp[0]<<"\tcandidate\tCDS\t"<<tmp[1]<<"\t"<<tmp[2]<<"\t.\t"<<tmp[4]<<"\t"<<tmp[3]<<"\t"<<"ID=group_num_"<<gro_num<<"_"<<num<<".mrna1.cds"<<n<<";Parent=group_num_"<<gro_num<<"_"<<num<<".mrna1\n";
    }
  }
}

void DP(const vector<pair<int,string>> &CDS_list,	\
	const int &gro_num,				\
	stringstream &ss,				\
	unordered_map<string,string> &hash,		\
	unordered_map<string,double> &hashcds,			\
	unordered_map<string,double> &hashintron,		\
	unordered_multimap<string,string> &hashinter,		\
	unordered_multimap<string,groupinter> &hashinter2,	\
	const string &w,					\
	vector<long long int> &score_list,			\
	vector<vector<string>> &correct_list,			\
	vector<vector<int>> &inter_list)
{
  //-----declare variable
  unordered_map <int,pair<long long int ,int > > maxdp;
  int inidp[CDS_list.size() + 1],intdp[CDS_list.size() + 1],lendp[CDS_list.size() + 1];

  //  long long int dp[CDS_list.size() + 1][CDS_list.size() + 1];
  //  long long int ldp[CDS_list.size() + 1][CDS_list.size() + 1];
  long long int** dp = new long long int*[CDS_list.size() + 1];
  long long int** ldp = new long long int*[CDS_list.size() + 1];
  for (int i = 0; i < CDS_list.size() + 1; ++i) {
    dp[i] = new long long int[CDS_list.size() + 1]; 
    ldp[i] = new long long int[CDS_list.size() + 1];
  }

  string key;
  vector<string> tmp,tmp2;
  int ter0=0,ter1=0,t=-1;
  long long int max=0;
  
  // Hash table for intergenic score
  
  // Starting point of groups
  tmp=Split(CDS_list[0].second,'\t');
  if(tmp[4]=="+"){ter0=stoi(tmp[1]);}else{ter1=stoi(tmp[2]);}
  tmp=Split(CDS_list[CDS_list.size()-1].second,'\t');
  if(tmp[4]=="+"){ter1=stoi(tmp[2]);}else{ter0=stoi(tmp[1]);}
  
  
  //----------------------Initialize DP table
  for(int i=0;i<(int)CDS_list.size();i++){

    //-----initialize table
    for(int j=0;j<i;j++){dp[j][i]=-1;ldp[j][i]=-1;}
    //-----initialize seq
    inidp[i]=-1;intdp[i]=-1;
    dp[CDS_list.size()][i]=-1;
    ldp[CDS_list.size()][i]=-1;
    
    //-----input exon score
    tmp=Split(CDS_list[i].second,'\t');
    key= ItoS(gro_num)+"_"+tmp[1]+"_"+tmp[2]+"_"+w+"_"+tmp[3];
    dp[i][i]=(long long int)hashcds[key];
    ldp[i][i]= (stoi(tmp[2])-stoi(tmp[1])+1);            
    lendp[i]=ldp[i][i];
    
    //-----Intergenic score at initial exon
    if(tmp[7]=="0" || tmp[7]=="3"){   //1 or single
      inidp[i]=1;
      if(tmp[4]=="+"){
        if(stoi(tmp[1])!=ter0){dp[CDS_list.size()][i]=(long long int)hashcds[key]+ (long long int)Score_intergenic(gro_num, ter0, stoi(tmp[1])-1,hashinter,hashinter2);}
        else{dp[CDS_list.size()][i]=(long long int)hashcds[key];}
      }else{
        if(stoi(tmp[2])!=ter1){dp[CDS_list.size()][i]=(long long int)hashcds[key]+ (long long int)Score_intergenic(gro_num, stoi(tmp[2])+1, ter1,hashinter,hashinter2);}
        else{dp[CDS_list.size()][i]=(long long int)hashcds[key];}
      }
    }
  }
  
  //single-single exon gene in group
  if(CDS_list.size()==1){t=0;max=dp[CDS_list.size()][0];}
  
  //----------------------Making DP table
  for(int i=0;i<(int)CDS_list.size();i++){
    
    long long int score=0;
    int path=-1;
    //initialization
    if(dp[CDS_list.size()][i] != -1){
      score=dp[CDS_list.size()][i];inidp[i] =1;intdp[i]=-1;
    }else{score=dp[i][i];}
    
    //Search
    for(int j=0;j<i;j++){
      int st=0,ed=0;

      tmp = Split(CDS_list[j].second,'\t');
      tmp2 = Split(CDS_list[i].second,'\t');
      if(tmp[4] == "+"){
	st = StoI(tmp[2]);
	ed = StoI(tmp2[1]);
      }else if(tmp[4] == "-"){
	st = stoi(tmp2[2]);
	ed = stoi(tmp[1]);
      }else{
	cout << "[WARNING] CDS with no strand exists in " << gro_num << " group!" << "\n";
      }
      key = ItoS(gro_num)+"_"+ ItoS(st+1)+"_"+ItoS(ed-1)+"_"+w;

      int flag = 0;
      if((tmp[7]=="2" || tmp[7]=="3") && (tmp2[7]=="0" || tmp2[7]=="3")){
	if(ed - st >= 20) flag=5;
      }
      
      if(flag!=5 && hashintron.find(key) == hashintron.end()) continue; //No intron
      if(flag!=5)flag = filter_fuc(CDS_list[j].second,CDS_list[i].second,hash,st,ed); //caluculate:1or5, not caluculate:0
      
      if(flag==1 || flag==5){
	//	key= ItoS(gro_num)+"_"+ ItoS(st+1)+"_"+ItoS(ed-1)+"_"+w;
	//if(hashintron.find(key) != hashintron.end() || flag == 5){
	  if(flag==1){
	    dp[j][i] = maxdp[j].first + dp[i][i] +(long long int)hashintron[key];
	    ldp[j][i] = lendp[j] + ldp[i][i];
	  }else if(flag==5){
	    dp[j][i] = maxdp[j].first + dp[i][i] +(long long int)Score_intergenic(gro_num,st,ed,hashinter,hashinter2);
	    ldp[j][i] = lendp[j] + ldp[i][i];
	  }else{cout <<"wrong flag.\n";}
	  
	  //when max score is updated
	  if(inidp[j]==1){
	    if(score < dp[j][i] ||(score == dp[j][i] && lendp[i] < ldp[j][i])){
	      path=j;score = dp[j][i];inidp[i] = inidp[j];lendp[i] =ldp[j][i];
	      if(flag==5){intdp[i]=1;}else{intdp[i]=-1;}
	    }}}
    }
    
    maxdp[i]=make_pair(score,path);
    
    //if an exon is terminal exon、others are set as intergenic
    tmp=Split(CDS_list[i].second,'\t');
    if(tmp[7]=="3"||tmp[7]=="2"){
      long long int predp = -1;
      if(tmp[4]=="+"){
	if(stoi(tmp[2]) != ter1){predp = maxdp[i].first +  (long long int)Score_intergenic(gro_num,stoi(tmp[2])+1,ter1,hashinter,hashinter2);
	}else{predp = maxdp[i].first;}
      }else{
	if(stoi(tmp[1]) != ter0){predp = maxdp[i].first +  (long long int)Score_intergenic(gro_num,ter0,stoi(tmp[1])-1,hashinter,hashinter2);
	}else{predp = maxdp[i].first;}
      }
      
      if(inidp[i]==1){
	if(t!=-1){
	  if(max < predp ||(max == predp && lendp[t] <  lendp[i])){
	    t=i;max = predp;
	  }
	}else{ //initial t
	  t=i;max = predp;
	}
      }
    }
    
  }

  //-----traceback
  //declare variable
  vector <string> correct;
  vector <int> inter;

  if(t!=-1){
    trace_back(CDS_list,maxdp,t,correct,inter,intdp);
    reverse(correct.begin(),correct.end());
    reverse(inter.begin(),inter.end());

    score_list.push_back(max);
    correct_list.push_back(correct);
    inter_list.push_back(inter);    
  }else{
    cout << "[WARNING] " << gro_num << " group has no path!" << "\n";
  }

  for (int i = 0; i < CDS_list.size() + 1; ++i) { // tt 0814
    delete[] dp[i];
    delete[] ldp[i];
  }  
  delete dp;
  delete ldp;
}

//trace back
void trace_back(const vector<pair<int,string>> &CDS_list,		\
		unordered_map <int,pair<long long int,int>> &maxdp,	\
		int t,							\
		vector <string> &correct,				\
		vector <int> &inter,					\
		int *intdp)
{
  vector<string> TMP;
  TMP=Split(CDS_list[t].second,'\t');
  correct.push_back(CDS_list[t].second);
  inter.push_back(intdp[t]);      
  if(maxdp[t].second != -1){trace_back(CDS_list,maxdp,maxdp[t].second,correct,inter,intdp);}
};
