//ver.2 intergenicのスコアの考え方を変更
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include "function.hpp"

using namespace std;

//構造体
struct tgroup
{
  double weight;
  int spos;
  int epos;
};

//intergenic用の構造体
struct groupinter
{
  int gspos;
  int gepos;
  double score;
};

struct groupinter2
{
  int gspos;
  int gepos;
  double score1;
  double score2;
};

//weight fileの読み込みhash
unordered_map<string, double> Weight_hash(const vector<string> &file1)
{
  string lin;
  vector<string> I;
  unordered_map<string, double> hash;
  for(int i=0;i<(int)file1.size();i++){
    lin = file1[i];
    I = Split(lin, '\t');
    hash[I[0]] = StoD(I[1]);
  };
  return hash;
}

//CDS用のscore hash
//file1:grouping file
//190911:グループファイルの9列目に記載されているスコアを使用するように実装を変更する
unordered_map<string, double> Score_cds_hash(const vector<string> &file1)
{
  
  string lin;
  vector<string> I,K,KeyVec;
  string key, strand, tool, key2;
  int cdsspos, cdsepos, groupnum;
  tgroup cdsinfo;
  unordered_multimap<int, tgroup> cdshash;
  unordered_map<string, double> cdshash2;
  
  //全groupに含まれる各ツールごとのexonを探索し、そのスコアを収納するハッシュを構築します
  for(int i=0;i<(int)file1.size();i++){
    lin = file1[i];
    I = Split(lin,'\t');
    if(I[5] == "mRNA"){
      groupnum = StoI(I[0]);
      strand = I[9];
      tool = I[4];
    }else if(I[5] == "CDS"){
      key2 = I[0] + "_" + tool + "_" + I[6] + "_" + I[7] + "_" + strand + "_" + I[10];
      KeyVec.push_back(key2);
      if(cdshash2.find(key2) == cdshash2.end()) cdshash2[key2] = StoD(I[8]);
      if(cdshash2.find(key2) != cdshash2.end() && cdshash2[key2] < StoD(I[8])) cdshash2[key2] = StoD(I[8]);
    }
  }
  
  sort(KeyVec.begin(),KeyVec.end());
  KeyVec.erase(unique(KeyVec.begin(),KeyVec.end()),KeyVec.end());
  
  for(int i=0;i<(int)KeyVec.size();i++){
    I = Split(KeyVec[i],'_');
    key = I[0] + "_" + I[2] + "_" + I[3] + "_" + I[4] + "_" + I[5];
    K.push_back(key);
    cdsinfo = {cdshash2[KeyVec[i]], StoI(I[2]), StoI(I[3])};
    cdshash.insert(pair<int, tgroup>(StoI(I[0]), cdsinfo));
  }
  
  //上記で構築したハッシュを元に、各exonの実際のスコアを計算、ハッシュを構築します。
  double tmpscore, score;
  unordered_map<string, double> hash;
  for(int i = 0; i < (int)K.size();i++){
    I = Split(K[i], '_');
    groupnum = StoI(I[0]);
    cdsspos = StoI(I[1]);
    cdsepos = StoI(I[2]);
    score = 0;
    
    auto range = cdshash.equal_range(groupnum);
    for (auto iterator = range.first; iterator != range.second; iterator++){
      auto target = *iterator;
      //case1	- < - > -
      if(target.second.spos <= cdsspos && cdsepos <= target.second.epos){
	tmpscore = (double)(cdsepos - cdsspos + 1) * target.second.weight;
	score += tmpscore;
      }
      //case2	 < - >
      else if(cdsspos <= target.second.spos && target.second.epos <= cdsepos){
	tmpscore = (double)(target.second.epos - target.second.spos + 1) * target.second.weight;
	score += tmpscore;
      }
      //case3 - < - >
      else if(target.second.spos <= cdsspos && cdsspos <= target.second.epos){
	tmpscore = (double)(target.second.epos - cdsspos + 1) * target.second.weight;
	score += tmpscore;
      }
      //case4 < - > -
      else if(target.second.spos <= cdsepos && cdsepos <= target.second.epos){
	tmpscore = (double)(cdsepos - target.second.spos + 1) * target.second.weight;
	score += tmpscore;
      }
    }
    
    hash[K[i]] = score;
  }
  return hash;
}

//intron用のscore hash
//file1:weight file, file2:grouping file
unordered_map<string, double> Score_intron_hash(const vector<string> &file1){
  
  string lin, tool;
  vector<string> I,K,KeyVec;
  bool check = false;
  double weight=0,weight2=0,min_weight=0;
  string strand, frame, key, key2, intronsposs, introneposs;
  int cdsspos=0, cdsepos=0, intronspos=0, intronepos=0, groupnum=0;
  tgroup introninfo;
  unordered_multimap<int,tgroup> intronhash;
  unordered_map<string,double> intronhash2;
  
  //全groupに含まれる各ツールごとのintronを探索し、そのスコアを収納するハッシュを構築します
  for(int i=0;i<(int)file1.size();i++){
    lin = file1[i];
    I = Split(lin, '\t');
    if(I[5] == "mRNA"){
      groupnum = StoI(I[0]);
      strand = I[9];
      tool = I[4];
      check = false;
    }
    else if(I[5] == "CDS"){
      //first exon
      if(check == false){
	weight = StoD(I[8]);
	cdsspos = StoI(I[6]);
	cdsepos = StoI(I[7]);
	if(strand == "+") intronspos = cdsepos + 1;
	else intronepos = cdsspos - 1;
	check = true;
      }
      //それ以外
      else{
	weight2 = StoD(I[8]);
	cdsspos = StoI(I[6]);
	cdsepos = StoI(I[7]);
	if(strand == "+") intronepos = cdsspos - 1;
	else intronspos = cdsepos + 1;
	intronsposs = ItoS(intronspos);
	introneposs = ItoS(intronepos);
	key2 = I[0] + "_" + tool + "_" + intronsposs + "_" + introneposs + "_" + strand;
	KeyVec.push_back(key2);
	min_weight = min(weight,weight2);
	if(intronhash2.find(key2) == intronhash2.end()) intronhash2[key2] = min_weight;
	if(intronhash2.find(key2) != intronhash2.end() && intronhash2[key2] < min_weight) intronhash2[key2]  = min_weight;
	if(strand == "+") intronspos = cdsepos + 1;
	else intronepos = cdsspos - 1;
	weight = StoD(I[8]);
      }
    }
  }

  //重複を削除
  sort(KeyVec.begin(),KeyVec.end());
  KeyVec.erase(unique(KeyVec.begin(),KeyVec.end()),KeyVec.end());

  for(int i=0;i<(int)KeyVec.size();i++){
    I = Split(KeyVec[i],'_');
    key = I[0] + '_' + I[2] + '_' + I[3] + '_' + I[4];
    K.push_back(key);
    introninfo = {intronhash2[KeyVec[i]], StoI(I[2]), StoI(I[3])};
    intronhash.insert(pair<int, tgroup>(StoI(I[0]), introninfo));
  }  
  
  //上記で構築したハッシュを元に、各intronの実際のスコアを計算、ハッシュを構築します。
  double tmpscore, score;
  unordered_map<string, double> hash2;
  for(int i = 0; i < (int)K.size(); i++){
    I = Split(K[i], '_');
    groupnum = StoI(I[0]);
    intronspos = StoI(I[1]);
    intronepos = StoI(I[2]);
    score = 0;
    
    auto range = intronhash.equal_range(groupnum);
    for (auto iterator = range.first; iterator != range.second; iterator++){
      auto target = *iterator;
      if(intronspos == target.second.spos && intronepos == target.second.epos){
	tmpscore = (double)(intronepos - intronspos + 1) * target.second.weight;
	score += tmpscore;
      }
    }
    hash2[K[i]] = score;
  }
  return hash2;
}

//intergenic領域計算用のhash2	file1:group
//key:groupnum_stat	value:group内に存在するtool名
unordered_multimap<string, string> Score_intergenic_hash(const vector<string> &file1)
{
  string lin;
  unordered_multimap<string, string> hash;
  vector<string> I;
  int checkp = 0;

  //group file
  for(int i=0;i<(int)file1.size();i++){
    lin = file1[i];
    I = Split(lin, '\t');
    size_t count = hash.count(I[0]);
    //keyが存在する場合
    if(count > 0){
      checkp = 0;
      auto range = hash.equal_range(I[0]);
      for (auto iterator = range.first; iterator != range.second; iterator++){
	auto target = *iterator;
	if(target.second == I[4]){
	  checkp = 1;
	}
      }
      if(checkp == 0){
	hash.insert(pair<string, string>(I[0], I[4]));
      }
    }
    //keyが存在しない場合
    else{
      hash.insert(pair<string, string>(I[0], I[4]));
    }
  }
  return hash;
}

//intergenic score計算用のhash3
//key:groupnum_tool  value:intergenic_score
unordered_multimap<string,groupinter> Score_intergenic_hash2(const vector<string> &file1){

  string lin,key,strand;
  int spos,epos,gro_num;
  double score,score2;
  groupinter groups;
  groupinter2 groups2;
  unordered_multimap<string,groupinter> output;
  unordered_multimap<string,groupinter2> hash;
  unordered_map<int,pair<int,int>> regionmap;
  int pre_group_num=0;
  
  regionmap = get_group_region(file1);

  //group gffからグループ番号区切りの二次元配列を構築する
  vector<string> tmp,tmp2;
  vector<vector<string>> group;
  
  for(int i=0;i<(int)file1.size();i++){
    lin = file1[i];
    tmp = Split(lin,'\t');
    if(stoi(tmp[0]) != pre_group_num){
      group.push_back(tmp2);
      pre_group_num = stoi(tmp[0]);
      tmp2.clear();
    }
    tmp2.push_back(lin);
  }
  group.push_back(tmp2);
  tmp2.clear();

  //各グループに対して、存在するツールを特定し、各ツール毎にintergenic領域とそのスコアを特定していく
  for(int i=0;i<(int)group.size();i++){    

    vector<string> processed_tool;

    for(int j=0;j<(int)group[i].size();j++){
      tmp = Split(group[i][j], '\t');
      key = tmp[0] + "_" + tmp[4];
      gro_num = StoI(tmp[0]);
      if(tmp[5] == "mRNA"){ //mRNA行を発見したら次の行を参照し、そのスコアとポジションを保存
	strand = tmp[9];
	processed_tool.push_back(tmp[4]);
	spos = 0;
	epos = 0;
	double score_p=0;
	double score_m=0;
	tmp2=Split(group[i][j+1],'\t');
	if(tmp2[5]=="CDS"){
	  if(strand=="+"){
	    score_p = StoD(tmp2[8]);
	    if(tmp2.back()=="P" || tmp2.back()=="PP" || tmp2.back()=="PC") score_p = 1 - score_p;
	    spos = StoI(tmp2[6]);
	  }else{
	    score_m = StoD(tmp2[8]);
	    if(tmp2.back()=="P" || tmp2.back()=="PP" || tmp2.back()=="CP") score_m = 1 - score_m;
	    epos = StoI(tmp2[7]);
	  }
	  while(1){//mRNAの末行に移動し、そのスコアとポジションを保存
	    j++;
	    if(j==(int)group[i].size()-1) break;
	    if(Split(group[i][j+1],'\t')[5]=="mRNA")break;
	  }
	  tmp2=Split(group[i][j],'\t');
	  if(strand=="+"){
	    score_m = StoD(tmp2[8]);
	    if(tmp2.back()=="P" || tmp2.back()=="PP" || tmp2.back()=="CP") score_m = 1 - score_m;
	    epos = StoI(tmp2[7]);
	  }else{
	    score_p = StoD(tmp2[8]);
	    if(tmp2.back()=="P" || tmp2.back()=="PP" || tmp2.back()=="PC") score_p = 1 - score_p;
	    spos = StoI(tmp2[6]);
	  }
	}
	groups2 = {spos, epos, score_p, score_m};
	hash.insert(pair<string, groupinter2>(key, groups2));
      }
    }
    //tool毎にintergenic領域とそのスコアを特定する
    sort(processed_tool.begin(),processed_tool.end());
    processed_tool.erase(unique(processed_tool.begin(),processed_tool.end()),processed_tool.end());
    
    for(int j=0;j<(int)processed_tool.size();j++){
      
      vector<groupinter2> hoge;      
      key = ItoS(gro_num)+"_"+processed_tool[j];
      auto range = hash.equal_range(key);
      for(auto iterator = range.first; iterator != range.second; iterator++){
	auto target = *iterator;
	hoge.push_back(target.second);
      }
      int group_stpos = regionmap[gro_num].first;
      int group_edpos = regionmap[gro_num].second;
      vector<pair<bool,pair<double,double>>> V;
      pair<bool,pair<double,double>> p = make_pair(false,make_pair(0,0));
      for(int k=group_stpos;k<=group_edpos;k++){
	V.push_back(p);
      }
      
      for(int k=0;k<(int)hoge.size();k++){
	spos = hoge[k].gspos - group_stpos;
	epos = hoge[k].gepos - group_stpos;
	score = hoge[k].score1;
	score2 = hoge[k].score2;
	for(int l=spos;l<=epos;l++){
	  V[l].first = true;
	  if(V[l].second.first < score)V[l].second.first = score;
	  if(V[l].second.second < score2)V[l].second.second = score2;
	}
      }
      
      int k=0;
      while(k<(int)V.size()-1){
	if(V[k].first == false){
	  int spos2,epos2;
	  double scr,scr2;
	  int temp = k;
	  bool flg1,flg2;
	  flg1 = false;
	  flg2 = false;
	  while(V[temp].first == false) {
	    if(temp==0) {flg1=true;break;}
	    temp--;
	  }
	  if(flg1==false){  
	    spos2 = temp + group_stpos + 1;
	    scr = V[temp].second.second;
	  }else{
	    spos2 = temp + group_stpos;
	    scr = 0;
	  }
	  temp = k;
	  while(V[temp].first == false){
	    if(temp==(int)V.size()-1){flg2=true;break;}
	    temp++;
	  }
	  if(flg2==false){
	    epos2 = temp + group_stpos - 1;
	    scr2 = V[temp].second.first;
	  }else{
	    epos2 = temp + group_stpos;
	    scr2 = 0;
	  }
	  
	  groups.gspos = spos2;
	  groups.gepos = epos2;
	  if(flg1==true || flg2==true)groups.score = scr + scr2;
	  else groups.score = (scr + scr2)/2;
	  output.insert(pair<string, groupinter>(key, groups));
	  k=temp;
	}else{
	  k++;
	}
      }
      
    }
  }

  return output;
}

//intergenic領域の計算(group内に存在するtoolでの計算)
double Score_intergenic(int groupnum, int spos, int epos, 
			 unordered_multimap<string, string> &hash,
			 unordered_multimap<string,groupinter> &hash2
			 )
{
  string groupnums, key;
  double score = 0;
  groupnums = ItoS(groupnum);
  auto range = hash.equal_range(groupnums);
  for (auto iterator = range.first; iterator != range.second; iterator++){
    auto target = *iterator;
    key = groupnums + "_" + target.second;

    auto range = hash2.equal_range(key);
    for(auto iterator = range.first; iterator != range.second; iterator++){
      auto target = *iterator;
      if(target.second.gspos <= spos &&  epos <= target.second.gepos){
        score += target.second.score * (double)(epos - spos + 1);
      }
      else if(spos <= target.second.gspos && target.second.gepos <= epos){
	score += target.second.score * (double)(target.second.gepos - target.second.gspos + 1);	
      }
      else if(target.second.gspos <= spos && spos <= target.second.gepos){
        score += target.second.score * (double)(target.second.gepos - spos + 1);
      }
      else if(target.second.gspos <= epos && epos <= target.second.gepos){
	score += target.second.score * (double)(epos - target.second.gspos + 1);
      }
    }
    
  }
  
  return score;
}
