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

//-----------------------------------  declare struct
struct group
{
      int len;
      int st;
      int ed;
      string str;
      vector<string> cds;
      bool operator<(const group &another) const
      {            
            return len == another.len ? st > another.st : len > another.len;     
      };

};
struct group2
{
      int st,ed;
      string chr,tool;
      vector<string> gff;
      bool operator<(const group2 &another) const
      { 
            return chr == another.chr ? st < another.st : chr < another.chr;
      };
};

//------------------------------------ declare function prototype
//vector<string> split(const string &str, char sep);
//vector<string> Split(const string &s, char delim);
//string ItoS(int number);      
//int StoI(string str);
//double StoD(string str);
//int frame_convert(int frame);
void Genome(ifstream &fin1,unordered_map<string,string>&hash);
void CDS_set(vector<pair<int , string> >&CDS_list,vector<pair<int,string> >&CDStmp,unordered_map<string,string>hash);
void CDS_filter(vector<pair<int , string> >&CDS_list);
int filter_fuc(string st1,string st2,unordered_map<string,string>&hash,int &st,int &ed);
vector<string> Grouping(vector<string> gffin);
void Grouping_output(vector<group2> &Box,int &gro_num,vector<string> &gffout);
void longest(vector<string> &gff,ofstream &fout);
void longest_output(vector<group> hash,ofstream &fout,int gst,int ged);

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

vector<string> Split(const string &s, char delim)
{
        vector<string> elems;
        stringstream ss(s);
        string item;
        while (getline(ss, item, delim))
        {
                if (!item.empty())
                {
                        elems.push_back(item);
                }
        }
        return elems;
}

//------------------------- ItoS
string ItoS(int number){
            stringstream ss;
            ss << number;
            return ss.str();
}

//------------------------- StoI
int StoI(string str)
{
        int number;
        istringstream iss(str);
        iss >> number;
        return number;
}

//------------------------- StoD
double StoD(string str)
{
        double number;
        istringstream iss(str);
        iss >> number;
        return number;
}

//------------------------- frame_convert
int frame_convert(int frame){
  if(frame==0) return 0;
  else if(frame==1) return 2;
  else if(frame==2) return 1;
  else return frame;
}

//---------------------------------------------------making genome hash table
void Genome(ifstream &fin1,unordered_map<string,string>&hash)
{
      string line;
      string g,h;

      getline(fin1,line);
      g=(line.substr(1));

      while(getline(fin1,line)){
            if(line.find('>') != string::npos){
                  //register scaffold
                  hash[g]=h;                  
                  //reassign scaffold
                  g=(line.substr(1));
                  h ="";
            }else{h += line;}
      }
      hash[g]=h;
}

void CDS_set(vector<pair<int , string> >&CDS_list,vector<pair<int,string> >&CDStmp,unordered_map<string,string>hash)
{
      vector<string> Vec;      
      sort(CDStmp.begin(),CDStmp.end());
      
      for(int i=0;i<CDStmp.size();i++){      
            Vec = split(CDStmp[i].second,'\t');
	    int frame = frame_convert(stoi(Vec[3]));
	    int r_frame = (stoi(Vec[2])-stoi(Vec[1])+1-frame)%3;
            if(CDStmp.size()==1){
	      string w = hash[Vec[0]].substr(stoi(Vec[1])-1+frame,3);
	      string w2 = hash[Vec[0]].substr(stoi(Vec[2])-3-r_frame,3);
                  if(Vec[4]=="+"){
		    if(w == "ATG" && (w2 == "TAA" || w2 == "TGA" || w2 == "TAG") && (stoi(Vec[2])-stoi(Vec[1])+1-frame-r_frame)%3==0 ){
                              CDStmp[i].second+="3";
                                     }
                        else{CDStmp[i].second+="4";}
                  }                  
                  else if(Vec[4]=="-"){     
                        if(w2 == "CAT" && (w == "TTA" || w == "TCA" || w == "CTA") && (stoi(Vec[2])-stoi(Vec[1])+1-frame-r_frame)%3==0 ){
                              CDStmp[i].second+="3";
                                     }
                        else{CDStmp[i].second+="4";}
                  }
                  else{CDStmp[i].second+="4";}
            }
            
            else if(Vec[4]=="+"){
                  if(i!= (CDStmp.size()-1)){
                        if(Vec[3]=="0"){
                              string w =hash[Vec[0]].substr(stoi(Vec[1])-1,3);
                              if(w == "ATG"){CDStmp[i].second+="0";}                              
                              else{CDStmp[i].second+="1";}
                        }
                        else{CDStmp[i].second+="1";}
                  }
                  else{      
                        string w =hash[Vec[0]].substr(stoi(Vec[2])-3-r_frame,3); 
                        if(w == "TAA" || w == "TGA" || w == "TAG"){CDStmp[i].second+="2";}
                        else{CDStmp[i].second+="1";}
                  }
            }
            else if(Vec[4]=="-"){
                  if(i==0){
                        string w =hash[Vec[0]].substr(stoi(Vec[1])-1+r_frame,3);
                        if(w == "TTA" || w == "TCA" || w == "CTA"){
                              CDStmp[i].second+="2";
                        }
                        else{CDStmp[i].second+="1";}
                  }
                  else{
                        if(Vec[3]=="0"){       
                              string w =hash[Vec[0]].substr(stoi(Vec[2])-3,3);
                              if(w == "CAT"){CDStmp[i].second+="0";}
                              else{CDStmp[i].second+="1";}
                        }
                        else{CDStmp[i].second+="1";}     
                  }     
            }
            CDS_list.push_back(CDStmp[i]);     
      }
      CDStmp.clear();
}

void CDS_filter(vector<pair<int , string> >&CDS_list)
{
      sort(CDS_list.begin(),CDS_list.end());
      CDS_list.erase(unique(CDS_list.begin(), CDS_list.end()), CDS_list.end());
}

int filter_fuc(string st1, string st2,unordered_map<string,string>&hash,int &st,int &ed){

      //delcare variable
      vector<string> st1_tmp =  split(st1,'\t');
      vector<string> st2_tmp =  split(st2,'\t');      

      if(st1_tmp[7] !=  "2" && st1_tmp[7] !=  "3") {
            
            if(st1_tmp[4] == st2_tmp[4]){
                  //2. ST[k]_ed < CDS[l]_st
                  if(st1_tmp[4] == "+"){st = stoi(st1_tmp[2]);ed = stoi(st2_tmp[1]);
                  }else if(st1_tmp[4] == "-"){st = stoi(st2_tmp[2]);ed = stoi(st1_tmp[1]);}
                  else{cout << "strandが記載されていません\n";return 10;}
                  
                  if(ed - st >= 1){
                        int t=(stoi(st1_tmp[2]) - stoi(st1_tmp[1]) +1 + stoi(st1_tmp[3])) % 3;
                        if( t == stoi(st2_tmp[3])){

                              string u="";
                              if(t != 0){
                                    if(st1_tmp[4] == "+"){u= hash[st1_tmp[0]].substr(stoi(st1_tmp[2])-t,t)+hash[st2_tmp[0]].substr(stoi(st2_tmp[1])-1,(3-t));}
                                    if(st1_tmp[4] == "-"){u= hash[st2_tmp[0]].substr(stoi(st2_tmp[2])-(3-t),(3-t))+hash[st1_tmp[0]].substr(stoi(st1_tmp[1])-1,t);}
                              }
                              
                              if(u == ""||(st1_tmp[4] =="+" && u != "TAG" && u != "TGA" && u != "TAA")||(st1_tmp[4] =="-" && u != "TTA" && u != "TCA" && u != "CTA")){     
                                    return 1;
                              }else{
                                    return 2;
                              }
                        }else{
                              return 3;
                        }      
                  }else{return 3;}
            }else{return 4;}
      }
      
      else if(st2_tmp[7] ==  "0" ||st2_tmp[7] ==  "3"){
            if(st1_tmp[4] == "+"){st = stoi(st1_tmp[2]);ed = stoi(st2_tmp[1]);
            }else if(st1_tmp[4] == "-"){st = stoi(st2_tmp[2]);ed = stoi(st1_tmp[1]);}
            else{cout << "strandが記載されていません\n";return 10;}
            if(ed - st >= 20){
                  return 5;
            }
            
      }
}

void longest(vector<string> &gff,ofstream &fout)
{
      string lin,strand,gro;
      vector<string> Vec,CDS,gff_repair;
      int st=0,ed=0,gst=2147483647,ged=0,CDSlen=0;
      group grape;
      vector<group> hash;

      gff_repair=gff;
//--------------------------------------Longest Selection
      cout <<"start Longest Selection\n";      
      
      Vec=split(gff_repair[0],'\t'); //Group
      gro=Vec[0],st=stoi(Vec[3+3]),ed=stoi(Vec[4+3]);strand=Vec[6+3];
      if(gst>st){gst=st;}if(ged<ed){ged=ed;}
      CDS.push_back(gff_repair[0]);

      for(int i=1;i<gff_repair.size();i++){
            Vec = split(gff_repair[i],'\t');
            //cout << gff_repair[i] <<"\n";
            if(gro == Vec[0]){
                  if(Vec[2+3] == "mRNA"){
                        grape = {CDSlen,st,ed,strand,CDS};
                        hash.push_back(grape); 
                        CDSlen=0;CDS.clear();
                        st=stoi(Vec[3+3]),ed=stoi(Vec[4+3]);strand=Vec[6+3];
                        if(gst>st){gst=st;}if(ged<ed){ged=ed;}
                        CDS.push_back(gff_repair[i]);
                        
                  }else if(Vec[2+3] == "CDS"){
                        CDSlen+=stoi(Vec[4+3])-stoi(Vec[3+3])+1;
                        CDS.push_back(gff_repair[i]);
                  }
            }else{
                  grape = {CDSlen,st,ed,strand,CDS};
                  hash.push_back(grape); 
                  sort(hash.begin(),hash.end());
                  longest_output(hash,fout,gst,ged);
                  CDSlen=0;CDS.clear();hash.clear();gst=2147483647,ged=0;
                  gro=Vec[0];st=stoi(Vec[3+3]),ed=stoi(Vec[4+3]);strand=Vec[6+3];
                  if(gst>st){gst=st;}if(ged<ed){ged=ed;}      
                  CDS.push_back(gff_repair[i]);
            }
      }

      grape = {CDSlen,st,ed,strand,CDS};
      hash.push_back(grape); 
      sort(hash.begin(),hash.end());
      longest_output(hash,fout,gst,ged);
      cout <<"Finish Longest Selection\n";
}


//Calculate longest
void longest_output(vector<group> hash,ofstream &fout,int gst,int ged)
{
      unordered_map<int,int> banana;
      vector<string> Vec;
      int st=0,ed=0,flag=0;
      
      for(int i=gst;i<=ged;i++){
            banana[i]=0;
      }

      for(int i=0;i<hash.size();i++){
            flag=0;
            for(int j=1;j<hash[i].cds.size();j++){
                  if(flag!=1){
                        Vec=split(hash[i].cds[j],'\t');
                        st=stoi(Vec[3+3]),ed=stoi(Vec[4+3]);
                        for(int k=st;k<=ed;k++){
                              if(banana[k]!=0){
                                    flag=1;break;
                              }
                        }
                  }
            }
            if(flag!=1){ 

                  //cout << hash[i].cds[0]<<"\n";
                  Vec=split(hash[i].cds[0],'\t');
                  for(int k=3;k<Vec.size();k++){fout << Vec[k] << "\t";}
                  fout << "\n";

                  for(int j=1;j<hash[i].cds.size();j++){
                        Vec=split(hash[i].cds[j],'\t');
                        st=stoi(Vec[3+3]),ed=stoi(Vec[4+3]);
                        for(int k=st;k<=ed;k++){banana[k]=1;}
                        //cout << hash[i].cds[j]<<"\n";
                        for(int k=3;k<Vec.size();k++){fout << Vec[k] << "\t";}
                        fout << "\n";
                  }
            }
      }
}
vector<string> Grouping(vector<string> gffin)
{
      group2 gfftmp;
      vector<string> I,gff_vec,gffout;
      vector<group2> Box,Box2;
      int st=0,ed=0,gro_num=1;
      string chr,tool,strand;

      gff_vec.push_back(gffin[0]);
      I = split(gffin[0],'\t');
      tool = I[1];st = stoi(I[3]);ed = stoi(I[4]);chr = I[0];strand = I[6];

      for(int i=1;i<gffin.size();i++){
            I = split(gffin[i],'\t');
            if(I[2] == "mRNA")
            {
                  gfftmp = {st,ed,chr,tool,gff_vec};
                  if(strand=="+"){Box.push_back(gfftmp);}else{Box2.push_back(gfftmp);}                  
                  tool = I[1];st = stoi(I[3]);ed = stoi(I[4]);chr = I[0];strand = I[6];
                  gff_vec.clear();gff_vec.push_back(gffin[i]);
            }    
            else if(I[2] == "CDS"){gff_vec.push_back(gffin[i]);}
      }
      gfftmp = {st,ed,chr,tool,gff_vec};
      if(strand=="+"){Box.push_back(gfftmp);}else{Box2.push_back(gfftmp);}                  
      sort(Box.begin(),Box.end());sort(Box2.begin(),Box2.end());
      Grouping_output(Box,gro_num,gffout);Grouping_output(Box2,gro_num,gffout);
      return gffout;
}
void Grouping_output(vector<group2> &Box,int &gro_num,vector<string> &gffout)
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
                        for(int j=i-ele_num;j<i;j++){
                              for(int k=0;k<Box[j].gff.size();k++){
                                    string tmp= ItoS(gro_num)+"\t"+ItoS(ele_num)+"\t"+ItoS(map_num)+"\t"+Box[j].gff[k];
                                    gffout.push_back(tmp);
                              }
                        }
                        gro_num++;ele_num=1;map_num=0;
                        tmpst = Box[i].st;
                        tmped = Box[i].ed;  
                        if(Box[i].tool == "mappingbase"){map_num ++;}
                  }
            }
            else{
                  //cout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num<<"\n";
                  for(int j=i-ele_num;j<i;j++){
                        for(int k=0;k<Box[j].gff.size();k++){
                              string tmp= ItoS(gro_num)+"\t"+ItoS(ele_num)+"\t"+ItoS(map_num)+"\t"+Box[j].gff[k];
                              gffout.push_back(tmp);
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
      //cout <<gro_num<<"\t"<<ele_num<<"\t"<<map_num<<"\n";
      
      for(int j=Box.size()-ele_num;j<Box.size();j++){
            for(int k=0;k<Box[j].gff.size();k++){
                  string tmp= ItoS(gro_num)+"\t"+ItoS(ele_num)+"\t"+ItoS(map_num)+"\t"+Box[j].gff[k];
                  gffout.push_back(tmp);
            }
      }
}
