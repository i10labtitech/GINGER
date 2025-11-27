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
#include "phase1_function.hpp"
using namespace std;

void repair(ifstream &ifs);
void calculate(ifstream &fin2,int num,int flag);
void calculate_fix(ifstream &fin1);
void longestc(ifstream &fin1);

int main(int argc,char**argv)
{
      ifstream fin; //mapping result
      int num = 1;
      int flag=0;
      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-f"){fin.open(argv[i+1]);}
            if(ss=="-num"){num = atoi(argv[i+1]);}
            if(ss=="-mlen"){num = atoi(argv[i+1]);flag=5;}
            if(ss=="-v"){flag=1;}
            if(ss=="-r"){flag=2;}
            if(ss=="-CDSfix"){flag=3;}
            if(ss=="-l"){flag=4;}
      }

      if(argc<=2){            
            cout << endl;
            cout <<"\t"<<"Gff_filter tool"<<"\n\n\n";
            cout <<"version 1.1" <<"\n";
            cout <<"updated 2018/11/09"<<"\t"<<"TAKAHIRO SHINODA"<<"\n\n";
            cout <<"Changed to cut at min length"<<"\n";
            cout <<"longest mode was corrected" <<"\n\n";

            cout <<"---How to use---"<<"\n";
            cout << "[Repair mode]"<<"\t\t\t"<< "0. gff_editor -f rna.gff -r " <<"\n";
            cout << "[Filter mode]"<<"\t\t\t"<<"1. gff_editor -f rna.gff -num 1 -v " <<"\n";
            cout << "[Filter mode]"<<"\t\t\t"<<"1. gff_editor -f rna.gff -mlen 500" <<"\n";
            cout << "[Gene region repair mode]" <<"\t"<<"2. gff_editor -f rna.gff -CDSfix" <<"\n";
            cout << "[Longest mode]" <<"\t\t\t"<<"3. gff_editor -f rna.gff -l" <<"\n";
            cout <<"------------"<<"\n\n\n";
            return 1;
      }
      
//-----------------------------------------------------------------------------------------------------
      if(flag==2){
            cout << "START Repairing." << endl;
            repair(fin);
            cout << "END Repairing." << endl;     
      }else if(flag==3){
            cout << "START Gene region repair." << endl;
            calculate_fix(fin);
            cout << "END Gene region repair." << endl;
      }else if(flag==4){
            cout << "START Longest mode." << endl;
            longestc(fin);
            cout << "END Longest mode." << endl;
      }else{
            cout << "START filtering." << endl;
            calculate(fin,num,flag);
            cout << "END filtering." << endl;
      }
      return 0;
}

void repair(ifstream &ifs)
{
      string zero;
      int checkp = 0;
      int checkp2 = 0;
      string lin, genelin;
      vector<string> I;
      ofstream fout("repaired.gff"); 

      while(getline(ifs,lin) )
      {
            I = Split(lin,'\t');
            if(I[2] == "gene")
            {
                  genelin = lin;checkp = 1;checkp2 = 1;     
            }
            else if(I[2] == "mRNA" && checkp == 1)
            {
                  if(checkp2 == 1)
                  {
                        fout << genelin << endl;
                        checkp2 = 0;     
                  }
                  fout << lin << endl;
                  checkp = 2;     
            }
            else if(I[2] == "mRNA" && checkp == 2)
            {fout << lin << endl;}
            else if(I[2] == "exon" && checkp == 2)
            {fout << lin << endl;}
            else if(I[2] == "CDS" && checkp == 2)
            {fout << lin << endl;}
            else if(checkp == 2 && I[2] != "exon" && I[2] != "CDS" && I[2] != "mRNA")
            {checkp = 1;}     
      }

}

void calculate(ifstream &fin,int num,int flag)
{
      string lin;
      int CDSlen=0;
      vector<string> ON,CDS;
      ofstream fout1("filtered.gff"); 
      
      getline(fin,lin);
      CDS.push_back(lin);
      
      while(getline(fin,lin)){
            if(lin != ""){
                  ON = Split(lin,'\t');
                  if(ON[2] == "mRNA" ||ON[2] == "gene"){
                        
                        if(flag==0){
                              if(CDS.size() == num+1){
                                    for(int i=0;i<CDS.size();i++){
                                          fout1 << CDS[i] << "\n";
                                    }
                              }
                        }
                        else if(flag==1){
                              if(CDS.size() != num+1){
                                    for(int i=0;i<CDS.size();i++){
                                          fout1 << CDS[i] << "\n";
                                    }
                              }
                        }
                        else if(flag==5){
                              if(CDSlen >= num){
                                    for(int i=0;i<CDS.size();i++){
                                          fout1 << CDS[i] << "\n";
                                    }
                              }
                        }
                        CDSlen=0;
                        CDS.clear();
                        CDS.push_back(lin);
                        //CDSlen+= stoi(ON[4])-stoi(ON[3])+1;
                        
                  }
                  else if(ON[2] == "CDS"){
                        CDS.push_back(lin);
                        CDSlen+= stoi(ON[4])-stoi(ON[3])+1;
                  }
            }
      }

      if(flag==0){
            if(CDS.size() == num+1){
                  for(int i=0;i<CDS.size();i++){
                        fout1 << CDS[i] << "\n";
                  }
            }
      }
      else if(flag==1){
            if(CDS.size() != num+1){
                  for(int i=0;i<CDS.size();i++){
                        fout1 << CDS[i] << "\n";
                  }
            }
      }
      else if(flag==5){
            if(CDSlen >= num){
                  for(int i=0;i<CDS.size();i++){
                        fout1 << CDS[i] << "\n";
                  }
            }
      }
}


void calculate_fix(ifstream &fin1)
{
      string lin,lin3;
      int st=2147483647,ed=0;
      vector<string> ON,ON2,TMP;
      ofstream fout("Gene_region_repair.gff");
      getline(fin1,lin);
      lin3=lin;
      while(getline(fin1,lin)){
            ON = Split(lin,'\t');
            if(ON[2] == "mRNA" ||ON[2] == "gene"){
                  TMP =  Split(lin3,'\t');            
                  fout << TMP[0]<<"\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<st<<"\t"<<ed<<"\t"<<TMP[5]<<"\t"<<TMP[6]<<"\t"<<TMP[7]<<"\t"<<TMP[8]<< "\n";
                  for(int i=0;i<ON2.size();i++){
                        fout << ON2[i] << "\n";     
                  }
                  lin3 = lin;
                  st=2147483647,ed=0;
                  ON2.clear();
            }
            else if(ON[2] == "CDS"){
                  if(stoi(ON[3]) <= st){st = stoi(ON[3]);}
                  if(stoi(ON[4]) <= st){st = stoi(ON[4]);}
                  if(ed <= stoi(ON[3])){ed = stoi(ON[3]);}
                  if(ed <= stoi(ON[4])){ed = stoi(ON[4]);}
                  //cout << st <<"\t" << ed <<"\n";
                  ON2.push_back(lin);     
            }     
      }
      TMP =  Split(lin3,'\t');
      fout << TMP[0]<<"\t"<<TMP[1]<<"\t"<<TMP[2]<<"\t"<<st<<"\t"<<ed<<"\t"<<TMP[5]<<"\t"<<TMP[6]<<"\t"<<TMP[7]<<"\t"<<TMP[8]<< "\n";
      for(int i=0;i<ON2.size();i++){
            fout << ON2[i] << "\n";
      }      
}

void longestc(ifstream &fin1)
{      
      ofstream fout("longest.gff");
      vector<string>Vec,gff,gfftmp;
      string lin;
      while(getline(fin1,lin)){
            Vec = Split(lin,'\t');
            if(Vec[2]=="mRNA"||Vec[2]=="CDS"){
                  gff.push_back(lin);
            }
      }
      gfftmp=Grouping(gff);      
      longest(gfftmp,fout);
}
