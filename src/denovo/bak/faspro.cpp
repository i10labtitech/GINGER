//長さが長い順から取り出す

#include <iostream>//標準データ入出力
#include <fstream>//ファイルの入出力
#include <vector>//vector(動的配列クラス)
#include <string>//string(文字列クラス)
#include <sstream>//文字列の入出力
#include <unordered_map>//unordered_map関数の導入
#include <map>//map関数の導入
#include <set>//set関数の導入
#include <stdlib.h>//絶対値を用いるための関数
#include <algorithm>    // std::copy
#include<iomanip>//少数点以下表示
#include<time.h>
using namespace std;

struct Getseq{
      string name;
      int start;
      int end;
      Getseq(string c="",int a=1,int b=0);
};
Getseq::Getseq(string c,int a,int b)
{
      name=c;
      start=a;
      end=b;
}

//bool関数宣言
bool is_ZZZZZ(string val){return (val == "ZZZZZ");}


//関数プロトタイプ宣言
vector<string> split(const string &str, char sep);
void output(vector<string>&a,vector<string>&b,int v);
void genome(vector<string>&a,vector<string>&b,ifstream &fin1);
void cutmin(vector<string>&a,vector<string>&b,int min);
void cutmax(vector<string>&a,vector<string>&b,int max);
void rmseq(vector<string>&a,vector<string>&b,string rm);
void getseq(vector<string>&a,vector<string>&b,Getseq get,int flag1);
void scalen(vector<string>&a,vector<string>&b);
void gfile(vector<string>&a,vector<string>&b,ifstream &fin2,int flag1);
void alngapfix(vector<string>&a,vector<string>&b,int gap);
void atoA(vector<string>&a,vector<string>&b);
void rfile(vector<string>&a,vector<string>&b,ifstream &fin2);
void chan(vector<string>&a,vector<string>&b,ifstream &fin2);
void lendis(vector<string>&a,vector<string>&b);
void rmN(vector<string>&a,vector<string>&b);
void rename(vector<string>&a,vector<string>&b,string v);


//split関数
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


//main関数
int main(int argc,char**argv)
{    
      ifstream fin1;
      ifstream fin2;
      int flag=0;
      int flag1=0;
      int min=0;
      int max=1000000000;
      int v=0;
      string rm;
      Getseq get;
      int gap=100;
      string line;
      string tmp="annotation";
//閾値設定

      for(int i=0;i<argc;i++){
            string ss=argv[i];
            if(ss=="-f"){
                  fin1.open(argv[i+1]);     
            }
            if(ss=="-min"){
                  min=atoi(argv[i+1]);     
            }
            if(ss=="-rmseq"){
                  line=argv[i+1];
                  if(line.find('.') == string::npos){
                        rm=(argv[i+1]);
                        flag=1;
                  }
                  else{
                        fin2.open(argv[i+1]);
                        flag=8;
                        if(!fin2){cout <<"\n"<< "Can't open the file. please check the file exist." <<"\n\n";return 1;}
                  }
            }
            if(ss=="-getseq"){
                  vector<string>a;
                  line=argv[i+1];
                  if(line.find(',') != string::npos){
                        a=split(argv[i+1],',');
                        get.name = a[0];
                        if(a.size()==3){get.start = stoi(a[1]);get.end = stoi(a[2]);}
                        flag=2;
                  }
                  else{
                        fin2.open(argv[i+1]);
                        flag=3;
                        if(!fin2){cout <<"\n"<< "Can't open the file. please check the file exist." <<"\n\n";return 1;}
                  }
            }
            if(ss=="-scalen"){
                  flag=4;
            }
            if(ss=="-aln_gap_fix"){
                  gap=atoi(argv[i+1]);
                  flag=5;
            }
            if(ss=="-max"){
                  max=atoi(argv[i+1]);     
                  flag=6;
            }
            if(ss=="-atoA"){
                  flag=7;
            }
            if(ss=="-chan"){
                  fin2.open(argv[i+1]);
                  flag=9;                  
                  if(!fin2){cout <<"\n"<< "Can't open the file. please check the file exist." <<"\n\n";return 1;}
            }
            if(ss=="-1"){
                  v=1;
            }
            
            if(ss=="-lendis"){
                  flag=10;
            }
            if(ss=="-rmN"){
                  flag=11;
            }
            if(ss=="-rename"){
                  flag=12;
                  tmp=argv[i+1];
            }

      }

//入力ファイルが存在しなかった時の出力
      if(!fin1 || argc<=2){
            
            cout << endl;
            cout <<"\t"<<"F   A   S   P   R   O"<<"\n\n\n";
            cout <<"version 1.7.0" <<"\n";
            cout <<"updated 2019/02/08"<<"\t"<<"TAKAHIRO SHINODA"<<"\n";
            cout <<"-scalenの機能を一部修正しました。"<<"\n";
            cout <<"scalenにNの本数を除いた長さを出力する機能を追加しました。"<<"\n\n\n";
            cout << "This is the tool to process fasta file." <<"\n\n\n";
            cout <<"example"<<"\n\n";
            cout <<"FASPRO -f aiueo.fa -min 500 > hoge.fa"<<"\n";
            cout <<"FASPRO -f aiueo.fa -max 500 > hoge.fa"<<"\n";
            cout <<"FASPRO -f aiueo.fa -rmseq itohlab > hoge.fa"<<"\n";
            cout <<"FASPRO -f aiueo.fa -rmseq itohlab.fa > hoge.fa"<<"\n";
            cout <<"FASPRO -f aiueo.fa -getseq itohlab,300,8000 > hoge.fa"<<"\n";
            cout <<"FASPRO -f aiueo.fa -scalen > hoge.fa"<<"\n";
            cout <<"FASPRO -f aiueo.fa -aln_gap_fix 90 > hoge.fa"<<"\n\n\n";
            cout << "-----------------------------------------------------------------------------------------" << endl;
            cout << "-f"<<"\t\t"<<" : input fastafile"<< endl;
            cout << "-min"<<"\t\t"<<" : min length of sequence"<<"\t\t\t\t"<<"[default 0]"<< endl;
            cout << "-max"<<"\t\t"<<" : max length of sequence"<<"\t\t\t\t"<<"[default 1000000000]"<< endl;
            cout << "-rmseq"<<"\t\t"<<" : remove specific sequence"<<"\t\t\t\t"<<"[default NONE]"<< endl;
            cout <<"\t\t   file format or only scaffold name"<<"\t\t\t"<< endl;
            cout << "-getseq"<<"\t\t"<<" : get specific sequence region"<<"\t\t\t\t"<< endl;
            cout <<"\t\t   name,start_position,end_position"<<"\t\t\t"<<"[default full length]"<< endl;
            cout <<"\t\t   file format or only scaffold name"<<"\t\t\t"<< endl;
            cout << "-scalen"<<"\t\t"<<" : output scaffold_name,scaffold_length"<<"\t\t\t\t"<<endl;
            cout << "-aln_gap_fix"<<"\t"<<" : input alignment result(fasta format) file"<<"\n"<<"\t\t   and output without -,x,N"<<"\t\t\t\t"<<"[default 0]"<<endl;
            cout << "-atoA"<<"\t\t"<<" : convert small letter to calital letter ra→A"<<"\t\t\t\t"<<endl;
            cout << "-1"<<"\t\t"<<" : output scaffold in 1 line"<<"\t\t\t\t"<<endl;
            cout << "-chan"<<"\t\t"<<" : convert change the scaffold name a→A"<<"\t\t\t\t"<<endl;
            cout << "-lendis"<<"\t\t"<<" : output length distributionA"<<"\t\t\t\t"<<endl;
            cout << "-rmN"<<"\t\t"<<" : remove contig including N at least 1"<<"\t\t\t\t"<<endl;
            cout << "-rename"<<"\t\t"<<" : rename scaffold name"<<"\t\t\t\t\t"<<"[default annotation]"<<endl;
            cout <<"\t\t   only file format separated by tab (before after)"<<"\t\t\t"<< endl;

            cout << "-----------------------------------------------------------------------------------------" << endl<<endl;
            cout << "ATTENTION! this tool output as standard output. please prepare the output file by yourself."<<endl<<endl<<endl;
      }

      //refgenomeの格納
      vector<string>a; //配列
      vector<string>b; //配列のID      
      
      genome(a,b,fin1);


      switch(flag){
      case 0:
            cutmin(a,b,min);
            output(a,b,v);
            break;
      case 1:
            if(rm==""){cout << "please input name of sequence you want to remove." <<endl;return 1;}
            rmseq(a,b,rm);
            output(a,b,v);
            break;
      case 2:
            if(get.name==""){cout << "please input name of sequence you want to get. [scaffold_name,start,end]" <<endl;return 1;}
            getseq(a,b,get,flag1);
            break;
      case 3:
            gfile(a,b,fin2,flag1);
            //output(a,b,v);
            break;
      case 4:
            scalen(a,b);
            break;
      case 5:
            alngapfix(a,b,gap);
            break;
      case 6:
            cutmax(a,b,max);
            output(a,b,v);
            break;
      case 7:
            atoA(a,b);
            output(a,b,v);
            break;
      case 8:
            rfile(a,b,fin2);
            output(a,b,v);
            break;
      case 9:
            chan(a,b,fin2);
            output(a,b,v);
            break;
      case 10:
            lendis(a,b);
            break;
      case 11:
            rmN(a,b);
            break;
      case 12:
            rename(a,b,tmp);
            break;
      default:
            cout << "flag don't perform correctly. please check program."<<endl;
            break;
      }
      
      return 0;
}

//output 出力
void output(vector<string>&a,vector<string>&b,int v)
{
      switch(v){
      case 0:
            for(int i=0;i<b.size();i++){
                  cout << b[i] << endl;
                  int p=0;
                  for (int j=0;j<a[i].size();j++){
                        cout << a[i][j];
                        p++;
                        if(p==80){p=0;cout << endl;}
                  }
                  cout << endl;
            }
            break;

      case 1:
            for(int i=0;i<b.size();i++){
                  cout << b[i] << endl;
                  cout << a[i] << endl;
            }
            break;

      default :
            cout << "output関数が正常に機能してません"<<endl;
            break;
            
      }
}

//genome配列切り出し関数
void genome(vector<string>&a,vector<string>&b,ifstream &fin1)
{
      string line;
      string h;
      double z=0;
      
      //最初の一行だけ先に配列に格納
      getline(fin1,line);
      b.push_back(line);  
      //残りのデータを格納していく
      while(getline(fin1,line)){
            if(line.find('>') != string::npos){
                  b.push_back(line);
                  z=1;
            }
            else{
                  if(z==1){
                        a.push_back(h); h =""; z =0; h += line;
                  }
                  else{
                        h += line;
                  }
            }
      }
      //最後に残っているhをaに格納する。
      a.push_back(h);
}

void cutmin(vector<string>&a,vector<string>&b,int min)
{
      if(a.size()!=b.size()){
            cout <<"not match the number of vector" <<endl;
            return ;
      }
      
      for(int i=0;i<a.size();i++){
            if(a[i].size()<min){
                  b[i]="ZZZZZ";                  
                  a[i]="ZZZZZ";
            }
      }
      a.erase(remove_if(a.begin(),a.end(),is_ZZZZZ),a.end());
      b.erase(remove_if(b.begin(),b.end(),is_ZZZZZ),b.end());
      //cout <<"hoge "<<endl;

}

void cutmax(vector<string>&a,vector<string>&b,int max)
{
      if(a.size()!=b.size()){
            cout <<"not match the number of vector" <<endl;
            return ;
      }
      
      for(int i=0;i<a.size();i++){
            if(a[i].size()>max){
                  b[i]="ZZZZZ";                  
                  a[i]="ZZZZZ";
            }
      }
      a.erase(remove_if(a.begin(),a.end(),is_ZZZZZ),a.end());
      b.erase(remove_if(b.begin(),b.end(),is_ZZZZZ),b.end());
      //cout <<"hoge "<<endl;

}

void rmseq(vector<string>&a,vector<string>&b,string rm)
{
      if(a.size()!=b.size()){
            cout <<"not match the number of vector" <<endl;
            return ;
      }
      
      for(int i=0;i<b.size();i++){
                  if(b[i].find(rm)!=string::npos){
                        b[i]="ZZZZZ";                  
                        a[i]="ZZZZZ";
                  }
      }
      a.erase(remove_if(a.begin(),a.end(),is_ZZZZZ),a.end());
      b.erase(remove_if(b.begin(),b.end(),is_ZZZZZ),b.end());
      //cout <<"hoge "<<endl;

}

void rfile(vector<string>&a,vector<string>&b,ifstream &fin2)
{
      if(a.size()!=b.size()){
            cout <<"not match the number of vector" <<endl;
            return ;
      }

      string line;
      while(getline(fin2,line)){     
            for(int i=0;i<b.size();i++){
                  
                  if(b[i].find(line)!=string::npos){
                        b[i]="ZZZZZ";                  
                        a[i]="ZZZZZ";
                        break;
                  }
            }
      }
      
      a.erase(remove_if(a.begin(),a.end(),is_ZZZZZ),a.end());
      b.erase(remove_if(b.begin(),b.end(),is_ZZZZZ),b.end());
      
}


void getseq(vector<string>&a,vector<string>&b,Getseq get,int flag1)
{
      if(a.size()!=b.size()){
            cout <<"not match the number of vector." <<endl;
            return ;
      }
      
      int s=0;
      int e=0;
      string seq;
      
      
      for(int i=0;i<b.size();i++){
            
            if(b[i].find(get.name)!=string::npos){
                  s=get.start-1;
                  e=get.end-1;
                  
                  if(get.end==0){
                        get.end=a[i].size();
                  }else{
                        e=(get.end)-1;
                  }
                  
                  int flag=0;
                  

                  if(s>e){
                        flag=1;
                        int tmp=s;
                        s=e;
                        e=tmp;
                  }
                  
                  seq=a[i].substr(s,e-s+1);
                  
                  if(flag==0){
                              cout << b[i] <<"\t"<<"["<<s+1<<","<<e+1<<"]"<< endl;                  
                              cout << seq << endl;
                  }
                  else if(flag==1){
                        cout << b[i] <<"\t"<<"["<<e+1<<","<<s+1<<"]"<< endl;
                        
                        for(int j=(seq.size()-1);j>=0;j--)
                        {
                              if(seq[j]=='A'){
                                    cout << 'T';
                              }else if(seq[j]=='T'){
                                    cout << 'A';
                              }else if(seq[j]=='G'){
                                    cout << 'C';
                              }else if(seq[j]=='C'){
                                    cout << 'G';
                              }else{
                                    cout << seq[j];
                              }
                        }
                        cout << endl;
                  }

                  break;
                  
            }
      }
}
void gfile(vector<string>&a,vector<string>&b,ifstream &fin2,int flag1)
{
      if(a.size()!=b.size()){
            cout <<"not match the number of vector" <<endl;
            return ;
      }
      
      string line;
      
      while(getline(fin2,line)){
            Getseq get;
            vector<string>q;
            q=split(line,',');
            get.name = q[0];
            if(q.size()==3){
                  get.start = stoi(q[1]);
                  get.end = stoi(q[2]);
            }
            int s=0;
            int e=0;
            string seq;
            int flag=0;            

            for(int i=0;i<b.size();i++){
                  if(b[i].find(get.name)!=string::npos){
                        s=get.start-1;
                        if(get.end==0){e=a[i].size()-1;
                        }else{e=(get.end)-1;
                        }if(s>e){
                              flag=1;
                              int tmp=s;
                              s=e;
                              e=tmp;
                        }

                        seq=a[i].substr(s,e-s+1);

                        if(flag1==1){     
                              if(flag==0){
                                    cout << b[i] <<"\t"<<"["<<s+1<<","<<e+1<<"]"<< endl;                  
                                    cout << seq << endl;
                              }
                              else if(flag==1){
                                    cout << b[i] <<"\t"<<"["<<e+1<<","<<s+1<<"]"<< endl;
                                    
                                    for(int j=(seq.size()-1);j>=0;j--)
                                    {
                                          if(seq[j]=='A'){
                                                cout << 'T';
                                          }else if(seq[j]=='T'){
                                                cout << 'A';
                                          }else if(seq[j]=='G'){
                                                cout << 'C';
                                          }else if(seq[j]=='C'){
                                          cout << 'G';
                                          }else{
                                                cout <<seq[j];
                                          }
                                    }
                                    cout << endl;
                              }
                        }

                        else if(flag1==0){
                              reverse(seq.begin(),seq.end());
                              if(flag==0){
                                    cout << b[i] <<"\t"<<"["<<s+1<<","<<e+1<<"]"<< endl;
                                    cout << seq << endl;
                              }
                              else if(flag==1){
                                    cout << b[i] <<"\t"<<"["<<e+1<<","<<s+1<<"]"<< endl;
                                    
                                    for(int j=(seq.size()-1);j>=0;j--)
                                    {
                                          if(seq[j]=='A'){
                                                cout << 'T';
                                          }else if(seq[j]=='T'){
                                                cout << 'A';
                                          }else if(seq[j]=='G'){
                                                cout << 'C';
                                          }else if(seq[j]=='C'){
                                                cout << 'G';
                                          }else{
                                                cout <<seq[j];
                                          }
                                    }
                                    cout << endl;
                              }
                        }
                        break;
                  }
            }
      }
}

void scalen(vector<string>&a,vector<string>&b)
{
      int N=0;
      for(int i=0;i<b.size();i++){
            for(int j=0;j<a[i].size();j++){
                  if(a[i][j]=='N'){
                        N++;
                  }     
            }
            cout << b[i].substr(1) <<"\t"<< a[i].size() << "\t" << (a[i].size() - N) << "\t" << endl;
            N=0;
      }     
}

void alngapfix(vector<string>&a,vector<string>&b,int gap)
{
      unordered_map<int,int>hash1;
      unordered_map<int,int>::iterator iterator1;
      

      for(int i=0;i<a[0].size();i++){
            hash1[i+1]=0;
      }

      for(int i=0;i<b.size();i++){
            for (int j=0;j<a[i].size();j++){
                  char x=a[i][j];                  
                  if(x!='-' && x!='X' && x!='x' && x!='n' && x!='N'){
                        hash1[j+1]++;
                  }            

            }
      }
      
      int num=(int)b.size();

      for(int i=0;i<b.size();i++){
            cout << b[i] << endl;
            int p=0;
            for (int j=0;j<a[i].size();j++){
                  if((double)hash1[j+1]/(double)num*100 >= gap){
                        cout << a[i][j];
                        p++;
                        if(p==80){p=0;cout << endl;}
                  }
            }
            cout << endl;
      }      
}

void atoA(vector<string>&a,vector<string>&b)
{
      for(int i=0;i<a.size();i++){
            for(int j=0;j<a[i].size();j++){
                  switch(a[i][j]){
                  case 'a':
                        a[i][j]='A';
                        break;
                  case 'g':
                        a[i][j]='G';
                        break;
                  case 'c':
                        a[i][j]='C';
                        break;
                  case 't':
                        a[i][j]='T';
                        break;
                  case 'n':
                        a[i][j]='N';
                        break;
                  case 'A':
                  case 'G':
                  case 'C':
                  case 'T':
                  case 'N':
                        break;
                  default:
                        cout <<"A,G,C,T以外の文字が使用されています"<<endl;
                        cout << b[i] << endl;
                        cout <<"position"<<"\t"<<j+1<<"\t"<<"使用された文字"<<"\t"<<a[i][j]<<endl;
                        break;
                  }
                  
            }      
      }
}

void chan(vector<string>&a,vector<string>&b,ifstream &fin2)
{

      if(a.size()!=b.size()){
            cout <<"not match the number of vector" <<endl;
            return ;
      }
      
      vector<string>p;
      string line;
      
      while(getline(fin2,line)){
            p=split(line,'\t');


            for(int i=0;i<b.size();i++){
                  if(b[i].find(p[0])!=string::npos){
                        //cout << b[i]<<endl;
                        b[i]=">"+p[1];                  
                        break;
                  }
            }
      }     
      
}

void lendis(vector<string>&a,vector<string>&b)
{
     unsigned int sum=0;
      
     for(int i=0;i<a.size();i++){
           sum += (unsigned int)a[i].size();
     }
     unsigned int ave=10000;
     //ave = sum/(unsigned int)a.size();
     
     unsigned int a0=0;
     unsigned int a1=0;
     unsigned int a2=0;
     unsigned int a3=0;
     unsigned int a4=0;
     unsigned int a5=0;
     unsigned int a6=0;
     unsigned int a7=0;
     unsigned int a8=0;
     unsigned int a9=0;
     unsigned int a10=0;
     unsigned int a11=0;
     unsigned int a12=0;
     unsigned int a13=0;
     unsigned int a14=0;
     unsigned int a15=0;
     unsigned int a16=0;
     unsigned int a17=0;
     unsigned int a18=0;
     unsigned int a19=0;
     unsigned int a20=0;
     unsigned int a21=0;
     unsigned int a22=0;

     unsigned int b0=0;
     unsigned int b1=0;
     unsigned int b2=0;
     unsigned int b3=0;
     unsigned int b4=0;
     unsigned int b5=0;
     unsigned int b6=0;
     unsigned int b7=0;
     unsigned int b8=0;
     unsigned int b9=0;
     unsigned int b10=0;
     unsigned int b11=0;
     unsigned int b12=0;
     unsigned int b13=0;
     unsigned int b14=0;
     unsigned int b15=0;
     unsigned int b16=0;
     unsigned int b17=0;
     unsigned int b18=0;
     unsigned int b19=0;
     unsigned int b20=0;
     unsigned int b21=0;
     unsigned int b22=0;
     


     for(int i=0;i<a.size();i++){
           int p=(int)a[i].size();
           
           if(p<500){
                 a0++;
                 b0+=p;
           }
           else if(p<1000){
                 a1++;
                 b1+=p;
           }
           else if(p<1500){
                 a2++;
                 b2+=p;
           }
           else if(p<2000){
                 a3++;
                 b3+=p;
           }          
           else if(p<2500){
                 a4++;
                 b4+=p;
           }       
           else if(p<3000){
                 a5++;
                 b5+=p;
           }       
           else if(p<3500){
                 a6++;
                 b6+=p;
           }
           else if(p<4000){
                 a7++;
                 b7+=p;
           }
           else if(p<4500){
                 a8++;
                 b8+=p;
           }
           else if(p<5000){
                 a9++;
                 b9+=p;
           }
           else if(p<5500){
                 a10++;
                 b10+=p;
           }
           else if(p<6000){
                 a11++;
                 b11+=p;
           }
           else if(p<6500){
                 a12++;
                 b12+=p;
           }
           else if(p<7000){
                 a13++;
                 b13+=p;
           }
           else if(p<7500){
                 a14++;
                 b14+=p;
           }
           else if(p<8000){
                 a15++;
                 b15+=p;
           }
           else if(p<8500){
                 a16++;
                 b16+=p;
           }
           else if(p<9000){
                 a17++;
                 b17+=p;
           }
           else if(p<9500){
                 a18++;
                 b18+=p;
           }
           else if(p<10000){
                 a19++;
                 b19+=p;
           }
           else if(p<10500){
                 a20++;
                 b20+=p;
           }
           else if(p<11000){
                 a21++;
                 b21+=p;
           }
           else {
                 a22++;
                 b22+=p;
           }
     }
     
     cout << "500" <<"\t" <<a0 <<"\t"<<b0<<endl;
     cout << "1000" <<"\t" <<a1<<"\t"<<b1<<endl;
     cout << "1500" <<"\t" <<a2<<"\t"<<b2<<endl;
     cout << "2000" <<"\t" <<a3<<"\t"<<b3<<endl;
     cout << "2500" <<"\t" <<a4<<"\t"<<b4<<endl;
     cout << "3000" <<"\t" <<a5<<"\t"<<b5<<endl;
     cout << "3500" <<"\t" <<a6<<"\t"<<b6<<endl;
     cout << "4000" <<"\t" <<a7<<"\t"<<b7<<endl;
     cout << "4500" <<"\t" <<a8<<"\t"<<b8<<endl;
     cout << "5000" <<"\t" <<a9<<"\t"<<b9<<endl;
     cout << "5500" <<"\t" <<a10<<"\t"<<b10<<endl;
     cout << "6000" <<"\t" <<a11<<"\t"<<b11<<endl;
     cout << "6500" <<"\t" <<a12<<"\t"<<b12<<endl;
     cout << "7000" <<"\t" <<a13<<"\t"<<b13<<endl;
     cout << "7500" <<"\t" <<a14<<"\t"<<b14<<endl;
     cout << "8000" <<"\t" <<a15<<"\t"<<b15<<endl;
     cout << "8500" <<"\t" <<a16<<"\t"<<b16<<endl;
     cout << "9000" <<"\t" <<a17<<"\t"<<b17<<endl;
     cout << "9500" <<"\t" <<a18<<"\t"<<b18<<endl;
     cout << "10000" <<"\t" <<a19<<"\t"<<b19<<endl;
     cout << "10500" <<"\t" <<a20<<"\t"<<b20<<endl;
     cout << "11000" <<"\t" <<a21<<"\t"<<b21<<endl;
     cout << "Inf" <<"\t" <<a22<<"\t"<<b22<<endl;
           
     
     
}

//output 出力
void rmN(vector<string>&a,vector<string>&b)
{
      vector<int> tmp;
      
      for(int i=0;i<a.size();i++){
            int flag=0;
            for (int j=0;j<a[i].size();j++){
                  if(a[i][j]=='N'){flag=1;continue;}
            }
            tmp.push_back(flag);
      }
      

      for(int i=0;i<b.size();i++){
            if(tmp[i]==0){
                  cout << b[i] << endl;
                  int p=0;
                  for (int j=0;j<a[i].size();j++){
                        cout << a[i][j];
                        p++;
                        if(p==80){p=0;cout << endl;}
                  }
                  cout << endl;
            }     
      }
}

//output 出力
void rename(vector<string>&a,vector<string>&b,string v)
{
      for(int i=0;i<b.size();i++){
            //cout <<v<<endl;
            cout << ">"<<v<<"_"<<i+1<< endl;
            int p=0;
            for (int j=0;j<a[i].size();j++){
                  cout << a[i][j];
                  p++;
                  if(p==80){p=0;cout << endl;}
            }
            cout << endl;
      }
}

