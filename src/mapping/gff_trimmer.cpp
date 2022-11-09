#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
using namespace std;
void split (const string &str, vector<string> &array){
    istringstream ss(str);
    string buffer;
    while(getline(ss, buffer, '\t')){
        array.push_back(buffer);
    }
}
void tag_trim (const string &str, const string tag, string &value){
    int start_pos = str.find(tag) + tag.length() + 1;
    int end_pos = min(str.find(";", start_pos), str.length());
    int length = end_pos - start_pos;
    value = str.substr(start_pos, length);
}
int main(int argc, char* argv[]){
    //manual
    if (argc != 4){
        cout << "---PURPOSE---" << endl;
        cout << " list上の配列に限定したgff fileを作成する" << endl;
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. .gff file" << endl;
        cout << " 2. seq name list" << endl;
        cout << " 3. output file name" << endl;
        return 0;
    }
    //data containor
    unordered_map <string, int> dict;
    stringstream ss;
    //input seq name list
    string templine;
    vector<string> temptab;
    ifstream fin1;
    fin1.open(argv[2]);
    while(getline(fin1, templine)){
        dict[templine] = 0;
    }
    fin1.close();
    //input gff file
    ifstream fin2;
    fin2.open(argv[1]);
    while(getline(fin2, templine)){
        if(templine[0] != '#'){
            split(templine, temptab);
            string seq_name;
            if(temptab[2] == "mRNA"){
                tag_trim(temptab[8], "ID", seq_name);
            }else if(temptab[2] == "exon"){
                tag_trim(temptab[8], "Parent", seq_name);
            }else if(temptab[2] == "CDS"){
                tag_trim(temptab[8], "Parent", seq_name);
            }
            auto itr = dict.find(seq_name);
            if(itr != dict.end()){
                ss << templine << '\n';
            }
            temptab.clear();
        }else{
            ss << templine << '\n';
        }
    }
    fin2.close();
    //output file
    ofstream fout;
    fout.open(argv[3]);
    fout << ss.str();
    fout.close();
    return 0;
}

    
        


