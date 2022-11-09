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
    if (argc != 3){
        cout << "---PURPOSE---" << endl;
        cout << " multi_orf_longest_nr.gff3のstrandに合致するようmulti.aln.gff3を変換する" << endl;
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. multi.aln.gff3" << endl;
        cout << " 2. multi_orf_longest_nr.gff3" << endl;
        return 0;
    }
    
    //data containor
    unordered_map<string, char> hash;

    //input argv[2]
    string tmpline;
    vector<string> tmptab;

    ifstream fin1;
    fin1.open(argv[2]);
    while(getline(fin1, tmpline)){
        split(tmpline, tmptab);
        if(tmptab.size() > 8){
            string tmpID = tmptab[0];
            char tmp_strand = tmptab[6][0];
            hash[tmpID] = tmp_strand;
        }
        tmptab.clear();
    }
    fin1.close();

    //input argv[1]
    ifstream fin2;
    fin2.open(argv[1]);
    while(getline(fin2, tmpline)){
        split(tmpline, tmptab);
        if(tmptab.size() > 8){
            string tmptag = tmptab[8];
            string tmpID;
            tag_trim(tmptag, "ID", tmpID);
            size_t trans_pos = tmpID.find("TRANS^");
            string realID = tmpID.substr(trans_pos + 6, realID.size() - trans_pos - 6);
            auto itr = hash.find(realID);
            if(itr != hash.end()){
                for(int i = 0; i < 6; i++){
                    cout << tmptab[i] << '\t';
                }
                cout << itr->second << '\t';
                for(int i = 7; i < tmptab.size()-1; i++){
                    cout << tmptab[i] << '\t';
                }
                cout << tmptab[tmptab.size()-1] << '\n';
            }else{
                cout << tmpline << '\n';
            }
        }
        tmptab.clear();
    }
    fin2.close();
    return 0;
}
