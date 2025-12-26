#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <set>
#include <unordered_map>
using namespace std;
struct transcript{
    string name;
    int len;
};
void split (const string &str, vector<string> &array){
    istringstream ss(str);
    string buffer;
    while(getline(ss, buffer, '\t')){
        array.push_back(buffer);
    }
}
int main(int argc, char* argv[]){
    //manual
    if (argc != 2){
        cout << "---PURPOSE---" << endl;
        cout << " 1gene-1transcriptになるようseq name listを作成" << endl;
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. input file [transcript:XXX] [gene=YYY] [length]" << endl;
        return 0;
    }
    //data containor
    unordered_multimap<string, transcript> hash;
    set<string> key_set;
    //input file
    ifstream fin1;
    fin1.open(argv[1]);
    string templine;
    vector<string> temptab;
    while(getline(fin1, templine)){
        split(templine, temptab);
        transcript tempset = {temptab[0], stoi(temptab[2])};
        hash.insert(pair<string, transcript>(temptab[1], tempset));
        key_set.insert(temptab[1]);
        temptab.clear();
    }
    fin1.close();

    //longest search
    for(auto itr = key_set.begin(); itr != key_set.end(); itr++){
        auto range = hash.equal_range(*itr);
        int max_len = 0;
        string max_name;
        for(auto itr2 = range.first; itr2 != range.second; itr2++){
            auto target = *itr2;
            if( max_len < target.second.len){
                max_name = target.second.name;
                max_len = target.second.len;
            }
        }
        cout << max_name << endl;
    }
    return 0;
}

                    


            


