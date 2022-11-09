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
void id_trim (const string &str, string &ans){
    int start_pos = str.find("transcript_id") + 15;
    int len = str.find("\"", start_pos) - start_pos;
    ans = str.substr(start_pos, len);
}
int main(int argc, char* argv[]){
    //manual
    if (argc != 5){
        cout << "---PURPOSE---" << endl;
        cout << " stringtie結果の.gtfファイルをexon数でfilteringする" << endl;
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. .gtf file" << endl;
        cout << " 2. output file name" << endl;
        cout << " 3. min exon# (integer)" << endl;
        cout << " 4. max exon# (integer)" << endl;
        return 0;
    }
    //data containor
    vector< vector<string> > gtf_data;
    int min = stoi(argv[3]);
    int max = stoi(argv[4]);
    stringstream ss;
    //input file
    string templine;
    vector<string> temptab;
    ifstream fin1;
    fin1.open(argv[1]);
    while (getline(fin1, templine)){
        if(templine[0] == '#'){
            ss << templine << '\n';
        }else{
            split(templine, temptab);
            gtf_data.push_back(temptab);
            temptab.clear();
        }
    }
    fin1.close();
    //make transcript hash
    unordered_map<string, int> t_hash;
    for(int i = 0; i < gtf_data.size(); i++){
        if( gtf_data[i][2] == "transcript"){
            string id;
            id_trim(gtf_data[i][8], id);
            t_hash[id] = 0;
        }
    }
    //count exon number
    for(int i = 0; i < gtf_data.size(); i++){
        if( gtf_data[i][2] == "exon"){
            string id;
            id_trim(gtf_data[i][8], id);
            auto itr = t_hash.find(id);
            if(itr != t_hash.end()){
                itr->second += 1;
            }
        }
    }
    //output
    for(int i = 0; i < gtf_data.size(); i++){
        string id;
        id_trim(gtf_data[i][8], id);
        auto itr = t_hash.find(id);
        if(itr != t_hash.end()){
            int count = itr->second;
            if(count >= min && count <= max){
                ss << gtf_data[i][0];
                for(int j = 1; j < gtf_data[i].size(); j++){
                    ss << '\t' << gtf_data[i][j];
                }
                ss << '\n';
            }
        }
    }
    ofstream fout;
    fout.open(argv[2]);
    fout << ss.str();
    fout.close();
    return 0;
}




