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
    int start_pos = str.find(tag) + tag.length() + 2;
    int end_pos = min(str.find(";", start_pos), str.length()) - 1;
    int length = end_pos - start_pos;
    value = str.substr(start_pos, length);
}

int main(int argc, char* argv[]){
    //manual
    if (argc != 4){
        cout << "---PURPOSE---" << endl;
        cout << " single exon gffにstrandを与える" << endl;
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. ORF_finder_output.gff3" << endl;
        cout << " 2. single_exon.gtf" << endl;
        cout << " 3. output file name" << endl;
        return 0;
    }
    //data containor
    unordered_map <string, string> hash;
    stringstream ss;

    //input gff
    string tmpline;
    vector<string> tmptab;
    ifstream fin1;
    fin1.open(argv[1]);
    while(getline(fin1, tmpline)){
        split(tmpline, tmptab);
        if(tmptab.size() > 8){
            if(tmptab[2] == "mRNA"){
                int e_pos = tmptab[0].find(".", 5);
                string genename = tmptab[0].substr(0, e_pos);
                hash[genename] = tmptab[6];
            }
        }
        tmptab.clear();
    }
    fin1.close();

    //input gtf
    ifstream fin2;
    fin2.open(argv[2]);
    while(getline(fin2, tmpline)){
        split(tmpline, tmptab);
        if(tmptab.size() > 8){
            string id;
            tag_trim(tmptab[8], "gene_id", id);
            auto itr = hash.find(id);
            if(itr != hash.end()){
                tmptab[6] = hash[id];
            }
            for(int i = 0; i < tmptab.size() - 1; i++){
                ss << tmptab[i] << '\t';
            }
            ss << tmptab[tmptab.size() - 1] << '\n';
        }else{
            ss << tmpline << '\n';
        }
        tmptab.clear();
    }
    fin2.close();

    //output
    ofstream fout;
    fout.open(argv[3]);
    fout << ss.str();
    fout.close();

    return 0;
}


