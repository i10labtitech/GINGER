#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
using namespace std;
struct Pos{
    int s_pos;
    int e_pos;
    int len;
};
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
        cout << " TransDecoderの結果をORF finderの結果を元に書き換える" << endl;
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. ORF_finder.gff3" << endl;
        cout << " 2. TransDecoder.gff3" << endl;
        cout << " 3. output file name" << endl;
        return 0;
    }
    
    //data containor
    unordered_map<string, Pos> hash;
    stringstream ss;

    //input ORF_finder file
    string tmpline;
    vector<string> tmptab;
    
    ifstream fin1;
    fin1.open(argv[1]);
    while(getline(fin1, tmpline)){
        split(tmpline, tmptab);
        if(tmptab.size() >= 9){
            if(tmptab[2] == "CDS"){
                Pos tmppos = {stoi(tmptab[3]), stoi(tmptab[4], 0)};
                hash[tmptab[0]] = tmppos;
            }
        }
        tmptab.clear();
    }
    fin1.close();
    fin1.open(argv[1]);
    while(getline(fin1, tmpline)){
        split(tmpline, tmptab);
        if(tmptab.size() >= 9){
            if(tmptab[2] == "mRNA"){
                int tmplen = stoi(tmptab[4]) - stoi(tmptab[3]);
                Pos tmppos = {hash[tmptab[0]].s_pos, hash[tmptab[0]].e_pos, tmplen};
                hash[tmptab[0]] = tmppos;
            }
        }
        tmptab.clear();
    }
    fin1.close();

    //replace TransDecoder file
    string tmpgene;
    ifstream fin2;
    fin2.open(argv[2]);
    while(getline(fin2, tmpline)){
        split(tmpline, tmptab);
        if(tmptab.size() >= 9){
            if(tmptab[2] == "gene"){
                tmpgene = tmpline;
            }else if(tmptab[2] == "mRNA"){
                string tmptag;
                tag_trim(tmptab[8], "ID", tmptag);
                auto itr = hash.find(tmptag);
                if(itr != hash.end()){
                    ss << tmpgene << '\n';
                    ss << tmpline << '\n';
                }
            }else if(tmptab[2] == "exon"){
                string tmptag;
                tag_trim(tmptab[8], "Parent", tmptag);
                auto itr = hash.find(tmptag);
                if(itr != hash.end()){
                    ss << tmpline << '\n';
                }
            }else if(tmptab[2] == "CDS"){
                string tmptag;
                tag_trim(tmptab[8], "Parent", tmptag);
                auto itr = hash.find(tmptag);
                if(itr != hash.end()){
                    int s = stoi(tmptab[3]);
                    int e = stoi(tmptab[4]);
                    int s_diff = itr->second.s_pos - 1;
                    int e_diff = itr->second.len - itr->second.e_pos + 1;
                    int re_s = s + s_diff;
                    int re_e = e - e_diff;
                    ss << tmptab[0] << '\t' << tmptab[1] << '\t' << tmptab[2] << '\t';
                    ss << re_s << '\t' << re_e << '\t';
                    ss << tmptab[5] << '\t' << tmptab[6] << '\t' << tmptab[7] << '\t' << tmptab[8] << '\n';
                }
            }else if(tmptab[2] == "five_prime_UTR"){
                string tmptag;
                tag_trim(tmptab[8], "Parent", tmptag);
                auto itr = hash.find(tmptag);
                if(itr != hash.end()){
                    int s = stoi(tmptab[3]);
                    int e = stoi(tmptab[4]);
                    int s_diff = itr->second.s_pos - 1;
                    int e_diff = itr->second.len - itr->second.e_pos + 1;
                    int re_s = s;
                    int re_e = e + s_diff;
                    ss << tmptab[0] << '\t' << tmptab[1] << '\t' << tmptab[2] << '\t';
                    ss << re_s << '\t' << re_e << '\t';
                    ss << tmptab[5] << '\t' << tmptab[6] << '\t' << tmptab[7] << '\t' << tmptab[8] << '\n';
                }
            }else if(tmptab[2] == "three_prime_UTR"){
                string tmptag;
                tag_trim(tmptab[8], "Parent", tmptag);
                auto itr = hash.find(tmptag);
                if(itr != hash.end()){
                    int s = stoi(tmptab[3]);
                    int e = stoi(tmptab[4]);
                    int s_diff = itr->second.s_pos - 1;
                    int e_diff = itr->second.len - itr->second.e_pos + 1;
                    int re_s = s - e_diff;
                    int re_e = e;
                    ss << tmptab[0] << '\t' << tmptab[1] << '\t' << tmptab[2] << '\t';
                    ss << re_s << '\t' << re_e << '\t';
                    ss << tmptab[5] << '\t' << tmptab[6] << '\t' << tmptab[7] << '\t' << tmptab[8] << '\n';
                }
            }
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



                    




