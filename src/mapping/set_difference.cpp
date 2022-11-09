#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <algorithm>
using namespace std;
int main(int argc, char* argv[]){
    //manual
    if (argc != 4){
        cout << "---PURPOSE---" << endl;
        cout << " 集合$1から集合$2に属する要素を間引いて差集合を得る" << endl;
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. set A" << endl;
        cout << " 2. set B" << endl;
        cout << " 3. A - B (output file name)" << endl;
        return 0;
    }
    // make dict
    unordered_set<string> hash;
    stringstream ss;

    // input B
    ifstream fin1;
    fin1.open(argv[2]);
    string templine;
    while(getline(fin1, templine)){
        hash.insert(templine);
    }
    fin1.close();
    
    //input A
    ifstream fin2;
    fin2.open(argv[1]);
    while(getline(fin2, templine)){
        auto itr = hash.find(templine);
        if(itr == hash.end()){
            ss << templine << "\n";
        }
    }
    fin2.close();
    
    //output
    ofstream fout;
    fout.open(argv[3]);
    fout << ss.str();
    fout.close();
    return 0;
}

