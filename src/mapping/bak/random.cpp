#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <stdio.h>
#include <random>
#include <unordered_set>
using namespace std;
int main(int argc, char* argv[]){
    //manual
    if (argc != 3 && argc != 4){
        cout << "---PURPOSE---" << endl;
        cout << " fileからランダムにN行を出力する" << endl;
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. input file" << endl;
        cout << " 2. the number of N" << endl;
        cout << " 3. seed (default=14)" << endl;
        return 0;
    }

    //data containor
    vector<string> data;
    int num = stoi(argv[2]);
    int seed = 14;
    if(argc == 4){
        seed = stoi(argv[3]);
    }
    unordered_set<int> id;

    //input file
    string tmpline;
    ifstream fin1;
    fin1.open(argv[1]);
    while(getline(fin1, tmpline)){
        data.push_back(tmpline);
    }
    fin1.close();

    //make rand
    int data_size = data.size();
    if(data_size < num){
        cout << "ERROR! N is larger than the number of lines in file." << endl;
        return 1;
    }
    mt19937 mt(seed);
    while(id.size() < num){
        int randn = mt() % data_size;
        id.insert(randn);
    }

    //output
    for(auto itr = id.begin(); itr != id.end(); itr++){
        cout << data[*itr] << '\n';
    }

    return 0;
}
        

