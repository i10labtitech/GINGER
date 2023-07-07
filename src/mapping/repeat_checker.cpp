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
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <set>
using namespace std;
void split (const string &str, vector<string> &array){
    istringstream ss(str);
    string buffer;
    while(getline(ss, buffer, '\t')){
        array.push_back(buffer);
    }
}
int main(int argc, char* argv[]){
    //manual
    if (argc != 5){
        cout << "---PURPOSE---" << endl;
        cout << " 配列がrepeat領域に被っているか確認する" << endl;
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. genome stats file ([seq name] [length])" << endl;
        cout << " 2. repeat stats file ([seq name] [start pos] [end pos] [repeat character]" << endl;
        cout << " 3. sequence stats file ([seq name] [start pos] [end pos] [ID])" << endl;
        cout << " 4. output file name" << endl;
        return 0;
    }
    //make database
    unordered_map<string, int> seq_name;
    vector< vector<int> > repeat;
    vector<string> repeat_char;

    string templine;
    vector<string> temptab;
    ifstream fin1;
    fin1.open(argv[1]);
    int counter = 0;
    while(getline(fin1, templine)){
        split(templine, temptab);
        seq_name[temptab[0]] = counter;
        vector<int> tempvec(stoi(temptab[1]), -1);
        repeat.push_back(tempvec);
        counter++;
        temptab.clear();
    }
    fin1.close();

    ifstream fin2;
    fin2.open(argv[2]);
    counter = 0;
    while(getline(fin2, templine)){
        split(templine, temptab);
        repeat_char.push_back(temptab[3]);
        int start = stoi(temptab[1]) - 1;
        int end = stoi(temptab[2]) - 1;
        for(int i = start; i <= end; i++){
            repeat[seq_name[temptab[0]]][i] = counter;
        }
        counter++;
        temptab.clear();
    }

    //search and output
    stringstream ss;
    ifstream fin3;
    fin3.open(argv[3]);
    while(getline(fin3, templine)){
        split(templine, temptab);
        int start = stoi(temptab[1]) - 1;
        int end = stoi(temptab[2]) - 1;
        set<int> stock;
        for(int i = start; i <= end; i++){
            if(repeat[seq_name[temptab[0]]][i] != -1){
                stock.insert(repeat[seq_name[temptab[0]]][i]);
            }
        }
        if(stock.empty()){
            ss << temptab[3] << "\t" << 0 << "\n";
        }else{
            ss << temptab[3] << "\t" << stock.size();
            for(auto itr = stock.begin(); itr != stock.end(); itr++){
                ss << "\t" << repeat_char[*itr];
            }
            ss << "\n";
        }
        temptab.clear();
    }
    fin3.close();

    //output
    ofstream fout;
    fout.open(argv[4]);
    fout << ss.str();
    fout.close();
    return 0;
}
