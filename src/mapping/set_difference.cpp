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
#include <unordered_set>
#include <algorithm>
using namespace std;
int main(int argc, char* argv[]){
    //manual
    if (argc != 4){
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

