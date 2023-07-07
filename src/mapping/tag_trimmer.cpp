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
using namespace std;
void split (const string &str, vector<string> &array){
    istringstream ss(str);
    string buffer;
    while(getline(ss, buffer, '\t')){
        array.push_back(buffer);
    }
}
void tag_trim (const string &str, const string tag, string &value){
    string tag_plus = tag + '=';
    int start_pos = str.find(tag_plus) + tag.length() + 1;
    int end_pos = min(str.find(";", start_pos), str.length());
    int length = end_pos - start_pos;
    value = str.substr(start_pos, length);
}
int main(int argc, char* argv[]){
    //manual
    if (argc != 3){
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. 8th row from .gff file" << endl;
        cout << " 2. property name" << endl;
        return 0;
    }
    string property = argv[2];
    string templine;
    ifstream fin1;
    fin1.open(argv[1]);
    while(getline(fin1, templine)){
        string tag;
        tag_trim(templine, property, tag);
        cout << tag << endl;
    }
    fin1.close();
    return 0;
}

    
        


