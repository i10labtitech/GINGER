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

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <map> 
using namespace std;

void split (const string &str, vector<string> &array){
    istringstream ss(str);
    string buffer;
    while(getline(ss, buffer, '\t')){
        array.push_back(buffer);
    }
}

void info(){
    cout << "Geneadd" << "\n\n";
    cout << "Version 1.0.0(19/11/15)" << "\n\n";
    cout << "Author: Fumiya Kobayashi" << "\n\n";
    cout << "Usage:" << '\n';
    cout << "  ./Geneadd [base gff] [additional gff]" << '\n';
    return;
}

int main(int argc, char* argv[]){
    //manual
    if (argc != 3){
        info();
        return 0;
    }

    //data structure
    unordered_map<string, map<int, char> > total_dict;
    stringstream ss;

    //input base gff
    string tmpline;
    vector<string> tmptab;
    
    ifstream fin1;
    fin1.open(argv[1]);
    while(getline(fin1, tmpline)){
        ss << tmpline << '\n';
        split(tmpline, tmptab);
        if(tmptab.size() > 8){
            if(tmptab[2] == "mRNA"){
                string tmp_scf = tmptab[0];
                int tmp_s = stoi(tmptab[3]);
                int tmp_e = stoi(tmptab[4]);
                string tmp_strand = tmptab[6];
                string tmp_key = tmp_scf + "_" + tmp_strand;
                auto itr = total_dict.find(tmp_key);
                if(itr == total_dict.end()){
                    map<int, char> empty_map;
                    total_dict[tmp_key] = empty_map;
                }
                total_dict[tmp_key][tmp_s] = 's';
                total_dict[tmp_key][tmp_e] = 'e';
            }
        }
        tmptab.clear();
    }
    fin1.close();

    //input additional gff
    ifstream fin2;
    fin2.open(argv[2]);
    bool flag = false;
    while(getline(fin2, tmpline)){
        split(tmpline, tmptab);
        if(tmptab.size() > 8){
            if(tmptab[2] == "mRNA"){
                string tmp_scf = tmptab[0];
                int tmp_s = stoi(tmptab[3]);
                int tmp_e = stoi(tmptab[4]);
                string tmp_strand = tmptab[6];
                string tmp_key = tmp_scf + "_" + tmp_strand;
                char outer_s, outer_e;
                if(total_dict[tmp_key].lower_bound(tmp_s) != total_dict[tmp_key].begin()){
                    outer_s = (--total_dict[tmp_key].lower_bound(tmp_s))->second;
                }else{
                    outer_s = 'n';
                }if(total_dict[tmp_key].upper_bound(tmp_e) != total_dict[tmp_key].end()){
                    outer_e = total_dict[tmp_key].lower_bound(tmp_e)->second;
                }else{
                    outer_e = 'n';
                }
                if((outer_s == 'e' || outer_s == 'n') && (outer_e == 's' || outer_e == 'n') && (total_dict[tmp_key].upper_bound(tmp_s)->first == total_dict[tmp_key].upper_bound(tmp_e)->first)){
                    flag = true;
                    ss << tmpline << '\n';
                }else{
                    flag = false;
                }
            }else{
                if(flag){
                    ss << tmpline << '\n';
                }
            }
        }
        tmptab.clear();
    }
    fin2.close();

    //output
    cout << ss.str();

    return 0;
}
