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
#include <stdexcept>
#include <algorithm>

using namespace std;

double threshold = 0;
double sum_of_weights = 0;

struct mRNA{
    int beg;
    int end;
    double total_score;
    double total_cds_len;
    string predict_source;
};

//split tab-separated values
void split (const string &str, vector<string> &array){
    istringstream ss(str);
    string buffer;
    while(getline(ss, buffer, '\t')){
        array.push_back(buffer);
    }
}

void split_tag (const string &str, vector<string> &array){
    istringstream ss(str);
    string buffer;
    while(getline(ss, buffer, ';')){
        array.push_back(buffer);
    }
}

//pick out "value" from attribute column
void tag_trim (const string &str, const string tag, string &value){
    int start_pos = str.find(tag) + tag.length() + 1;
    int end_pos = min(str.find(";", start_pos), str.length());
    int length = end_pos - start_pos;
    value = str.substr(start_pos, length);
}

void info(){
    cout << "subgroup" << "\n\n";
    cout << "Version 2.2.0(19/11/15)" << "\n\n";
    cout << "Author: Fumiya Kobayashi" << "\n\n";
    cout << "Usage:" << '\n';
    cout << "  ./subgroup <Group.gff> <threshold> <sum of method weights>" << '\n';
    return;
}

vector<string> id_rename(vector<string> record, const int break_counter){
    if (record[2] == "mRNA"){
        vector<string> attribute;
        string new_attribute = "";
        split_tag(record[8], attribute);
        for (auto &item : attribute){
            if(item.substr(0,3) == "ID="){
                new_attribute += (item + "_" + to_string(break_counter) + ";");
            }else{
                new_attribute += (item + ";");
            }
        }
        new_attribute.pop_back();
        record[8] = new_attribute;
    }else if (record[2] == "CDS"){
        vector<string> attribute;
        string new_attribute = "";
        split_tag(record[8], attribute);
        for (auto &item : attribute){
            if(item.substr(0,7) == "Parent="){
                new_attribute += (item + "_" + to_string(break_counter) + ";");
            }else{
                new_attribute += (item + ";");
            }
        }
        new_attribute.pop_back();
        record[8] = new_attribute;
    }
    return record;
}

//input Group.gff file
void input_gff(ifstream &fin, vector< vector< vector<string> > > &total_data){
    size_t index = 0;
    string tmpline;
    vector<string> tmptab;
    while(getline(fin, tmpline)){
        split(tmpline, tmptab);
        if(tmptab.size() > 11){
            vector<string> tmp_data;
            tmp_data.assign(tmptab.begin()+3, tmptab.end());
            size_t tmp_index = stoi(tmptab[0]);
            if(index != tmp_index){
                vector< vector<string> > tmp_capsule;
                total_data.push_back(tmp_capsule);
                index++;
            }
            total_data[tmp_index-1].push_back(tmp_data);
        }
        tmptab.clear();
    }
}

//making dictionary about {mrna:list[cds]}
void make_mrna_dict(const vector< vector<string> > &group, unordered_map<string, mRNA> &dict, int &max_pos, int &min_pos){
    for (auto &record : group){
        if (record[2] == "mRNA"){
            int mrna_beg = stoi(record[3]);
            int mrna_end = stoi(record[4]);
            if (min_pos == 0 || min_pos > mrna_beg){
                min_pos = mrna_beg;
            }
            if (max_pos == 0 || max_pos < mrna_end){
                max_pos = mrna_end;
            }
            string mrna_id;
            tag_trim(record[8], "ID", mrna_id);
            mRNA tmp_mrna = {mrna_beg, mrna_end, 0, 0, record[1]};
            dict[mrna_id] = tmp_mrna;
        }else if (record[2] == "CDS"){
            int cds_beg = stoi(record[3]);
            int cds_end = stoi(record[4]);
            double cds_score = stod(record[5]);
            string cds_id;
            tag_trim(record[8], "Parent", cds_id);
            double cds_len = cds_end - cds_beg + 1;
            dict[cds_id].total_score += cds_len * cds_score;
            dict[cds_id].total_cds_len += cds_len;
        }
    }
}

void lookup_break_points(const unordered_map<string, mRNA> &dict, vector<int> &break_pos_list, const int max_pos, const int min_pos){
    vector<double> pos_list(max_pos - min_pos + 1, 0);
    double max_score = 0.0;
    //merge mrna score that was predicted in same method
    unordered_map<string, vector<double> > tmp_pos_dict;
    for (auto &item : dict){
        auto itr = tmp_pos_dict.find(item.second.predict_source);
        if(itr == tmp_pos_dict.end()){
            tmp_pos_dict[item.second.predict_source] = pos_list;
        }
        int tmp_beg = item.second.beg - min_pos;
        int tmp_end = item.second.end - min_pos;
        double tmp_score = item.second.total_score / item.second.total_cds_len;
        for (int i = 0; i < (tmp_end - tmp_beg + 1); ++i){
            if (tmp_pos_dict[item.second.predict_source][i + tmp_beg] < tmp_score){
                tmp_pos_dict[item.second.predict_source][i + tmp_beg] = tmp_score;
            }
        }
    }
    for (auto &method_score : tmp_pos_dict){
        for(int i = 0; i < (int)method_score.second.size(); ++i){
            pos_list[i] += method_score.second[i];
            if(pos_list[i] > max_score){
                max_score = pos_list[i];
            }
        }
    }
    //lookup break points
    bool flag = false;
    int tmp_s_pos = 0;
    for(int i = 0; i < (int)pos_list.size(); ++i){
        if (pos_list[i] <= (min(max_score, sum_of_weights) * threshold) && !flag){
            tmp_s_pos = i;
            flag = true;
        }else if(pos_list[i] > (min(max_score, sum_of_weights) * threshold) && flag){
            flag = false;
            if (tmp_s_pos != 0){
                int tmp_break_pos = ( (i - 1 + min_pos ) + (tmp_s_pos + min_pos) ) / 2;
                break_pos_list.push_back(tmp_break_pos);
            }
        }
    }
}

void classify_into_new_groups(const unordered_map<string, mRNA> &dict, const vector<int> &break_pos_list, unordered_map<string, int> &new_group_dict){
    for(auto &item : dict){
        bool flag = false;
        for(int i = 1; i < (int)break_pos_list.size()+1; ++i){ // when i = 0, the mRNA is separated into more than 1 new groups
            if(item.second.beg < break_pos_list[i-1] && item.second.end < break_pos_list[i-1]){
                flag = true;
                new_group_dict[item.first] = i;
                break;
            }else if(item.second.beg < break_pos_list[i-1] && item.second.end >= break_pos_list[i-1]){ //changed in ver2.1
                flag = true;
                new_group_dict[item.first] = 0;
                break;
            }
        }
        if(!flag){
            new_group_dict[item.first] = break_pos_list.size() + 1;
        }
    }
}

void make_new_gff(vector< vector<string> > &group, const unordered_map<string, int> &new_group_dict, vector< vector< vector<string> > > &new_group, vector<int> &new_groups_num, const vector<int> &break_pos_list){
    vector<string> tmp_break_mrna;
    int break_counter = 0;
    int empty_flag = -1;
    for (auto record : group){
        string tmp_ID;
        if (record[2] == "mRNA"){
            tag_trim(record[8], "ID", tmp_ID);
        }else if(record[2] == "CDS"){
            tag_trim(record[8], "Parent", tmp_ID);
        }
        int tmp_index = new_group_dict.at(tmp_ID);
        if(tmp_index != 0){
            new_group[tmp_index-1].push_back(record);
            if(record[2] == "mRNA"){
                new_groups_num[tmp_index-1] += 1;
            }
        }else{
            //specifying CDSs in break mrna to new groups
            if(record[2] == "mRNA"){
                tmp_break_mrna = record;
                break_counter = 0;
                empty_flag = -1;
            }else{
                int cds_beg = stoi(record[3]);
                int cds_end = stoi(record[4]);
                bool flag = false;
                for(int i = 0 ; i < (int)break_pos_list.size(); ++i){
                    if(cds_beg < break_pos_list[i] && cds_end < break_pos_list[i]){
                        flag = true;
                        if(empty_flag != i){
                            ++break_counter;
                            new_group[i].push_back(id_rename(tmp_break_mrna, break_counter));
                            new_groups_num[i] += 1;
                            empty_flag = i;
                        }
                        new_group[i].push_back(id_rename(record, break_counter));
                        break;
                    }else if(cds_beg < break_pos_list[i] && cds_end >= break_pos_list[i]){
                        flag = true;
                        break;
                    }
                }
                if(!flag){
                    if(cds_beg > break_pos_list.back()){
                        if(empty_flag != (int)break_pos_list.size()){
                            ++break_counter;
                            new_group.back().push_back(id_rename(tmp_break_mrna, break_counter));
                            new_groups_num.back() += 1;
                            empty_flag = (int)break_pos_list.size();
                        }
                        new_group.back().push_back(id_rename(record, break_counter));
                    }
                }
            }
        }
    }
}

void subgroup(vector< vector<string> > &group, stringstream &ss, int &counter){
    unordered_map<string, mRNA> record_dict;
    int max_pos = 0;
    int min_pos = 0;
    //making dictionary about {mrna:list[cds]}
    make_mrna_dict(group, record_dict, max_pos, min_pos);

    vector<int> break_pos_list;
    lookup_break_points(record_dict, break_pos_list, max_pos, min_pos);
    
    unordered_map<string, int> new_group_dict;
    classify_into_new_groups(record_dict, break_pos_list, new_group_dict);
    
    vector< vector< vector<string> > > new_group(break_pos_list.size() + 1);
    vector<int> new_groups_num((int)break_pos_list.size() + 1, 0);
    make_new_gff(group, new_group_dict, new_group, new_groups_num, break_pos_list);
    
    /*---debug---
    ss << '#' << counter << ":";
    for(auto x : break_pos_list){
        ss << x << ",";
    }
    ss << endl;
    */

    //output result
    for(int i = 0; i < (int)new_group.size(); ++i){
        if(new_groups_num[i] != 0){
            for(auto &item : new_group[i]){
                ss << counter << '\t' << new_groups_num[i] << "\t0";
                for(auto &value : item){
                    ss << '\t' << value;
                }
                ss << '\n';
            }
            ++counter;
        }
    }
}

int main(int argc, char* argv[]){
    //manual
    if (argc != 4){
        info();
        return 0;
    }
    //set parameter
    try{
        threshold = stod(argv[2]);
    }catch(const std::invalid_argument &e){
        cerr << "Error: invalid threshold" << endl;
        exit(0);
    }catch(const std::out_of_range &e){
        cerr << "Error: out of range threshold" << endl;
        exit(0);
    }

    try{
        sum_of_weights = stod(argv[3]);
    }catch(const std::invalid_argument &e){
        cerr << "Error: invalid sum-of-weights" << endl;
        exit(0);
    }catch(const std::out_of_range &e){
        cerr << "Error: out of range sum-of-weights" << endl;
        exit(0);
    }

    //prepair data containor
    vector< vector <vector<string> > > total_data;

    //input Group.gff file
    ifstream fin1;
    fin1.open(argv[1]);
    input_gff(fin1, total_data);
    
    stringstream ss;
    int counter = 1;
    for (auto &i : total_data){
        subgroup(i, ss, counter);
    }
    
    cout << ss.str();
    return 0;
}
