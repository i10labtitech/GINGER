#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
using namespace std;
//strand 変換関数
void reverse_st(const string &seq, string &rev_seq){
    for(int i = seq.size()-1; i >= 0; i--){
        switch(seq[i]){
            case 'A': rev_seq += 'T'; break;
            case 'C': rev_seq += 'G'; break;
            case 'G': rev_seq += 'C'; break;
            case 'T': rev_seq += 'A'; break;
            case 'N': rev_seq += 'N'; break;
            default : break;
        }
    }
}
//+ strand 最長ORF探索関数
void p_finder(const string &seq, int &s_pos, int &e_pos, int &len, const int &frame){
    int temp_len = 0;
    int temp_s = 0;
    int temp_e = 0;
    for(int i = frame; i < seq.size() - 2; i += 3){
        if (seq.substr(i, 3) == "ATG"){
            int now_s = i;
            //start codon
            for(int j = i + 3; j < seq.size() - 2; j += 3){
                if( seq.substr(j, 3) == "TAA" || seq.substr(j, 3) == "TAG" || seq.substr(j, 3) == "TGA"){
                    if( (j - now_s + 3 ) > temp_len){
                        temp_s = now_s;
                        temp_e = j+2;
                        temp_len = j-now_s+3;
                    }
                    break;
                }
                i = j + 3;
            }
        }
    }
    s_pos = temp_s;
    e_pos = temp_e;
    len = temp_len;
}
//- strand 最長ORF探索関数
void m_finder(const string &seq, int &s_pos, int &e_pos, int &len, const int &frame){
    int temp_len = 0;
    int temp_s = 0;
    int temp_e = 0;
    for(int i = seq.size() - frame - 1; i > 1; i -= 3){
        if (seq.substr(i-2, 3) == "CAT"){
            int now_s = i;
            //start codon
            for(int j = i - 3; j > 1; j -= 3){
                if( seq.substr(j-2, 3) == "TTA" || seq.substr(j-2, 3) == "CTA" || seq.substr(j-2, 3) == "TCA"){
                    if( (now_s - j + 3 ) > temp_len){
                        temp_s = now_s;
                        temp_e = j-2;
                        temp_len = now_s-j+3;
                    }
                    break;
                }
                i = j - 3;
            }
        }
    }
    s_pos = temp_s;
    e_pos = temp_e;
    len = temp_len;
}
//最長ORF探索関数
void orf_finder(const string &seq, int &start, int &end, int &len, bool &strand){
    len = 0;
    vector<int> s_set (6, 0);
    vector<int> e_set (6, 0);
    vector<int> len_set (6, 0);
    p_finder(seq, s_set[0], e_set[0], len_set[0], 0);
    p_finder(seq, s_set[1], e_set[1], len_set[1], 1);
    p_finder(seq, s_set[2], e_set[2], len_set[2], 2);
    if (len_set[0] > len_set[1] && len_set[0] > len_set[2]){
        start = s_set[0];
        end = e_set[0];
        len = len_set[0];
    }else if(len_set[1] > len_set[2]){
        start = s_set[1];
        end = e_set[1];
        len = len_set[1];
    }else{
        start = s_set[2];
        end = e_set[2];
        len = len_set[2];
    }
    if (strand == false){
        m_finder(seq, s_set[3], e_set[3], len_set[3], 0);
        m_finder(seq, s_set[4], e_set[4], len_set[4], 1);
        m_finder(seq, s_set[5], e_set[5], len_set[5], 2);
        if (len_set[3] > len_set[4] && len_set[3] > len_set[5] && len_set[3] > len){
            start = s_set[3];
            end = e_set[3];
            len = len_set[3];
        }else if(len_set[4] > len_set[5] && len_set[4] > len){
            start = s_set[4];
            end = e_set[4];
            len = len_set[4];
        }else if(len_set[5] > len){
            start = s_set[5];
            end = e_set[5];
            len = len_set[5];
        }
    }
}
//塩基配列 -> アミノ酸配列置換関数
/*
void translation (const string &str, string &pep_str){
    string amino; //アミノ酸配列
    int i = 0;
    string codon = "";
    bool stop = true;
    while(stop){
        codon = str.substr(i,3);
        if (codon[0] == 'T'){
            if (codon[1] == 'T'){
                if (codon[2] == 'T' || codon[2] == 'C'){
                    amino += "F";
                }else{
                    amino += "L";
                }
            }else if (codon[1] == 'C'){
                amino += "S";
            }else if (codon[1] == 'A'){
                if (codon[2] == 'T' || codon[2] == 'C'){
                    amino += "Y";
                }else{
                    stop = false;
                }
            }else{
                if (codon[2] == 'T' || codon[2] == 'C'){
                    amino += "C";
                }else if (codon[2] == 'A'){
                    stop = false;
                }else{
                    amino += "W";
                }
            }
        }else if (codon[0] == 'C'){
            if (codon[1] == 'T'){
                amino += "L";
            }else if (codon[1] == 'C'){
                amino += "P";
            }else if (codon[1] == 'A'){
                if (codon[2] == 'T' || codon[2] == 'C'){
                    amino += 'H';
                }else{
                    amino += 'Q';
                }
            }else{
                amino += 'R';
            }
        }else if (codon[0] == 'A'){
            if (codon[1] == 'T'){
                if (codon[2] == 'T' || codon[2] == 'C' || codon[2] == 'A'){
                    amino += "I";
                }else{
                    amino += "M";
                }
            }else if (codon[1] == 'C'){
                amino += "T";
            }else if (codon[1] == 'A'){
                if (codon[2] == 'T' || codon[2] == 'C'){
                    amino += "N";
                }else{
                    amino += "K";
                }
            }else{
                if (codon[2] == 'T' || codon[2] == 'C'){
                    amino += "S";
                }else{
                    amino += "R";
                }
            }
        }else{
            if (codon[1] == 'T'){
                amino += "V";
            }else if (codon[1] == 'C'){
                amino += "A";
            }else if (codon[1] == 'A'){
                if (codon[2] == 'T' || codon[2] == 'C'){
                    amino += "D";
                }else{
                    amino += "E";
                }
            }else{
                amino += "G";
            }
        }
        i += 3;
    }
    pep_str = amino;
}
*/
//gff3作成関数
void make_gff(const string &seq_name, const int &s, const int &e, const int &all_len, stringstream &ss){
    int space = seq_name.find(" ");
    string seq_parts = seq_name.substr(1, space-1);
    if(s < e){
        ss << seq_parts << "\tORF_finder\tgene\t1\t" << all_len << "\t.\t+\t.\tID=GENE." << seq_parts << "\n";
        ss << seq_parts << "\tORF_finder\tmRNA\t1\t" << all_len << "\t.\t+\t.\tID=" << seq_parts << ";Parent=GENE." << seq_parts << "\n";
        if(s != 0){
            ss << seq_parts << "\tORF_finder\tfive_prime_UTR\t1\t" << s << "\t.\t+\t.\tID=" << seq_parts << ".utr5;Parent=" << seq_parts << "\n";
        }
        ss << seq_parts << "\tORF_finder\texon\t1\t" << all_len << "\t.\t+\t.\tID=" << seq_parts << ".exon;Parent=" << seq_parts << "\n";
        ss << seq_parts << "\tORF_finder\tCDS\t" << s + 1 << "\t" << e + 1 << "\t.\t+\t.\tID=cds." << seq_parts << ";Parent=" << seq_parts << "\n";
        if(e != all_len - 1){
            ss << seq_parts << "\tORF_finder\tthree_prime_UTR\t" << e + 2 << "\t" << all_len << "\t.\t+\t.\tID=" << seq_parts << ".utr3;Parent=" << seq_parts << "\n";
        }
    }else{
        ss << seq_parts << "\tORF_finder\tgene\t1\t" << all_len << "\t.\t-\t.\tID=GENE." << seq_parts << "\n";
        ss << seq_parts << "\tORF_finder\tmRNA\t1\t" << all_len << "\t.\t-\t.\tID=" << seq_parts << ";Parent=GENE." << seq_parts << "\n";
        if(e != 0){
            ss << seq_parts << "\tORF_finder\tfive_prime_UTR\t" << s + 2 << "\t" << all_len << "\t.\t-\t.\tID=" << seq_parts << ".utr5;Parent=" << seq_parts << "\n";
        }
        ss << seq_parts << "\tORF_finder\texon\t1\t" << all_len << "\t.\t-\t.\tID=" << seq_parts << ".exon;Parent=" << seq_parts << "\n";
        ss << seq_parts << "\tORF_finder\tCDS\t" << e + 1 << "\t" << s + 1 << "\t.\t-\t.\tID=cds." << seq_parts << ";Parent=" << seq_parts << "\n";
        if(s != all_len - 1){
            ss << seq_parts << "\tORF_finder\tthree_prime_UTR\t1\t" << e << "\t.\t-\t.\tID=" << seq_parts << ".utr3;Parent=" << seq_parts << "\n";
        }
    }
}
int main(int argc, char* argv[]){
    //manual
    if (argc <= 2){
        cout << "---PURPOSE---" << endl;
        cout << " fasta file から各seqの最長ORFを検出する" << endl;
        cout << " .cdsと.gff3を出力する" << endl;
        cout << "---ARGUMENTs---" << endl;
        cout << " 1. fasta file" << endl;
        cout << " 2. output file PREFIX" << endl;
        cout << " 3. minimum length (base)" << endl;
        cout << " 4. strand specific (true or false) (default = false)" << endl;
        cout << "---!ATTENTION!---" << endl;
        cout << " minimum lengthは300が推奨. あまり小さいと擬陽性が増えるだけでなく, バグが発生するかも..." << endl;
        cout << " fasta file中の塩基は全て大文字で表記する(seqkit seq -u で変換可能)" << endl;
        return 0;
    }
    //min len
    int min_len;
    if(argc > 2){
        min_len = stoi(argv[3]);
    }
    //strand specific
    bool strand = false;
    string strand_state = argv[4];
    if(argc > 3 && strand_state == "true"){
        strand = true;
    }
    //input fasta file
    stringstream ss_cds, ss_pep, ss_gff;
    string templine;
    string tempseqname;
    string tempseq = "";
    ifstream fin;
    fin.open(argv[1]);
    while(getline(fin, templine)){
        if(templine[0] != '>'){
            tempseq += templine;
        }else if(templine[0] == '>' && tempseq == ""){
            tempseqname = templine;
        }else{
            int start, end, len;
            orf_finder(tempseq, start, end, len, strand);
            if(start < end && len >= min_len){
                ss_cds << tempseqname << " len:" << len << " " << start+1 << "-" << end+1 << "(+)" << endl;
                ss_cds << tempseq.substr(start, len) << endl;
                string pep;
                //translation(tempseq.substr(start, len), pep);
                //ss_pep << tempseqname << " len:" << len/3 << " " << start+1 << "-" << end+1 << "(+)" << endl;
                //ss_pep << pep << "*" <<endl;
                make_gff(tempseqname, start, end, tempseq.length(), ss_gff);
            }else if(len >= min_len){
                ss_cds << tempseqname << " len:" << len << " " << start+1 << "-" << end+1 << "(-)" << endl;
                string tempsubseq = tempseq.substr(end, len);
                string rev_seq = "";
                reverse_st(tempsubseq, rev_seq);
                ss_cds << rev_seq << endl;
                string pep;
                //translation(rev_seq, pep);
                //ss_pep << tempseqname << " len:" << len/3 << " " << start+1 << "-" << end+1 << "(-)" << endl;
                //ss_pep << pep << "*" << endl;
                make_gff(tempseqname, start, end, tempseq.length(), ss_gff);
            }
            tempseq = "";
            tempseqname = templine;
        }
    }
    int start, end, len;
    orf_finder(tempseq, start, end, len, strand);
    if(start < end && len >= min_len){
        ss_cds << tempseqname << " len:" << len << " " << start+1 << "-" << end+1 << "(+)" << endl;
        ss_cds << tempseq.substr(start, len) << endl;
        string pep;
        //translation(tempseq.substr(start, len), pep);
        //ss_pep << tempseqname << " len:" << len/3 << " " << start+1 << "-" << end+1 << "(+)" << endl;
        //ss_pep << pep << "*" <<endl;
        make_gff(tempseqname, start, end, tempseq.length(), ss_gff);
    }else if(len >= min_len){
        ss_cds << tempseqname << " len:" << len << " " << start+1 << "-" << end+1 << "(-)" << endl;
        string tempsubseq = tempseq.substr(end, len);
        string rev_seq = "";
        reverse_st(tempsubseq, rev_seq);
        ss_cds << rev_seq << endl;
        string pep;
        //translation(rev_seq, pep);
        //ss_pep << tempseqname << " len:" << len/3 << " " << start+1 << "-" << end+1 << "(-)" << endl;
        //ss_pep << pep << "*" << endl;
        make_gff(tempseqname, start, end, tempseq.length(), ss_gff);
    }
    fin.close();
    //output file
    string cds_file = argv[2];
    cds_file += ".cds";
    //string pep_file = argv[2];
    //pep_file += ".pep";
    string gff_file = argv[2];
    gff_file += ".gff3";
    ofstream fout;
    fout.open(cds_file);
    fout << ss_cds.str();
    fout.close();
    ofstream fout2;
    //fout2.open(pep_file);
    //fout2 << ss_pep.str();
    //fout2.close();
    ofstream fout3;
    fout3.open(gff_file);
    fout3 << ss_gff.str();
    fout3.close();
    
    return 0;
}
