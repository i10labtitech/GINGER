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

//
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <algorithm> 

using namespace std;

vector<string> Split(const string &s, char delim)
{
	vector<string> elems;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim))
	{
		if (!item.empty())
		{
			elems.push_back(item);
		}
	}
	return elems;
}

int StoI(string str)
{
	int number;
	istringstream iss(str);
	iss >> number;
	return number;
}

double StoD(string str)
{
	double number;
	istringstream iss(str);
	iss >> number;
	return number;
}

string ItoS(int number)
{
	stringstream ss;
	ss << number;
	return ss.str();
}

string Erace(string testString, string tmp)
{
	for(int c = testString.find_first_of(tmp); c != string::npos; c = c = testString.find_first_of(tmp))
	{
		testString.erase(c,1);
	}
	return testString;
}

string PChange(string nucl)
{
	string nucl2, prot;
	int protnum;
	nucl2 = Erace(nucl,"-");
	nucl = nucl2;
	protnum = nucl.size()/3;
	for(int i = 0; i < protnum; i++)
	{
		if(nucl[3*i] == 'C' && nucl[3*i+1] == 'T')
		{
			prot += 'L';
		}
		else if(nucl[3*i] == 'G' && nucl[3*i+1] == 'T')
		{
			prot += 'V';
		}
		else if(nucl[3*i] == 'T' && nucl[3*i+1] == 'C')
		{
			prot +=  'S';
		}
		else if(nucl[3*i] == 'C' && nucl[3*i+1] == 'C')
		{
			prot +=  'P';
		}
		else if(nucl[3*i] == 'A' && nucl[3*i+1] == 'C')
		{
			prot +=  'T';
		}
		else if(nucl[3*i] == 'G' && nucl[3*i+1] == 'C')
		{
			prot +=  'A';
		}
		else if(nucl[3*i] == 'C' && nucl[3*i+1] == 'G')
		{
			prot +=  'R';
		}
		else if(nucl[3*i] == 'G' && nucl[3*i+1] == 'G')
		{
			prot +=  'G';
		}
		else if(nucl[3*i] == 'T' && nucl[3*i+1] == 'T' && (nucl[3*i+2] =='T' || nucl[3*i+2] == 'C') )
		{
			prot += 'F';
		}
		else if(nucl[3*i] == 'T' && nucl[3*i+1] == 'T' && (nucl[3*i+2] =='A' || nucl[3*i+2] == 'G') )
		{
			prot += 'L';
		}
		else if(nucl[3*i] == 'T' && nucl[3*i+1] == 'A' && (nucl[3*i+2] =='T' || nucl[3*i+2] == 'C') )
		{
			prot += 'Y';
		}
		else if(nucl[3*i] == 'C' && nucl[3*i+1] == 'A' && (nucl[3*i+2] =='T' || nucl[3*i+2] == 'C') )
		{
			prot += 'H';
		}
		else if(nucl[3*i] == 'C' && nucl[3*i+1] == 'A' && (nucl[3*i+2] =='A' || nucl[3*i+2] == 'G') )
		{
			prot += 'Q';
		}
		else if(nucl[3*i] == 'A' && nucl[3*i+1] == 'A' && (nucl[3*i+2] =='T' || nucl[3*i+2] == 'C') )
		{
			prot += 'N';
		}
		else if(nucl[3*i] == 'A' && nucl[3*i+1] == 'A' && (nucl[3*i+2] =='A' || nucl[3*i+2] == 'G') )
		{
			prot += 'K';
		}
		else if(nucl[3*i] == 'G' && nucl[3*i+1] == 'A' && (nucl[3*i+2] =='T' || nucl[3*i+2] == 'C') )
		{
			prot += 'D';
		}
		else if(nucl[3*i] == 'G' && nucl[3*i+1] == 'A' && (nucl[3*i+2] =='A' || nucl[3*i+2] == 'G') )
		{
			prot += 'E';
		}
		else if(nucl[3*i] == 'T' && nucl[3*i+1] == 'G' && (nucl[3*i+2] =='T' || nucl[3*i+2] == 'C') )
		{
			prot += 'C';
		}
		else if(nucl[3*i] == 'A' && nucl[3*i+1] == 'G' && (nucl[3*i+2] =='T' || nucl[3*i+2] == 'C') )
		{
			prot += 'S';
		}
		else if(nucl[3*i] == 'A' && nucl[3*i+1] == 'G' && (nucl[3*i+2] =='A' || nucl[3*i+2] == 'G') )
		{
			prot += 'R';
		}
		else if(nucl[3*i] == 'A' && nucl[3*i+1] == 'T' && (nucl[3*i+2] =='T' || nucl[3*i+2] == 'C' || nucl[3*i+2] == 'A') )
		{
			prot += 'I';
		}
		else if(nucl[3*i] == 'A' && nucl[3*i+1] == 'T' && nucl[3*i+2] =='G' )
		{
			prot += 'M';
		}
		else if(nucl[3*i] == 'T' && nucl[3*i+1] == 'G' && nucl[3*i+2] =='G' )
		{
			prot += 'W';
		}
		else if(nucl[3*i] == 'T' && nucl[3*i+1] == 'A' && nucl[3*i+2] =='A' )
		{
			prot += '*';
		}
		else if(nucl[3*i] == 'T' && nucl[3*i+1] == 'A' && nucl[3*i+2] =='G' )
		{
			prot += '*';
		}
		else if(nucl[3*i] == 'T' && nucl[3*i+1] == 'G' && nucl[3*i+2] =='A' )
		{
			prot += '*';
		}
		else
		{
			prot += '#';
		}
	}
	return prot;
}


vector<int> Histmake(vector<double> V, double stepsize, int max)
{
	vector<int> V2(max,0);
	int vsize = V.size();
	int tmp;
	for(int i = 0; i < vsize; i++)
	{
		for(int j = 0; j < max; j++)
		{
			if(j == max-1)
			{
				if( (double)j * stepsize <= V[i] && (double)(j+1) * stepsize >= V[i] )
				{
					tmp = V2[j];
					tmp++;
					V2[j] = tmp;
				}
			}
			else
			{
				if( (double)j * stepsize <= V[i] && (double)(j+1) * stepsize > V[i] )
				{
					tmp = V2[j];
					tmp++;
					V2[j] = tmp;
				}
			}
		}
	}
	return V2;
}


string Reverseorder(string pre)
{
	string rev;
	for(int i = 0; i < pre.length(); i++)
	{
		rev += pre[pre.length()-1-i];
	}
	return(rev);
}

string Greverse(string frev)
{
	string rev;
	for (int rv = 0; rv < frev.length(); ++rv) 
	{
		if (frev[frev.length()-1 - rv] == 'A')
    	{	
        	rev += "T";
		}	
		else if (frev[frev.length()-1 - rv] == 'T')
    	{
			rev += "A";
    	}
		else if (frev[frev.length()-1 - rv] == 'C')
    	{
       		rev += "G";
    	}
		else if (frev[frev.length()-1 - rv] == 'G')
    	{
        	rev += "C";
    	}
		else
		{
			rev += "N";
		}
	}
	return(rev);
}

string stoS(string bef)
{
	transform(bef.begin(), bef.end(), bef.begin(), ::toupper);
	return(bef);
}


string Stos(string bef)
{
	transform(bef.begin(), bef.end(), bef.begin(), ::tolower);
	return(bef);
}

unordered_map<string,string> Fasta_Hash(istream &file)
{
	unordered_map<string,string> hash;
	string lin, scfn, fas, zero;
	while(getline(file,lin) )
	{
		string::size_type a1 = lin.find(">");
		if(a1 != string::npos && fas == zero)
		{
			scfn = lin.substr(1);
		}
		else if(a1 != string::npos)
		{
			hash[scfn] = fas;
			fas = zero;
			scfn = lin.substr(1);
		}
		else
		{
			fas += lin;
		}
	}
	hash[scfn] = fas;
	return hash;
}

vector<vector<string> > Fasta_V(istream &file)
{
	vector<string> V;
	vector<vector<string> > VV;

	string lin, scfn, fas, zero;
	while(getline(file,lin) )
	{
		string::size_type a1 = lin.find(">");
		if(a1 != string::npos && fas == zero)
		{
			scfn = lin.substr(1);
			V.push_back(scfn);
		}
		else if(a1 != string::npos)
		{
			V.push_back(fas);
			VV.push_back(V);
			V.clear();
			fas = zero;
			scfn = lin.substr(1);
			V.push_back(scfn);
		}
		else
		{
			fas += lin;
		}
	}
	V.push_back(fas);
	VV.push_back(V);
	return VV;
}

string Replaced_V(string str, string replacestr, int replace)
{
	string output;
	vector<string> I;
	I = Split(str, '\t');
	for(int i = 0; i < replace - 1; i++)
	{
		output += I[i];
		output += '\t';
	}
	output += replacestr;
	for(int i = replace; i < I.size(); i++)
	{
		output+= '\t';
		output+= I[i];
	}
	return(output);
}

string Replace(string str, string target, string replacement)
{
	string output;
	output = str;
	if (!target.empty()) 
	{
		string::size_type pos = 0;
		while ((pos = output.find(target, pos)) != string::npos)
		{
			output.replace(pos, target.length(), replacement);
			pos += replacement.length();
		}
	}
	return output;
}
