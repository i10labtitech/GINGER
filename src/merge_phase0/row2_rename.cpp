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

using namespace std;

vector<string> Split(const string &s, char delim);
string Replaced_V(string str, string replacestr, int replace);

int main(int argc, char*argv[])
{

	if(argc != 4)
	  {
	    cout << "error" << endl;
	    cout << "./a.out [input.gff3] [output.gff3] [2列目に入れたい文字列]" << endl;
	    return 0;
	  }

	ifstream ifs(argv[1]);
	if(ifs.fail())
	  {
		cout << "error:fin file not open" << endl;
		return 0;
	  }

	ofstream fout(argv[2]);
	if(fout.fail())
	  {
	    cout << "error:fout file not open" << endl;
	    return 0;
	  }
	
	string row2;
	row2 = argv[3];
	
	string zero;
	
	string lin, tmplin;
	

	while(getline(ifs,lin) )
	  {
	    if(lin=="") continue;
	    tmplin = Replaced_V(lin, row2, 2);
	    fout << tmplin << endl;
	  }
	return 0;
}


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
