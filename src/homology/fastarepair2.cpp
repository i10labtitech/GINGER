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
#include "function.hpp"

using namespace std;

int main(int argc, char*argv[])
{

	if(argc != 3)
	{
		cout << "error" << endl;
		cout << "./fastarepair2 [prerepair.faa] [output]" << endl;
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

	string zero;

	string lin;

	string tmpstr, fas, cutfas, cutfas2;
	int length;

	while(getline(ifs,lin) )
	{
		string::size_type p1 = lin.find(">");
		if(p1 != string::npos)
		{
			if(fas != zero)
			{

				cutfas = Erace(fas," ");
				cutfas2 = Erace(cutfas,"*");
				cutfas = Erace(cutfas2,"B");
				cutfas2 = Erace(cutfas,"J");
				cutfas = Erace(cutfas2,"O");
				cutfas2 = Erace(cutfas,"U");
				cutfas = Erace(cutfas2,"X");
				cutfas2 = Erace(cutfas,"Z");
				cutfas = Erace(cutfas2,"#");
				fout << cutfas << endl;
			}
			string::size_type p2 = lin.find(" ");
			if(p2 != string::npos)
			{
				tmpstr = lin.substr(0, p2);
				fout << tmpstr << endl;
			}
			else
			{
				fout << lin << endl;
			}
			fas = zero;
		}
		else
		{
			fas += lin;
		}
	}
	cutfas = Erace(fas," ");
	cutfas2 = Erace(cutfas,"*");
	cutfas = Erace(cutfas2,"B");
	cutfas2 = Erace(cutfas,"J");
	cutfas = Erace(cutfas2,"O");
	cutfas2 = Erace(cutfas,"U");
	cutfas = Erace(cutfas2,"X");
	cutfas2 = Erace(cutfas,"Z");
	cutfas = Erace(cutfas2,"#");
	fout << cutfas << endl;
	return 0;
}
