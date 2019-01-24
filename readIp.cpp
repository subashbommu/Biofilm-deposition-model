#include "stdafx.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "var_surf.h"
using namespace std;
void readIp()
{
	vector <vector <string> > data;
	ifstream infile("txt1.dat");

	while (infile)
	{
		string s;
		if (!getline(infile, s)) break;

		istringstream ss(s);
		vector <string> record;

		while (ss)
		{
			string s;
			if (!getline(ss, s, '\t')) break;
			record.push_back(s);
		}

		data.push_back(record);
	}
	if (!infile.eof())
	{
		cerr << "Fooey!\n";
	}
//	double t[45][45];
	//t = stod(data[0][0]);
	//cout << t << '\n';
	for (int i = 0; i < x_grid_no; i++)
	{
		for (int j = 0; j < y_grid_no; j++)
		{
			z[i][j] = stod(data[i][j]);
			//cout << z[i][j] << '\t';
		}
		//cout << '\n';
	}
}
