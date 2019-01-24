#include "stdafx.h"
#include<iostream>
#include<string>
#include<string>
#include<fstream>

#include "var_surf.h"
void writeOpd(double *data1, int size, std::string& d_type)
{
	char t1;
	


		std::string name=d_type;
		//std::
		//	cout << "Enter the output file name for " << d_type << "\t";
		//std::
		//	cin >> name;
		std::ofstream outfile(name);
		if (outfile.is_open())
		{

			for (int i = 0; i < size; i++)
			{

				outfile << data1[i] << "\n";

			}
			t1 = 'N';
		}
		else
		{
			std::
				cout << "\n\n\t\t\t" << "FILE NOT OPEN!!!" << "\n";
		}
		outfile.close();
	
	
}
