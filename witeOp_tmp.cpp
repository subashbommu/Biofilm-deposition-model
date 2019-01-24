/* FUNCTION OF MAIN_CONTINUUM.CPP FOR WRITING THE OUTPUT 2D ARRAYS IN A FILE --  FILE SHOULD BE ENTERED WITH EXTENSION AND LOCATION IF NEEDED*/

#include "stdafx.h"
#include<iostream>
#include<string>
#include<string>
#include<fstream>

#include "var_surf.h"
void writeOp(double **data, int x_dimension, int y_dimension, std::string& d_type)
{

	std::string name;
	std::
		cout << "Enter the output file name for " << d_type << "\t";
	std::
		cin >> name;
	std::ofstream outfile(name);
	if (outfile.is_open())
	{

		for (int i = 0; i < y_grid_no; i++)
		{
			for (int j = 0; j < x_grid_no; j++)
			{
				outfile << data[i][j] << "\t";

			}
			outfile << "\n";

		}

	}
	else
		std::
		cout << "\n\n\t\t\t" << "FILE NOT OPEN!!!" << "\n";
	outfile.close();
}
