#include "stdafx.h"
#include "var_surf.h"
#include<iostream>
#include<string>
#include<string>
#include<fstream>


void writeDat()
{
	
		std::ofstream outfile("final.dat", std::ios::app|std::ios::out);
		if (outfile.is_open())
		{

			
				
outfile<<"xi"<<"\t" <<xi <<"\n";
outfile<<"yi"<<"\t" <<yi <<"\n";
outfile<<"x_c"<<"\t" <<x_c <<"\n";
outfile<<"y_c"<<"\t" <<y_c <<"\n";
outfile<<"Ap"<<"\t" <<Ap <<"\n";
outfile<<"ra"<<"\t" <<ra <<"\n";
outfile<<"Amax"<<"\t" << Amax<<"\n";	 
outfile<<"xs"<<"\t" <<xs <<"\n";
outfile<<"ys"<<"\t" <<ys <<"\n";
outfile<<"D"<<"\t" << D<<"\n";	
outfile<<"Ap2"<<"\t" << Ap2<<"\n"; 
outfile<<"Amin"<<"\t" <<Amin <<"\n";
outfile<<"ri"<<"\t" << ri<<"\n";
outfile<<"h_pro"<<"\t" <<h_pro <<"\n";
outfile<<"h_dep"<<"\t" <<h_dep <<"\n";
outfile<<"h_no"<<"\t" << h_no<<"\n";
outfile<<"h2"<<"\t" <<h2 <<"\n";

			
		}
		else
		{
			std::
				cout << "\n\n\t\t\t" << "FILE NOT OPEN!!!" << "\n";
		}
		outfile.close();
	
}
