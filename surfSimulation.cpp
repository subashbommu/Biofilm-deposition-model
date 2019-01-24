
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>
#include <new>
#include "time.h"
#include "var_surf.h"
#define pi 3.1416
double normRand(double mean, double stdDev);

int surfSimulation()
{
	srand(time(0));
	nrmldist = new double[Na]; surface = new double[Na];
	while (sad_diff > desired_sad_diff)
	{
		count = 0.0;
		surf_total = 0.0;
		for (int i = 0; i < Na; i++)
		{
			nrmldist[i] = Rq*normRand(0, 1);
			surface[i] =  pi*pow(abs(nrmldist[i]), 2);  //2*
			surf_total = surf_total + surface[i];
			count += 1;
		}

		sim_sad_tmp = (surf_total / xy);

		Rq_tmp = 0.0; Ra_tmp = 0.0;
		for (int i = 0; i < Na; i++)
		{
			Ra_tmp = Ra_tmp + abs(nrmldist[i]);
			Rq_tmp = Rq_tmp + pow(abs(nrmldist[i]), 2);
		}

		Ra_tmp = Ra_tmp / Na;
		Rq_tmp = Rq_tmp / Na;
		Rq_tmp = sqrt(Rq_tmp);
		sa_tmp =  pi*pow(Ra_tmp, 2);  //2*
		if (sim_sad_tmp>sad)                 //CHECK WITH JUST NA_TMP=NA-1;
		{
			Na_tmp = Na - 1;
		}
		else if (sim_sad_tmp < sad)
		{
			Na_tmp = Na + 1;
		}


		if (Na_tmp < 0)
		{

			std::
				cout << "Program exit --> ASPERITY NUMBER NEGATIVE";

			return EXIT_FAILURE;
		}
		sim_sad = sim_sad_tmp;
		sad_diff = abs((sad - sim_sad) / sad);
		//if (sad_diff >= desired_sad_diff)     //CHANGE AND CHECK
		//{
			Na = Na_tmp;
			Ra=Ra_tmp;
			Rq=Rq_tmp;
		//}

	}

}
