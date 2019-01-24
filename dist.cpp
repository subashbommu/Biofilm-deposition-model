/*NORMRAND.CPP --- FUNCTION OF MAIN_CONTINUUM.CPP TO CALCUATE NORMALLY DISTRIBUTED RANDOM VARIABLES FOR SURFACE RECREATION*/
/*INPUT ARGUMENTS --- MEAN AND STANDARD DEVIATION*/
/*OUTPUT VARIABLES --- NORMAL RANDOM NUMBER*/

#include "stdafx.h"
#include<iostream>
#include<string>
#include "time.h"
#include<math.h>
double normRand(double mean, double stdDev) {
	double u, v, s;
	
	do {
		u = ((double)rand() / (double)RAND_MAX) * 2.0 - 1.0;
		v = ((double)rand() / (double)RAND_MAX) * 2.0 - 1.0;
		s = u * u + v * v;
	} while (s >= 1 || s == 0);
	double mul = sqrt(-2.0 * log(s) / s);

	return mean + stdDev * u * mul;
}
