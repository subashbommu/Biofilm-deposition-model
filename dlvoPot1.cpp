//edit change fion1 to 10.0 ----> Thu Jan 21 16:41:52 2016 
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <ctgmath>

#include <new>
#include "time.h"
#include "var_surf.h"
#define pi 3.141
double dlvoPot(double h)
{
	double tempf = 298.15;
	double epsvid1 = 8.85e-12;
	double epseaul = 78.5;
	double fion1 = 20.0;
	double phiel1 = -0.030;
	double phiel2 = -0.030;
	double LW_w = pow(0.0218, 0.5);
	double LW_s = pow(0.01293, 0.5);
	double LW_p = pow(0.02526, 0.5);
	double h0 = 0.158e-9;
	double D_LW = abs(2 * (LW_w - LW_s)*(LW_p - LW_w));
D_LW=0.000766031;
	double hamaker = 12 * pi*pow(h0, 2)*D_LW;
	double kb = 1.38065e-23;
	double ee = 1.6e-19;
	double rhop = 1050;
	double rhof = 998;
	double debye = pow(((epseaul*epsvid1*kb*tempf) / (2 * pow(ee, 2)*6.0221213e23*fion1)), 0.5);
	double kappa = 1 / debye;
	double C1 = hamaker / (12 * pi);
	double asp_vdw = -C1 / (pow(h, 2));
	double C2 = 0.5*epseaul*epsvid1*kappa;
	double C3 = pow(phiel1, 2) + pow(phiel2, 2);
	double C4 = 2 * phiel1*phiel2;
	double xt = h*kappa;
	double asp_el = C2*(C3*(1 - (1/tanh(xt))) + C4 / sinh(xt));
	double lambda_AB = 6e-10;
	double Rw_p = pow(0.0255, 0.5);
	double Rw_n = pow(0.0255, 0.5);
	double Rp_p = pow(0.000430, 0.5);
	double Rp_n = pow(0.004686, 0.5);
	double Rm_p = pow(0.0000563, 0.5);
	double Rm_n = pow(0.00307, 0.5);
	double delta_G = 2 * Rw_p*(Rm_n + Rp_n - Rw_n) + 2 * Rw_n*(Rm_p + Rp_p - Rw_p);
	delta_G = delta_G - 2 * (Rm_p*Rp_n + Rm_n*Rp_p);
	double asp_AB = delta_G*exp((h0 - h) / lambda_AB);
	double y = asp_el + asp_vdw + asp_AB;
	double z = y / (kb*tempf);
	return z;
}
