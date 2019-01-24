/* MAIN_CONTINUUM.CPP  ---  MAIN PROGRAM FOR THE CONTINUUM MODEL WHERE ALL THE PROCESS ARE MENTIONED IN DIFFERENT SUBROUTINES*/


#include "stdafx.h"
#include<iostream>
#include<string>
#include<fstream>
#include"variables.h"
#include<string>

//#define nx 200
//#define ny 80

double l = 0.0002, w = 0.00008, dt = 0.1;
double t = 0.0, t_lim = 1800.0;
double c[nx][ny], tmp_c[nx][ny], u[nx][ny], v[nx][ny];
double D = 4e-13;
int e = 0;
double dx = l / (nx - 1);
double dy = w / (ny - 1);
double beta = 0.9;
double rate[nx];
double mass[nx],mass1[40],rate1[40];
#include"modules.h"
std::string velocity = "velocity field";
std::string diff = "diffusion field";
int c1 = 0;
int main()
{
	
	double avg1 = 0.0;
	double avg = 0.0;
	char write;
	for (int i = 0; i < nx; i++)
	{
		rate[i] = 0.0;
		mass[i] = 0.0;
		for (int j = 0; j < ny; j++)
		{
			c[i][j] = 0.0001;
			tmp_c[i][j] = 1;
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
		c1 = i;
	}          
	std::cout << "avg1" << c1 << "\n";
	while (t < t_lim)
	{
		diffEqn_50();
		t = t + dt;


		std::cout << c[0][0] << "\t" << t << "\t" << mass[nx - 1] << "\n\n";
		for (int i = 0; i < nx; i++)
		{
			avg = mass[i] + avg;
		}
		//avg = avg / nx;
		for (int i = 0; i < 40; i++)
		{
			avg1 = mass1[i] + avg1;
		}
		//avg1 = avg1 / 40;
		//avg = (avg + avg1) / 2;
	}
	//std::cout <<"avg_50"<< avg << "\n";
	std::cout << "avg1" << avg1 << "\n";
	std::cout << "avg_50" << avg << "\n";
}

