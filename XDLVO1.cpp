// XDLVO.cpp : Defines the entry point for the console application.
//1st edit  ----->  Wed Jan 20 10:22:54 2016 
//2nd edit  ----->  Thu Jan 21 10:48:59 2016 

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>
#include <new>
#include "time.h"
#include "var_surf.h"
#include<fstream>
#define pi 3.1416

double dlvoPot(double h);
double normRand(double mean, double stdDev);
int surfSimulation();
void readIp();

//double Rq;
double sad_diff = 1.0;
double surf_total = 0.0;
double xmax = 1500, ymax = 1500, sa_tmp, Na_tmp, sad = 0.498, Ra =195.0 , Rq = 245.0;
double xy = 4 * xmax*ymax, sa = pi*pow(Ra, 2);  //
int Na = sad*xy / sa;
double *nrmldist, *surface, sim_sad_tmp, sim_sad;
double Ra_tmp = 0.0, Rq_tmp = 0.0;
double desired_sad_diff = 0.00001;   // 0.01
int count = 0;
double surface_delta_x, surface_delta_y, **x, **y, **r_grid, **sgrid, **xgrid, **ygrid, **tmp_r, **z, **xc, **yc, ***U, ***U0;
void writeOp(double **data, int x_dimension, int y_dimension, std::string& d_type);
int x_grid_no;
int y_grid_no;
int peak = 0, valley = 0;
double Rmax_loc, Rmin_loc, Rmax = 0.0, Rmin = 0.0, Rm;
int Np = 1000;
std::string xp = "x plot";
std::string yp = "y plot";
std::string yq = "y original";
std::string xq = "x original";
std::string r = "r points";
std::string Udata = "U points";
std::string Uzero = "U_0 points";
std::string height = "height points-";
void writeOpd(double *data1, int size, std::string& d_type);

int main()
{

	
	surfSimulation();
	Rq = Rq_tmp;
	Ra = Ra_tmp;
	x_grid_no = (sqrt(Na - 1) + 1);
	y_grid_no = (sqrt(Na - 1) + 1);
	x = new double*[x_grid_no];
	y = new double*[x_grid_no];
	xc = new double*[x_grid_no];
	yc = new double*[x_grid_no];
	z = new double*[x_grid_no];
	U = new double**[x_grid_no];
	U0 = new double**[x_grid_no];
	r_grid = new double *[x_grid_no];
	tmp_r = new double *[x_grid_no];                  ////double U[45][45][1000];
	sgrid = new double *[x_grid_no];
	ygrid = new double *[x_grid_no];
	xgrid = new double *[x_grid_no];
	for (int i = 0; i < x_grid_no; i++)
	{
		x[i] = new double[y_grid_no];
		y[i] = new double[y_grid_no];
		xc[i] = new double[y_grid_no];
		yc[i] = new double[y_grid_no];
		U[i] = new double *[y_grid_no];
		U0[i] = new double *[y_grid_no];
		r_grid[i] = new double[y_grid_no];
		tmp_r[i] = new double[y_grid_no];
		sgrid[i] = new double[y_grid_no];
		xgrid[i] = new double[y_grid_no];
		ygrid[i] = new double[y_grid_no];
		z[i] = new double[y_grid_no];
	}
	for (int i = 0; i < x_grid_no; i++)
	{
		for (int j = 0; j < y_grid_no; j++)
		{
			U[i][j] = new double[Np];
			U0[i][j] = new double[Np];
		}
	}
	std::
		cout << Na << '\n' << x_grid_no << '\n' << y_grid_no;
	std::
		cout << "Ra" << Ra<<'\n' << "Sim sad" <<sim_sad<< '\n' << "Rq"<<Rq<<"\t";
	surface_delta_x = (1 / sqrt(Na - 1)) * 2 * xmax;
	surface_delta_y = (1 / sqrt(Na - 1)) * 2 * ymax;
	for (int i = 0; i < (y_grid_no); i++)
	{
		for (int j = 0; j < x_grid_no; j++)
		{
			x[i][j] = -xmax + (surface_delta_x*j);
			y[i][j] = -ymax + (surface_delta_y*i);

		}
	}


	/******************************************************************************************************************/
	for (int i = 0; i < x_grid_no; i++)
	{
		for (int j = 0; j < y_grid_no; j++)
		{
			r_grid[i][j] = Rq*normRand(0, 1);
			if (r_grid[i][j] >0.0)
			{
				peak = peak + 1;
				Rmax = r_grid[i][j];
				if (Rmax>Rmax_loc)
				{
					Rmax_loc = Rmax;
				}
			}
			if (r_grid[i][j] <0.0)
			{
				valley = valley + 1;
				Rmin = r_grid[i][j];
				if (Rmin<Rmin_loc)
				{
					Rmin_loc = Rmin;
				}
			}
			sgrid[i][j] = pi *pow((abs(r_grid[i][j])), 2);  //2*
			xgrid[i][j] = xmax*normRand(0, 1);
			ygrid[i][j] = ymax*normRand(0, 1);
			tmp_r[i][j] = r_grid[i][j];
		}
	}

	writeOp(xgrid, x_grid_no, y_grid_no, xp);
	writeOp(ygrid, x_grid_no, y_grid_no, yp);
	writeOp(tmp_r, x_grid_no, y_grid_no, r);
	writeOp(x, x_grid_no, y_grid_no, xq);
	writeOp(y, x_grid_no, y_grid_no, yq);


	std::cin.ignore();
	std::
		cout << "Press Enter to Continue";
	getchar();
	readIp();


	Rm = Rmax_loc - Rmin_loc;
	double delta = 15.0;
	double start = 0.05;
	double imax = (delta / start) + 1;
	double rp = 500.0;
	double rapp = 0.0;// rp / 10.0;
	double lapp = 0.0;// rp / 5.0;
	int bact_x = x_grid_no, bact_y = y_grid_no;
	for (int i = 0; i < (y_grid_no); i++)
	{
		for (int j = 0; j < x_grid_no; j++)
		{
			xc[i][j] = -xmax + (surface_delta_x*j);
			yc[i][j] = -ymax + (surface_delta_y*i);

		}
	}

	double Zc[1000];
	double mon_height = 15.0;
	for (int i = 0; i < Np; i++)
	{
		Zc[i] = (mon_height / (double)Np)*i;
	}
	int qw = Np;
	 writeOpd(Zc, qw, height);
	/*--------------------------------------------------------------------------------------------------------*/
	/*DISCRETISATION IN CYLINDRICAL COORDINATES*/
	int NI = 30.0;
	double rmax = rp;
	double rmin = 0.0;
	double dr = rmax / NI;
	double r[30];
	for (int i = 0; i < NI; i++)
	{
		r[i] =(rmax/(NI-1))*i;
	}
	double theta_min = 0.0;
	double theta_max = 2 * pi;
	double dtheta1 = (theta_max - theta_min) / (NI-1);
	double dtheta = (theta_max - theta_min) / (NI );
	double theta[30];

	double Uc[1000];
	double h2;
	double Uav[1000];
	double Uav_0[1000];
	double countD[1000];
	double countD_0[1000];
	int neg = 0;
	double peak = 0, valley = 0, flat = 0;
	for (int i = 0; i < NI; i++)
	{
		theta[i] = dtheta1*i;
	}
	double xi, yi, x_c, y_c, Ap, ra, Amax, h;
	double xs, ys, D, Ap2, h0, Amin, ri;

	for (int p = 0; p < Np; p++)
	{
		Uc[p] = 0.0;
		Uav[p] = 0.0;
		Uav_0[p] = 0.0;
		countD[p] = 0.0;
		countD_0[p] = 0.0;
		for (int m = 0; m < x_grid_no; m++)
		{
			for (int n = 0; n < y_grid_no; n++)
			{
				U[m][n][p] = 0.0;
				U0[m][n][p] = 0.0;

			}
		}
	}

	for (int p = 0; p < Np; p++)
	{

		std::cout << p << "\n";

		for (int n = 0; n < bact_x; n++)
		{


			for (int m = 0; m < bact_y; m++)
			{


				for (int i = 0; i < x_grid_no; i++)
				{
					for (int j = 0; j < y_grid_no; j++)
					{
						xi = x[i][j];
						yi = y[i][j];
						x_c = xc[n][m];
						y_c = yc[n][m];
						Ap = pow((xi - x_c), 2) + pow((yi - y_c), 2);
						
						ra = abs(z[i][j]);
						Amax = pow((rp + ra), 2);
					//	std::cout << "Ap" << "\t" << Ap << "\t" << "Amax" << Amax << "\n";
						h = 0.0;
						if (Ap < Amax)
						{
							
							//std::cout << "p" << "\t" << p << "\t" << xi << "\n" << yi << "\n" << x_c << "\n" << y_c << "\n";
							//getchar();
							for (int q = 0; q < NI; q++)
							{
								for (int s = 0; s < NI; s++)
								{
									xs = r[q] * cos(theta[s]);
									ys = r[q] * sin(theta[s]);
									D = Zc[p];
									Ap2 = pow((xi - xs - x_c), 2) + pow((yi - ys - y_c), 2);
									//Ap2 = pow((x_c + xs - xi), 2) + pow((y_c + ys - yi), 2);
									//h = -1.0;
									h0 = 0.158;
									Amin = pow((ra+h0 ), 2);
									ri = abs(pow((pow(yi, 2) + pow(xi, 2)), 0.5) - pow((pow((y_c + ys), 2) + pow((x_c + xs), 2)), 0.5));
									if (Ap2>Amin)
									{
										//std::cout << Ap2 << "\n";
										if (ri < ra)
										
										{
											
											if (z[i][j]>0.0)
											{
												peak = peak + 1;
												double h_pro = D + rp - pow((abs(pow(rp, 2) - pow(r[q], 2))), 0.5) - pow((abs(pow(ra, 2) - pow(ri, 2))), 0.5);
												h = h_pro;
												if (h < 0.0)
												{
													h = -1;

												}
											}
											else 
											{
												valley = valley + 1;
												double h_dep = D + rp - pow((abs(pow(rp, 2) - pow(r[q], 2))), 0.5) + pow((abs(pow(ra, 2) - pow(ri, 2))), 0.5);
												h = h_dep;
												if (h < 0.0)
												{
													h = -1;
												}
											}
											
										}
										
										else //if (ri<ra)
										{
											flat = flat + 1;
											double h_no = D + rp - pow((abs(pow(rp, 2) - pow(r[q], 2))), 0.5);
											h = h_no;
											if (h < 0)
											{
												h = -1;
											}
										
									}
									
										//	std::cout << h << "\n";
									}
									//std::cout << "--------------------------" << "\n";
									//std::cout <<h << "\n";
									//std::cout << ra << "\n";
									//std::cout << r[q] << "\n";
									//std::cout << D << "\n";
									//std::cout << rp << "\n";
									double ht = h;
									if (h > 0.0 )
									{
										
										h = h*1e-9;
										
										U[m][n][p] = U[m][n][p] + dlvoPot(h)*r[q] * 1e-9*dr*1e-9*dtheta;
										Uc[p] = Uc[p] + dlvoPot(h)*r[q] * 1e-9*dr*1e-9*dtheta;
										//std::cout << U[m][n][p] << "\n";
										/*if (p>1 &&h > 30e-10 && h < 200e-10)
										{

											std::cout << U[m][n][p] << "\n";
											getchar();
										}*/
										
									}
									/*if (h<1e-20)
									{
										getchar();
									}*/
									
									h2 = D + rp - pow((abs(pow(rp, 2) - pow(r[q], 2))), 0.5);
									if (h2 > 0.0)
									{

										h2 = h2*1e-9;

										U0[m][n][p] = U0[m][n][p] + dlvoPot(h2)*r[q] * 1e-9*dr*1e-9*dtheta;

									}
									//else
									//{
									//	 neg=neg+1;
									// }
								}
							}
						}
					}
				}
			}
		}
	}

	/*-------------------------------------------------------------------------------------------------------------------------------------------------*/
	//**x, **y, **r_grid, **sgrid, **xgrid, **ygrid, **tmp_r, **z, **xc, **yc
	delete[] x, y, r_grid, sgrid, xgrid, ygrid, tmp_r, z, xc, yc;

	for (int p = 0; p < Np; p++)
	{
		for (int n = 0; n < bact_x; n++)
		{
			for (int m = 0; m < bact_y; m++)
			{
				//std::cout << "1   " << U[m][n][p] << "\n";
				if (abs(U[n][m][p])>1e-15 && U[n][m][p]>-20)
				{
					Uav[p] = Uav[p] + U[n][m][p];

					countD[p] = countD[p] + 1;
					//std::cout <<"11   "<< Uav[p] << "\n";
				}
				//std::cout << "1   " << U0[m][n][p] << "\n";
				if (abs(U0[n][m][p])>1e-15 && U0[n][m][p]>-20)
				{
					Uav_0[p] = Uav_0[p] + U0[n][m][p];
					countD_0[p] = countD_0[p] + 1;
					if (p == 930 || p == 950)
					{
						std::cout << U0[n][m][p] << "\t" << p << "\t" << Uav_0[p] << "\n";
					}

				}

			}
		}
	}
	for (int p = 0; p < Np; p++)
	{
		if (countD[p]>0.0)
		{
			Uav[p] = Uav[p] / countD[p];
			//std::cout << Uav[p] << "\n";
		}
		if (countD_0[p] > 0.0)
		{
			Uav_0[p] = Uav_0[p] / countD_0[p];
		}

	}
	std::cout << neg;
	std::cout << "peak"<<peak << "\n";
	std::cout << "valley"<<valley << "\n";
	std::cout << "flat"<<flat << "\n";
	writeOpd(Uav, qw, Udata);
	writeOpd(Uav_0, qw, Uzero);




}
