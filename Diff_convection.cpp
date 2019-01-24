/* DiffEqn.cpp  ---- CODE FOR SOLVING DIFFUSION CONVECTION EQUATION*/
#include "stdafx.h"
#include<iostream>
#include<string>
#include"modules.h"
#include"variables.h"

void diffEqn_50()
{
	
		for (int i = 0; i < nx; i++)
		{
			c[i][ny - 1] = c[i][ny - 2];
			c[i][0] = beta*c[i][1];
		}
		for (int i = 0; i < 20; i++)
		{
			c[i][20] = beta*c[i][21];
			c[i][0] = 0.0;
		}
		for (int i = 0; i < 21; i++)
		{
			c[20][i] = beta*c[21][i];
		}
		for (int i = 0; i < 21; i++)
		{
			c[60][i] = beta*c[59][i];
		}
		for (int i = 60; i < nx; i++)
		{
			c[i][20] = beta*c[i][21];
			c[i][0] = 0.0;
		}
		for (int i = 1; i < ny - 2; i++)
		{
			c[0][i] = c[1][i];
			c[nx - 1][i] = c[nx - 2][i];
		}
		for (int i = 1; i < nx - 1; i++)
		{
			for (int j = 1; j < ny - 1; j++)
			{
				if ((i >= 0 && j <= 20 && i < 20) || (i >= 60 && j <= 20 && i < nx))
				{
					tmp_c[i][j] = 0.0;
				}
				else
				{

					tmp_c[i][j] = ((dy*dy*(c[i + 1][j] + c[i - 1][j])) + (dx*dx*(c[i][j + 1] + c[i][j - 1])) - (2 * c[i][j] * ((dx*dx) + (dy*dy))))*D / (dx*dx*dy*dy)
						- ((u[i][j] * (c[i + 1][j] - c[i - 1][j]) / (2 * dx)) + v[i][j] * ((c[i][j + 1] - c[i][j - 1]) / (2 * dy)));
				}
				c[i][j] = c[i][j] + dt*tmp_c[i][j];
				for (int k = 20; k <= 60; k++)
				{

					
						rate[k] = (D*(c[k][1] - c[k][0]) / dy)*dx;
						mass[k] = mass[k] + dt*rate[k];
				
				}
				for (int k =0; k < 20; k++)
				{


					rate[k] = (D*(c[k][21] - c[k][20]) / dy)*dx;
					mass[k] = mass[k] + dt*rate[k];

				}
				for (int k = 61; k < nx; k++)
				{


					rate[k] = (D*(c[k][21] - c[k][20]) / dy)*dx;
					mass[k] = mass[k] + dt*rate[k];

				}
				for (int k = 0; k < 20; k++)
				{


					rate1[k] = (D*(c[21][k] - c[20][k]) / dy)*dx;
					mass1[k] = mass1[k] + dt*rate1[k];

				}
				for (int k = 20; k <40 ; k++)
				{


					rate1[k] = (D*(c[59][k] - c[60][k]) / dy)*dx;
					mass1[k] = mass1[k] + dt*rate1[k];

				}
				
			}
		}
	
}

