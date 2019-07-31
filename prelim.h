#ifndef _prelim_h_included_
#define _prelim_h_included_

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <stdlib.h>

using namespace std;

	// *** Processing parameters *** 	//
	// material constants needed to update the betha variable 
	// units are length[m] ; viscosity[Pa.S]; permeability[m^2];
	const double size	= 2.0;			// Size of the laminate; unit: [m]
	const double d_time = 240.0;		// debulk time [min]
	const double ac_press = 1.0;		// Autoclave pressure [atm] NOTE: this pressue is 
	// Vacuum pressure is set to 1.0 if AC pressure > 1.0 // applied after debulk time (d_time)
	double p_v			= 0.0;			// pressrue at vacuum vavle or right boundary
	const double p_0	= 1.0;			// initial pressure
	const double L		= 1.0;			// dimentional length (NOTE: do not change)
	const double K_yy	= 1.0e-16;		// through-thickness permeability; resin transport
	const double K_xx	= 3.2e-15;		// in-plane permeability; resin transport
	const double mu_g	= 1.81e-05;		// gas (air) viscosity
	const double h		= 0.1e-03;		//thickness of the fiber-bed
	const double p0		= 101320;		// normal atmospheric pressure [Pa]
	const double phi0	= 0.46;			// fiber-bed porosity (phi0 = 1.0 - fiber_volume_fraction)
	const double p_c	= 0.00;			// capillary pressure
	const double init_imp	= 0.803;	// initial impregnation
	const double klin	= 0.1283/p_0;	// this should be nondimentionalized
										// Klinkenberg constant suggested by James Kay for MTM45-1
	// impregnation coeff or betta = betta_cst / (mu_r)
	const double betta_cst	= mu_g*K_yy*(L*L)/K_xx/(h*h)/phi0; 
	// time-constant; dimensional_time = time-constant * tau
	const double time_cst	= mu_g*(L*L)/K_xx/p0; 
	// *******************************	 //

	// ****** Model parameters ****** 	//
	const double dx		= 0.01;			// set space discritization
	const int nx		= size/dx;		// cell number
	const double t_max	= 1000.;		// time limit
	const double dt		= 0.0001;		// time step
	const double p_a_0	= 1.0;			// atmospheric pressure
	// ****************************** 	//

	void update_f(double aux_f[nx+2],  const double ap_a[nx+1], const double aomega[nx+1],	// update field variable
				  const double a_beta[nx+1], double time);
	void update_fi(const double auxf[nx+2], double auxfi[nx+1], const double ap_a[nx+1],	// update flux at faces
   				  const double aomega[nx+1], const double abetta[nx+1]);			 
	void update_omega (const double aux_f[nx+2], const double ap_a[nx+1],		// update cell impregnation
   					  double omega_aux[nx+1], const double abetta[nx+1]); 	
	void update_beta(const double time, const double debulk_time, double abetta[nx+1], 
					 double data [1000][3], const int auxnumb);// update impregnation coeff accor. viscosity data
	void set_gcells (double aux[nx+2], double time);// update ghost cells
	void SolveThomas(double LHS[nx+2][3], double RHS[nx+2],			// Thomas solver for a tri-diagonal system
					const int iSize);
	void generate_csv(const double auxf[nx+2], const double aux_omega[nx+1],
					  const double aux_res[nx+1],int s);
	int read_data(double aux_data[][3]);	// reads data from a file and stores in array;
											// aux_data[][0]:Temp; aux_data[][1]:DOC; aux_data[][2]:viscosity						
	double error_anly(const double aux1[nx+2],const double aux2[nx+2]);	// returns l2 error of solution convegence

#endif
