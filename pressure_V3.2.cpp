/** Please cite: Mohseni M, Zobeiry N, Fernlund G. Journal of Reinforced Plastics
 and Composites. 2019 Jul 19:0731684419865783.**/
/** coupled gas/resin flow simulation using FVM and Implicit Euler **/
/** V3: Klinkenberg effect has been added	**/

/** defaul format for material data file:

	- every line represents 1 minute time progress
	- evey row contains 4 columns:
		"Time(min)"	"Temperature[c]"	"Degree of Cure" "Resin Viscosity (Pas)"
 **/

#include "prelim.h"			// definition of all constants, global variables and functions
 
 ///////////////***** M  A  I  N *****///////////////
int main ()
{
	FILE *out2, *out3, *out4;
	out2 = fopen ("Report.txt", "wt");
	fprintf (out2,"%s\t\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n\n","time", \
	"Di-time(min)","Porosity(%)", "Average press.", "Vac.Pres(Pa)" , "Atm.Pres(Pa)"\
	, "Betha_coeff", "Convergence");
 
	out3 = fopen ("Pressure.txt", "wt");
	out4 = fopen ("Porosity.txt", "wt"); 
 
	double t = 0.;					// time variable 
	// assign and initialize global variables 
 	double f[nx+2];		std::fill_n(f, nx+2, p_0);			// set domain variables+2 ghost cells 
 	double p_a[nx+1];	std::fill_n(p_a, nx+1, p_a_0);		// set atmospheric pressure field 
 	double omega[nx+1];	std::fill_n(omega, nx+1, init_imp);	// impregnation coefficient 
															//## 1 cell is added just to count from 1 to nx
	double betaa[nx+1];										// source term (squeeze flow) coefficient
															// defined as vector in case spacial variation exists  

	//*// comment this part if material file available
	//*// otherwise viscosity is set to a constant
	double const_viscosity = 10000.0;
	std::fill_n(betaa, nx+1, betta_cst/const_viscosity);


	double data [1000][3];
	int lnumb;			    // records the number of lines in the input data file (equivalent to time in min)
	int q =0;				// counter of step number
	lnumb = read_data(data);

	while (t < t_max)				/////**** IMPLICIT EULER ****/////
	{
		if (int((t*time_cst)/60.) > d_time)
		{
			if (ac_press > 1.0)	// if autoclave pressure (hgiher that 1 atm)
				p_v = 1.0;		// disconnect vacuum, i.e. Pv = 1 atm
			std::fill_n(p_a, nx+1, ac_press);
		}
		//*// uncomment if material file available
		//*// otherwise visocisty and betha paramater remain constant
		//update_beta (t,d_time,betaa,data,lnumb);
		double pre_f [nx+2] = { };
		for (int i=1; i < nx+1; i++) pre_f [i] = f[i];	// record previous solution for error analysis
		update_f (f,p_a,omega,betaa,t);					// update field variable - advance 1 time step		
		update_omega (f,p_a,omega,betaa);				// updates impregnation of the cell at the current times step
														// inputs are resin viscosity and pressure distribution
		double phi = 0.0;								// stores average impregnation
		double ave_pressure = 0.0 ;						// stores average pressure
		for (int i=1; i < nx+1; i++){
			phi += omega[i]/(nx); // calcualtes the average poroisyt throughout the length of the laminate
			ave_pressure += f[i]/nx;
		}
		
		fprintf (out2,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %2.5e \n",t, \
		(t*time_cst)/60., (1.-phi)*100, ave_pressure, p_v, p_a[20], betaa[20], error_anly (pre_f,f));		// Print report
 
		t += dt ;   // update time
 
		q++;
		if (q%100 == 0)
		{
			fprintf (out3,"%f\t ", ((t*time_cst)/60.));
   			fprintf (out4,"%f\t ", ((t*time_cst)/60.));
			double gas_resid [nx+1];
   			for (int i=1; i < nx+1; i++)
			   {
			   		double aux  = f[i];
			   		double aux1 = (1.0- omega[i])*100;
   					fprintf (out3,"%f\t ", aux);
   					fprintf (out4,"%f\t ", aux1);
			   }
			fprintf (out3,"\n ");
   			fprintf (out4,"\n ");  
		}
 
		if ((int((t * time_cst)/60.) - d_time) > lnumb){
			if (error_anly (pre_f,f) < 1.0e-8) cout << "Convergence reached ..." << endl;
			phi=0;
			for (int i=1; i < nx+1; i++) phi += omega[i]/(nx);
			cout << (1-phi)*100 << endl;
			break;
		}
	}
    double gas_resid [nx+1];
    for (int i=1; i < nx+1; i++) gas_resid[i] = (1.- omega[i])*f[i];
	generate_csv(f,omega,gas_resid,t);

    fclose (out2);
    fclose (out3);
    fclose (out4);
	return 0;
}

///////////////***** UPDATE VARIABLES IN TIME ******///////////////CCCCCC
void update_f (double aux_f[nx+2],  const double ap_a[nx+1], const double aomega[nx+1], 
			   const double a_beta[nx+1], double time)
{
	double fi[nx+1] = { }; 						// variable to record f flux integral
	set_gcells(aux_f,time);						// update ghost cells 
	update_fi (aux_f,fi,ap_a, aomega,a_beta);	// update flux intergal							
	double mrhs [nx+2] ;						// b (in Ax=b) as a vector
   	mrhs [0] = mrhs [nx+1] = 0.;				// impose bdc (ghost cells)
   	double mlhs [nx+2][3];						// update tridiagonal matrix A (in Ax=b) as 3 vetors
   	mlhs[0][1] = -1. ; mlhs[0][2] = 1.; 		// impose Neumaan bdc on left bdr
	mlhs[nx+1][0] = 1.; mlhs[nx+1][1] = 1.;		// impose Dirichlet bdc on right bdr

	for (int i=1; i < nx+1; i++)
	{	
		double n_term = aux_f[i]*(1.0+(klin/aux_f[i]));	// non-linear term including the Klinkenberg term
		mrhs [i] = fi[i]; //### source term is added when updating flux integral
   		mlhs [i][0] = -(n_term/dx/dx);
   		mlhs [i][1] = ((1. - aomega[i])/dt) + (2.*n_term/dx/dx) - (a_beta[i]*(ap_a[i] - (2.*aux_f[i])+p_c)/aomega[i]);
   		mlhs [i][2] = -(n_term/dx/dx);	
	}

	SolveThomas(mlhs, mrhs,nx+2); 				// Thomas solver
	for (int i=1; i < nx+1; i++){
		aux_f[i] += mrhs [i];
	}
}

 ////////////////****** UPDATE GHOST CELLS *****////////////////
void set_gcells(double aux[nx+2], double time)
{
	aux[0] = aux[1];					// left bdc = Newmann - fully-developed
	aux[nx+1] = (2.*p_v) - aux[nx];		// right bdc = Dirichlt - constant
}

 ////////////////****** UPDATE FLUX *****////////////////
void update_fi (const double auxf[nx+2], double auxfi[nx+1], const double ap_a[nx+1], 
				   const double aomega[nx+1], const double abetta[nx+1])
{
	for (int i=1; i < nx+1; i++)
	{
		double upper = auxf [i+1];						// auxf is defined as constant in this function
		double lower = auxf [i-1];
		if (aomega[i+1] == 1.0) upper = auxf[i];		// impose Newmann bdc if right cell is fully-impregnated
		if (aomega[i-1] == 1.0) lower = auxf[i];		// impose Newmann bdc if left cell is fully-impregnated
		double n_term = auxf[i]*(1.0+(klin/auxf[i]));	// non-linear term including the Klinkenberg term
		auxfi[i] = (n_term*((upper - 2.*auxf[i] + lower)/dx/dx))
				  + (abetta[i]*auxf[i]*(ap_a[i] - auxf[i])/aomega[i]);

		if (aomega[i] == 1.0)	auxfi[i] = 0.0;			// flux is equal to zero if the cell is fully-impregnated
	}
}

////////////////****** UPDATE CELL IMPREGNATION *****////////////////
void update_omega (const double aux_f[nx+2], const double ap_a[nx+1], 
					double omega_aux[nx+1], const double abetta[nx+1])
{
	for (int i=1; i < nx+1; i++) 
	{
		omega_aux[i] += dt*abetta[i]*(ap_a[i] - aux_f[i]+p_c)/omega_aux[i];
		if (omega_aux[i] > 1.) omega_aux[i] = 1.;
	}
}

////////////////****** UPDATE IMPREGNATION COEFFICIENT *****////////////////
void update_beta(const double time, const double debulk_time, double abetta[nx+1], 
				  double data [1000][3], const int auxnumb)
{
	// Assumption: data in material data file
	// is recorded one line per minute, starting with the time = 0 min

	double di_t = (time * time_cst) / 60.;		// dimentional time -  unit : minute
	double aux_beta = 0.0;
	
	if (di_t < debulk_time){  // assign room_temp viscosity (first line of the material file) if debulk is running
	   aux_beta = betta_cst / data [0][2];
	}
	else{
		double aux_time = di_t - debulk_time;
		aux_beta = betta_cst / data [int(aux_time)][2];
	} 
	// if simulation time exceeds recorded time in data_input file
	if (int(di_t-debulk_time) > auxnumb-1) aux_beta = betta_cst / data [auxnumb-1][2]; 
	std::fill_n(abetta, nx+1, aux_beta);		// assuming squeeze flow coeff is spatially constant
}

 ////////////////***** THOMAS SOLVER ******** /////////////////
void SolveThomas(double LHS[nx+2][3], double RHS[nx+2], const int iSize)
{
  int i;
  // This next line actually has no effect, but it -does- make clear that
  //  the values in those locations have no impact. 
  LHS[0][0] = LHS[iSize-1][2] = 0;
  // Forward elimination //
  for (i = 0; i < iSize-1; i++)
  {
    LHS[i][2] /= LHS[i][1];
    RHS[i] /= LHS[i][1];
    LHS[i+1][1] -= LHS[i][2]*LHS[i+1][0];
    RHS[i+1] -= LHS[i+1][0]*RHS[i];
  }
  RHS[iSize-1] /= LHS[iSize-1][1];// Last line of elimination //
  // Back-substitution //
  for (i = iSize-2; i >= 0; i--)    RHS[i] -= RHS[i+1]*LHS[i][2];
}

 ////////////////***** OUTPU CSV FILE ********/////////////////
void generate_csv(const double auxf[nx+2], const double aux_omega[nx+1],const double aux_res[nx+1],int s)
{
	FILE *out1;
    char name[50];
	int m = 0;
	m = (int) (s);
 
    sprintf (name, "Results_%d.dat",m);
    out1 = fopen(name, "wt");
	for (int i = 1; i < nx+1; i++){
	    	fprintf (out1,"%d\t %2.10e\t %2.10e\n",i, auxf[i], aux_omega[i]);
	    }
    fclose (out1);
}

 ////////////////***** READ DATA FORM FILE, STORE IN ARRAY ********/////////////////
int read_data(double aux_data[1000][3])
{
	ifstream data;
    data.open("data_nodebulk.txt");
    int x = 0;
	std::string line;
	while (std::getline(data, line))
	{
       istringstream ss(line);
   		double a, b, c, d; 	
		while( ss >> a >> b >> c >>d)
		{
			aux_data [x][0] = b; 
			aux_data [x][1] = c;
			aux_data [x][2] = d;
			++x;
		}
	}
	return x;
}

 ////////////////***** SOLUTION CONVERGENCE - L2 ERROR ********/////////////////
double error_anly(const double aux1[nx+2],const double aux2[nx+2])
{
	double l2 = 0.;
	for (int i = 1; i<nx+1; ++i) 
		l2 += (aux1[i] - aux2[i])*(aux1[i] - aux2[i]);
	return sqrt(l2);
}
