#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>







// DEFINITIONS
// -----------


#define		RES_DT		0.02
#define		RES_FINT	20.0
#define		C		2.997956376932163e+08
#define		R_E		6370000.0
#define		H_IONO		1E5
#define		PI		3.141592653589793115997963468544185
#define		D2R		PI/180.0
#define		M_EL		9.1E-31
#define		E_EL		5.105396765648739E5
#define		NAPe		2.71828182845905

#define		DE_EXP		0.003
#define		E_EXP_BOT	1.477
#define		E_EXP_TOP	7.477
#define		NUM_E		2200
#define		E_MIN		0 // min E for N,Q calc [in keV]
#define		SQUARE		0 // pitch angle distribution

#define		LONG_BOT	0
#define		DLONG		0.5
#define		LONG_TOP	30






// GLOBAL VARIABLES:
// -----------------

double		L_TARG;
int		NUM_LONGS;
int		NUM_TIMES;






// Constants to perform Gauss quadrature integration

double t5[] =	{	-0.9061798459, 
			-0.5384693101, 
			0, 
			0.9061798459, 
			0.5384693101	};

double beta5[] = {	0.2369268851,
			0.4786286745, 
			0.56888888889, 
			0.2369268851,
			0.4786286745	};







// PROTOTYPES
// ----------

float *getArr(int len1, int len2);
void compFlux(float *arr, double L, int k, char *dir);
void readJ(float J[][100]);
double getJdiff(float J[][100], double E, double alpha_lc);
double reduceFactor(double dlong);






/*
 * FUNCTION: main
 * --------------
 * This is an upgrade of calcFluxE.  It now reads AE8 
 * (L-dependent) flux values converts into differential 
 * flux, and assumes "sin" shaped loss-cone.
 *
 */

int main(int argc, char *argv[])
{
  FILE *inPtr, *alphaPtr;
  int i, k, m, nin, ei, ti;
  char *dir, sysCmd[512], filename[64], *NS, alphaFile[128];
  float L, *arr, *arrarr[128]; // <- array of pointers to precip arrays
  float J[100][100];
  double Jext;


  // Get directory name from command line
  if(argc != 3) {
    printf("\n\aWrong number of input arguments\n");
    exit(0);
  }
  dir	= argv[1];
  L_TARG = atof(argv[2]);


  NUM_TIMES = (int)ceil( RES_FINT / RES_DT ); 
  sprintf( alphaFile, "%s/alpha_%g_%s", dir, L_TARG, "N" );
  arr = getArr(NUM_E, NUM_TIMES);

  if(  (alphaPtr=fopen(alphaFile, "r"))!=NULL  ) {
    
    nin = fread(arr,sizeof(float),(NUM_E*NUM_TIMES),alphaPtr);
    fclose(alphaPtr);
    arrarr[0] = arr;
    
    
    sprintf( alphaFile,   "%s/alpha_%g_%s", dir, L_TARG, "S" );
    if((alphaPtr=fopen(alphaFile, "r"))==NULL) 
      printf("\n\aprob opening %s\n", alphaFile);  
    arr = getArr(NUM_E, NUM_TIMES);
    nin = fread(arr,sizeof(float),(NUM_E*NUM_TIMES),alphaPtr);
    fclose(alphaPtr);
    arrarr[1] = arr;
    
    
  } else {

    printf("\n%s cannot be opened\n", alphaFile);

  }	// if alphaFile exists
  
  
  for(k=0; k<2; k++) {
    compFlux(arrarr[k], L_TARG, k, dir);
  } // N/S - hemisphere



  return 0;
}










/* 
 * FUNCTION: getArr
 * ----------------
 * This function simply allocates dynamically a block of memory the 
 * size of NUM_E * NUM_TIMES of type float, and initializes it.
 * It returns the pointer to the memory.
 *
 */
float *getArr(int len1, int len2)
{
  float *arr;
  int ei, ti;

  arr = (float *) malloc( len1 * len2 * sizeof(float) );
  if(arr == NULL) {
    printf("\nProb assigning mem in calcFlux\n");
    exit(0);
  }

  for(ei=0; ei<len1; ei++) {
    for(ti=0; ti<len2; ti++) {
      arr[ei*len2+ti] = 0.0;
    }
  }

  return arr;

}










/*
 * FUNCTION: compFlux
 * ------------------
 * This function will open a file with the appropriate hemisphere and
 * L-shell inserted into the name, calculate the precipitated flux  
 * due to each cell in the dAlpha_RMS matrix and write it into the 
 * file.
 *
 */

void compFlux(float *arr, double L, int k, char *dir)
{
  FILE *QPtr, *NPtr;
  int ei, ti, i, li, nout, iofEbot;
  double mag, P, Q, epsm, alpha_eq, v, crunch, gamma;
  double I=0.0, x, g, field, Phi_p, b=1.0, vcm, fcgs;
  double v_tot_arr[NUM_E], E_tot_arr[NUM_E], dE_arr[NUM_E];
  double Jdiff[NUM_E], long_reduce;
  float alpha, J[100][100], *arrtmp, *Qarr, *Narr;
  char *NS, QFile[128], NFile[128];
  


  // Open up Phi, Q, N file for writing
  if(k==0) {NS = "N";} else {NS="S";}  

  sprintf( QFile,   "%s/QLong_%g_%s", dir, L, NS );
  sprintf( NFile,   "%s/NLong_%g_%s", dir, L, NS );

  printf("Calculating Q, N for L = %g, E > %g", L, E_MIN);
  NUM_LONGS = (int) (LONG_TOP - LONG_BOT)/DLONG;
  if(E_MIN < 1e-3) {
    iofEbot = 0;
  } else {	// remember E_MIN is in KeV
    iofEbot = floor((log10(E_MIN*1000)-E_EXP_BOT)/DE_EXP);
  }
  arrtmp = getArr(NUM_E, NUM_TIMES);  
  Qarr = getArr( NUM_LONGS, NUM_TIMES );
  Narr = getArr( NUM_LONGS, NUM_TIMES );
  

  if( (QPtr=fopen(QFile, "w"))==NULL ) {
    printf("\nProblem opening %s\n", QFile);
    exit(0);
  }
  
  if( (NPtr=fopen(NFile, "w"))==NULL ) {
    printf("\nProblem opening %s\n", NFile);
    exit(0);
  }
  
  
  epsm = (1/L)*(R_E+H_IONO)/R_E;

  crunch	= sqrt(1+3*(1-epsm))/pow(epsm,3) ;
  alpha_eq	= asin(sqrt( 1/crunch ));
  readJ(J);


  // Precalculate energy and velocity values
  for(i=0; i<NUM_E; i++) {

    //Energy in eV
    E_tot_arr[i] = pow(10, (E_EXP_BOT+(DE_EXP/2)+DE_EXP*i) );

    // Energy differential dE in keV
    dE_arr[i] = 1e-3 * pow(10, (E_EXP_BOT+ (DE_EXP/2))) * 
      exp(DE_EXP*i / log10(NAPe)) * DE_EXP / log10(NAPe);

    // Differential flux @ this Energy
    Jdiff[i] = getJdiff( J, E_tot_arr[i], alpha_eq );

    // Velocity corresponding to this Energy
    v_tot_arr[i] = C*sqrt(1 - pow( (E_EL/(E_EL+E_tot_arr[i])) ,2) );
  }


  printf("\n");
  for(li=0; li < NUM_LONGS; li++) {  // loop over longitudes

    long_reduce = reduceFactor( li*DLONG  );
    printf("Now doing long: %g\n", li*DLONG);

    for(ei=0; ei<NUM_E; ei++) {
      
      if(SQUARE) {
	v = v_tot_arr[ei];
	vcm = v*100;	// v in cm for distrib fn calculation
	gamma = 1/sqrt( 1 - v*v/(C*C) );
	fcgs =	4.9e5/pow( (vcm*gamma) ,4) - 
	  8.3e14/pow( (vcm*gamma) ,5) + 
	  5.4e23/pow( (vcm*gamma) ,6);
     
	// fcgs = 7.034e26 / pow(vcm,6); //10^8/E^2 distribution
	
	b = (v*v/M_EL)*pow( sqrt(1 - (v*v)/(C*C)), 3) * 1.6e-8 * fcgs;
	
	// b = 1e8 / pow(E_tot_arr[i],2);
      } else {
	b = Jdiff[ei]*1000;
      }
      
      
      for(ti=0; ti<NUM_TIMES; ti++) {
	
	alpha = long_reduce * arr[ei*NUM_TIMES+ti] ;

	mag = 1.4142135623731*alpha;	// sqrt(2)*alpha_RMS = peak
	
	P = mag/2;		//[ alpha_lc - (alpha_lc-mag) ] /2
	Q = alpha_eq - mag/2;	//[ alpha_lc + (alpha_lc-mag) ] /2
	
	I = 0.0;

	if(mag > 1e-8 && mag < 1) {
	  for(i=0; i<5; i++) {
	    x = P*t5[i] + Q ;
	    
	    if(SQUARE) {
	      g = (P/PI)*sin(2*x)*(  asin((x-alpha_eq)/mag)+ (PI/2) );
	    } else {
	      g = (P/PI)*sin(2*x)*((x - alpha_eq) *
				   ( asin((x-alpha_eq)/mag)+ (PI/2) ) +
				   sqrt(mag*mag-pow((x-alpha_eq),2)));
	    }
	    I += ( beta5[i]*g );
	    
	  } // for(i ... ) -> Gauss quad integration
	} // if mag != 0
	Phi_p = PI*crunch*b*I;
		
	if(Phi_p<0) {
	  printf("\nPhi_p<0, at ti: %d, ei: %d, li: %d, NS: %s, L: %g\n",
		 ti,ei,li, NS, L_TARG);
	  printf("PI: %g, crunch: %g, b: %g, I: %g, alpha: %g\n",
		 PI, crunch, b, I, alpha);
	  Phi_p = 0;
	}

	if(  isnan(Phi_p)  ) {
	  printf("\nPhi_p isnan, at ti: %d, ei: %d, li: %d, NS: %s, L: %g\n",
		 ti,ei,li, NS, L_TARG);
	  printf("PI: %g, crunch: %g, b: %g, I: %g, alpha: %g\n",
		 PI, crunch, b, I, alpha);
	  Phi_p = 0;
	}


	// Now do the integration to get Q and N
	if(ei >= iofEbot) {
	  Qarr[li*NUM_TIMES+ti] += (float)( Phi_p * 
					    E_tot_arr[ei] * 
					    dE_arr[ei] * 
					    1.602e-9);
	  Narr[li*NUM_TIMES+ti] += (float)( Phi_p * 
					    dE_arr[ei]); //eV->keV
	}
	
      } // for(ti ... )
    } // for(ei ... )        
  }  //li


  nout=fwrite(Qarr, sizeof(float), (NUM_LONGS*NUM_TIMES), QPtr);
  if(nout!=NUM_LONGS*NUM_TIMES) printf("\n\aProblem writing Q\n");
  nout=fwrite(Narr, sizeof(float), (NUM_LONGS*NUM_TIMES), NPtr);
  if(nout!=NUM_LONGS*NUM_TIMES) printf("\n\aProblem writing N\n");

  fclose(QPtr);
  fclose(NPtr);
}





/*
 * FUNCTION: readJ
 * ---------------
 * This function simply looks for a file EQFLUXMA.dat and reads 
 * it in.  The columns are energies and the rows are L-shells.
 * The first column is just a list of L-shells and the first row 
 * is just a list of energies (in MeV).
 *
 */
void readJ(float J[][100])
{
  char *filename;
  FILE *filePtr;
  int i,j;

  filename = "/home/TIPER/1_WIPP_code_BB/runfiles/EQFLUXMA.dat";

  if( (filePtr = fopen( filename ,"r")) == NULL ) {
    printf("Hey buddy, I can't open the file!\n");
    exit(1);
  }
  

  // INITIALIZE
  for(i=0; i<100; i++) {
    for(j=0; j<100; j++) {
      J[i][j] = 0.0;
    }
  }

  // READ IN VALUES
  for(i=0; i<47; i++) {
    for(j=0; j<31; j++) {
      fscanf(filePtr, "%e", &(J[i][j]));
    }
  }

  fclose(filePtr);
}






/*
 * FUNCTION: getJdiff
 * ------------------
 * Using the AE8 data stored in J, calculate the differential flux 
 * by taking the (energy) derivative of the integral flux, dividing
 * by the total solid angle and extrapolating down in energy if 
 * need be.
 *
 */
double  getJdiff(float J[][100], double E, double alpha_lc)
{
  int row, i, topCol, botE;
  double J1, J2, I, x1, x2, y1, y2, m, c, x_ext, y_ext, J_ext;

  row = (int)floor((L_TARG+0.11 - J[1][0])/0.1); // to make sure! 
  if(  fabs((double)J[row][0]-L_TARG) > 1e-3   ) 
    printf("\nL-shell not matching data\n\a");

  I = PI * cos(alpha_lc) * (PI - 2*alpha_lc);

  // Find column corresponding to highest energy value
  for(i=0; i<100; i++) {
    if(J[0][i+1] < 0.01) { 
      topCol = i; 
      break; 
    }
  }



  // Case 1. E < 100 keV
  // -------------------

  if( E <= 1e5 ) {
 
    // diff flux @ 100 keV and 200 keV
    J1 = 1e-6*fabs(J[row][2] - J[row][1]) / (J[0][2] - J[0][1]); 
    J2 = ((1e-6*fabs(J[row][3] - J[row][2]) / (J[0][3] - J[0][2])) 
	  + J1 )/2; // central difference

    // do extrapolation in log-log space for best fit 
    x1 = log10( J[0][1]*1e6 );
    x2 = log10( J[0][2]*1e6 );
    y1 = log10( J1 );
    y2 = log10( J2 );

    m = (y2-y1)/(x2-x1);		// gradient of line
    c = (y1*x2 - y2*x1)/(x2-x1)	;	// offset of line, i.e.
					// y = m*x + c
    x_ext = log10( E );
    y_ext = m*x_ext + c;
    J_ext = pow(10, y_ext);

    return (J_ext/I);

  }


  

  // Case 2. E > 7 MeV
  // -----------------

  if( E >= 7e6 ) {
  
    // If flux at 7 Mev = 0, flux above it is zero too
    if( J[row][topCol]==0 )  return 0;

    // Otherwise need to extrapolate as in case 1.
    // diff flux @ 6.5 MeV and 7 MeV
    J2 = 1e-6*fabs( J[row][topCol] - J[row][topCol-1] ) 
      / (J[0][topCol] - J[0][topCol-1]); 

    J1 = ((1e-6*fabs( J[row][topCol-1] - J[row][topCol-2]) / 
	   (J[0][topCol-1] - J[0][topCol-2]) ) + J2 )/2; // cdiff

    // do extrapolation in log-log space for best fit 
    x1 = log10( J[0][topCol-1]*1e6 );
    x2 = log10( J[0][topCol]*1e6 );
    y1 = log10( J1 );
    y2 = log10( J2 );

    m = (y2-y1)/(x2-x1);		// gradient of line
    c = (y1*x2 - y2*x1)/(x2-x1)	;	// offset of line, i.e.
					// y = m*x + c
    x_ext = log10( E );
    y_ext = m*x_ext + c;
    J_ext = pow(10, y_ext);

    if(J_ext < 1e-10 ) J_ext = 0.0;

    return (J_ext/I);
  }


  // Case 3. 100 keV < E < 7 MeV
  if( E<7e6 && E>1e5 ) {


    // Find column corresponding lower energy value
    for(i=1; i<100; i++) {
      if( (J[0][i+1]*1e6) > E ) { 
	botE = i; 
	break; 
      }
    }


    // central diff flux @ lower and higher energies
    J1 = ( (1e-6 * fabs( J[row][botE] - J[row][botE-1] )
	    / ( J[0][botE] - J[0][botE-1] ) ) + 
	   (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
	    / ( J[0][botE+1] - J[0][botE] ) )  ) / 2;

    J2 = ( (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
	    / ( J[0][botE+1] - J[0][botE] ) ) + 
	   (1e-6 * fabs( J[row][botE+2] - J[row][botE+1] )
	    / ( J[0][botE+2] - J[0][botE+1] ) )  ) / 2;

    if(botE == 1)
      J1 =  (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
	      / ( J[0][botE+1] - J[0][botE] ) );
    
    if(botE == (topCol-1))
      J2 = (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
	    / ( J[0][botE+1] - J[0][botE] ) );
    



    // If J1 = J2 = 0, interpolated value also 0
    if( J1==0 && J2==0 ) return 0;



    // If only J2 = 0, do linear interpolation
    if( J2 == 0 ) {
      J_ext = J1*( ( J[0][botE+1]-(E*1e-6) )/
		   ( J[0][botE+1] - J[0][botE] ) );
      return (J_ext/I);
    }



    // Otherwise interpolate as in case 1 (log-log space)

    x1 = log10( J[0][botE]*1e6 );
    x2 = log10( J[0][botE+1]*1e6 );
    y1 = log10( J1 );
    y2 = log10( J2 );

    m = (y2-y1)/(x2-x1);		// gradient of line
    c = (y1*x2 - y2*x1)/(x2-x1)	;	// offset of line, i.e.
					// y = m*x + c
    x_ext = log10( E );
    y_ext = m*x_ext + c;
    J_ext = pow(10, y_ext);

    return (J_ext/I);
  }

}











/*
 * FUNCTION reduceFactor
 * ---------------------
 * This function calculates the factor by which the electric 
 * field is reduced when we move dlong degrees off-center in 
 * longitude away from the source.
 *
 */
double reduceFactor(double dlong)
{
  double d0, d1, Rl0, Rl1;

  d0  = 0.9*D2R*(R_E + H_IONO/2);
  Rl0 = sqrt(d0*d0 + H_IONO*H_IONO);
  
  // maybe use spherical geometry?

  d1 =  sqrt( d0*d0 + pow( (dlong*D2R*(R_E+H_IONO)) , 2) );
  Rl1 = sqrt(d1*d1 + H_IONO*H_IONO);

  return ((d1*Rl0*Rl0)/(Rl1*Rl1*d0));
}
