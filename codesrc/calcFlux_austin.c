/*---------------------------------------------------------------
* 8/10/06:
  Already made some changes.  Unfortunately, just started documenting
  them now.

  1) Instead of using large number of frequencies (necessary for broad-band
  whistlers), reading in 2 frequency components from command line.  The lower is
  the transmitter fkc, and the upper is a small bandwidth added to the signal.

  Other than this modification, the code is basically unchanged.
------------------------------------------------------------------*/



#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include "consts.h"

// long freqs[] = {  200,
//       210,  
//       220,    
//       230,
//       240,
//       250,
//       260,
//       270,
//       280,
//       290,
      
//       300,
//       310, 
//             320,
//       330,
//       340,
//       350,
//       360,
//       370,
//       380,
//       390,
      
//       400,
//       410,
//       420,
//       430,
//       440,
//       460,
//       480,
//       500,
//       520,
//       540,
      
//       560,
//       580,
//       600,
//       620,
//       640,
//       660,
//       690,
//       720,
//       750,
//       780, 
      
//       810,
//       840,
//       880,
//       920,
//       960, 
//       1000,  
//       1040,
//       1080,
//       1120,
//       1160,
      
//       1200,
//       1250,
//       1300,
//       1360,
//       1420,
//       1480,
//       1540,
//       1600,
//       1670,
//       1740,
        
//       1800,
//       1880,
//       1960,
//       2050,
//       2140,
//       2230,
//       2320,
//       2410,
//       2500,
//       2600,
      
//       2700,
//       2800,
//       2900, 
//       3000,
//       3100,
//       3200,
//       3300,
//       3450,
//       3600,
//       3750,
      
//       3900,
//       4050,
//       4200,
//       4400,
//       4600, 
//       4800,
//       5000,
//       5200,
//       5400,
//       5600,
      
//       5800,
//       6000,
//       6200,
//       6400,
//       6600,
//       6900,
//       7200,
//       7500,
//       7800,
//       8100,
      
//       8400,
//       8800,
//       9200,
//       9600,
//       10000,
//       10500,
//       11000,
//       11500,
//       12000,
//       12500,
      
//       13000,
//       13500,
//       14000,
//       14500,
//       15000,
//       16000,
//       17000,
//       18000,
//       19000,
//       20000,
      
//       22000,
//       24000,
//       26000,
//       28000,
//       30000,
//       35000,
//       40000,
//       45000,
//       50000,
//       60000     };




// DEFINITIONS
// -----------


// #define		RES_DT		0.02
// #define		RES_FINT	5.0
// #define		C		2.997956376932163e+08
// #define		R_E		6370000.0
// #define		H_IONO		1E5
// #define		PI		3.141592653589793115997963468544185
// #define		D2R		PI/180.0
// #define		M_EL		9.1E-31
// #define		E_EL		5.105396765648739E5

// #define		DE_EXP		0.003
// #define		E_EXP_BOT	1.477
// #define		E_EXP_TOP	7.477
// #define		NUM_E		2200
// #define		SQUARE		1





// GLOBAL VARIABLES:
// -----------------

double		L_TARG;
int		NUM_TARGS	=	1;
int		NUM_TIMES;

int nancounter;


/*
------MAKING CHANGES----------
Making changes on 9/14/05
Praj Kulkarni

Since we're doing narrowband frequencies, just
going to pass it in from command line. So i'm 
going to comment out the 'freqs' array. 

*/
// long freqs[] = {	200,
//       210,	
// 			220,		
// 			230,
// 			240,
// 			250,
// 			260,
// 			270,
// 			280,
// 			290,
			
// 			300,
// 			310, 
// 		       	320,
// 			330,
// 			340,
// 			350,
// 			360,
// 			370,
// 			380,
// 			390,
			
// 			400,
// 			410,
// 			420,
// 			430,
// 			440,
// 			460,
// 			480,
// 			500,
// 			520,
// 			540,
			
// 			560,
// 			580,
// 			600,
// 			620,
// 			640,
// 			660,
// 			690,
// 			720,
// 			750,
// 			780, 
			
// 			810,
// 			840,
// 			880,
// 			920,
// 			960, 
// 			1000,  
// 			1040,
// 			1080,
// 			1120,
// 			1160,
			
// 			1200,
// 			1250,
// 			1300,
// 			1360,
// 			1420,
// 			1480,
// 			1540,
// 			1600,
// 			1670,
// 			1740,
				
// 			1800,
// 			1880,
// 			1960,
// 			2050,
// 			2140,
// 			2230,
// 			2320,
// 			2410,
// 			2500,
// 			2600,
			
// 			2700,
// 			2800,
// 			2900, 
// 			3000,
// 			3100,
// 			3200,
// 			3300,
// 			3450,
// 			3600,
// 			3750,
			
// 			3900,
// 			4050,
// 			4200,
// 			4400,
// 			4600, 
// 			4800,
// 			5000,
// 			5200,
// 			5400,
// 			5600,
			
// 			5800,
// 			6000,
// 			6200,
// 			6400,
// 			6600,
// 			6900,
// 			7200,
// 			7500,
// 			7800,
// 			8100,
			
// 			8400,
// 			8800,
// 			9200,
// 			9600,
// 			10000,
// 			10500,
// 			11000,
// 			11500,
// 			12000,
// 			12500,
			
// 			13000,
// 			13500,
// 			14000,
// 			14500,
// 			15000,
// 			16000,
// 			17000,
// 			18000,
// 			19000,
// 			20000,
			
// 			22000,
// 			24000,
// 			26000,
// 			28000,
// 			30000,
// 			35000,
// 			40000,
// 			45000,
// 			50000,
// 			60000     };


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

float *getArr(void);
void updateArr(float *arr1, float *arr2);
void compFlux(float *arr, double L, int k, char *dir);
void readJ(float J[][100]);
double getJdiff(float J[][100], double E, double alpha_lc);




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
  int numFreqs, i, k, m, nin, ei, ti;
  char *dir, sysCmd[512], filename[64], *NS, alphaFile[128];
  float L, *arr, *arrarr[128]; // <- array of pointers to precip arrays
  float J[100][100];
  double Jext;

  nancounter = 0;
  // Get directory name from command line
  dir   = argv[1];
  L_TARG = atof(argv[2]);

  
  printf("DIRECTORY:%s\n",dir);
  printf("L_TARG: %f\n",L_TARG);
  

  NUM_TIMES = (int)ceil( RES_FINT / RES_DT );
  printf("\n\aWill do %d time steps\n", NUM_TIMES);

  // 64 bits now
  numFreqs = (sizeof freqs)/sizeof(long) - 1; // don't want to do f=60000!

  printf("\n\aWill do %d frequencies\n", numFreqs);
  //numFreqs = 1;
  sprintf( alphaFile, "%s/alpha_%g_%s", dir, L_TARG, "N" );
  printf("alphaFile:%s\n",alphaFile);
  arr = getArr();

  if(  (alphaPtr=fopen(alphaFile, "r"))!=NULL  ) {

    nin = fread(arr,sizeof(float),(NUM_E*NUM_TIMES),alphaPtr);
    fclose(alphaPtr);
    for(ei=0; ei<NUM_E; ei++) {
      for(ti=0; ti<NUM_TIMES; ti++) {
	      arr[ei*NUM_TIMES+ti] = pow(arr[ei*NUM_TIMES+ti],2);
      }
    }
    arrarr[0] = arr;


    sprintf( alphaFile,   "%s/alpha_%g_%s", dir, L_TARG, "S" );
    if((alphaPtr=fopen(alphaFile, "r"))==NULL) 
      printf("\n\aprob opening %s\n", alphaFile);  
    arr = getArr();
    nin = fread(arr,sizeof(float),(NUM_E*NUM_TIMES),alphaPtr);
    fclose(alphaPtr);
    for(ei=0; ei<NUM_E; ei++) {
      for(ti=0; ti<NUM_TIMES; ti++) {
	       arr[ei*NUM_TIMES+ti] = pow(arr[ei*NUM_TIMES+ti],2);
      }
    }
    arrarr[1] = arr;



  } else {


    for(i=0; i<numFreqs; i++) {
      
      // unpack the archive & unzip
      //printf("i=%d , FILE: p%d_%g\n",i,freqs[i],L_TARG);
      // sprintf(sysCmd,"cd %s; tar xf p%d_%g; gunzip *.gz",
	     //  dir, freqs[i], L_TARG );
      // printf("%s\n",sysCmd);
      // system(sysCmd);
  
            
      for(k=0; k<2; k++) {
	if(k==0) {NS = "N";} else {NS="S";}	
	if(i==0)  arrarr[k]=getArr();
	sprintf(filename,"%s/p%s%d_%g.dat",dir,NS,freqs[i],L_TARG);
	//printf("i: %d, k: %d, filename: %s\n", i, k, filename);
	inPtr = fopen(filename, "r");
	nin = fread(arr, sizeof(float), (NUM_E*NUM_TIMES), inPtr);
	fclose(inPtr);
	
	updateArr( arrarr[k] , arr );
      } // for(k ... )  N/S - hemisphere
      
      // need to rm the unzipped file
      // sprintf(sysCmd, "cd %s; rm -f pN%d* pS%d*",
	     //  dir, freqs[i], freqs[i]);
      // system(sysCmd);
      
    } // freqs
    
  } // if alpha array exists
  
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
float *getArr(void)
{
  float *arr;
  int ei, ti;

  arr = (float *) malloc( NUM_E * NUM_TIMES * sizeof(float) );
  if(arr == NULL) {
    printf("\nProb assigning mem in calcFlux\n");
    exit(0);
  }

  for(ei=0; ei<NUM_E; ei++) {
    for(ti=0; ti<NUM_TIMES; ti++) {
      arr[ei*NUM_TIMES+ti] = 0.0;
    }
  }

  return arr;

}




/* 
 * FUNCTION: updateArr
 * -------------------
 * This function updates the values of arr1 with those of arr2.
 *
 */
void updateArr(float *arr1, float *arr2)
{
  int ei, ti;

  for(ei=0; ei<NUM_E; ei++) {
    for(ti=0; ti<NUM_TIMES; ti++) {
      if(arr2[ei*NUM_TIMES+ti]>0.0)
	arr1[ei*NUM_TIMES+ti] += arr2[ei*NUM_TIMES+ti];
    }
  }
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
  FILE *phiPtr, *QPtr, *NPtr, *alphaPtr;
  int ei, ti, i, nout;
  double mag, P, Q, epsm, alpha_eq, v, crunch;
  double I=0.0, x, g, field, Phi_p, b=1.0, vcm, fcgs;
  double v_tot_arr[NUM_E], E_tot_arr[NUM_E], gamma;
  double Jdiff[NUM_E];
  float Phi_float, alpha, J[100][100];
  char *NS, PhiFile[128], QFile[128], NFile[128], alphaFile[128];

  float arg_in1;
  float arg_in2;
  float arg_in3;
  float arg1;
  float arg2;
  float arg3;

  // Open up Phi, Q, N file for writing
  if(k==0) {NS = "N";} else {NS="S";}  

  sprintf( PhiFile, "%s/phi_%g_%s", dir, L, NS );
  sprintf( QFile,   "%s/Q_%g_%s", dir, L, NS );
  sprintf( NFile,   "%s/N_%g_%s", dir, L, NS );
  sprintf( alphaFile,   "%s/alpha_%g_%s", dir, L, NS );  

  printf("writing %s\n", PhiFile);

  if( (phiPtr=fopen(PhiFile, "w"))==NULL ) {
    printf("\nProblem opening %s\n", PhiFile);
   exit(0);
   }
  /*
  if( (QPtr=fopen(QFile, "w"))==NULL ) {
    printf("\nProblem opening %s\n", QFile);
   exit(0);
   }

  if( (NPtr=fopen(NFile, "w"))==NULL ) {
    printf("\nProblem opening %s\n", NFile);
   exit(0);
   }
  */  
  if( (alphaPtr=fopen(alphaFile, "w"))==NULL ) {
   exit(0);
  }


  epsm = (1/L)*(R_E+H_IONO)/R_E;

  crunch	= sqrt(1+3*(1-epsm))/pow(epsm,3) ;

  alpha_eq	= asin(sqrt( 1/crunch ));

  readJ(J);

  // Precalculate energy and velocity values
  for(i=0; i<NUM_E; i++) {
    E_tot_arr[i] = pow(10, (E_EXP_BOT+(DE_EXP/2)+DE_EXP*i) ); //E in eV
    Jdiff[i] = getJdiff( J, E_tot_arr[i], alpha_eq );
    v_tot_arr[i] = C*sqrt(1 - pow( (E_EL/(E_EL+E_tot_arr[i])) ,2) );
  }

  for(ei=0; ei<NUM_E; ei++) {

    if(SQUARE) {
      v = v_tot_arr[ei];
      vcm = v*100;	// v in cm for distrib fn calculation
      gamma = 1/sqrt( 1 - v*v/(C*C) );
      fcgs =	4.9e5/pow( (vcm*gamma) ,4) - 
		8.3e14/pow( (vcm*gamma) ,5) + 
		5.4e23/pow( (vcm*gamma) ,6);
     



      //fcgs = fcgs*50;
      // fcgs = 7.034e26 / pow(vcm,6); //10^8/E^2 distribution

      b = (v*v/M_EL)*pow( sqrt(1 - (v*v)/(C*C)), 3) * 1.6e-8 * fcgs;

      // b = 1e8 / pow(E_tot_arr[i],2);
    } else {
      b = Jdiff[ei]*1000;
    }

    //display values of b and E
    // USEFUL BUT ANNOYING RIGHT NOW - 5/2015
  //  printf("b = %f\tE = %f\n",b,E_tot_arr[ei]);

    for(ti=0; ti<NUM_TIMES; ti++) {
      
      alpha = sqrt( arr[ei*NUM_TIMES+ti] );

      

      nout=fwrite(&alpha, sizeof(float), 1, alphaPtr);      

      mag = 1.4142135623731*alpha;	// sqrt(2)*alpha_RMS = peak
      
      P = mag/2;		//[ alpha_lc - (alpha_lc-mag) ] /2
      Q = alpha_eq - mag/2;	//[ alpha_lc + (alpha_lc-mag) ] /2
      
      I = 0.0;
      if(mag != 0) {
	     for(i=0; i<5; i++) {
	       x = P*t5[i] + Q ;

         // ------------- Austin's edits: Trying to fix NaNs being
         //               thrown by asin() function
	 

      	  if(SQUARE) {

            // --------- Austin's careful re
            arg1 = (P/PI)*sin(2*x);
            arg_in1 = (x - alpha_eq)/mag;
            arg_in2 = 0;
            arg_in2 = asin(abs(arg_in1));
            // if (arg_in1 < 0) { arg_in2 = -1.0*asin(abs(arg_in1)); }
            // else { arg_in2 = asin(arg_in1); }
          
            arg2 = arg_in2 + (PI/2);

            g = arg1*arg2;
      	   // g = (P/PI)*sin(2*x)*(  asin((x-alpha_eq)/mag)+ (PI/2) );

      	    } else {

            arg1 = (P/PI)*sin(2*x);
            arg_in1 = (x - alpha_eq)/mag;
            arg_in2 = 0;
            arg_in2 = asin(abs(arg_in1));
            // if (arg_in1 < 0) { arg_in2 = -1.0*asin(abs(arg_in1)); }
            // else { arg_in2 = asin(arg_in1); }
          
            arg2 = arg_in2 + (PI/2);
            arg_in3 = (float)(mag*mag-pow((x-alpha_eq),2));

            g = arg1*(abs(arg_in1*arg2 + sqrt(arg_in3)));

      	    // g = (P/PI)*sin(2*x)*((x - alpha_eq)*( asin((x-alpha_eq)/mag)+ (PI/2) ) +
      				 // sqrt(mag*mag-pow((x-alpha_eq),2)));
	        }; // Square

          if isnan(g) { printf("G ISNAN\n"); };
          if isnan(arg1) { printf("Fuckin arg1\n");};
          if isnan(arg2) { printf("ergh, arg2\n");};

          if isnan(arg_in1) { printf("jeez, arg_in1\n");};
          if isnan(arg_in2) { printf("jeez, arg_in2: input was: %f\n",arg_in1);};
          if isnan(arg_in3) { printf("jeez, arg_in3");};
          
	     I += ( beta5[i]*g );

	     } // for(i ... ) -> Gauss quad integration
      } // if mag != 0
      Phi_p = PI*crunch*b*I;
      Phi_float = ( float ) Phi_p;

      //if isnan(I) { printf("I ISNAN\n"); };


      if isnan(Phi_p) { 
        nancounter=nancounter + 1;
        //printf("Total NaNs: %i\n",nancounter);
      };

      nout=fwrite(&Phi_float, sizeof(float), 1, phiPtr);
      
    } // for(ti ... )
  } // for(ei ... )
  
  fclose(phiPtr);
  fclose(alphaPtr);

  printf("Total NaNs: %i\n",nancounter);

  // Now calculate Q and N
  //
  // Need to integrate over E
  //for(ti=0; ti<NUM_TIMES; ti++) {
  //  for(ei=0; ei<NUM_E; ei++) {
  //    } // ei
  // }  // ti

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

  filename = "EQFLUXMA.dat";

  if( (filePtr = fopen( filename ,"r")) == NULL ) {
    printf("Hey buddy, I can't open the FLUX file!\n");
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


