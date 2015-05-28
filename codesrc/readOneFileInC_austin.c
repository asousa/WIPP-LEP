#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
# include <stdlib.h>






// Global variables:
//------------------

float		LAM_S;
double		MAX_POWER;
double		DIV_LAT_NUM;
double		DIV_FREQ_NUM;
int		NUM_TIMES;

int		EA_NUMROWS;
int		EA_NUMCOLS;

double		L_TARG;
double		L_MARGIN[]	=	{ 6E-4};
double		LAT_TARG[]	=	{ 0.0 };
double		LAT_MARGIN[]	=	{ 80.0 };
int		NUM_TARGS	=	1;







// Definitions:
//-------------

#define		PI		3.141592653589793115997963468544185161590576
#define		D2R		PI/180.0
#define		R2D		180.0/PI

#define		NUM_STEPS	15000
#define		T_STEP		0.0004

#define		Q_EL		1.602E-19
#define		M_EL		9.1E-31
#define		E_EL		5.105396765648739E5
#define		MU0		PI*4E-7		
#define		EPS0		8.854E-12
#define		C		2.997956376932163e+08
#define		Z0		377.0
#define		R_E		6370000.0
#define		H_MAGNETO	1E6
#define		H_IONO		1E5

#define		A		5E3	
#define		B		1E5
#define		H_E		5000.0
#define		I0		-100000.0
// peak current defined here, in Amps!

#define		P_DIST		0.0
#define		Q_DIST		2.0
#define		AN_CM_DIST	2E5
#define		V0_DIST		1.0
#define		M_RES		0

#define		EALimS		-40.0
#define		EALimN		40.0
#define		EAIncr		1.0
#define		dL0		6E-4
#define		DF		50.0
#define		DT		0.01
#define		NUMLATS		((EALimN - EALimS)/EAIncr + 1)

// #define		DV_TOT		1e5
#define		DE_EXP		0.003
#define		E_EXP_BOT	1.477
#define		E_EXP_TOP	7.477	
#define		RES_DT		0.02
#define		RES_FINT	5.0
#define		EA_SPLIT	1
#define		MULT		2.0
#define		NUM_E		2200




// Data structures:
//-----------------

// Store each ray in rayT
typedef struct {
  float	tg[NUM_STEPS];		// group time
  float	distre[NUM_STEPS];     	// geocentric dist in Re
  float	lat[NUM_STEPS];		// latitude
  float	delta[NUM_STEPS];      	// angle wrt vert
  float	t[NUM_STEPS];		// phase time	
  float	l_sh[NUM_STEPS];       	// local L-shell
  float	psi[NUM_STEPS];		// angle wrt B-field
  float	psiray[NUM_STEPS];     	// ray angle wrt B
  float	psires[NUM_STEPS];     	// res. cone angle wrt B
  float	mu[NUM_STEPS];		// refr. index
  float	dens[NUM_STEPS];       	// electron density /cm^3
  float	anH[NUM_STEPS];		// concentrations of: H+
  float	anHe[NUM_STEPS];       	//      ,,            He+
  float	anO[NUM_STEPS];		//      ,,            O+
  float	fH[NUM_STEPS];		// gyrofreq. in kHz
  float	stixP[NUM_STEPS];		
  float	stixR[NUM_STEPS];
  float	stixL[NUM_STEPS];	
  int	numSteps;		// num of steps in this ray
  float lam_s;			// latitude of lightning source
  long	f;			// ray frequency
  float pwr[NUM_STEPS];	       	// pwr/m^2 @ every point
} rayT;


// Store each cell in the spectrogram in cellT
typedef struct cellT {
  double		t;
  double		f;
  double		pwr;
  double		psi;
  double		mu;
  double		stixP;
  double		stixR;
  double		stixL;
  int			num_rays;	// that added into this cell
  struct cellT		*link;
} cellT;










// Prototypes:
//------------

FILE *readFile(char *fileName, FILE *filePtr, rayT *ray, long f);
FILE *readDamp(char *filename, FILE *filePtr, rayT *ray, long f);
void initRay(rayT *ray );
void initStep(rayT *ray, int i);
float interpPt(float *xI, float *yI, int n, float xO);
void doInterp(rayT *BL, rayT *TL, rayT *BR, rayT *TR, 
	      cellT *lat_arr[][NUM_TARGS], double *EAarrPtr);
float smallest(float n1, float n2, float n3, float n4);
FILE *writeFile(long lower_freq, char *pref, char *suff);
void addToLFile(double BL_fact, double TL_fact, double BR_fact, 
		double TR_fact, double t, FILE *outPtr, 
		rayT *BL, rayT *TL, rayT *BR, rayT *TR, 
		double div_lat, double div_freq, 
		double l_int, double lat_int);
void checkCross(double BL_fact, double TL_fact, double BR_fact,
		double TR_fact, double t, cellT *lat_arr[][NUM_TARGS], 
		rayT *BL, rayT *TL, rayT *BR, rayT *TR, 
		double div_lat, double div_freq, 
		double l_int, double lat_int, double *EAarrPtr,
		double l_int_prev, double lat_int_prev, int iTarg);
double ltgPwr(float lam_s, float dlat, float dfreq, long f, 
	      float lat, float dlong );
double ionoAbsorp(float lat, long f);
void initPwr(rayT *BL, rayT *TL, rayT *BR, rayT *TR);
int outsideL(double BL_l, double TL_l, double BR_l, double TR_l, int *iTarg);
int outsideLat(double BL_lat, double TL_lat, double BR_lat, 
	       double TR_lat, int *iTarg);
double *initEA_Arr(void);
void writeEntry(int row, int col, double entry, double *arrPtr);
double getEntry(int row, int col, double *arrPtr);
void initLatArr( cellT *lat_arr[][NUM_TARGS] );
void addToArr(cellT *lat_arr[][NUM_TARGS],double t, double f, 
	      double pwr, double psi, double mu, int lami, 
	      int iTarg, double stixP, double stixR, double stixL);
int compare(double t_try, double f_try, double t_next, double f_next);
void dispLatArr(cellT *lat_arr[][NUM_TARGS], FILE *outPtr);
void freeArr(cellT *lat_arr[][NUM_TARGS]);
void calcRes(cellT *lat_arr[][NUM_TARGS], long lower_freq);
void getFltConst(double L, double lat, double alpha_eq, 
		 double *flt_const_N, double *flt_const_S);
void Fresnel(double x0, double *FS, double *FC);
void addToAlphaArr(float **arr_N, float **arr_S, 
		   double iono_time, int e_toti, double dalpha);













/*
 * FUNCTION: main
 * --------------
 * This function reads in 2 different ray files of different 
 * frequencies, and isolates 4 rays at a time: 
 * BL - Bottom Left ray, lower frequency, lower launch latitude
 * TL - Top Left ray, lower frequency, higher launch latitude
 * BR - Bottom Right ray, higher frequency, lower launch lat
 * TR - Top Right ray, higher freq. higher launch lat
 *
 * main accepts the following arguments at the command line:
 * low_freq - lower frequency in Hz 
 * high_freq - higher frequency in Hz
 * lam_s - center of lightning strike
 * div_lat_num - number of subdivisions between adjacent latitudes
 * freq_step - frequency step of interpolation, i.e. every 1 Hz
 * L_TARG - Target L-shell where we want spectra, precip, etc.
 *
 */

int main(int argc, char *argv[])
{
  FILE	*outPtr=NULL, *ptrL=NULL, *ptrR=NULL;
  FILE *dampL=NULL, *dampR=NULL;
  char sysCmd[128];
  char fileL[128], fileR[128], dampFileL[128], dampFileR[128];
  char *freqPref="f", freqSuff[32];
  float yO, step; 
  double AN=0.0, v0_sq, freq_step, *EA_Arr;
  rayT BL, TL, BR, TR, *tmp;
  int i, j;
  long lower_freq, upper_freq;
  cellT *lat_arr[ (int)NUMLATS ][ NUM_TARGS ];





  // --------------  get arguments from command line ---------------
  if(argc != 7) 
    { printf("Wrong number of input arguments\n");  exit(0); }

  lower_freq = atol(argv[1]);
  upper_freq = atol(argv[2]);
  LAM_S = atof(argv[3]);
  DIV_LAT_NUM = atof(argv[4]);
  freq_step = atof(argv[5]);
  L_TARG = atof(argv[6]);
  DIV_FREQ_NUM = (upper_freq - lower_freq) / freq_step ;

  sprintf(fileL, "./newray%d.dat",lower_freq);
  sprintf(fileR, "./newray%d.dat",upper_freq);
  sprintf(dampFileL, "./d%d.dat",lower_freq);
  sprintf(dampFileR, "./d%d.dat",upper_freq);
  printf("\nlower_freq: %d\nupper_freq: %d\n",lower_freq,upper_freq);
  printf("lam_s: %g,\ndiv_lat_num: %g \n\n",LAM_S, DIV_LAT_NUM);

  // ----------------------------------------------------------------




  // ----------------  Run initialization scripts -------------------
  EA_Arr = (double *)initEA_Arr();
  sprintf(freqSuff, "_%g", L_TARG);
  outPtr = (FILE *)writeFile(lower_freq, freqPref, freqSuff);
  initLatArr(lat_arr);
  // ----------------------------------------------------------------




  // ---------- Step through all rays in one rayfile ----------------
 
  ptrL  = (FILE *)readFile( fileL, ptrL, &BL, lower_freq);
  dampL = (FILE *)readDamp( dampFileL, dampL, &BL, lower_freq);
  ptrR  = (FILE *)readFile( fileR, ptrR, &BR, upper_freq );
  dampR = (FILE *)readDamp( dampFileR, dampR, &BR, upper_freq);

  while( 1 ) {
    ptrL = (FILE *)readFile( fileL, ptrL, &TL, lower_freq);
    dampL = (FILE *)readDamp( dampFileL, dampL, &TL, lower_freq);
    ptrR = (FILE *)readFile( fileR, ptrR, &TR, upper_freq);
    dampR = (FILE *)readDamp( dampFileR, dampR, &TR, upper_freq);
    

    if(ptrL==NULL || ptrR==NULL) break;
    
    initPwr(&BL, &TL, &BR, &TR);
   
    doInterp( &BL, &TL, &BR, &TR, lat_arr, EA_Arr );
    
    tmp=&BL; BL=TL; TL=*tmp;
    tmp=&BR; BR=TR; TR=*tmp;
  }
  printf("made it past the while loop\n\n");

  dispLatArr(lat_arr, outPtr);
  fclose(outPtr);
  if(ptrL != NULL) fclose(ptrL);
  if(ptrR != NULL) fclose(ptrR);
  free(EA_Arr);
  // ------------------------------------------------------------------





  // ---------------- Do the resonance calculation --------------------
  calcRes(lat_arr, lower_freq);
  // zip and archive files for network transfer 
  // printf("\nOK before gzip, tar and rm\n");
  // sprintf(sysCmd, "gzip pN* pS*; tar cf p%d_%g ./*.gz; rm -f ./*.gz",
	 //  lower_freq, L_TARG); 
  // system(sysCmd);
  // printf("OK after gzip, tar and rm\n");
  // ------------------------------------------------------------------





  // ------------------ Close all files and return ---------------------
  freeArr(lat_arr);
  printf("OK after free lat_arr\n\n");
  return 0;
  
  // --------------------------------------------------------------------

}









/*
 * FUNCTION: readFile
 * ------------------
 * This will take in the name of a file and return the next ray as 
 * a rayT* until the end of the file, when the file will be closed.
 * 
 */

FILE  *readFile( char *fileName, FILE *filePtr, rayT *ray, long f)
{
  int i=0;
  

  //  Open the file for reading  //
  if(filePtr == NULL) {
    if( (filePtr = fopen( fileName ,"r")) == NULL ) {
      printf("Hey buddy, I can't open the ray file!\n");
      exit(1);
    } else {
      printf("\nFile %s\nOpened successfully! \n\n", fileName);
    }
    // Place the pointer after the header for reading in data
    fseek(filePtr, 146*sizeof(float), SEEK_SET);
    //fseek(filePtr, 578 , SEEK_SET);
  }


  // Read in each ray //
  initRay( ray );
  while( 1 ) { 
    if(fscanf(filePtr, 
	      "%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e", 
	      &(ray->tg[i]), &(ray->distre[i]), &(ray->lat[i]), 
	      &(ray->delta[i]), &(ray->t[i]), &(ray->l_sh[i]),
	      &(ray->psi[i]), &(ray->psiray[i]), &(ray->psires[i]), 
	      &(ray->mu[i]), &(ray->dens[i]), &(ray->anH[i]), 
	      &(ray->anHe[i]), &(ray->anO[i]), &(ray->fH[i]), 
	      &(ray->stixP[i]), &(ray->stixR[i]), 
	      &(ray->stixL[i]) ) == EOF ) {
      //fclose( filePtr );
      ray->numSteps=0;
      return NULL;
    }
    i++;
    
    if( ray->tg[(i-1)]>=99998.0 ) {
      initStep( ray, (i-1) );
      break; 
    }

    if(i==NUM_STEPS) {
      i--;	// back up a step, record next step over this one until
		// eventually it is rewritten by last line & initialised
      printf("\nnumber of NUM_STEPS reached\n");
    }

  }
  ray->numSteps = i-1;

  // Set the ray frequency 
  ray->f = f;

  // Determine the center of the lightning, i.e. source lat
  ray->lam_s = LAM_S;


  printf("\nNow reading, f: %d, lat: %g\n",ray->f, ray->lat[0]);

  return filePtr; 
}

  









/*
 * FUNCTION: readDamp
 * ------------------
 * This function is very similar to the readFile function and takes the
 * name of a "damping" file, and reads it.  It places the data into   
 * ray->pwr to be used in later parts othe program.  Note that this
 * damping is only Landau damping + geometrical focusing/defocusing, and 
 * is normalized to 1 at the beginning of the ray.
 *
 */

FILE  *readDamp( char *fileName, FILE *filePtr, rayT *ray, long f)
{
  int i=0;
  

  //  Open the file for reading  //
  if(filePtr == NULL) {
    if( (filePtr = fopen( fileName ,"r")) == NULL ) {
      printf("Hey buddy, I can't open the damping file!\n");
      exit(1);
    } else {
      printf("File %s\nOpened successfully! \n", fileName);
    }
  }


  // Read in each ray //
  while( 1 ) { 
    if(fscanf( filePtr, "%g", &(ray->pwr[i]) ) == EOF ) {
      fclose( filePtr );  // should I close this here?
      return NULL;
    }
    i++;
    
    if( ray->pwr[(i-1)]>=99998.0 ) {
      ray->pwr[(i-1)]=0.0;
      break; 
    }

    if(i==NUM_STEPS) {
      i--;	// back up a step, record next step over this one until
		// eventually it is rewritten by last line & initialised
      printf("\nnumber of NUM_STEPS reached\n");
    }

  }
  return filePtr; 
}

  






/* 
 * FUNCTION initRay
 * ----------------
 * Takes in a ray structure and initialises all the variables
 *
 */

void initRay(rayT *ray )
{ 
  int i;

  for( i=0 ; i<NUM_STEPS ; i++)
    initStep(ray,i);
      
  ray->numSteps = 0;
  ray->f = 0;
  ray->lam_s = 0;
}


void initStep(rayT *ray, int i)
{
   ray->tg[i]=0.0;	     	// group time
   ray->distre[i]=0.0;     	// geocentric dist in Re
   ray->lat[i]=0.0;		// latitude
   ray->delta[i]=0.0;      	// angle wrt vert
   ray->t[i]=0.0;      		// phase time	
   ray->l_sh[i]=0.0;       	// local L-shell
   ray->psi[i]=0.0;		// angle wrt B-field
   ray->psiray[i]=0.0;     	// ray angle wrt B
   ray->psires[i]=0.0;     	// res. cone angle wrt B
   ray->mu[i]=0.0;		// refr. index
   ray->dens[i]=0.0;       	// electron density /cm^3
   ray->anH[i]=0.0;		// concentrations of: H+
   ray->anHe[i]=0.0;       	//      ,,            He+
   ray->anO[i]=0.0;		//      ,,            O+
   ray->fH[i]=0.0;		// gyrofreq. in kHz
   ray->stixP[i]=0.0;		
   ray->stixR[i]=0.0;
   ray->stixL[i]=0.0;
   ray->pwr[i]=0.0;
}




/*
 * FUNCTION: interpPt
 * ------------------
 * This function is designed to work pretty much the same as matlab's
 * interp1 function.  Specify the input x and y vectors, and the output
 * x point.
 *
 */

float interpPt(float *xI, float *yI, int n, float xO)
{
  int i, iHigh, iLow, iMid;
  float yO;

  // Check that xO is within bounds
  if( (xO < xI[0]) || xO > xI[n-1] ) {
    printf("\nPoint is out of bounds!\a\n");
    return 0.0;
  }
  
  // Do a binary search for the correct index 
  iHigh = n-1;
  iLow = 0;  
  while(1) {
    iMid = (iHigh+iLow)/2;
    if( (xO >= xI[iLow]) && (xO < xI[iMid]) ) {
      iHigh = iMid;
    } else {
      iLow = iMid;
    }
    if(xO==xI[n-1]){printf("\nin interpPt\n"); return(yI[n-1]);}
    if(iHigh==iLow) {
      printf("\nexiting from binary search in 1st condtion\n");
      break;
    }
    if( (xO>=xI[iMid]) && (xO<xI[iMid+1]) ) break;
    if( (xO>=xI[iMid-1]) && (xO<xI[iMid]) ) {
      iMid--;
      break;
    }
  }     
  yO = (  ( yI[iMid+1]-yI[iMid] ) / ( xI[iMid+1]-xI[iMid] )  )
    *( xO-xI[iMid] ) + yI[iMid];
  return(yO);
}





/*
 * doInterp
 * --------
 * This function takes in 4 rays, designated by T, B, L, and R,  
 * which are Top, Bottom, Left, and Right respectively.  When 
 * propagating, the lower frequency is usually (grapichally) on
 * the Left, top freq on Right, lower lat is Bottom, upper lat is
 * Top. If the configuration changes ,the interpolation will still 
 * work.
 *
 */

void doInterp(rayT *BL, rayT *TL, rayT *BR, rayT *TR, 
	      cellT *lat_arr[][NUM_TARGS], double *EAarrPtr)
{
  double t_end, t, BL_l, TL_l, BR_l, TR_l, this_f, this_ll;
  double BL_l_prev, TL_l_prev, BR_l_prev, TR_l_prev;
  double BL_lat, TL_lat, BR_lat, TR_lat, lat_int, t_prev;
  double BL_lat_prev, TL_lat_prev, BR_lat_prev, TR_lat_prev;
  double div_lat_step, div_freq_step, div_lat, div_freq;  
  double BL_fact, TL_fact, BR_fact, TR_fact, l_int, lat_time=-1.0;
  double l_int_prev, lat_int_prev;
  int iTarg;

  
  t_end = smallest( BL->tg[ BL->numSteps-1 ],
		    TL->tg[ TL->numSteps-1 ],
		    BR->tg[ BR->numSteps-1 ],
		    TR->tg[ TR->numSteps-1 ] );

  div_lat_step = 1/DIV_LAT_NUM;
  div_freq_step = 1/DIV_FREQ_NUM;

  for(t=T_STEP ; t<=t_end ; t+= T_STEP ) {

    //printf("now doing t= %g\n", t);

    if(t==T_STEP) {	
      // this is first time step, calc prev interpolated L
      BL_l_prev = interpPt( BL->tg, BL->l_sh, BL->numSteps, 0);
      TL_l_prev = interpPt( TL->tg, TL->l_sh, TL->numSteps, 0);
      BR_l_prev = interpPt( BR->tg, BR->l_sh, BR->numSteps, 0);
      TR_l_prev = interpPt( TR->tg, TR->l_sh, TR->numSteps, 0);
    } else {
      // assign prev L values from previous loop
      BL_l_prev = BL_l; 
      TL_l_prev = TL_l; 
      BR_l_prev = BR_l; 
      TR_l_prev = TR_l; 
    }

    BL_l = interpPt( BL->tg, BL->l_sh, BL->numSteps, t);
    TL_l = interpPt( TL->tg, TL->l_sh, TL->numSteps, t);
    BR_l = interpPt( BR->tg, BR->l_sh, BR->numSteps, t);
    TR_l = interpPt( TR->tg, TR->l_sh, TR->numSteps, t);

    if( outsideL(BL_l, TL_l, BR_l, TR_l, &iTarg)  &&
	outsideL(BL_l_prev, TL_l_prev, BR_l_prev, TR_l_prev, &iTarg) ) { 
      
      // then, they're all outside L-shell of interest
    
    } else { 
      // calculate interpolated latitudes of guiding points

      t_prev = t- T_STEP;

      if( lat_time == t_prev ) {	
	// then lats already worked out at previous time step
	BL_lat_prev = BL_lat;
	TL_lat_prev = TL_lat;
	BR_lat_prev = BR_lat;
	TR_lat_prev = TR_lat;
      } else {
	BL_lat_prev = interpPt( BL->tg, BL->lat, BL->numSteps, t_prev); 
	TL_lat_prev = interpPt( TL->tg, TL->lat, TL->numSteps, t_prev);
	BR_lat_prev = interpPt( BR->tg, BR->lat, BR->numSteps, t_prev);
	TR_lat_prev = interpPt( TR->tg, TR->lat, TR->numSteps, t_prev);	
      }

      BL_lat = interpPt( BL->tg, BL->lat, BL->numSteps, t); 
      TL_lat = interpPt( TL->tg, TL->lat, TL->numSteps, t);
      BR_lat = interpPt( BR->tg, BR->lat, BR->numSteps, t);
      TR_lat = interpPt( TR->tg, TR->lat, TR->numSteps, t);
      lat_time = t;
      
      // We don't need to check if it's outside latitude of interest 
      // since we're doing the whole L-shell.  Can put this part of the 
      // code back in if we later restrict our allowable latitudes.
      //
      //if( outsideLat(BL_lat, TL_lat, BR_lat, TR_lat, &iTarg) && 
      //  outsideLat(BL_lat_prev, TL_lat_prev, 
      //	     BR_lat_prev, TR_lat_prev, &iTarg)  ) {
	
	// outside lat-range of interest - don't do anything

      //      } else {
	// do interpolations within the 4 rays!

	for(div_lat=0; div_lat<.99999; div_lat += div_lat_step) {
	  for(div_freq=0; div_freq<.99999; div_freq+= div_freq_step ) {
	    BL_fact = 1-div_lat-div_freq+div_lat*div_freq;
	    TL_fact = div_lat - div_lat*div_freq;
	    BR_fact = div_freq - div_lat*div_freq;
	    TR_fact = div_lat*div_freq;
	  
	    l_int = BL_l*BL_fact + TL_l*TL_fact + 
	      BR_l*BR_fact +	TR_l*TR_fact ;
	    
	    lat_int = BL_lat*BL_fact + TL_lat*TL_fact + 
	      BR_lat*BR_fact +	TR_lat*TR_fact ;

	    l_int_prev = BL_l_prev*BL_fact + TL_l_prev*TL_fact + 
	      BR_l_prev*BR_fact + TR_l_prev*TR_fact ;
	    
	    lat_int_prev = BL_lat_prev*BL_fact + TL_lat_prev*TL_fact + 
	      BR_lat_prev*BR_fact + TR_lat_prev*TR_fact ;
	    
	    checkCross( BL_fact, TL_fact, BR_fact, TR_fact,
			t, lat_arr, BL, TL, BR, TR, 
			div_lat, div_freq, l_int, lat_int, EAarrPtr,
			l_int_prev, lat_int_prev, iTarg );

	    /*addToLFile( BL_fact, TL_fact, BR_fact, TR_fact,
			t, outPtr, BL, TL, BR, TR, 
			div_lat, div_freq, l_int, lat_int); */
	    
	  } // for div_freq
	} // for div_lat
	//  } // else (not outside lat of interest)
    } // else (not outside L of interest)
  } // for TIME_STEP
} // END function





/*
 * FUNCTION: smallest
 * ------------------
 * Takes in 4 floating point numbers and returns the smallest
 *
 */

float smallest(float n1, float n2, float n3, float n4)
{
  float small1, small2;
  
  small1 = (n1<=n2)? n1 : n2 ;
  small2 = (n3<=n4)? n3 : n4 ;

  return ((small1<=small2)? small1 : small2);
}







/*
 * FUNCTION: writeFile
 * -------------------
 * This function just gets an appropriate name for the output file 
 * and opens it, checking that there is no problem.  Returns a 
 * pointer to the open file.
 *
 */

FILE *writeFile(long lower_freq, char *pref, char *suff)
{
  char outFileName[64];
  FILE *outPtr;

  // make some kind of filename convention, return pointer
  sprintf( outFileName, "%s%d%s.dat", pref, lower_freq, suff);  

  if((outPtr=fopen(outFileName,"w"))== NULL) {
    printf("Problem opening output file!\n");
    exit(1);
  } else {
    printf("File %s \nopened successfully for writing!\n", outFileName);
  }
  return outPtr;
}






/*
 * FUNCTION: addToLFile
 * --------------------
 * This function is called from doInterp when a certain interpolated  
 * ray is within the right bounds and the point needs to be recorded
 * in the output file. 
 *
 * Parameters to record are:
 * 1. tg	- group time
 * 2. lat	- local latitude on this L shell
 * 3. f		- ray frequency
 * 4. ll	- launch latitude
 * 5. pwr	- landau damped ??? Initial power computed ???
 * 6. psi	- wn wrt +ve B
 * 7. l_int	- L value of ray
 *
 */

void addToLFile(double BL_fact, double TL_fact, double BR_fact,
		double TR_fact, double t, FILE *outPtr, 
		rayT *BL, rayT *TL, rayT *BR, rayT *TR, 
		double div_lat, double div_freq, 
		double l_int, double lat)
{
  long f;
  double ll, psi, pwr, dlat, dfreq, avePwr;
  
  // ray frequency
  f = TL->f + div_freq*( TR->f - TL->f );
  
  // ray launch latitude
  ll = BL->lat[0] + div_lat*( TL->lat[0] - BL->lat[0] );

    
  // ray wn angle wrt B0
  psi = ( interpPt( BL->tg, BL->psi, BL->numSteps, t)*BL_fact + 
	  interpPt( TL->tg, TL->psi, TL->numSteps, t)*TL_fact +
	  interpPt( BR->tg, BR->psi, BR->numSteps, t)*BR_fact +
	  interpPt( TR->tg, TR->psi, TR->numSteps, t)*TR_fact  );

  // ray pwr with Landau damping
  pwr = ( interpPt( BL->tg, BL->pwr, BL->numSteps, t)*BL_fact + 
	  interpPt( TL->tg, TL->pwr, TL->numSteps, t)*TL_fact +
	  interpPt( BR->tg, BR->pwr, BR->numSteps, t)*BR_fact +
	  interpPt( TR->tg, TR->pwr, TR->numSteps, t)*TR_fact  );
  
  // print information to file if greater than 1e-6 of initial pwr
  avePwr = ( TL->pwr[0] + BL->pwr[0] + TR->pwr[0] + BR->pwr[0] )/4.0;
  if( pwr > (1e-6 * avePwr) )
  fprintf(outPtr,"%g %g %d %g %g %g %g\n",t,lat,f,ll,pwr,psi,l_int);
  
}






/*
 * FUNCTION: ltgPwr
 * ----------------
 * This function calculates the initial power of the ray at 1000km
 * taking into account ionospheric absorption. Inputs are:
 * lam_s:	lightning source latitude
 * dlat:	latitudinal separation between adjacent rays
 * dfreq:	frequency separation between adjacent rays
 * f:		freq of this particular ray*
 * lat:		latitude of this particular ray 
 * dlong:	long offset from center to get rid of central null
 *
 */

double ltgPwr(float lam_s, float dlat, float dfreq, long f, 
	      float lat, float dlong )
{
 
  double w, w_sq, dist_lat, dist_long, dist_iono;
  double S, S_vert, dist_tot, xi, attn_factor;  


  // calc dist on 1 long, then on 1 lat, then total ground dist,
  // then tot distance to sub-iono point
  dist_lat = (R_E+H_IONO/2) * fabs(lam_s-lat) * D2R ;
  dist_long = (R_E+H_IONO/2) * dlong * D2R ;

  dist_iono = hypot(dist_lat, dist_long) ;
  dist_tot = hypot(dist_iono, H_IONO) ;
  
  xi = atan2(dist_iono, H_IONO);

  // Pwr at sub-iono point
  w = 2*PI*f;
  w_sq =  pow( w , 2 );
  S = ( (1/Z0) * pow( (H_E*I0*2E-7*(sin(xi)/dist_tot)*w*(A-B)) , 2 ) 
	/  (  (w_sq+pow(A,2))*(w_sq+pow(B,2))  )      ) ;
  S_vert = S * cos(xi) ;	// factor for vert prop.
  
  // now DO THE ATTENUATION
  attn_factor = pow(10,-(ionoAbsorp(lat,f)/10)  );
  S_vert = S_vert * attn_factor ;


  // Integrate wrt space, assume slice-width=1m.  The correction 
  // factor 0.877... comes about so that its a 1m slice at 1000km
  // which means that at 100km its only 0.877m wide
  return ( S_vert * dlat*D2R*(R_E+H_IONO) * dfreq * 0.87788331); 
}






/*
 * FUNCTION: ionoAbsorp
 * --------------------
 * This function just calculates the ionospheric absorption of a 
 * ray in moving through the ionosphere from ~100km to ~1000km.
 * Uses data from Helliwell 1965, fig. 3-35, and interpolates for 
 * latitudes between 0-90 degrees, then uses data at 2 kHz, and 
 * 20 kHz to interpolate/extrapolate for other frequencies.
 *
 */

double ionoAbsorp(float lat, long f)
{
  float db2i, db20i;
  double  db2iLog, db20iLog, m, c; 

  float lats[46] = {	0.0,	2.0,	4.0,	6.0,	8.0,
			10.0,	12.0,	14.0,	16.0,	18.0,
			20.0,	22.0,	24.0,	26.0,	28.0,
			30.0,	32.0,	34.0,	36.0,	38.0,	
			40.0,	42.0,	44.0,	46.0,	48.0,	
			50.0,	52.0,	54.0,	56.0,	58.0,
			60.0,	62.0,	64.0,	66.0,	68.0,
			70.0,	72.0,	74.0,	76.0,	78.0,	
			80.0,	82.0,	84.0,	86.0,	88.0,	
			90.0	};

  float db20kHzLog[46] = {
    2.72890330135884,
    2.61413531724660,
    2.49936733313435,
    2.38459934902211,
    2.26983136490987,
    2.15506338079763,
    2.04029539668538,
    1.92531644206607,
    1.79924052160795,
    1.65746834723986,
    1.49802217224091,
    1.36186709269283,
    1.25598388505325,
    1.15652386401487,
    1.07820410788996,
    1.01300984614119,
    0.95644081081742,
    0.90646976883951,
    0.85864979711507,
    0.81518987704864,
    0.77367088822273,
    0.73407595275903,
    0.70105485532855,
    0.67236287394032,
    0.64568469779523,
    0.62025316330165,
    0.59707886563282,
    0.57435687349111,
    0.55630869316853,
    0.53826051284596,
    0.52280663628542,
    0.50796595891253,
    0.49540084969300,
    0.48392405611603,
    0.47327480980675,
    0.46355655505289,
    0.45432489421971,
    0.44715189958245,
    0.43997890494519,
    0.43289667529356,
    0.42589049202632,
    0.41888430875909,
    0.41279384916331,
    0.40818263512987,
    0.40428972648146,
    0.40189872848716  };

  float db2kHzLog[46] = {
    2.20354296279991,
    2.10139102288977,
    1.99923908297962,
    1.89708714306948,
    1.79493520315934,
    1.69278326324919,
    1.59070536198610,
    1.47122984416933,
    1.32946932452408,
    1.19113861945499,
    1.04312210740998,
    0.92225656425647,
    0.81909999940519,
    0.73129314437689,
    0.65802549162668,
    0.59490689076827,
    0.54352690196227,
    0.49737249399863,
    0.45587327008181,
    0.41437404616500,
    0.38243818008768,
    0.35131375279102,
    0.32279494989472,
    0.29651211576080,
    0.27148893501543,
    0.24935600431553,
    0.22747724000215,
    0.20903314606947,
    0.19058905213679,
    0.17380919954012,
    0.15966172931239,
    0.14551425908466,
    0.13161845036475,
    0.11976152945062,
    0.10790460853649,
    0.09604768762236,
    0.08688693663841,
    0.07777734644718,
    0.06866775625595,
    0.05958751374305,
    0.05050956953555,
    0.04345615362238,
    0.03810140109787,
    0.03270477356885,
    0.02717155386097,
    0.02163833415309   };
    


  db2iLog  = interpPt(lats, db2kHzLog  , 46, lat);
  db20iLog = interpPt(lats, db20kHzLog , 46, lat);
  
  m = (db20iLog - db2iLog);	// officially should be 
				// m'=m/( log10(20)-log10(2) )
				// but denominator = 1

  c = (db2iLog+db20iLog)/2 - m * 0.8010; // this is just
					 // (log10(2)+log10(20))/2

  // now extrapolate to desired f value in log10 domain
  // get 10^result, and return it
  return ( pow(10.0 , ( m*log10( (f/1000.0) ) + c ) ) );
}










/*
 * FUNCTION: initPwr
 * -----------------
 * This function is responsible for initialising the power of the 
 * 4 "outer" rays, within which we will be interpolating.  It first 
 * finds the normalising constant AN, then the initial (integrated) 
 * lightning power at the top of the ionosphere, and finally cal-
 * culates the damnping as the ray traverses the magnetosphere.
 *
 */

void initPwr(rayT *BL, rayT *TL, rayT *BR, rayT *TR)
{
  double dlat, dfreq;
  int i, j;


  // calculate ray power and interpolate
  dlat = (TL->lat[0] - BL->lat[0])/ DIV_LAT_NUM ;
  dfreq = (TR->f - TL->f) / DIV_FREQ_NUM ;
 
  // Some potential debugging
  if(dlat<0 || dfreq<0) {
    printf("\n\aProblem: dlat or dfreq are negative!!!???\n");
    printf("TL->lat[0]: %g, BL->lat[0]: %g \n",TL->lat[0],BL->lat[0]);
    printf("TR->f: %d, TL->f: %d \n",TR->f,TL->f);
    printf("DIV_LAT_NUM: %g, DIV_FREQ_NUM: %g \n",DIV_LAT_NUM,DIV_FREQ_NUM);
    printf("dlat: %g, dfreq: %g\n\n", dlat, dfreq);
    printf("Perhaps NUM_STEPS is too small ... \n");
    for(i=0;i<1;i++) 
    printf( "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n", 
	    TL->tg[i], TL->distre[i], TL->lat[i], 
	    TL->delta[i], TL->t[i], TL->l_sh[i],
	    TL->psi[i], TL->psiray[i], TL->psires[i], 
	    TL->mu[i], TL->dens[i], TL->anH[i], 
	    TL->anHe[i], TL->anO[i], TL->fH[i], 
	    TL->stixP[i], TL->stixR[i], 
	    TL->stixL[i] );
  }


  // determine maximum power in the simulation run (assume 4kHz)
  MAX_POWER = ltgPwr(BL->lam_s,dlat,dfreq,4000,BL->lam_s,0.7);


  if(BL->pwr[0]>=.9) {
    BL->pwr[0]=ltgPwr(BL->lam_s,dlat,dfreq,BL->f,BL->lat[0],0.7);
    for(j=1;j<BL->numSteps;j++) BL->pwr[j]=BL->pwr[j]*BL->pwr[0];
  }

  if(TL->pwr[0]>=.9) {
    TL->pwr[0]=ltgPwr(TL->lam_s,dlat,dfreq,TL->f,TL->lat[0],0.7);
    for(j=1;j<TL->numSteps;j++) TL->pwr[j]=TL->pwr[j]*TL->pwr[0];
  }

  if(BR->pwr[0]>=.9) {
    BR->pwr[0]=ltgPwr(BR->lam_s,dlat,dfreq,BR->f,BR->lat[0],0.7);
    for(j=1;j<BR->numSteps;j++) BR->pwr[j]=BR->pwr[j]*BR->pwr[0];
  }

  if(TR->pwr[0]>=.9) {
    TR->pwr[0]=ltgPwr(TR->lam_s,dlat,dfreq,TR->f,TR->lat[0],0.7);
    for(j=1;j<TR->numSteps;j++) TR->pwr[j]=TR->pwr[j]*TR->pwr[0];
  }

}






/*
 * FUNCTION: outsideL
 * ------------------
 * This function simply evaluates if the current 4 surrounding ray
 * points are possibly within the L region of interest, i.e. if any one 
 * of the 4 rays is inside in any way.  It returns 1 if all are outside
 * or 0 if something is potentially inside.
 *
 */
int outsideL(double BL_l, double TL_l, double BR_l, double TR_l, int *iTarg)
{
  int i;
  double lowBound, upBound;

  for(i=0; i<NUM_TARGS; i++) {
    lowBound = L_TARG - L_MARGIN[i];
    upBound = L_TARG + L_MARGIN[i];    
    if(  ( (BL_l<=lowBound) &&		
	    (TL_l<=lowBound) &&		// all too low
	    (BR_l<=lowBound) &&
	    (TR_l<=lowBound)  ) ||
	  ( (BL_l>=upBound) &&		// all too high
	    (TL_l>=upBound) &&
	    (BR_l>=upBound) &&
	    (TR_l>=upBound) )   ) {
      // yes - all rays are OUTSIDE L of interest
    } else { 
      *iTarg = i; // tell index of L-shell we're in
      return 0; // FALSE - one point not outside 
    }
  }  // for i ...
  return 1; // TRUE - all points are outside
}







/*
 * FUNCTION: outsideLat
 * --------------------
 * This function simply evaluates if the current 4 surrounding ray
 * points are possibly within the latitudinal region of interest, 
 * i.e. if any one of the 4 rays is inside in any way.  It returns 
 * 1 if all are outside or 0 if something is potentially inside. 
 * Note, this function checks all LAT_TARG's of a given L_TARG, so 
 * if we have L_TARG=2 with LAT_TARG=0 and LAT_TARG=30, it will 
 * check both and return iTarg which matches LAT and L targets. 
 *
 */
int outsideLat(double BL_lat, double TL_lat, double BR_lat, 
	       double TR_lat, int *iTarg)
{
  double lowBound, upBound;
  int i;

  for(i=0; i<NUM_TARGS; i++) {
    
    lowBound = LAT_TARG[i] - LAT_MARGIN[i];
    upBound = LAT_TARG[i] + LAT_MARGIN[i];    
    
    if(  ( (BL_lat<=lowBound) &&		
	   (TL_lat<=lowBound) &&	// all too low
	   (BR_lat<=lowBound) &&
	   (TR_lat<=lowBound)  ) ||
	 ( (BL_lat>=upBound) &&	// all too high
	   (TL_lat>=upBound) &&
	   (BR_lat>=upBound) &&
	   (TR_lat>=upBound) )   ) {
      // All rays are outside our lat-range of interest
    } else { 
      *iTarg = i;
      return 0; // if one is inside, return FALSE (i.e. NOT outside)
      }	// if BL_lat, TL_lat ... outside lat range
  }	// for i ...
  
  return 1; // TRUE - all points outside lat-range of interest
  
}







/*
 * FUNCTION: initEA_Arr
 * --------------------
 * This function creates an array of effective area (EA) coordinates.  It
 * divides up the space between EALimN and EALimS into EAIncr-sized pieces 
 * to create the latitude list (lam).  At each lam, the function computes 
 * the (x,y) corodinates of the ends of a line-segment which is perpen-
 * dicular to the field line, and dL0 L-shells wide.  The final array has
 * columns:	
 *		0	1	2	3	4     
 *		lam	x1	x2	y1	y2  	
 *
 *			5	6	7	8	       
 *			EA_a	EA_b	EA_c	EA_length	
 *								
 * The parameters EA_a, EA_b_, EA_c describe the straight line passing 
 * through (x1,y1) (x2,y2), from: EA_a*y + EA_b*x + EA_c = 0, and EA_length
 * is just the length of the aperture.
 *
 */
double *initEA_Arr(void)
{
  double *arrPtr;
  int i, k;
  double lam, slam, clam, dL_lam, ptR, ptX, ptY, x1, x2, y1, y2;
  double slam2, clam2, rootTerm, x_unit_vect, y_unit_vect, ptL;
  double EA_a, EA_b, EA_c, EA_length;
  
  EA_NUMROWS = (EALimN - EALimS) / EAIncr + 1 ;
  EA_NUMCOLS = 1 + 8*NUM_TARGS;

  printf("rows: %d, columns: %d \n", EA_NUMROWS, EA_NUMCOLS);

  if( !(arrPtr = (double *)malloc( EA_NUMROWS * 
				   EA_NUMCOLS * 
				   sizeof(double) ) ) ){
    printf("Problem allocating memory in initEA_Arr");
  }
  
  
  ptL = L_TARG;
  
  for( i=0 ; i<EA_NUMROWS ; i++ ) {
    
    lam = EALimS + i*EAIncr;
    
    clam = cos(lam*D2R);
    slam = sin(lam*D2R);
    clam2 = pow(clam,2);
    slam2 = pow(slam,2);
    rootTerm = sqrt(1+3*slam2);
      
    dL_lam = clam2*clam / rootTerm * dL0 ;
    
    x_unit_vect = (3*clam2 - 2) / rootTerm ;
    y_unit_vect = (3*slam*clam) / rootTerm ;
    
    ptR = ptL*clam2;
    ptX = ptR*clam;
    ptY = ptR*slam;
      
    
    x1 = ptX - x_unit_vect*dL_lam ;
    x2 = ptX + x_unit_vect*dL_lam ;
    
    y1 = ptY - y_unit_vect*dL_lam ;
    y2 = ptY + y_unit_vect*dL_lam ;
    
    EA_a = y1 - y2;
    EA_b = x2 - x1;
    EA_c = x1*y2 - y1*x2;
    
    EA_length = sqrt( (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2) )*R_E;  
    
    writeEntry(i, 0, lam, arrPtr);
    writeEntry(i, 1, x1, arrPtr);
    writeEntry(i, 2, x2, arrPtr);
    writeEntry(i, 3, y1, arrPtr);
    writeEntry(i, 4, y2, arrPtr);
    writeEntry(i, 5, EA_a, arrPtr);
    writeEntry(i, 6, EA_b, arrPtr);
    writeEntry(i, 7, EA_c, arrPtr);
    writeEntry(i, 8, EA_length, arrPtr);
    
  }
  return arrPtr;
}







/*
 * FUNCTION: writeEntry
 * --------------------
 * This is a simple function that writes an entry into a dynamically 
 * defined 'double' array, using regular [row, column] notation.
 *
 */
void writeEntry(int row, int col, double entry, double *arrPtr)
{
  *(arrPtr + col + row*EA_NUMCOLS)	= entry;
}







/*
 * FUNCTION: getEntry
 * --------------------
 * This is a simple function that gets an entry from a dynamically 
 * defined 'double' array, using regular [row, column] notation.
 *
 */
double getEntry(int row, int col, double *arrPtr)
{
  return *(arrPtr + col + row*EA_NUMCOLS);
}







/*
 * FUNCTION: checkCross
 * --------------------
 * This function checks if there is one or more EA's between the start
 * and end latitudes of the line segment.  If there is, then the function
 * actually checks to see if there is a real EA crossing and if there is,
 * it records this entry into the lat_arr data structure, by calling the
 * addToArr function.
 *
 */
void checkCross(double BL_fact, double TL_fact, double BR_fact,
		double TR_fact, double t, cellT *lat_arr[][NUM_TARGS], 
		rayT *BL, rayT *TL, rayT *BR, rayT *TR, 
		double div_lat, double div_freq, 
		double l_int, double lat_int, double *EAarrPtr,
		double l_int_prev, double lat_int_prev, int iTarg)
{
  double lat_high, lat_low;
  double x1ray, x2ray, y1ray, y2ray, r1ray, r2ray, Aray, Bray, Cray;
  double x1EA, x2EA, y1EA, y2EA, A_EA, B_EA, C_EA, EA_length; 
  double val1, val2, val3, val4;
  int EA_iLow, EA_iHigh, EA_i;
  long f;
  double ll, psi, pwr, mu, avePwr, stixP, stixR, stixL;

  // First, check if there is an EA between latitude limits
  // Order the limits
  if(lat_int_prev <= lat_int) {
    lat_low = lat_int_prev;
    lat_high = lat_int;
  } else {
   lat_high = lat_int_prev;
   lat_low = lat_int;
  }

  // get index numbers into EA array
  EA_iLow  = (int)ceil(  (lat_low - EALimS) / EAIncr );
  EA_iHigh = (int)floor( (lat_high - EALimS) / EAIncr );

  // if Low > High then there are no EA between lat limits
  //printf("latLow: %f, latHigh: %f\n",lat_low, lat_high);
  //printf("EA_iLow: %d, EA_iHigh: %d\n",EA_iLow, EA_iHigh);
  if(	( EA_iLow > EA_iHigh )	|| 
	( EA_iLow < 0 )		|| 
	( EA_iHigh > (EALimN-EALimS)/EAIncr )  ) return;

  //  printf("made it past check\n");
  // Go through all EA's between the limits and check crossings
  for( EA_i=EA_iLow ; EA_i <= EA_iHigh ; EA_i++ ) {

    r1ray = l_int_prev * pow( cos(lat_int_prev*D2R) , 2 ) ;
    x1ray = r1ray * cos(lat_int_prev*D2R) ;
    y1ray = r1ray * sin(lat_int_prev*D2R) ;

    r2ray = l_int * pow( cos(lat_int*D2R) , 2 ) ;
    x2ray = r1ray * cos(lat_int*D2R) ;
    y2ray = r1ray * sin(lat_int*D2R) ;

    Aray = y1ray - y2ray;
    Bray = x2ray - x1ray;
    Cray = x1ray*y2ray - y1ray*x2ray;
  
    x1EA = getEntry(EA_i, iTarg*8+1, EAarrPtr);
    x2EA = getEntry(EA_i, iTarg*8+2, EAarrPtr);
    y1EA = getEntry(EA_i, iTarg*8+3, EAarrPtr);
    y2EA = getEntry(EA_i, iTarg*8+4, EAarrPtr);
    A_EA = getEntry(EA_i, iTarg*8+5, EAarrPtr);
    B_EA = getEntry(EA_i, iTarg*8+6, EAarrPtr);
    C_EA = getEntry(EA_i, iTarg*8+7, EAarrPtr);    
    EA_length = getEntry(EA_i, iTarg*8+8, EAarrPtr);

    val1 = Aray*x1EA + Bray*y1EA + Cray;
    val2 = Aray*x2EA + Bray*y1EA + Cray;
    val3 = A_EA*x1ray + B_EA*y1ray + C_EA;
    val4 = A_EA*x2ray + B_EA*y2ray + C_EA;   

    if( (val1*val2 <= 0.0) && (val3*val4 <= 0.0) ) {
      // OKAY! This is a genuine crossing - figure out what to do next!
      //printf("yes!  It's a crossing\n");
      //if(val3*val4 ==0) printf("ray is ON the EA\n");

      // ray frequency
      f = TL->f + div_freq*( TR->f - TL->f );
      
      // ray launch latitude
      ll = BL->lat[0] + div_lat*( TL->lat[0] - BL->lat[0] );
      
      
      // ray wn angle wrt B0
      psi = ( interpPt( BL->tg, BL->psi, BL->numSteps, t)*BL_fact + 
	      interpPt( TL->tg, TL->psi, TL->numSteps, t)*TL_fact +
	      interpPt( BR->tg, BR->psi, BR->numSteps, t)*BR_fact +
	      interpPt( TR->tg, TR->psi, TR->numSteps, t)*TR_fact  );
      
      // ray pwr with Landau damping - divided by length of EA
      pwr = ( interpPt( BL->tg, BL->pwr, BL->numSteps, t)*BL_fact + 
	      interpPt( TL->tg, TL->pwr, TL->numSteps, t)*TL_fact +
	      interpPt( BR->tg, BR->pwr, BR->numSteps, t)*BR_fact +
	      interpPt( TR->tg, TR->pwr, TR->numSteps, t)*TR_fact  )
              / EA_length ;  // NB - we divided by EA length!
      
      // refractive index
      mu =  ( interpPt( BL->tg, BL->mu, BL->numSteps, t)*BL_fact + 
	      interpPt( TL->tg, TL->mu, TL->numSteps, t)*TL_fact +
	      interpPt( BR->tg, BR->mu, BR->numSteps, t)*BR_fact +
	      interpPt( TR->tg, TR->mu, TR->numSteps, t)*TR_fact  ); 
      
      // Stix P
      stixP=( interpPt( BL->tg, BL->stixP, BL->numSteps, t)*BL_fact + 
	      interpPt( TL->tg, TL->stixP, TL->numSteps, t)*TL_fact +
	      interpPt( BR->tg, BR->stixP, BR->numSteps, t)*BR_fact +
	      interpPt( TR->tg, TR->stixP, TR->numSteps, t)*TR_fact  ); 
      // Stix R
      stixR=( interpPt( BL->tg, BL->stixR, BL->numSteps, t)*BL_fact + 
	      interpPt( TL->tg, TL->stixR, TL->numSteps, t)*TL_fact +
	      interpPt( BR->tg, BR->stixR, BR->numSteps, t)*BR_fact +
	      interpPt( TR->tg, TR->stixR, TR->numSteps, t)*TR_fact  ); 
      
      // Stix L
      stixL=( interpPt( BL->tg, BL->stixL, BL->numSteps, t)*BL_fact + 
	      interpPt( TL->tg, TL->stixL, TL->numSteps, t)*TL_fact +
	      interpPt( BR->tg, BR->stixL, BR->numSteps, t)*BR_fact +
	      interpPt( TR->tg, TR->stixL, TR->numSteps, t)*TR_fact  ); 
       
      // print information to file if greater than 1e-6 of initial pwr
      avePwr = ( TL->pwr[0] + BL->pwr[0] + TR->pwr[0] + BR->pwr[0] )/4.0;
      //if( pwr > (1e-6 * avePwr) )
      //fprintf(outPtr,"%g %g %d %g %g %g %g\n",t,lat,f,ll,pwr,psi,l_int); 
      addToArr(lat_arr, t, f, pwr, psi, mu, EA_i, iTarg, stixP, stixR, stixL);
    }
  }
}




/* 
 * FUNCTION: initLatArr
 * --------------------
 * This function takes the newly created array of cellT pointers, and 
 * initialises them all to point to NULL.
 *
 */
void initLatArr(cellT *lat_arr[][NUM_TARGS])
{
  int i, j;

  for(i=0 ; i < NUMLATS ; i++) {
    for(j=0; j < NUM_TARGS ; j++) {
      lat_arr[i][j] = NULL;
    }
  }
}




/*
 * FUNCTION: addToArr
 * ------------------
 * Given the parameters fora certain entry, i.e., f, t, power, etc. this 
 * function finds the right place in the data array for the new entry and 
 * either creates a new cell in the correct place or updates the correct
 * cell with the new entry.
 *
 */
void addToArr(cellT *lat_arr[][NUM_TARGS],double t, double f, 
	      double pwr, double psi, double mu, int lami, 
	      int iTarg, double stixP, double stixR, double stixL)
{

  cellT	*cell, *next;
  double t_grid, f_grid, lat;
  int comp, comp_next;


  t_grid = floor(t/DT)*DT ;
  f_grid = floor(f/DF)*DF ;


  // CASE 1: if NO entries exist at lami
  
  if(lat_arr[lami][iTarg] == NULL) {

    cell = (cellT *)malloc(sizeof(cellT));
    if(cell==NULL) printf("Problem in addToArr in allocating memory\n");

    // do all the pointer "hook-ups" and enter in values
    cell->link		= NULL;
    lat_arr[lami][iTarg]= cell; 
    cell->t		= t_grid;
    cell->f		= f_grid;
    cell->pwr		= pwr;
    cell->psi		= psi;
    cell->mu		= mu ;
    cell->stixP		= stixP;
    cell->stixR		= stixR;
    cell->stixL		= stixL;
    cell->num_rays	= 1;

    return;
  }


  // CASE 2: if there are entries at lami 
  //	but this entry must be inserted first
  
  next = lat_arr[lami][iTarg];
  if( compare(t_grid, f_grid, next->t, next->f ) == -1 ) {

    cell = (cellT *)malloc(sizeof(cellT));
    if(cell==NULL) printf("Problem in addToArr in allocating memory\n");

    // do all the pointer "hook-ups" and enter in values
    cell->link		= next;
    lat_arr[lami][iTarg]= cell; 
    cell->t		= t_grid;
    cell->f		= f_grid;
    cell->pwr		= pwr;
    cell->psi		= psi;
    cell->mu		= mu; 
    cell->stixP		= stixP;
    cell->stixR		= stixR;
    cell->stixL		= stixL;
    cell->num_rays	= 1;

    return;
  }


  // CASE 3: if there are entries at lami
  //	and this entry is not first or exists

  next = lat_arr[lami][iTarg];

  while( next != NULL ) {

    comp = compare(t_grid,f_grid, next->t, next->f);


    // Cell already exists, we must update it
    if(comp==0) {

      next->pwr		+= pwr;
      next->psi		+= psi;
      next->mu		+= mu ;
      next->stixP	+= stixP;      
      next->stixR	+= stixR;
      next->stixL	+= stixL;
      next->num_rays	+= 1;  
      break;
    }	//   if( comp == 0)


    // we're checking the last cell and cell does not exist
    if(comp==1 && (next->link == NULL) ) {

      cell = (cellT *)malloc(sizeof(cellT));
      if(cell==NULL) printf("Problem in addToArr in allocating memory\n");

      // do all the pointer "hook-ups" and enter in values
      cell->link	= NULL;
      next->link	= cell; 
      cell->t		= t_grid;
      cell->f		= f_grid;
      cell->pwr		= pwr;
      cell->psi		= psi;
      cell->mu		= mu ;
      cell->stixP	= stixP;
      cell->stixR	= stixR;
      cell->stixL	= stixL;
      cell->num_rays	= 1;      
      break;
    }	//   if( comp ==1 && next->link==NULL )


    // the current entry is between this cell and the next
    comp_next = compare(t_grid, f_grid, next->link->t, next->link->f);
    if(comp==1 && comp_next==-1 ) {

      cell = (cellT *)malloc(sizeof(cellT));
      if(cell==NULL) printf("Problem in addToArr in allocating memory\n");
      
      // do all the pointer "hook-ups" and enter in values
      cell->link	= next->link;
      next->link	= cell; 
      cell->t		= t_grid;
      cell->f		= f_grid;
      cell->pwr		= pwr;
      cell->psi		= psi;
      cell->mu		= mu ;
      cell->stixP	= stixP;
      cell->stixR	= stixR;
      cell->stixL	= stixL;
      cell->num_rays	= 1;      
      break;
    }	//   if( comp ==1 && comp_next == -1 )

    next = next->link;

  }	// while loop
}


 




/*
 * FUNCTION: compare
 * -----------------
 * This function takes in a (t,f) pair and compares it to some given 
 * (f_next, t_next) pair.  This is necessary when searching for the 
 * right place to put a cell into the array.
 *
 */
int compare(double t_try, double f_try, double t_next, double f_next)
{

  if(t_try < t_next) return -1;

  if(t_try > t_next) return 1;

  if(t_try == t_next) {

    if(f_try < f_next) return -1;    
    if(f_try > f_next) return 1;
    if(f_try == f_next) return 0;
  }

}







/*
 * FUNCTION: dispLatArr
 * --------------------
 * This function just goes through the whole array and prints out what we 
 * have entered into it.
 *
 */
void dispLatArr(cellT *lat_arr[][NUM_TARGS],  FILE *outPtr)
{
  int i;
  cellT *next;
  double lat;

  
  for(i=0; i<NUMLATS; i++) {
    next = lat_arr[i][0];
    lat = EALimS + EAIncr*i;
   
    if( fabs(lat)<0.01 ) {
      while(next != NULL) {
	
	fprintf(outPtr, "%g %g %g %g %g %g %g %g %g %g\n", 
		L_TARG,
		lat,
		next->t,
		next->f,
		((next->pwr)/(next->num_rays)/DF/DT),
		((next->psi)/(next->num_rays)),
		((next->mu)/(next->num_rays)),
		((next->stixP)/(next->num_rays)),
		((next->stixR)/(next->num_rays)),
		((next->stixL)/(next->num_rays))   );
      /*
	
      fwrite( &(L_TARG[j]), sizeof(double), 1, outPtr );
      fwrite( &(lat), sizeof(double), 1, outPtr );
      fwrite( &(next->t), sizeof(double), 1, outPtr );
      fwrite( &(next->f), sizeof(double), 1, outPtr );
      fwrite( &(next->pwr), sizeof(double), 1, outPtr );
      fwrite( &(next->psi), sizeof(double), 1, outPtr );
      fwrite( &(next->mu), sizeof(double), 1, outPtr );
      fwrite( &(next->num_rays), sizeof(int), 1, outPtr );
      */
	
	
	next = next->link;
      }  // while loop
    }  // ONLY DO EQUATORIAL CROSSING
  }  // for i
}






/*
 * FUNCTION: calcRes
 * -----------------
 * This function will calculate the resonant pitch angle change of the 
 * particle using the lat_arr data  
 *
 */
void calcRes(cellT *lat_arr[][NUM_TARGS], long lower_freq)
{
  FILE *resPtrN, *resPtrS;
  int i, j=0, kk, mres, noutN, noutS, ei, ti, e_toti;
  cellT *next;
  char *prefN="pN", *prefS="pS", suff[64];
  double lat, L, t, f, pwr, psi, mu, stixP, stixR, stixL, latk;
  double Bxw, Byw, Bzw, Exw, Eyw, Ezw, stixD, stixS, stixA;
  double stixB, stixX, n_x, n_z, k, kx, kz, rho1, rho2, Byw_sq;
  double flt_const_N[EA_SPLIT], flt_const_S[EA_SPLIT], flt_time, eta_dot;
  double wh, dwh_ds, gamma, alpha1, alpha2, beta, v_para, v_perp;
  double spsi, cpsi, spsi_sq, cpsi_sq, mu_sq, w, R1, R2, w1, w2;
  double alpha_lc, alpha_eq, epsm, slat, clat, slat_term, ds;
  double t1, t2, t3, direction, v_para_res, v_tot, v_tot_res;
  double salph, calph, wtau_sq, Y, dv_para_ds, AA, BB, T1;
  double Farg, Farg0, Fs, Fc, Fs0, Fc0, dFs_sq, dFc_sq, dalpha;
  double alpha_eq_p, dalpha_eq, e_starti, e_endi, E_res;
  double v_tot_arr[NUM_E], E_tot_arr[NUM_E], v_para_star, v_para_star_sq;
  float *arr_N, *arr_S;
  time_t start, end;

 
  start = time(NULL);
  arr_N = NULL;
  arr_S = NULL;
  
  L = L_TARG;
  epsm = (1/L)*(R_E+H_IONO)/R_E;
  alpha_eq = asin(sqrt( pow(epsm,3)/sqrt(1+3*(1-epsm)) ));
  sprintf(suff, "_%g", L);
  resPtrN = writeFile(lower_freq, prefN, suff);
  resPtrS = writeFile(lower_freq, prefS, suff);
  printf("\nNow doing L: %g\n", L);
  
  //initialize the velocity and energy arrays
  for(i=0; i<NUM_E; i++) {
    E_tot_arr[i] = pow(10, (E_EXP_BOT+ DE_EXP*i) ); // energy in eV
    v_tot_arr[i] = C*sqrt(1 - pow( (E_EL/(E_EL+E_tot_arr[i])) ,2) );
  }


  // Go through all latitudes
  for(i=0; i<NUMLATS; i++) {
    next = lat_arr[i][j];
    
    // CALCULATE:
    // wh
    // dwh/ds
    // flight-time constant
    // alpha_lc
    
    lat = EALimS + EAIncr*i ;
    slat = sin( lat*D2R );
    clat = cos( lat*D2R );
    slat_term = sqrt(1+3*slat*slat);
    wh = 2*PI*880000/pow(L,3)*slat_term/pow(clat,6);
    dwh_ds = 3*wh/(L*R_E)*slat/slat_term*
      (1/(slat_term*slat_term) + 2/(clat*clat));

    for(kk=0; kk < EA_SPLIT; kk++) {
      latk = (lat-.5*EAIncr+.5*EAIncr/EA_SPLIT)+kk*(EAIncr/EA_SPLIT);
      getFltConst(L,latk,alpha_eq,&(flt_const_N[kk]),&(flt_const_S[kk]));
    }

    alpha_lc = asin(sqrt( slat_term/pow(clat,6) )*sin(alpha_eq));
    salph = sin(alpha_lc);
    calph = cos(alpha_lc);
    // ds = L*R_E* slat_term*clat*EAIncr/180*PI;    
    ds = L*R_E* slat_term*clat*EAIncr/180*PI * MULT ; 
    dv_para_ds = -0.5*pow(salph,2)/calph/wh*dwh_ds;
    
    printf("\nLat: %g", lat);
    //printf("lat: %g, wh: %g, dwh_ds: %g, flt_const_N: %g \n",
    //     lat, wh, dwh_ds, flt_const_N);
    //printf("flt_const_S: %g, alpha_lc: %g, \nds: %g, dv_para_ds: %g\n", 
    //     flt_const_S, alpha_lc, ds, dv_para_ds);
    
    
    // Go through the cells in each latitude
    while(next != NULL) {
      
      t = next->t + DT/2;	// We want the time and freq to be in the 
      f = next->f + DF/2;	// center of the cell, so add DT/2 or DF/2
      pwr = (next->pwr)/(next->num_rays)/DT; // <-I know this looks wrong 
      psi = (next->psi)/(next->num_rays)*D2R;   // but it's right!
						// Grrr... back off!


      if(pwr > 1.0e-50) {


      mu  = (next->mu)/(next->num_rays);
      stixP =	(next->stixP)/(next->num_rays);
      stixR = (next->stixR)/(next->num_rays);
      stixL =	(next->stixL)/(next->num_rays);
      
      //printf("t: %g, f: %g, pwr: %g, psi: %g,\nmu: %g, stixP: %g, stixR:
      //       %g, stixL: %g \n\n", 
      //     t, f, pwr, psi, mu, stixP, stixR, stixL);
      
      
      spsi = sin(psi);
      cpsi = cos(psi);
      spsi_sq = pow(spsi,2);
      cpsi_sq = pow(cpsi,2);
      n_x = mu*fabs(spsi);
      n_z = mu*cpsi;
      mu_sq = mu*mu;
      w = 2.0*PI*f;
      k = w*mu/C;
      kx = w*n_x/C;
      kz = w*n_z/C;
      Y = wh / w ;
      stixS = ( stixR + stixL ) /2.0;
      stixD = ( stixR - stixL ) /2.0;
      stixA = stixS + (stixP-stixS)*cpsi_sq;
      stixB = stixP*stixS+stixR*stixL+(stixP*stixS-stixR*stixL)*cpsi_sq;
      stixX = stixP/(stixP- mu_sq*spsi_sq);
      
      rho1=((mu_sq-stixS)*mu_sq*spsi*cpsi)/(stixD*(mu_sq*spsi_sq-stixP));
      rho2 = (mu_sq - stixS) / stixD ;
      
      // OLD - INCORRECT WAY!
      //Byw_sq = 2.0*MU0/C* pwr *stixX*stixX *mu* fabs(cpsi)*
      // sqrt( pow( ((spsi/cpsi)-rho1*rho2) ,2) + 
      //     pow( (1+rho2*rho2*stixX)  ,2) );
      Byw_sq =  2.0*MU0/C*pwr*stixX*stixX*rho2*rho2*mu*fabs(cpsi)/
	sqrt(  pow((tan(psi)-rho1*rho2*stixX),2) + 
	       pow( (1+rho2*rho2*stixX), 2 ) );




      //printf("\n2*mu0/c: %g, pwr*stixX*stixX: %g, \n mu*fabs(cosi): %g\n",
      //       (2.0*MU0/C),(pwr*stixX*stixX),(mu*fabs(cpsi)));
      
      //printf("pow( ((spsi/cpsi)-rho1*rho2) ,2): %g,
      //        \npow( (1+rho2*rho2*stixX)  ,2): %g, sqrt: %g\n",
      //     (pow( ((spsi/cpsi)-rho1*rho2) ,2)), 
      //     pow((1+rho2*rho2*stixX), 2),  
      //      sqrt(pow(((spsi/cpsi)-rho1*rho2),2)+pow((1+rho2*rho2*stixX),2)) );
      
      //printf("\nByw_sq: %g, rho1: %g, rho2: %g, \nstixS: %g, stixB: %g\n", 
      //       Byw_sq, rho1, rho2, stixS, stixB);
      
      //printf("\nn_x: %g, n_z: %g, stixX: %g, stixD: %g, stixA: %g,
      //mu_sq: %g\n", 
      //       n_x, n_z, stixX, stixD, stixA, mu_sq);
      
      // get all RMS wave components
      
      Byw = sqrt(Byw_sq);
      Exw = fabs(C*Byw * (stixP - n_x*n_x)/(stixP*n_z)); 
      Eyw = fabs(Exw * stixD/(stixS-mu_sq));
      Ezw = fabs(Exw *n_x*n_z / (n_x*n_x - stixP));
      Bxw = fabs(Exw *stixD*n_z /C/ (stixS - mu_sq));
      Bzw = fabs((Exw *stixD *n_x) /(C*(stixX - mu_sq)));
      
      /* printf("\nExw: %g, stixD: %g, n_x: %g, prod: %g\n", 
	 Exw, stixD, n_x, (Exw*stixD*n_x) );
	 printf("\nC: %g, stixX: %g, mu_sq: %g, prod: %g\n",
	 C, stixX, mu_sq, (C*(stixX-mu_sq)) );
	 printf("div: %g, fabs: %g\n",
	 ((Exw *stixD *n_x) /(C*(stixX - mu_sq))), 
	 fabs((Exw *stixD *n_x) /(C*(stixX - mu_sq))));
	 printf("Bxw: %g, Byw: %g, Bzw: %g, \nExw: %g, Eyw: %g, Ezw: %g\n", 
	 Bxw, Byw, Bzw, Exw, Eyw, Ezw);
      */
      // Oblique integration quantities
      R1 = (Exw + Eyw)/(Bxw+Byw);
      R2 = (Exw - Eyw)/(Bxw-Byw);
      w1 = Q_EL/(2*M_EL)*(Bxw+Byw);
      w2 = Q_EL/(2*M_EL)*(Bxw-Byw);
      alpha1 = w2/w1;
      
      // printf("R1: %g, R2: %g, w1: %g, w2: %g\n", R1, R2, w1, w2);
      
      //begin MRES loop here
      for(mres=-5; mres <6; mres++) {
	
	// get parallel resonance velocity
	t1 = w*w*kz*kz;
	t2 = pow((mres*wh),2)-w*w;
	t3 = kz*kz + pow((mres*wh),2)/(pow(C*cos(alpha_lc),2));
	if(mres==0) {
	  direction = -kz/fabs(kz);
	} else {
	  direction = kz/fabs(kz) * mres/fabs(mres) ;
	}
	v_para_res = ( direction*sqrt(t1 + t2*t3) - w*kz ) / t3;
	v_tot_res = v_para_res / cos(alpha_lc); 
	E_res = E_EL*( 1/sqrt( 1-(v_tot_res*v_tot_res/(C*C)) ) -1 );

	// get starting and ending indices, +-20% energy band  
	e_starti = floor(( log10(E_res)-E_EXP_BOT-0.3)/DE_EXP);
	e_endi = ceil(( log10(E_res)-E_EXP_BOT+0.3)/DE_EXP);
	if(e_endi>NUM_E) e_endi=NUM_E;
	if(e_starti>NUM_E) e_starti=NUM_E;
	if(e_endi<0) e_endi=0;
	if(e_starti<0) e_starti=0;
	

	// begin V_TOT loop here
	for(e_toti=e_starti; e_toti < e_endi; e_toti++) {
	  v_tot = direction*v_tot_arr[e_toti];
	  v_para = v_tot * calph;
	  v_perp = fabs(v_tot * salph);
	  
	  gamma = 1 / sqrt(1 - pow((v_tot/C),2)); 
	  alpha2 = Q_EL*Ezw /(M_EL*gamma*w1*v_perp);
	  beta = kx*v_perp / wh ;
	  wtau_sq = pow((-1),(mres-1)) * w1/gamma * 
	    ( jn( (mres-1), beta ) - 
	      alpha1*jn( (mres+1) , beta ) +
	      gamma*alpha2*jn( mres , beta ) ); 
	  T1 = -wtau_sq*(1+ ( (calph*calph) / (mres*Y-1) ) );
	  
	  // Now - start analytical evaluation!!!
	  
	  if( fabs(lat)< 1e-3) {
	    
	    eta_dot = mres*wh/gamma - w - kz*v_para;

	    if(fabs(eta_dot)<10) {
	      dalpha_eq = fabs(T1/v_para)*ds/sqrt(2); 
	    } else {
	      dalpha_eq = fabs(T1/eta_dot)*sqrt(1-cos(ds*eta_dot/v_para)); 
	    }

	  } else {  
	    
	    // AA = (mres/(2.0*v_para*gamma))*dwh_ds - (kz/2.0)*dv_para_ds;
	    // doesn't the second term need  to get divided by v_para?


	    //BB = mres/v_para/gamma*(wh-dwh_ds*ds/2) - 
	    //  w/v_para - kz + kz/v_para*dv_para_ds*ds/2;

	    v_para_star = v_para - dv_para_ds*ds/2.0;
	    v_para_star_sq = v_para_star * v_para_star;

	    AA = (mres/(2.0*v_para_star*gamma))*dwh_ds* 
	      (1 + ds/(2.0*v_para_star)*dv_para_ds) - 
	      mres/(2.0*v_para_star_sq*gamma)*wh*dv_para_ds + 
	      w/(2.0*v_para_star_sq)*dv_para_ds ;


	    BB = mres/(gamma*v_para_star)*wh - 
	      mres/(gamma*v_para_star)*dwh_ds*(ds/2.0) -
	      w/v_para_star - kz;

	    //	    printf("\nAA: %g, BB: %g, v_para_star: %g\n", 
	    //	   AA, BB, v_para_star);

	    // MATLAB code
	    //aa1 = mres/2./v_para_star/gammaiA.*dwhi_ds.* ...
	    //  (1+dsiA./2./v_para_star.*dv_para_ds) - ...
	    //  mres/2./(v_para_star.^2)./gammaiA.*whiA.*dv_para_ds + ...
	    //  wi./2./(v_para_star.^2).*dv_para_ds ;
  
	    //  bb1 = mres./gammaiA./v_para_star.*whiA - ...
	    // mres./gammaiA./v_para_star.*dwhi_ds.*dsiA/2 - ...
	    // wi./v_para_star - kzi;



	    
	    Farg = (BB + 2*AA*ds) / sqrt(2*PI*fabs(AA));
	    Farg0 = BB / sqrt(2*PI*fabs(AA));  
	    
	    Fresnel(Farg, &Fs, &Fc);
	    Fresnel(Farg0, &Fs0, &Fc0);
	    
	    dFs_sq = pow((Fs - Fs0),2);
	    dFc_sq = pow((Fc - Fc0),2);
	    
	    dalpha = sqrt(PI/4/fabs(AA))*fabs(T1/v_para)*sqrt(dFs_sq+dFc_sq);
	    
	    alpha_eq_p = asin( sin(alpha_lc+dalpha)*pow(clat,3) / 
			       sqrt(slat_term) );
	    dalpha_eq = alpha_eq_p - alpha_eq;
	    
	  }
	  for(kk=0; kk<EA_SPLIT; kk++) {
	    if(direction>0) {
	      flt_time = fabs(flt_const_N[kk]/v_para);
	    } else {
	      flt_time = fabs(flt_const_S[kk]/v_para);
	    }
	    addToAlphaArr(&arr_N, &arr_S, (t+flt_time), 
			  direction*e_toti, dalpha_eq);
	  } // kk loop

	 
	} // v_para
	
      } // mres loop
      

      } // if pwr > 0.0

      end = time(NULL);
      
      // exit(0); return;
      
      next = next->link;
    }  // while loop
  }  // for i
  printf("\nL= %g precipitation took: %g sec\n",L,difftime(end,start));  
  // Make fwrite work in completely unbuffered mode to avoid errors
  setbuf(resPtrN, NULL);
  setbuf(resPtrS, NULL);
  
  for(ei=0; ei<NUM_E; ei++) {
    for(ti=0; ti<NUM_TIMES; ti++) {
      noutN=fwrite(&(arr_N[ei*NUM_TIMES+ti]), sizeof(float), 1, resPtrN);
      noutS=fwrite(&(arr_S[ei*NUM_TIMES+ti]), sizeof(float), 1, resPtrS);
      arr_N[ei*NUM_TIMES+ti] = 0.0;
      arr_S[ei*NUM_TIMES+ti] = 0.0;    
    }
  }
  
  printf("Successfully wrote elements for L=%g\n", L);
  //noutN=fwrite(arr_N, sizeof(float), (NUM_E*NUM_TIMES), resPtrN);
  //noutS=fwrite(arr_S, sizeof(float), (NUM_E*NUM_TIMES), resPtrS);
  //printf("Wrote %d elements N, and %d elements S\n", noutN, noutS);
  free(arr_N);
  free(arr_S);
  arr_N = NULL;
  arr_S = NULL;  
  

}







/*
 * FUNCTION: getFltConst
 * ---------------------
 * This function calculates the 'flight-time-const' for a particular 
 * latitude,  and returns the constant for a particle flying adiabatically
 * from that latitude, to the Northern hemisphere, and a similar constant 
 * to the Southern hemisphere.  This constant is then multiplied by 
 * 1/v_tot for the particular particle, which gives the total flight
 * time in seconds to the appropriate hemisphere.
 *
 */
void getFltConst(double L, double lat, double alpha_eq, 
		 double *flt_const_N, double *flt_const_S)
{
  double sin_alpha_eq_sq, x=0.0, dx, endx, walt_tau, xterm, I=0.0;
  int num_div=10000, n;

  sin_alpha_eq_sq = pow(sin(alpha_eq),2);
  endx = sin(lat*D2R);
  dx = endx/num_div;
  
  // Evaluate Walt's flight-time constant, from equator to mirror pt.
  walt_tau = 0.02925/R_E*C*(1 - 0.4635*pow(sin_alpha_eq_sq,0.375));

  // Evaluate flight-time integral from equator to lat  
  //for(x=0; x/endx < 1; x+=dx) {
  for(n=0; n<num_div; n++) {
    xterm = sqrt(1+3*x*x);
    I += dx*xterm/sqrt( 1 - (sin_alpha_eq_sq/pow((1-x*x),3))*xterm);
    x+=dx;
  }
  
  *flt_const_N = L*R_E*(walt_tau - I);
  *flt_const_S = L*R_E*(walt_tau + I);

}








/*
 * FUNCTION: Fresnel
 * -----------------
 * This function calculates the sine and cose Fresnel integrals and 
 * returns the values of either.  The integrals are defined as:
 *
 *    /x			     /x
 *    | sin( pi/2*t^2 ) dt           | cos( pi/2*t^2 ) dt
 *    /0			     /0
 *
 * This uses the methodology of  Klaus D. Mielenz, "Computation of
 * Fresnel Integrals II" Journal of Research of the National 
 * Institute of Standards and Technology, 105, 589 (2000)
 *
 */
void Fresnel(double x0, double *FS, double *FC)
{

  double fn[12] = { 0.318309844,
		    9.34626e-8 , 
		    -0.09676631, 
		    0.000606222, 
		    0.325539361, 
		    0.325206461, 
		    -7.450551455,
		    32.20380908, 
		    -78.8035274, 
		    118.5343352, 
		    -102.4339798,
		    39.06207702 } ;

  double gn[12] = { 0.0 ,
		    0.101321519, 
		    -4.07292e-5, 
		    -0.152068115, 
		    -0.046292605, 
		    1.622793598, 
		    -5.199186089, 
		    7.477942354, 
		    -0.695291507, 
		    -15.10996796,
		    22.28401942, 
		    -10.89968491 };
  
  double cn[12] = { 1.0 ,
		    -0.24674011002723,
		    0.02818550087789 ,
		    -0.00160488313564 ,
		    5.407413381408390e-05 ,
		    -1.200097255860028e-06,
		    1.884349911527268e-08,
		    -2.202276925445466e-10, 
		    1.989685792418021e-12,
		    -1.430918973171519e-14,
		    8.384729705118549e-17,
		    -4.079981449233875e-19 } ;

  double sn[12] = {    0.52359877559830,
		       -0.09228058535804,
		       0.00724478420420,
		       -3.121169423545791e-04,
		       8.444272883545251e-06,
		       -1.564714450092211e-07,
		       2.108212193321454e-09,
		       -2.157430680584343e-11,
		       1.733410208887483e-13,
		       -1.122324478798395e-15,
		       5.980053239210401e-18,
		       -2.667871362841397e-20 };

  double xpow, x_sq, fx=0, gx=0, x;
  int n;

  x = fabs(x0);
  *FS = 0.0;
  *FC = 0.0;
  x_sq = x*x;



  if(x<=1.6) {

    *FS =   sn[0]*pow(x,3) +	// it takes longer to write this out
	    sn[1]*pow(x,7) +	// but we save valuable CPU cycles!
	    sn[2]*pow(x,11) + 
	    sn[3]*pow(x,15) + 
	    sn[4]*pow(x,19) + 
	    sn[5]*pow(x,23) + 
	    sn[6]*pow(x,27) + 
	    sn[7]*pow(x,31) + 
	    sn[8]*pow(x,35) + 
	    sn[9]*pow(x,39) + 
	    sn[10]*pow(x,43) + 
	    sn[11]*pow(x,47)  ; 

    *FC =   cn[0]*x +	
	    cn[1]*pow(x,5) +	
	    cn[2]*pow(x,9) + 
	    cn[3]*pow(x,13) + 
	    cn[4]*pow(x,17) + 
	    cn[5]*pow(x,21) + 
	    cn[6]*pow(x,25) + 
	    cn[7]*pow(x,29) + 
	    cn[8]*pow(x,33) + 
	    cn[9]*pow(x,37) + 
	    cn[10]*pow(x,41) + 
	    cn[11]*pow(x,45)  ; 

  } else {
      
    
    for(n=0; n<=11; n++) {
      xpow = pow(x, (-2*n-1) );
      fx += fn[n]*xpow;
      gx += gn[n]*xpow;
    }     
    *FC = 0.5 + fx*sin(PI/2*x_sq) - gx*cos(PI/2*x_sq);
    *FS = 0.5 - gx*sin(PI/2*x_sq) - fx*cos(PI/2*x_sq);      
  }
 
  if(x0<0) {
    *FC = -(*FC);
    *FS = -(*FS);
  }
}











/*
 * FUNCTION: addToAlphaArr
 * -----------------------
 * This function takes a certain arrival-time-at-the-ionoshpere, and 
 * particle velocity pair, and enters the corresponding change in 
 * equatorial pitch angle in the correct place in the table.  If the 
 * array has no yet been defined, this functino defines and initializes
 * the array.
 *
 */
void addToAlphaArr(float **arr_N, float **arr_S, 
		   double iono_time, int e_toti, double dalpha)
{
  int veli, timei, ei, ti, direction; 


  // If arrays don't exist, create + initialise
  if(*arr_N==NULL || *arr_S==NULL) {
    // NUM_E = ceil(C / DV_TOT); // used to determine NUM_VELS
    NUM_TIMES = ceil( RES_FINT / RES_DT );   
    
    *arr_N=(float *)malloc(NUM_E*NUM_TIMES*sizeof(float));
    *arr_S=(float *)malloc(NUM_E*NUM_TIMES*sizeof(float));

    printf("\nNUM_E: %d, NUM_TIMES: %d\n", NUM_E, NUM_TIMES);

    if(*arr_N==NULL || *arr_S==NULL) {
      printf("\nProblem allocating arr_N or arr_S\a\n");
      exit(0);
    }
    
    
    for(ei=0; ei<NUM_E; ei++) {
      for(ti=0; ti<NUM_TIMES; ti++) {
	*(*arr_N + (ei*NUM_TIMES+ti)) = 0.0;
	*(*arr_S + (ei*NUM_TIMES+ti)) = 0.0;    
      }
    }  
    printf("\nMade it past initialization loop\n"); 
  } // end of creating + initializing
  
  
  // get Correct index into array
  timei = floor(iono_time / RES_DT);
  
  // Enter dalpha^2 into array
  if(timei < NUM_TIMES) {  // check that arrival time < final time
    if(e_toti>= 0 && e_toti < NUM_E ) {
      *(*arr_N+(e_toti*NUM_TIMES+timei)) += ((dalpha*dalpha)/EA_SPLIT)/MULT;
    }
    if(e_toti<= 0 && e_toti > -NUM_E ) {    
      *(*arr_S + (abs(e_toti)*NUM_TIMES + timei)) += 
	((dalpha*dalpha)/EA_SPLIT)/MULT;
    }
  }
}









/*
 * FUNCTION: freeArr
 * -----------------
 * This function goes through the array of f-t cells and simply frees
 * the memory sequentially.
 *
 */
void freeArr( cellT *lat_arr[][NUM_TARGS] )
{
  int i, j;
  cellT *this, *next;
  double lat;

  for(j=0; j<NUM_TARGS; j++) {
    for(i=0; i<NUMLATS; i++) {
      this = lat_arr[i][j];
      lat = EALimS + EAIncr*i;
      
      while(this != NULL) {
	
	if(this->link != NULL) {
	  next = this->link;
	  free(this);
	  this = next;
	} else {
	  free(this);
	  this = NULL;
	} // if - else

      }  // while  
    }  // for i<NUMLATS
  }  // for j<NUM_TARGS
}


















